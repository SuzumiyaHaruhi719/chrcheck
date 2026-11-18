"""Streaming parsers for BAM, SAM, VCF, BED, GTF/GFF3, FAI, DICT, and FASTA."""

import gzip
import os
import re
import struct
import zlib
from typing import Dict, IO, List, Optional, Tuple

from chrcheck.assemblies import (
    classify_contig,
    infer_assemblies_from_coordinates,
    infer_assemblies_from_header,
    is_alt_or_decoy,
)
from chrcheck.models import ContigInfo, FileProfile

_VCF_CONTIG_RE = re.compile(r"^##contig=<(.*)>$")
_SAM_SQ_TAG_RE = re.compile(r"^([A-Za-z0-9]{2}):(.*)$")
_KEY_VAL_RE = re.compile(r'([A-Za-z0-9_]+)=("[^"]*"|[^,]*)')


def detect_format(path: str) -> str:
    """Infer genomics file format from extension (and magic bytes for BAM)."""
    lower = path.lower()
    if lower.endswith(".bam"):
        return "bam"
    if lower.endswith(".sam"):
        return "sam"
    if lower.endswith((".fa.fai", ".fasta.fai", ".fna.fai", ".fai")):
        return "fai"
    if lower.endswith(".dict"):
        return "dict"
    if lower.endswith(".chrom.sizes"):
        return "chrom_sizes"

    stem = lower[:-3] if lower.endswith(".gz") else lower
    if stem.endswith((".vcf", ".bcf")):
        return "vcf"
    if stem.endswith((".bed", ".narrowpeak", ".broadpeak", ".bedgraph")):
        return "bed"
    if stem.endswith((".gtf", ".gff", ".gff3")):
        return "gtf"
    if stem.endswith((".fa", ".fasta", ".fna")):
        return "fasta"

    # Fallback magic inspection for binary BAM
    try:
        with open(path, "rb") as fh:
            head = fh.read(4)
        if head[:2] == b"\x1f\x8b" or head == b"BAM\x01":
            return "bam"
    except OSError:
        pass
    raise ValueError(f"Unsupported or unrecognized genomics file format: {path}")


class _BgzfHeaderReader:
    """Read just enough BGZF blocks to extract the uncompressed BAM header."""

    def __init__(self, fh: IO[bytes]) -> None:
        self._fh = fh
        self._buf = bytearray()
        self._pos = 0
        self._eof = False

    def _fill_one_block(self) -> bool:
        if self._eof:
            return False
        hdr = self._fh.read(12)
        if not hdr:
            self._eof = True
            return False
        if len(hdr) < 12:
            raise ValueError("Truncated gzip/BGZF header in BAM file")

        if hdr[:4] == b"BAM\x01" and not self._buf:
            # Uncompressed BAM stream (rare, e.g. samtools -u)
            self._buf.extend(hdr)
            self._buf.extend(self._fh.read(65536))
            return True

        if hdr[0] != 0x1F or hdr[1] != 0x8B or hdr[2] != 8:
            raise ValueError("Invalid BAM/BGZF header magic")

        flg = hdr[3]
        if flg & 0x04:
            xlen = struct.unpack("<H", hdr[10:12])[0]
            extra = self._fh.read(xlen)
            if len(extra) < xlen:
                raise ValueError("Truncated BGZF extra subfields")
            bsize: Optional[int] = None
            idx = 0
            while idx + 4 <= len(extra):
                si1, si2, slen = struct.unpack("<BBH", extra[idx : idx + 4])
                if si1 == ord("B") and si2 == ord("C") and slen == 2 and idx + 6 <= len(extra):
                    bsize = struct.unpack("<H", extra[idx + 4 : idx + 6])[0]
                    break
                idx += 4 + slen
            if bsize is not None:
                # Compressed deflate payload ends 8 bytes before block end (CRC32 + ISIZE)
                cdata_len = bsize - xlen - 19
                if cdata_len < 0:
                    raise ValueError("Corrupt BGZF block size")
                cdata = self._fh.read(cdata_len)
                trailer = self._fh.read(8)
                if len(cdata) < cdata_len or len(trailer) < 8:
                    raise ValueError("Truncated BGZF block payload")
                chunk = zlib.decompress(cdata, -zlib.MAX_WBITS)
                self._buf.extend(chunk)
                return True

        # Non-BGZF plain gzipped BAM fallback: decompress prefix via streaming zlib
        self._fh.seek(0)
        decomp = zlib.decompressobj(16 + zlib.MAX_WBITS)
        raw = self._fh.read(65536)
        while raw:
            out = decomp.decompress(raw)
            if out:
                self._buf.extend(out)
                return True
            if decomp.eof:
                break
            raw = self._fh.read(65536)
        self._eof = True
        return False

    def read_exact(self, n: int) -> bytes:
        while len(self._buf) - self._pos < n:
            if not self._fill_one_block():
                raise ValueError("Unexpected EOF while reading BAM header")
        out = bytes(self._buf[self._pos : self._pos + n])
        self._pos += n
        return out


def _open_text_stream(path: str) -> IO[str]:
    with open(path, "rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _record_observed(
    profile: FileProfile,
    chrom: str,
    end_pos: Optional[int],
) -> None:
    if chrom not in profile.observed_counts:
        profile.observed_contigs.append(chrom)
        profile.observed_counts[chrom] = 0
    profile.observed_counts[chrom] += 1
    if end_pos is not None:
        prev = profile.observed_max_pos.get(chrom, 0)
        if end_pos > prev:
            profile.observed_max_pos[chrom] = end_pos


def _parse_sam_sq_line(line: str, profile: FileProfile) -> None:
    parts = line.rstrip("\r\n").split("\t")
    if not parts:
        return
    if parts[0] == "@SQ":
        name: Optional[str] = None
        length: Optional[int] = None
        asm: Optional[str] = None
        md5: Optional[str] = None
        for tag in parts[1:]:
            m = _SAM_SQ_TAG_RE.match(tag)
            if not m:
                continue
            k, v = m.group(1), m.group(2)
            if k == "SN":
                name = v
            elif k == "LN":
                try:
                    length = int(v)
                except ValueError:
                    length = None
            elif k == "AS":
                asm = v
                if v:
                    profile.reference_hints.append(v)
            elif k == "M5":
                md5 = v
        if name:
            canon, style, _ = classify_contig(name)
            info = ContigInfo(
                name=name,
                length=length,
                canonical=canon,
                style=style,
                assembly_hint=asm,
                md5=md5,
            )
            profile.header_contigs.append(info)
            profile.header_order.append(name)
            if length is not None:
                profile.header_lengths[name] = length
            if is_alt_or_decoy(name):
                profile.has_alt_or_decoy = True
    elif parts[0] in ("@HD", "@PG", "@CO"):
        for field_txt in parts[1:]:
            if "assembly" in field_txt.lower() or "grch" in field_txt.lower() or "hg" in field_txt.lower():
                profile.reference_hints.append(field_txt)


def _finalize_profile(profile: FileProfile) -> FileProfile:
    styles_seen: List[str] = []
    mito_candidates: List[str] = []

    names_for_style = profile.header_order + [
        c for c in profile.observed_contigs if c not in profile.header_lengths and c not in profile.header_order
    ]
    for contig in names_for_style:
        canon, style, _ = classify_contig(contig)
        if canon is not None and style in ("ucsc", "ensembl", "refseq"):
            if style not in styles_seen:
                styles_seen.append(style)
        elif style == "broken_mito":
            if style not in styles_seen:
                styles_seen.append(style)
        if canon == "MT" and contig not in mito_candidates:
            mito_candidates.append(contig)

    # Also track observed styles separately in case records disagree with headers
    for contig in profile.observed_contigs:
        canon, style, _ = classify_contig(contig)
        if canon is not None and style in ("ucsc", "ensembl", "refseq", "broken_mito"):
            if style not in profile.observed_styles:
                profile.observed_styles.append(style)
        if canon == "MT" and contig not in mito_candidates:
            mito_candidates.append(contig)

    if len(styles_seen) == 1:
        profile.naming_style = styles_seen[0]
    elif len(styles_seen) > 1:
        profile.naming_style = "mixed"
    else:
        profile.naming_style = "unknown"

    if len(mito_candidates) == 1:
        profile.mito_style = mito_candidates[0]
    elif len(mito_candidates) > 1:
        profile.mito_style = "/".join(mito_candidates)

    if profile.header_lengths:
        profile.inferred_assemblies = infer_assemblies_from_header(
            profile.header_lengths,
            profile.reference_hints,
        )
    else:
        profile.inferred_assemblies = infer_assemblies_from_coordinates(
            profile.observed_max_pos,
            profile.reference_hints,
        )
    return profile


def _parse_bam(path: str, profile: FileProfile) -> FileProfile:
    with open(path, "rb") as fh:
        reader = _BgzfHeaderReader(fh)
        magic = reader.read_exact(4)
        if magic != b"BAM\x01":
            raise ValueError(f"Not a valid BAM file (bad magic {magic!r}): {path}")
        l_text = struct.unpack("<I", reader.read_exact(4))[0]
        header_bytes = reader.read_exact(l_text) if l_text > 0 else b""
        header_text = header_bytes.decode("utf-8", errors="replace")
        for line in header_text.splitlines():
            if line.startswith("@"):
                _parse_sam_sq_line(line, profile)

        n_ref = struct.unpack("<I", reader.read_exact(4))[0]
        binary_names: List[Tuple[str, int]] = []
        for _ in range(n_ref):
            l_name = struct.unpack("<I", reader.read_exact(4))[0]
            raw_name = reader.read_exact(l_name)
            name = raw_name.rstrip(b"\x00").decode("utf-8", errors="replace")
            l_ref = struct.unpack("<I", reader.read_exact(4))[0]
            binary_names.append((name, l_ref))

        if not profile.header_order and binary_names:
            for name, length in binary_names:
                canon, style, _ = classify_contig(name)
                profile.header_contigs.append(
                    ContigInfo(name=name, length=length, canonical=canon, style=style)
                )
                profile.header_order.append(name)
                profile.header_lengths[name] = length
                if is_alt_or_decoy(name):
                    profile.has_alt_or_decoy = True
        elif binary_names:
            # Verify @SQ headers match binary dictionary
            sq_pairs = [(c.name, c.length) for c in profile.header_contigs]
            if [(n, l) for n, l in binary_names] != sq_pairs:
                profile.parse_warnings.append(
                    "BAM text @SQ header differs from binary reference dictionary"
                )
    return _finalize_profile(profile)


def _parse_sam(path: str, profile: FileProfile, sample_records: Optional[int]) -> FileProfile:
    with _open_text_stream(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("@"):
                _parse_sam_sq_line(line, profile)
                continue
            if sample_records is not None and profile.total_records_scanned >= sample_records:
                profile.truncated_scan = True
                break
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 4:
                continue
            rname = cols[2]
            if rname == "*":
                continue
            profile.total_records_scanned += 1
            try:
                pos = int(cols[3])
            except ValueError:
                profile.malformed_intervals.append(f"SAM invalid POS on {rname}: {cols[3]}")
                pos = 0
            _record_observed(profile, rname, pos if pos > 0 else None)
    return _finalize_profile(profile)


def _parse_fai_or_sizes(path: str, profile: FileProfile) -> FileProfile:
    with _open_text_stream(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            cols = stripped.split("\t") if "\t" in stripped else stripped.split()
            if len(cols) < 2:
                profile.parse_warnings.append(f"Line {lineno}: expected at least 2 columns")
                continue
            name = cols[0]
            try:
                length = int(cols[1])
            except ValueError:
                profile.parse_warnings.append(f"Line {lineno}: non-integer contig length {cols[1]!r}")
                continue
            canon, style, _ = classify_contig(name)
            profile.header_contigs.append(
                ContigInfo(name=name, length=length, canonical=canon, style=style)
            )
            profile.header_order.append(name)
            profile.header_lengths[name] = length
            if is_alt_or_decoy(name):
                profile.has_alt_or_decoy = True
    return _finalize_profile(profile)


def _parse_dict(path: str, profile: FileProfile) -> FileProfile:
    with _open_text_stream(path) as fh:
        for line in fh:
            if line.startswith("@"):
                _parse_sam_sq_line(line, profile)
    return _finalize_profile(profile)


def _parse_vcf(path: str, profile: FileProfile, sample_records: Optional[int]) -> FileProfile:
    with _open_text_stream(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\r\n")
            if not raw:
                continue
            if raw.startswith("##"):
                if raw.startswith("##contig=<"):
                    m = _VCF_CONTIG_RE.match(raw)
                    if m:
                        attrs: Dict[str, str] = {}
                        for k, v in _KEY_VAL_RE.findall(m.group(1)):
                            attrs[k] = v.strip('"')
                        cid = attrs.get("ID")
                        if cid:
                            length: Optional[int] = None
                            if "length" in attrs:
                                try:
                                    length = int(attrs["length"])
                                except ValueError:
                                    length = None
                            asm = attrs.get("assembly")
                            if asm:
                                profile.reference_hints.append(asm)
                            canon, style, _ = classify_contig(cid)
                            profile.header_contigs.append(
                                ContigInfo(
                                    name=cid,
                                    length=length,
                                    canonical=canon,
                                    style=style,
                                    assembly_hint=asm,
                                )
                            )
                            profile.header_order.append(cid)
                            if length is not None:
                                profile.header_lengths[cid] = length
                            if is_alt_or_decoy(cid):
                                profile.has_alt_or_decoy = True
                elif raw.startswith("##reference=") or raw.startswith("##assembly="):
                    profile.reference_hints.append(raw.split("=", 1)[1])
                continue
            if raw.startswith("#CHROM"):
                continue
            if raw.startswith("#"):
                continue

            if sample_records is not None and profile.total_records_scanned >= sample_records:
                profile.truncated_scan = True
                break

            cols = raw.split("\t")
            if len(cols) < 4:
                profile.parse_warnings.append(f"VCF line {lineno}: fewer than 4 columns")
                continue
            chrom = cols[0]
            profile.total_records_scanned += 1
            try:
                pos = int(cols[1])
            except ValueError:
                profile.malformed_intervals.append(f"line {lineno}: invalid POS {cols[1]!r}")
                _record_observed(profile, chrom, None)
                continue
            if pos < 1:
                profile.malformed_intervals.append(f"line {lineno}: VCF POS {pos} < 1 on {chrom}")
            ref_len = len(cols[3]) if cols[3] and cols[3] != "." else 1
            end_pos = pos + max(ref_len - 1, 0)
            _record_observed(profile, chrom, end_pos)
    return _finalize_profile(profile)


def _parse_bed(path: str, profile: FileProfile, sample_records: Optional[int]) -> FileProfile:
    with _open_text_stream(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\r\n")
            stripped = raw.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith(("track ", "browser ")):
                continue
            if sample_records is not None and profile.total_records_scanned >= sample_records:
                profile.truncated_scan = True
                break
            cols = raw.split("\t") if "\t" in raw else raw.split()
            if len(cols) < 3:
                profile.parse_warnings.append(f"BED line {lineno}: fewer than 3 columns")
                continue
            chrom = cols[0]
            profile.total_records_scanned += 1
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                profile.malformed_intervals.append(f"line {lineno}: non-integer BED interval on {chrom}")
                _record_observed(profile, chrom, None)
                continue
            if start < 0 or start > end:
                profile.malformed_intervals.append(
                    f"line {lineno}: invalid BED interval {chrom}:{start}-{end}"
                )
            elif start == end:
                profile.zero_length_bed_count += 1
                if len(profile.zero_length_examples) < 5:
                    profile.zero_length_examples.append(f"{chrom}:{start}-{end} (line {lineno})")
            if is_alt_or_decoy(chrom):
                profile.has_alt_or_decoy = True
            _record_observed(profile, chrom, end)
    return _finalize_profile(profile)


def _parse_gtf(path: str, profile: FileProfile, sample_records: Optional[int]) -> FileProfile:
    with _open_text_stream(path) as fh:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\r\n")
            stripped = raw.strip()
            if not stripped:
                continue
            if stripped == "##FASTA":
                break
            if stripped.startswith("##sequence-region"):
                parts = stripped.split()
                if len(parts) >= 4:
                    seqid = parts[1]
                    try:
                        end = int(parts[3])
                    except ValueError:
                        end = None
                    canon, style, _ = classify_contig(seqid)
                    profile.header_contigs.append(
                        ContigInfo(name=seqid, length=end, canonical=canon, style=style)
                    )
                    profile.header_order.append(seqid)
                    if end is not None:
                        profile.header_lengths[seqid] = end
                    if is_alt_or_decoy(seqid):
                        profile.has_alt_or_decoy = True
                continue
            if stripped.startswith("#"):
                if "genome-build" in stripped.lower():
                    profile.reference_hints.append(stripped)
                continue
            if sample_records is not None and profile.total_records_scanned >= sample_records:
                profile.truncated_scan = True
                break
            cols = raw.split("\t")
            if len(cols) < 5:
                profile.parse_warnings.append(f"GTF/GFF line {lineno}: fewer than 5 columns")
                continue
            seqid = cols[0]
            profile.total_records_scanned += 1
            try:
                start = int(cols[3])
                end = int(cols[4])
            except ValueError:
                profile.malformed_intervals.append(f"line {lineno}: non-integer GTF/GFF coords on {seqid}")
                _record_observed(profile, seqid, None)
                continue
            if start < 1 or start > end:
                profile.malformed_intervals.append(
                    f"line {lineno}: invalid GTF/GFF interval {seqid}:{start}-{end}"
                )
            if is_alt_or_decoy(seqid):
                profile.has_alt_or_decoy = True
            _record_observed(profile, seqid, end)
    return _finalize_profile(profile)


def _parse_fasta(path: str, profile: FileProfile) -> FileProfile:
    current_name: Optional[str] = None
    current_len = 0
    with _open_text_stream(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current_name is not None:
                    canon, style, _ = classify_contig(current_name)
                    profile.header_contigs.append(
                        ContigInfo(name=current_name, length=current_len, canonical=canon, style=style)
                    )
                    profile.header_order.append(current_name)
                    profile.header_lengths[current_name] = current_len
                    if is_alt_or_decoy(current_name):
                        profile.has_alt_or_decoy = True
                current_name = stripped[1:].split()[0]
                current_len = 0
            else:
                current_len += len(stripped)
        if current_name is not None:
            canon, style, _ = classify_contig(current_name)
            profile.header_contigs.append(
                ContigInfo(name=current_name, length=current_len, canonical=canon, style=style)
            )
            profile.header_order.append(current_name)
            profile.header_lengths[current_name] = current_len
            if is_alt_or_decoy(current_name):
                profile.has_alt_or_decoy = True
    return _finalize_profile(profile)


def parse_file(path: str, sample_records: Optional[int] = 50000) -> FileProfile:
    """Parse a genomics file into a FileProfile."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")
    fmt = detect_format(path)
    profile = FileProfile(path=path, fmt=fmt)

    if fmt == "bam":
        return _parse_bam(path, profile)
    if fmt == "sam":
        return _parse_sam(path, profile, sample_records)
    if fmt in ("fai", "chrom_sizes"):
        return _parse_fai_or_sizes(path, profile)
    if fmt == "dict":
        return _parse_dict(path, profile)
    if fmt == "vcf":
        return _parse_vcf(path, profile, sample_records)
    if fmt == "bed":
        return _parse_bed(path, profile, sample_records)
    if fmt == "gtf":
        return _parse_gtf(path, profile, sample_records)
    if fmt == "fasta":
        return _parse_fasta(path, profile)
    raise ValueError(f"Unsupported format {fmt!r} for {path}")
