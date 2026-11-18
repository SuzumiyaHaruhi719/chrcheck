"""Safe contig rename-map generation and column-1 BED/GTF/GFF3 rewriting."""

import gzip
from typing import List, Optional, Tuple

from chrcheck.assemblies import (
    builds_compatible,
    classify_contig,
    map_contig_to_style,
)
from chrcheck.checker import check_profiles
from chrcheck.models import FileProfile
from chrcheck.parsers import _open_text_stream


class UnsafeBuildMappingError(ValueError):
    """Raised when attempting to emit a rename map across incompatible builds."""


def _canonical_to_contig_map(profile: FileProfile) -> dict:
    mapping = {}
    for contig in profile.all_contigs_ordered():
        canon, _, _ = classify_contig(contig)
        if canon is not None and canon not in mapping:
            mapping[canon] = contig
    return mapping


def _infer_target_style(target: FileProfile) -> str:
    if target.naming_style in ("ucsc", "ensembl", "refseq"):
        return target.naming_style
    for contig in target.all_contigs_ordered():
        canon, style, _ = classify_contig(contig)
        if canon is not None and style in ("ucsc", "ensembl", "refseq"):
            return style
    raise ValueError(f"Could not determine chromosome naming style of target {target.path}")


def build_chr_map(
    source: FileProfile,
    target: Optional[FileProfile] = None,
    target_style: Optional[str] = None,
    force: bool = False,
) -> List[Tuple[str, str]]:
    """Return list of (old_contig, new_contig) pairs suitable for bcftools --rename-chrs."""
    if target is not None:
        report = check_profiles([target, source])
        blocking = [f for f in report.findings if f.code in ("CHR-001", "CHR-002")]
        if blocking and not force:
            raise UnsafeBuildMappingError(
                f"{blocking[0].code}: {blocking[0].message} Use --force to override."
            )
        if not builds_compatible(target.inferred_assemblies, source.inferred_assemblies) and not force:
            raise UnsafeBuildMappingError(
                "Incompatible genome builds detected; coordinate liftover is required."
            )
        resolved_style = target_style or _infer_target_style(target)
        target_build = target.inferred_assemblies[0] if len(target.inferred_assemblies) == 1 else None
        target_mito = target.mito_style
        target_contigs = _canonical_to_contig_map(target)
    else:
        if not target_style:
            raise ValueError("Either target FileProfile or target_style must be provided")
        resolved_style = target_style
        target_build = source.inferred_assemblies[0] if len(source.inferred_assemblies) == 1 else None
        target_mito = "chrM" if resolved_style == "ucsc" else ("MT" if resolved_style == "ensembl" else None)
        target_contigs = None

    pairs: List[Tuple[str, str]] = []
    seen = set()
    for contig in source.all_contigs_ordered():
        if contig in seen:
            continue
        seen.add(contig)
        mapped = map_contig_to_style(
            contig,
            target_style=resolved_style,
            target_build=target_build,
            target_mito=target_mito,
            target_contigs=target_contigs,
        )
        if mapped and mapped != contig:
            pairs.append((contig, mapped))
    return pairs


def format_chr_map(pairs: List[Tuple[str, str]]) -> str:
    """Format (old, new) pairs as tab-separated lines for `bcftools annotate --rename-chrs`."""
    return "".join(f"{old}\t{new}\n" for old, new in pairs)


def rewrite_bed_or_gff(
    input_path: str,
    output_path: str,
    source_profile: FileProfile,
    target: Optional[FileProfile] = None,
    target_style: Optional[str] = None,
    force: bool = False,
) -> Tuple[int, int]:
    """Safely rewrite column 1 (and ##sequence-region) in BED/GTF/GFF3 without altering gene names.

    Returns (rows_rewritten, rows_total).
    """
    pairs = dict(build_chr_map(source_profile, target=target, target_style=target_style, force=force))
    rewritten = 0
    total = 0

    out_ctx = (
        gzip.open(output_path, "wt", encoding="utf-8")
        if output_path.lower().endswith(".gz")
        else open(output_path, "wt", encoding="utf-8")
    )
    with _open_text_stream(input_path) as fin, out_ctx as fout:
        in_fasta = False
        for line in fin:
            if in_fasta:
                fout.write(line)
                continue
            stripped = line.strip()
            if stripped == "##FASTA":
                in_fasta = True
                fout.write(line)
                continue
            if stripped.startswith("##sequence-region"):
                parts = line.rstrip("\r\n").split()
                if len(parts) >= 4 and parts[1] in pairs:
                    parts[1] = pairs[parts[1]]
                    fout.write(" ".join(parts) + "\n")
                else:
                    fout.write(line)
                continue
            if not stripped or stripped.startswith("#") or stripped.startswith(("track ", "browser ")):
                fout.write(line)
                continue

            total += 1
            if "\t" in line:
                col0, rest = line.split("\t", 1)
                if col0 in pairs:
                    rewritten += 1
                    fout.write(f"{pairs[col0]}\t{rest}")
                else:
                    fout.write(line)
            else:
                parts = line.rstrip("\r\n").split(None, 1)
                if parts and parts[0] in pairs:
                    rewritten += 1
                    suffix = f"\t{parts[1]}" if len(parts) > 1 else ""
                    fout.write(f"{pairs[parts[0]]}{suffix}\n")
                else:
                    fout.write(line)
    return rewritten, total
