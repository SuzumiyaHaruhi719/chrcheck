"""Unit tests for chrcheck."""

import gzip
import io
import json
import os
import struct
import tempfile
import unittest
import zlib
from contextlib import redirect_stderr, redirect_stdout

from chrcheck.assemblies import ASSEMBLIES, classify_contig
from chrcheck.checker import check_profiles
from chrcheck.cli import main
from chrcheck.harmonize import UnsafeBuildMappingError, build_chr_map
from chrcheck.parsers import parse_file


def _pack_bgzf_block(raw_chunk: bytes) -> bytes:
    """Wrap one uncompressed chunk (<= 65536 bytes) in a valid BGZF block."""
    comp = zlib.compressobj(zlib.Z_DEFAULT_COMPRESSION, zlib.DEFLATED, -zlib.MAX_WBITS)
    cdata = comp.compress(raw_chunk) + comp.flush()
    bsize = 18 + 8 + len(cdata) - 1  # 12 hdr + 6 BC subfield + cdata + 8 trailer - 1
    header = struct.pack(
        "<BBBBIBBHBBHH",
        0x1F,
        0x8B,
        8,
        0x04,
        0,
        0,
        255,
        6,
        ord("B"),
        ord("C"),
        2,
        bsize,
    )
    crc = zlib.crc32(raw_chunk) & 0xFFFFFFFF
    trailer = struct.pack("<II", crc, len(raw_chunk) & 0xFFFFFFFF)
    return header + cdata + trailer


def _write_test_bam(path: str, contigs, as_tag: str = ""):
    """Write a minimal valid BGZF-compressed binary BAM file with SAM @SQ + n_ref."""
    sq_lines = []
    for name, length in contigs:
        tag = f"\tAS:{as_tag}" if as_tag else ""
        sq_lines.append(f"@SQ\tSN:{name}\tLN:{length}{tag}")
    header_text = ("@HD\tVN:1.6\tSO:coordinate\n" + "\n".join(sq_lines) + "\n").encode("utf-8")

    payload = bytearray()
    payload.extend(b"BAM\x01")
    payload.extend(struct.pack("<I", len(header_text)))
    payload.extend(header_text)
    payload.extend(struct.pack("<I", len(contigs)))
    for name, length in contigs:
        nb = name.encode("utf-8") + b"\x00"
        payload.extend(struct.pack("<I", len(nb)))
        payload.extend(nb)
        payload.extend(struct.pack("<I", length))

    eof_block = _pack_bgzf_block(b"")
    with open(path, "wb") as fh:
        # Split into 2 blocks if large, or 1 block + EOF block
        raw_bytes = bytes(payload)
        for idx in range(0, len(raw_bytes), 32768):
            fh.write(_pack_bgzf_block(raw_bytes[idx : idx + 32768]))
        fh.write(eof_block)


def _human_contigs(build: str, style: str = "ucsc", subset=None):
    spec = ASSEMBLIES[build]
    order = subset or ([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
    out = []
    for canon in order:
        if style == "ucsc":
            name = spec.ucsc_names[canon]
        elif style == "ensembl":
            name = spec.ensembl_names[canon]
        elif style == "refseq":
            name = spec.refseq_names[canon]
        else:
            raise ValueError(style)
        out.append((name, spec.lengths[canon]))
    return out


class ChrcheckTests(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.root = self._tmp.name

    def tearDown(self) -> None:
        self._tmp.cleanup()

    def _p(self, name: str) -> str:
        return os.path.join(self.root, name)

    def test_classify_contig_styles(self) -> None:
        self.assertEqual(classify_contig("chr1"), ("1", "ucsc", None))
        self.assertEqual(classify_contig("1"), ("1", "ensembl", None))
        self.assertEqual(classify_contig("chrM"), ("MT", "ucsc", None))
        self.assertEqual(classify_contig("MT"), ("MT", "ensembl", None))
        self.assertEqual(classify_contig("M"), ("MT", "broken_mito", None))
        self.assertEqual(classify_contig("NC_000001.11"), ("1", "refseq", "GRCh38"))
        self.assertEqual(classify_contig("NC_060925.1"), ("1", "refseq", "T2T-CHM13v2.0"))
        self.assertEqual(classify_contig("NC_000067.7"), ("1", "refseq", "GRCm39"))

    def test_bam_bgzf_header_and_compatible_check(self) -> None:
        bam = self._p("sample.bam")
        _write_test_bam(bam, _human_contigs("GRCh38", "ucsc"))
        bed = self._p("targets.bed")
        with open(bed, "wt", encoding="utf-8") as fh:
            fh.write("chr1\t1000\t2000\tEGFR\nchr2\t500\t1500\tTP53\nchrM\t10\t200\tMT-ND1\n")

        p_bam = parse_file(bam)
        p_bed = parse_file(bed)
        self.assertEqual(p_bam.fmt, "bam")
        self.assertEqual(p_bam.inferred_assemblies, ["GRCh38"])
        self.assertEqual(p_bam.naming_style, "ucsc")
        self.assertEqual(p_bam.mito_style, "chrM")

        rep = check_profiles([p_bam, p_bed])
        self.assertTrue(rep.compatible)
        self.assertEqual(rep.exit_code, 0)
        self.assertEqual([f.code for f in rep.findings], ["CHR-012"])

    def test_chr001_grch38_vs_grch37_blocks_rename(self) -> None:
        bam = self._p("grch38.bam")
        _write_test_bam(bam, _human_contigs("GRCh38", "ucsc"))
        fai = self._p("grch37.fa.fai")
        with open(fai, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh37", "ensembl"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        p38 = parse_file(bam)
        p37 = parse_file(fai)
        rep = check_profiles([p38, p37])
        codes = {f.code for f in rep.findings}
        self.assertIn("CHR-001", codes)
        self.assertFalse(rep.compatible)
        self.assertEqual(rep.exit_code, 1)

        with self.assertRaises(UnsafeBuildMappingError):
            build_chr_map(p37, target=p38)

        # --force overrides safety refusal
        forced = build_chr_map(p37, target=p38, force=True)
        self.assertIn(("1", "chr1"), forced)
        self.assertIn(("MT", "chrM"), forced)

    def test_chr001_grch38_vs_chm13_detected(self) -> None:
        d38 = self._p("hg38.dict")
        with open(d38, "wt", encoding="utf-8") as fh:
            fh.write("@HD\tVN:1.6\n")
            for name, length in _human_contigs("GRCh38", "ucsc"):
                fh.write(f"@SQ\tSN:{name}\tLN:{length}\n")

        chm13 = self._p("chm13.chrom.sizes")
        with open(chm13, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("T2T-CHM13v2.0", "ucsc"):
                fh.write(f"{name}\t{length}\n")

        rep = check_profiles([parse_file(d38), parse_file(chm13)])
        codes = {f.code for f in rep.findings}
        self.assertIn("CHR-001", codes)
        self.assertEqual(rep.exit_code, 1)

    def test_chr002_hg19_16571_vs_grch37_16569(self) -> None:
        fai_hg19 = self._p("hg19.fa.fai")
        with open(fai_hg19, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("hg19", "ucsc"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        fai_b37 = self._p("b37.fa.fai")
        with open(fai_b37, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh37", "ensembl"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        p19 = parse_file(fai_hg19)
        p37 = parse_file(fai_b37)
        self.assertEqual(p19.inferred_assemblies, ["hg19"])
        self.assertEqual(p37.inferred_assemblies, ["GRCh37"])
        rep = check_profiles([p19, p37])
        codes = {f.code for f in rep.findings}
        self.assertIn("CHR-002", codes)
        self.assertEqual(rep.exit_code, 1)

    def test_chr003_and_chr004_ucsc_vs_ensembl_and_broken_sed_m(self) -> None:
        bam = self._p("hg38.bam")
        _write_test_bam(bam, _human_contigs("GRCh38", "ucsc"))
        vcf = self._p("clinvar.vcf.gz")
        with gzip.open(vcf, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            for name, length in _human_contigs("GRCh38", "ensembl"):
                fh.write(f"##contig=<ID={name},length={length}>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("1\t1000\trs1\tA\tG\t.\tPASS\t.\n")
            fh.write("MT\t250\trs2\tC\tT\t.\tPASS\t.\n")

        bed_bad_m = self._p("sed_broken.bed")
        with open(bed_bad_m, "wt", encoding="utf-8") as fh:
            fh.write("1\t10\t20\nM\t1\t50\n")

        rep = check_profiles([parse_file(bam), parse_file(vcf), parse_file(bed_bad_m)])
        codes = {f.code for f in rep.findings}
        self.assertIn("CHR-003", codes)
        self.assertIn("CHR-004", codes)
        self.assertNotIn("CHR-001", codes)

    def test_chr005_mixed_styles_within_single_file(self) -> None:
        bed = self._p("mixed.bed")
        with open(bed, "wt", encoding="utf-8") as fh:
            fh.write("chr1\t10\t20\n2\t30\t40\n")
        rep = check_profiles([parse_file(bed)])
        self.assertIn("CHR-005", {f.code for f in rep.findings})

        vcf = self._p("header_vs_body.vcf")
        with open(vcf, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("1\t100\trs1\tA\tG\t.\tPASS\t.\n")
        rep2 = check_profiles([parse_file(vcf)])
        self.assertIn("CHR-005", {f.code for f in rep2.findings})

    def test_chr006_coordinate_exceeds_contig_length(self) -> None:
        fai = self._p("hg38.fa.fai")
        with open(fai, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh38", "ucsc"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        bed = self._p("overflow.bed")
        with open(bed, "wt", encoding="utf-8") as fh:
            # 249500000 exceeds all human assemblies; 249000000 also exceeds GRCh38 (248956422)
            fh.write("chr1\t249400000\t249500000\tover\n")

        rep = check_profiles([parse_file(fai), parse_file(bed)])
        codes = {f.code for f in rep.findings}
        self.assertIn("CHR-006", codes)

        vcf_own = self._p("own_overflow.vcf")
        with open(vcf_own, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("chr1\t1500\t.\tA\tG\t.\tPASS\t.\n")
        rep_own = check_profiles([parse_file(vcf_own)])
        self.assertIn("CHR-006", {f.code for f in rep_own.findings})

    def test_chr007_and_chr008_malformed_and_zero_length_bed(self) -> None:
        bed_zero = self._p("snps.bed")
        with open(bed_zero, "wt", encoding="utf-8") as fh:
            fh.write("chr1\t100\t100\trs123\n")
        rep_warn = check_profiles([parse_file(bed_zero)], fail_on="error")
        self.assertIn("CHR-008", {f.code for f in rep_warn.findings})
        self.assertEqual(rep_warn.exit_code, 0)
        rep_strict = check_profiles([parse_file(bed_zero)], fail_on="warning")
        self.assertEqual(rep_strict.exit_code, 1)

        bed_bad = self._p("bad.bed")
        with open(bed_bad, "wt", encoding="utf-8") as fh:
            fh.write("chr1\t200\t100\tbackwards\nchr1\t-5\t10\tneg\n")
        rep_bad = check_profiles([parse_file(bed_bad)])
        self.assertIn("CHR-007", {f.code for f in rep_bad.findings})
        self.assertEqual(rep_bad.exit_code, 1)

    def test_chr009_sort_order_mismatch(self) -> None:
        fai = self._p("ref.fa.fai")
        with open(fai, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh38", "ucsc"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        bed_lex = self._p("lex.bed")
        with open(bed_lex, "wt", encoding="utf-8") as fh:
            fh.write("chr1\t10\t20\nchr10\t10\t20\nchr2\t10\t20\n")
        rep = check_profiles([parse_file(fai), parse_file(bed_lex)])
        self.assertIn("CHR-009", {f.code for f in rep.findings})

        # Single file check also catches lexicographical chr1, chr10, chr2
        rep_single = check_profiles([parse_file(bed_lex)])
        self.assertIn("CHR-009", {f.code for f in rep_single.findings})

    def test_chr010_unlisted_contig_and_chr011_alt_and_sex(self) -> None:
        vcf = self._p("unlisted.vcf")
        with open(vcf, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("chrUn_gl000220\t100\t.\tA\tG\t.\tPASS\t.\n")
        rep = check_profiles([parse_file(vcf)])
        self.assertIn("CHR-010", {f.code for f in rep.findings})

        d_primary = self._p("primary.dict")
        with open(d_primary, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh38", "ucsc", subset=[str(i) for i in range(1, 23)] + ["X", "MT"]):
                fh.write(f"@SQ\tSN:{name}\tLN:{length}\n")

        d_alt_y = self._p("with_alt_y.dict")
        with open(d_alt_y, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh38", "ucsc"):
                fh.write(f"@SQ\tSN:{name}\tLN:{length}\n")
            fh.write("@SQ\tSN:chr1_KI270706v1_random\tLN:175055\n")

        rep_11 = check_profiles([parse_file(d_primary), parse_file(d_alt_y)])
        self.assertIn("CHR-011", {f.code for f in rep_11.findings})

    def test_refseq_to_ucsc_and_mouse_assemblies(self) -> None:
        mm39 = self._p("mm39.chrom.sizes")
        with open(mm39, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCm39", "ucsc", subset=[str(i) for i in range(1, 20)] + ["X", "Y", "MT"]):
                fh.write(f"{name}\t{length}\n")
        mm10 = self._p("mm10.chrom.sizes")
        with open(mm10, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCm38", "ucsc", subset=[str(i) for i in range(1, 20)] + ["X", "Y", "MT"]):
                fh.write(f"{name}\t{length}\n")

        p39 = parse_file(mm39)
        p10 = parse_file(mm10)
        self.assertEqual(p39.inferred_assemblies, ["GRCm39"])
        self.assertEqual(p10.inferred_assemblies, ["GRCm38"])
        self.assertIn("CHR-001", {f.code for f in check_profiles([p39, p10]).findings})

        vcf_rs = self._p("refseq.vcf")
        with open(vcf_rs, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            for name, length in _human_contigs("GRCh38", "refseq", subset=["1", "2", "MT"]):
                fh.write(f"##contig=<ID={name},length={length}>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        bam_ucsc = self._p("ucsc38.bam")
        _write_test_bam(bam_ucsc, _human_contigs("GRCh38", "ucsc", subset=["1", "2", "MT"]))

        pairs = build_chr_map(parse_file(vcf_rs), target=parse_file(bam_ucsc))
        self.assertEqual(
            pairs,
            [("NC_000001.11", "chr1"), ("NC_000002.12", "chr2"), ("NC_012920.1", "chrM")],
        )

    def test_fix_bed_preserves_chr_gene_names_and_gtf_directives(self) -> None:
        bam = self._p("hg38.bam")
        _write_test_bam(bam, _human_contigs("GRCh38", "ucsc"))

        bed_in = self._p("panel.bed")
        bed_out = self._p("panel.ucsc.bed")
        with open(bed_in, "wt", encoding="utf-8") as fh:
            fh.write("#header comment chr1\n1\t100\t200\tCHRNA1\nMT\t10\t90\tchrRNA_control\n")

        rc = main(["fix-bed", bed_in, "--to", bam, "-o", bed_out])
        self.assertEqual(rc, 0)
        with open(bed_out, "rt", encoding="utf-8") as fh:
            out_lines = fh.read().splitlines()
        self.assertEqual(out_lines[0], "#header comment chr1")
        self.assertEqual(out_lines[1], "chr1\t100\t200\tCHRNA1")
        self.assertEqual(out_lines[2], "chrM\t10\t90\tchrRNA_control")

        gff_in = self._p("genes.gff3")
        gff_out = self._p("genes.ensembl.gff3")
        with open(gff_in, "wt", encoding="utf-8") as fh:
            fh.write("##gff-version 3\n##sequence-region chr1 1 248956422\n")
            fh.write("chr1\tENSEMBL\texon\t100\t500\t.\t+\t.\tID=exon1;Name=CHRNA1\n")
            fh.write("chrM\tENSEMBL\tgene\t1\t500\t.\t+\t.\tID=mt1\n")
        rc2 = main(["fix-bed", gff_in, "--to", "ensembl", "-o", gff_out])
        self.assertEqual(rc2, 0)
        with open(gff_out, "rt", encoding="utf-8") as fh:
            gff_txt = fh.read()
        self.assertIn("##sequence-region 1 1 248956422", gff_txt)
        self.assertIn("1\tENSEMBL\texon\t100\t500\t.\t+\t.\tID=exon1;Name=CHRNA1", gff_txt)
        self.assertIn("MT\tENSEMBL\tgene\t1\t500", gff_txt)

    def test_sam_fasta_and_cli_outputs_json_markdown_and_emit_map(self) -> None:
        sam = self._p("aln.sam")
        with open(sam, "wt", encoding="utf-8") as fh:
            fh.write("@HD\tVN:1.6\n")
            for name, length in _human_contigs("GRCh38", "ucsc", subset=["1", "2", "MT"]):
                fh.write(f"@SQ\tSN:{name}\tLN:{length}\n")
            fh.write("r1\t0\tchr1\t500\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n")

        fa = self._p("mini.fa.gz")
        with gzip.open(fa, "wt", encoding="utf-8") as fh:
            fh.write(">1\nACGTACGT\n>MT\nGGGG\n")

        vcf = self._p("vars.vcf")
        with open(vcf, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            for name, length in _human_contigs("GRCh38", "ensembl", subset=["1", "2", "MT"]):
                fh.write(f"##contig=<ID={name},length={length}>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            fh.write("1\t100\t.\tA\tT\t.\tPASS\t.\n")

        map_path = self._p("map.tsv")
        out_buf = io.StringIO()
        err_buf = io.StringIO()
        with redirect_stdout(out_buf), redirect_stderr(err_buf):
            rc = main([sam, vcf, "--format", "json", "--emit-chr-map", map_path])
        self.assertEqual(rc, 1)
        data = json.loads(out_buf.getvalue())
        self.assertFalse(data["compatible"])
        self.assertTrue(os.path.exists(map_path))
        with open(map_path, "rt", encoding="utf-8") as fh:
            self.assertIn("1\tchr1\n", fh.read())

        md_buf = io.StringIO()
        with redirect_stdout(md_buf):
            rc_md = main(["check", sam, "--format", "markdown"])
        self.assertEqual(rc_md, 0)
        self.assertIn("## chrcheck Preflight Report", md_buf.getvalue())
        self.assertIn("CHR-012", md_buf.getvalue())

        # inspect and fasta smoke
        ins_buf = io.StringIO()
        with redirect_stdout(ins_buf):
            rc_ins = main(["inspect", fa, "--format", "json"])
        self.assertEqual(rc_ins, 0)
        ins_data = json.loads(ins_buf.getvalue())
        self.assertEqual(ins_data["files"][0]["naming_style"], "ensembl")

    def test_rename_map_cli_refuses_incompatible_builds(self) -> None:
        fai38 = self._p("hg38.fa.fai")
        with open(fai38, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh38", "ucsc"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")
        fai37 = self._p("b37.fa.fai")
        with open(fai37, "wt", encoding="utf-8") as fh:
            for name, length in _human_contigs("GRCh37", "ensembl"):
                fh.write(f"{name}\t{length}\t0\t80\t81\n")

        err_buf = io.StringIO()
        with redirect_stderr(err_buf):
            rc = main(["rename-map", fai37, "--to", fai38])
        self.assertEqual(rc, 1)
        self.assertIn("CHR-001", err_buf.getvalue())

    def test_custom_organism_contigs_and_gz_fix_bed(self) -> None:
        fai = self._p("yeast.fa.fai")
        with open(fai, "wt", encoding="utf-8") as fh:
            fh.write("I\t230218\t0\t80\t81\nII\t813184\t0\t80\t81\n")
        bed = self._p("yeast_genes.bed.gz")
        with gzip.open(bed, "wt", encoding="utf-8") as fh:
            fh.write("I\t100\t500\tYAL068C\n")
        rep = check_profiles([parse_file(fai), parse_file(bed)])
        self.assertEqual([f.code for f in rep.findings], ["CHR-012"])

        bed_in = self._p("ensembl.bed.gz")
        bed_out = self._p("ucsc.bed.gz")
        with gzip.open(bed_in, "wt", encoding="utf-8") as fh:
            fh.write("1\t10\t20\tCHRNA1\nMT\t5\t15\tND1\n")
        rc = main(["fix-bed", bed_in, "--to", "ucsc", "-o", bed_out])
        self.assertEqual(rc, 0)
        with gzip.open(bed_out, "rt", encoding="utf-8") as fh:
            txt = fh.read()
        self.assertIn("chr1\t10\t20\tCHRNA1", txt)
        self.assertIn("chrM\t5\t15\tND1", txt)

    def test_sample_records_truncation_and_table_output(self) -> None:
        bed = self._p("many.bed")
        with open(bed, "wt", encoding="utf-8") as fh:
            for i in range(25):
                fh.write(f"chr1\t{i * 10}\t{i * 10 + 5}\n")
        p_trunc = parse_file(bed, sample_records=10)
        self.assertTrue(p_trunc.truncated_scan)
        self.assertEqual(p_trunc.total_records_scanned, 10)

        p_full = parse_file(bed, sample_records=None)
        self.assertFalse(p_full.truncated_scan)
        self.assertEqual(p_full.total_records_scanned, 25)

        out_buf = io.StringIO()
        with redirect_stdout(out_buf):
            rc = main([bed, "--no-color"])
        self.assertEqual(rc, 0)
        self.assertIn("CHR-012 OK", out_buf.getvalue())

    def test_rename_map_output_file_and_cli_error_codes(self) -> None:
        vcf = self._p("ens.vcf")
        with open(vcf, "wt", encoding="utf-8") as fh:
            fh.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n##contig=<ID=MT,length=16569>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out_tsv = self._p("rename.tsv")
        err_buf = io.StringIO()
        with redirect_stderr(err_buf):
            rc = main(["rename-map", vcf, "--to", "ucsc", "-o", out_tsv])
        self.assertEqual(rc, 0)
        with open(out_tsv, "rt", encoding="utf-8") as fh:
            self.assertEqual(fh.read(), "1\tchr1\nMT\tchrM\n")

        with redirect_stderr(io.StringIO()):
            rc_missing = main([self._p("does_not_exist.bed")])
        self.assertEqual(rc_missing, 2)


if __name__ == "__main__":
    unittest.main()
