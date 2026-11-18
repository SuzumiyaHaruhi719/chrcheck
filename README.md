<div align="center">

<img src="assets/banner.svg" alt="chrcheck banner" width="840" />

### Catch genome build conflicts, `chr1` vs `1` naming drift, `chrM` vs `MT` traps, and contig order mismatches before your genomics pipeline silently emits empty files.

[![CI](https://github.com/SuzumiyaHaruhi719/chrcheck/actions/workflows/ci.yml/badge.svg)](https://github.com/SuzumiyaHaruhi719/chrcheck/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/SuzumiyaHaruhi719/chrcheck?color=0ea5e9&logo=github)](https://github.com/SuzumiyaHaruhi719/chrcheck/releases)
[![License: MIT](https://img.shields.io/badge/License-MIT-f59e0b.svg)](LICENSE)
[![Downloads](https://img.shields.io/github/downloads/SuzumiyaHaruhi719/chrcheck/total?color=8b5cf6&logo=python&logoColor=white)](https://github.com/SuzumiyaHaruhi719/chrcheck/releases)

<img src="assets/demo.gif" alt="chrcheck terminal demo" width="820" />

</div>

---

## Why

Bioinformatics tools fail quietly when chromosome names or genome builds disagree:

- `bedtools intersect` and `bcftools annotate` often exit with status `0` and write header-only empty outputs when one input uses `chr1` (UCSC) and another uses `1` (Ensembl) or `NC_000001.11` (RefSeq).
- Mitochondrial annotations on `chrM` vanish against `MT` references, and running `sed 's/^chr//'` turns `chrM` into an invalid `M` instead of `MT` while corrupting gene symbols that start with `chr` (such as `CHRNA1`).
- Mixing `hg19` (`chrM` length `16571`) with `GRCh37` / `b37` (`MT` length `16569`), or `GRCh38` (`chr1` length `248956422`) with `T2T-CHM13v2.0` (`chr1` length `248387328`), passes naive prefix renaming and produces shifted coordinates.
- Lexicographic sort order (`chr1`, `chr10`, `chr2`) breaks `bedtools intersect -sorted` and `GATK` sequence dictionary checks hours into an HPC job.

`chrcheck` inspects `.bam`, `.sam`, `.vcf(.gz)`, `.bed(.gz)`, `.gtf`/`.gff3(.gz)`, `.fa.fai`, and `.dict` headers together in milliseconds without installing `samtools`, `bcftools`, `pysam`, or Java.

## Features

- **Cross-Format Preflight in One Command**: Compare binary `.bam` headers (BGZF-streamed without reading alignment records), `.sam`, `.vcf`/`.vcf.gz`, `.bed`/`.narrowPeak`/`.broadPeak`, `.gtf`/`.gff3`, `.fa.fai`, `.dict`, and `.chrom.sizes` in a single run.
- **Assembly Build Fingerprinting**: Identifies `GRCh38`/`hg38`, `GRCh37`/`b37`, `hg19`, `T2T-CHM13v2.0`/`hs1`, `GRCm39`/`mm39`, and `GRCm38`/`mm10` by contig length signatures and coordinate upper bounds, distinguishing `hg19` (`chrM=16571`) from `GRCh37` (`MT=16569`).
- **Build-Safe Contig Mapping (`rename-map` & `--emit-chr-map`)**: Generates `old<TAB>new` TSV files ready for `bcftools annotate --rename-chrs`, and refuses to emit a rename map if inputs come from different genome assemblies (`CrossMap`/`liftOver` required).
- **Safe Column-1 BED/GTF/GFF3 Rewrite (`fix-bed`)**: Converts column 1 and `##sequence-region` headers to match a target `.bam`, `.fai`, `.dict`, `.vcf.gz`, or style (`ucsc`, `ensembl`, `refseq`) without touching gene names or INFO attributes.
- **12 Concrete Diagnostic Rules (`CHR-001`–`CHR-012`)**: Flags build mismatches, mixed styles within a single file, coordinate overflows (`pos > contig_len`), 1-based `start == end` zero-length BED traps, sort-order drift, and missing `##contig` headers.
- **Zero Runtime Dependencies**: Pure Python 3.9+ standard library (`zlib`, `gzip`, `struct`), running natively on Linux, macOS, and Windows.

## Quickstart

Install directly from GitHub (or run without global install via `pipx` / `uvx`):

```bash
pip install "git+https://github.com/SuzumiyaHaruhi719/chrcheck.git"
```

Run a preflight check across mixed pipeline inputs:

```bash
# 1. Check BAM, BED, and VCF compatibility before running downstream tools
chrcheck alignments.bam capture_panel.bed clinvar.vcf.gz --emit-chr-map chr_map.tsv

# 2. Rewrite BED column 1 to match the BAM reference dictionary safely
chrcheck fix-bed capture_panel.bed --to alignments.bam -o capture_panel.ucsc.bed

# 3. Re-check everything before submitting the cluster job
chrcheck alignments.bam capture_panel.ucsc.bed
```

## Usage

### 1. Check compatibility across files (`check`)

`check` is the default subcommand when passing file paths directly:

```bash
chrcheck reference.fa.fai sample.bam targets.bed variants.vcf.gz
```

Useful flags:

- `--format table|json|markdown`: output format (`table` default)
- `--fail-on error|warning`: exit `1` on errors only (`error`, default) or also on warnings (`warning`)
- `--emit-chr-map path.tsv`: write a `bcftools annotate --rename-chrs` mapping from the last mismatched input to the first/reference input when assembly builds are compatible
- `--sample-records 50000` / `--full`: cap sampled records on huge plain/gz VCF/BED/GTF files (`50000` default) or scan every record
- `--no-color`: disable ANSI colors (also honors `NO_COLOR`)

### 2. Inspect a single file (`inspect`)

```bash
chrcheck inspect clinvar.vcf.gz
chrcheck inspect alignments.bam --format json
```

Reads only the BGZF header blocks for `.bam` files, so inspecting a 50 GB WGS BAM completes in a few milliseconds.

### 3. Generate a `bcftools --rename-chrs` map (`rename-map`)

```bash
chrcheck rename-map clinvar.vcf.gz --to alignments.bam -o chr_map.tsv
bcftools annotate --rename-chrs chr_map.tsv -Oz -o clinvar.ucsc.vcf.gz clinvar.vcf.gz
```

If `clinvar.vcf.gz` is `GRCh37` and `alignments.bam` is `GRCh38`, `chrcheck` stops with exit code `1` and advises liftover instead of generating a coordinate-corrupting rename map (override only if intentional with `--force`).

### 4. Fix a BED, narrowPeak, GTF, or GFF3 file (`fix-bed`)

```bash
chrcheck fix-bed targets.bed --to alignments.bam -o targets.ucsc.bed
chrcheck fix-bed genes.gtf --to ensembl -o genes.ensembl.gtf
```

Only column 1 (and GFF3 `##sequence-region` directives) is remapped, including `chrM` <-> `MT` and RefSeq `NC_*` accessions, leaving gene symbols (`CHRNA1`, `CHRNG`) and attributes untouched.

### 5. Pipeline & CI integration

Use `chrcheck` as a fast preflight stage in Nextflow, Snakemake, or GitHub Actions:

```bash
chrcheck --fail-on warning "${bam}" "${bed}" "${vcf}"
```

```yaml
- name: Preflight reference compatibility
  run: |
    chrcheck --format markdown --fail-on error data/ref.fa.fai data/panel.bed data/variants.vcf.gz >> "$GITHUB_STEP_SUMMARY"
```

## Comparison with Alternatives

| Capability | `chrcheck` | `dpryan79/ChromosomeMappings` | `vcfverifier` | `vep-preflight` | `picard CompareSAMs` | `cut` / `sed` one-liners |
|---|---|---|---|---|---|---|
| Cross-format (`BAM + VCF + BED + GTF + FAI + DICT`) | **Yes** | Static TSV tables only | VCF + FASTA only | Single VEP VCF only | SAM/BAM only | Manual per format |
| Zero runtime deps (`samtools`/`pysam`/Java not needed) | **Pure Python stdlib** | Manual curl/grep | Rust binary | Rust binary | Requires JVM + Picard | POSIX shell |
| Reads 50 GB `.bam` header in milliseconds | **Yes (BGZF header only)** | N/A | N/A | N/A | Reads full header | Fails on binary BAM |
| Detects `GRCh38` vs `GRCh37` vs `hg19` (`16571` vs `16569`) | **Yes** | Static files per pair | Ref-base check only | Single-file headers | Length match only | No |
| Refuses rename map across incompatible builds | **Yes (`CHR-001`/`002`)** | No | N/A | N/A | N/A | Silently corrupts coords |
| Catches `chrM`/`MT`, order drift & `start==end` BED traps | **Yes (`CHR-004/008/009`)** | No | No | `chr` mix only | Order only | `sed` turns `chrM` -> `M` |

## Diagnostic Rules

| Code | Severity | Meaning |
|---|---|---|
| `CHR-001` | `ERROR` | Genome assembly build mismatch across inputs (`GRCh38` vs `GRCh37`, `CHM13v2.0`, `mm39` vs `mm10`) |
| `CHR-002` | `ERROR` | `hg19` (`chrM=16571`) vs `GRCh37`/`b37` (`MT=16569`) mitochondrial length mismatch |
| `CHR-003` | `ERROR` | Chromosome naming convention mismatch (`UCSC chr1` vs `Ensembl 1` vs `RefSeq NC_*`) |
| `CHR-004` | `ERROR` | Mitochondrial contig naming mismatch (`chrM` vs `MT` vs `M`) |
| `CHR-005` | `ERROR` | Mixed chromosome naming styles inside one file or between header and records |
| `CHR-006` | `ERROR` | Record coordinate exceeds contig length in header or reference dictionary |
| `CHR-007` | `ERROR` | Malformed interval (`start < 0`, `start > end`, or 1-based `pos < 1`) |
| `CHR-008` | `WARNING` | Zero-length BED intervals (`start == end`) suggesting unshifted 1-based coordinates |
| `CHR-009` | `WARNING` | Contig sort-order mismatch across inputs (`chr1, chr10, chr2` vs karyotypic order) |
| `CHR-010` | `WARNING` | Data records reference contigs absent from the file's own header or reference input |
| `CHR-011` | `WARNING` | Sex chromosome or alt/decoy scaffold mismatch across inputs |
| `CHR-012` | `INFO` | All checked inputs are compatible |

## FAQ

**Does `chrcheck` load an entire `.bam` or `.vcf.gz` into memory?**
No. `.bam` inspection stops immediately after parsing the compressed BGZF header blocks and binary `n_ref` dictionary. Text and `.gz` inputs stream line by line up to `--sample-records` (default `50000`, or all lines with `--full`).

**Does `chrcheck` work on non-human organisms?**
Yes. Built-in assembly signatures cover human (`GRCh38`/`hg38`, `GRCh37`/`b37`, `hg19`, `T2T-CHM13v2.0`), mouse (`GRCm39`/`mm39`, `GRCm38`/`mm10`), zebrafish (`GRCz11`/`danRer11`), fruit fly (`BDGP6`/`dm6`), and *C. elegans* (`WBcel235`/`ce11`). For any other organism or custom reference, `chrcheck` leaves the assembly label as `custom/unknown` and still validates cross-file contig names, lengths, sort order, mitochondrial style, and coordinate bounds.

**Why not just run `sed 's/^chr//' input.bed`?**
Three reasons: `sed` turns `chrM` into `M` (not `MT`), alters feature names like `CHRNA1` or `chr22` in later columns if not column-anchored, and gives zero warning when your BED file is actually `hg19` while your BAM is `GRCh38`.

## Roadmap

- CRAM binary container header parsing with `@SQ` `M5` checksum comparison ([#1](https://github.com/SuzumiyaHaruhi719/chrcheck/issues/1))
- SARIF 2.1.0 output (`--format sarif`) for GitHub Code Scanning annotations ([#2](https://github.com/SuzumiyaHaruhi719/chrcheck/issues/2))
- Binary `.bigWig` / `.bigBed` chromosome B+ tree header inspection ([#4](https://github.com/SuzumiyaHaruhi719/chrcheck/issues/4))

## Contributing

Contributions are welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) and our [Code of Conduct](CODE_OF_CONDUCT.md) before opening a pull request.

## License

Released under the [MIT License](LICENSE).
