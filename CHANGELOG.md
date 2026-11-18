# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-11-18

### Added

- Streaming parsers for binary `.bam` (BGZF header-only decompression), `.sam`, `.vcf`/`.vcf.gz`, `.bed`/`.bed.gz`/`.narrowPeak`/`.broadPeak`, `.gtf`/`.gff`/`.gff3` (`.gz`), `.fa.fai`, `.dict`, `.chrom.sizes`, and `.fa`/`.fasta` (`.gz`).
- Built-in contig length fingerprints and UCSC/Ensembl/RefSeq equivalence maps for `GRCh38`/`hg38`, `GRCh37`/`b37`, `hg19`, `T2T-CHM13v2.0`/`hs1`, `GRCm39`/`mm39`, and `GRCm38`/`mm10`.
- 12 compatibility diagnostic rules (`CHR-001` through `CHR-012`), including `hg19` (`16571`) vs `GRCh37` (`16569`) mitochondrial length checking, `sed 's/^chr//'` `M` detection, sort-order drift, coordinate overflow, and zero-length `start == end` BED warnings.
- Safe `rename-map` TSV generation for `bcftools annotate --rename-chrs` with automatic refusal when assembly builds conflict.
- Column-1 `fix-bed` rewriting for BED, GTF, and GFF3 (`##sequence-region` included) without modifying gene symbols or attributes.
- Table, JSON, and GitHub-flavored Markdown report outputs (`--format table|json|markdown`) and `--fail-on error|warning` CI exit gating.
