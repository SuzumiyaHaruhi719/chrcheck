"""Cross-file genome build, naming convention, and coordinate compatibility checks."""

from typing import Dict, List, Optional, Tuple

from chrcheck.assemblies import (
    builds_compatible,
    classify_contig,
    is_alt_or_decoy,
    karyotypic_rank,
)
from chrcheck.models import CheckReport, FileProfile, Finding


def _mito_length_for(profile: FileProfile) -> Optional[int]:
    for contig, length in profile.header_lengths.items():
        canon, _, _ = classify_contig(contig)
        if canon == "MT":
            return length
    return None


def _canonical_order(profile: FileProfile) -> List[str]:
    seen = set()
    order: List[str] = []
    for contig in profile.all_contigs_ordered():
        canon, _, _ = classify_contig(contig)
        if canon is not None and canon not in seen:
            seen.add(canon)
            order.append(canon)
    return order


def _reference_canonical_lengths(profiles: List[FileProfile]) -> Dict[str, Tuple[int, str]]:
    """Collect known canonical contig lengths across reference-backed inputs."""
    ref_lengths: Dict[str, Tuple[int, str]] = {}
    for p in profiles:
        for contig, length in p.header_lengths.items():
            canon, _, _ = classify_contig(contig)
            if canon is not None and canon not in ref_lengths:
                ref_lengths[canon] = (length, p.path)
    return ref_lengths


def check_profiles(
    profiles: List[FileProfile],
    fail_on: str = "error",
) -> CheckReport:
    """Evaluate CHR-001 through CHR-012 across one or more FileProfile objects."""
    findings: List[Finding] = []
    file_paths = [p.path for p in profiles]

    # CHR-001 & CHR-002: Assembly build compatibility and hg19 (16571) vs GRCh37 (16569)
    build_conflict = False
    mito_len_conflict = False
    for i in range(len(profiles)):
        for j in range(i + 1, len(profiles)):
            a, b = profiles[i], profiles[j]
            mt_a, mt_b = _mito_length_for(a), _mito_length_for(b)
            if mt_a is not None and mt_b is not None and mt_a != mt_b:
                if {mt_a, mt_b} == {16569, 16571} and (
                    set(a.inferred_assemblies) <= {"GRCh37", "hg19"}
                    and set(b.inferred_assemblies) <= {"GRCh37", "hg19"}
                ):
                    mito_len_conflict = True
                    findings.append(
                        Finding(
                            code="CHR-002",
                            severity="ERROR",
                            title="hg19 vs GRCh37/b37 mitochondrial length mismatch (16571 vs 16569)",
                            message=(
                                f"{a.path} has mitochondrial contig length {mt_a} bp while "
                                f"{b.path} has {mt_b} bp. Autosomal coordinates match, but MT/chrM "
                                "variants cannot be merged without realigning or dropping MT."
                            ),
                            files=[a.path, b.path],
                            remediation=(
                                "Realign mitochondrial reads to rCRS (16,569 bp) or filter MT/chrM "
                                "before combining hg19 and GRCh37/b37 inputs."
                            ),
                            details={"lengths": {a.path: mt_a, b.path: mt_b}},
                        )
                    )
                    continue

            if not builds_compatible(a.inferred_assemblies, b.inferred_assemblies):
                build_conflict = True
                findings.append(
                    Finding(
                        code="CHR-001",
                        severity="ERROR",
                        title="Genome assembly build mismatch",
                        message=(
                            f"{a.path} matches {', '.join(a.inferred_assemblies)} while "
                            f"{b.path} matches {', '.join(b.inferred_assemblies)}. "
                            "Renaming contigs without coordinate liftover will silently map intervals to wrong loci."
                        ),
                        files=[a.path, b.path],
                        remediation=(
                            "Convert coordinates to a shared assembly using CrossMap or UCSC liftOver "
                            "with the appropriate chain file before intersecting or annotating."
                        ),
                        details={
                            a.path: a.inferred_assemblies,
                            b.path: b.inferred_assemblies,
                        },
                    )
                )

    # CHR-005: Mixed naming styles inside one file or between header and data records
    for p in profiles:
        header_styles = set()
        for c in p.header_contigs:
            if c.canonical is not None and c.style in ("ucsc", "ensembl", "refseq", "broken_mito"):
                header_styles.add(c.style)
        record_styles = set(p.observed_styles)
        if p.naming_style == "mixed" or (
            header_styles and record_styles and not (header_styles & record_styles)
        ) or len(record_styles) > 1:
            findings.append(
                Finding(
                    code="CHR-005",
                    severity="ERROR",
                    title="Mixed chromosome naming styles within a single file",
                    message=(
                        f"{p.path} mixes naming conventions "
                        f"(header={sorted(header_styles) or ['none']}, "
                        f"records={sorted(record_styles) or ['none']})."
                    ),
                    files=[p.path],
                    remediation=(
                        "Normalize all headers and records in this file to one convention "
                        "(e.g. `chrcheck fix-bed` for BED/GTF or `bcftools annotate --rename-chrs` for VCF)."
                    ),
                )
            )

    # CHR-003: Cross-file chromosome naming convention mismatch
    styles_by_file: Dict[str, List[str]] = {}
    for p in profiles:
        if p.naming_style in ("ucsc", "ensembl", "refseq"):
            styles_by_file.setdefault(p.naming_style, []).append(p.path)
    if len(styles_by_file) > 1:
        summary = ", ".join(f"{style} ({', '.join(paths)})" for style, paths in styles_by_file.items())
        findings.append(
            Finding(
                code="CHR-003",
                severity="ERROR",
                title="Chromosome naming convention mismatch across inputs",
                message=f"Inputs use incompatible naming styles: {summary}.",
                files=file_paths,
                remediation=(
                    "Harmonize contig names to match your reference/BAM using "
                    "`chrcheck rename-map` (for bcftools) or `chrcheck fix-bed` (for BED/GTF)."
                ),
                details={"styles": styles_by_file},
            )
        )

    # CHR-004: Mitochondrial contig naming mismatch or broken 'M' from sed
    broken_mito_files = [
        p.path for p in profiles if p.mito_style and "M" in p.mito_style.split("/")
    ]
    mito_styles = {p.mito_style for p in profiles if p.mito_style}
    if broken_mito_files or len(mito_styles) > 1:
        detail_parts = [f"{p.path}: {p.mito_style}" for p in profiles if p.mito_style]
        remediation = (
            "Map mitochondrial contigs explicitly (`chrM` <-> `MT`, never `chrM` -> `M`)."
        )
        findings.append(
            Finding(
                code="CHR-004",
                severity="ERROR",
                title="Mitochondrial contig name mismatch",
                message=(
                    "Mitochondrial contig names disagree or use bare 'M' (a common `sed 's/^chr//'` bug): "
                    + "; ".join(detail_parts)
                ),
                files=[p.path for p in profiles if p.mito_style],
                remediation=remediation,
            )
        )

    # CHR-006: Record coordinate exceeds contig length
    ref_lengths = _reference_canonical_lengths(profiles)
    for p in profiles:
        overflows: List[str] = []
        for contig, max_pos in p.observed_max_pos.items():
            own_len = p.contig_length(contig)
            if own_len is not None and max_pos > own_len:
                overflows.append(f"{contig}:{max_pos} > {own_len} (header)")
                continue
            if not p.header_lengths or (not build_conflict and not mito_len_conflict):
                canon, _, _ = classify_contig(contig)
                if canon and canon in ref_lengths:
                    rlen, rpath = ref_lengths[canon]
                    if max_pos > rlen:
                        overflows.append(f"{contig}:{max_pos} > {rlen} ({rpath})")
        if overflows:
            findings.append(
                Finding(
                    code="CHR-006",
                    severity="ERROR",
                    title="Record coordinate exceeds contig length",
                    message=f"{p.path} contains coordinates past contig end: {', '.join(overflows[:5])}.",
                    files=[p.path],
                    remediation="Verify genome build alignment and trim/filter out-of-bounds intervals.",
                )
            )

    # CHR-007: Malformed intervals (start < 0, start > end, pos < 1)
    for p in profiles:
        if p.malformed_intervals:
            findings.append(
                Finding(
                    code="CHR-007",
                    severity="ERROR",
                    title="Malformed genomic coordinates",
                    message=(
                        f"{p.path} has {len(p.malformed_intervals)} malformed interval(s), "
                        f"e.g. {'; '.join(p.malformed_intervals[:3])}."
                    ),
                    files=[p.path],
                    remediation="Fix negative, inverted (start > end), or zero-based VCF positions.",
                )
            )

    # CHR-008: Zero-length BED intervals (start == end) suggesting 1-based VCF mistake
    for p in profiles:
        if p.zero_length_bed_count > 0:
            findings.append(
                Finding(
                    code="CHR-008",
                    severity="WARNING",
                    title="Zero-length BED intervals (start == end)",
                    message=(
                        f"{p.path} contains {p.zero_length_bed_count} zero-length BED interval(s) "
                        f"(e.g. {', '.join(p.zero_length_examples[:3])}), commonly caused by writing "
                        "1-based VCF POS into BED without `start = POS - 1`."
                    ),
                    files=[p.path],
                    remediation="Convert 1-based SNP positions to 0-based half-open BED intervals (`POS - 1`, `POS`).",
                )
            )

    # CHR-009: Contig sort-order mismatch across inputs (or lexicographical vs karyotypic)
    if len(profiles) >= 2:
        ref_order = _canonical_order(profiles[0])
        ref_pos = {c: idx for idx, c in enumerate(ref_order)}
        for p in profiles[1:]:
            other_order = _canonical_order(p)
            shared_in_other = [c for c in other_order if c in ref_pos]
            shared_in_ref = [c for c in ref_order if c in set(shared_in_other)]
            if len(shared_in_other) >= 2 and shared_in_other != shared_in_ref:
                findings.append(
                    Finding(
                        code="CHR-009",
                        severity="WARNING",
                        title="Contig sort-order mismatch across inputs",
                        message=(
                            f"{p.path} contig order ({', '.join(shared_in_other[:5])}...) differs from "
                            f"{profiles[0].path} ({', '.join(shared_in_ref[:5])}...), which breaks "
                            "`bedtools intersect -sorted` and GATK sequence dictionary checks."
                        ),
                        files=[profiles[0].path, p.path],
                        remediation=(
                            f"Reorder intervals to match {profiles[0].path} (e.g. "
                            "`bedtools sort -faidx ref.fa.fai -i input.bed` or `bcftools sort`)."
                        ),
                    )
                )
    elif len(profiles) == 1:
        contigs = profiles[0].all_contigs_ordered()
        primary_only = [c for c in contigs if classify_contig(c)[0] is not None]
        if len(primary_only) >= 3 and primary_only != sorted(primary_only, key=karyotypic_rank):
            findings.append(
                Finding(
                    code="CHR-009",
                    severity="WARNING",
                    title="Non-karyotypic (e.g. lexicographical) contig order",
                    message=(
                        f"{profiles[0].path} orders primary contigs as "
                        f"{', '.join(primary_only[:6])}... instead of karyotypic order."
                    ),
                    files=[profiles[0].path],
                    remediation="Sort according to your reference `.fai` dictionary before using `-sorted` tools.",
                )
            )

    # CHR-010: Observed data contigs missing from own header or reference dictionary
    ref_canon_sets = [
        {classify_contig(c)[0] for c in p.header_order if classify_contig(c)[0] is not None}
        for p in profiles
        if p.header_order
    ]
    ref_exact_sets = [set(p.header_order) for p in profiles if p.header_order]
    for p in profiles:
        if p.header_order and p.observed_contigs:
            missing_own = [c for c in p.observed_contigs if c not in p.header_lengths and c not in p.header_order]
            if missing_own:
                findings.append(
                    Finding(
                        code="CHR-010",
                        severity="WARNING",
                        title="Record contigs missing from file header",
                        message=(
                            f"{p.path} has data records on contig(s) absent from its header: "
                            f"{', '.join(missing_own[:5])}."
                        ),
                        files=[p.path],
                        remediation="Regenerate header `##contig` lines (`bcftools reheader --fai`) before indexing.",
                    )
                )
        elif ref_canon_sets and p.observed_contigs:
            # Headerless file (BED/GTF without ##sequence-region) checked against primary ref
            ref_p = next((x for x in profiles if x.header_order), None)
            if ref_p is not None:
                ref_exact = ref_exact_sets[0]
                ref_canon = ref_canon_sets[0]
                unmatched = []
                for c in p.observed_contigs:
                    canon, _, _ = classify_contig(c)
                    if c not in ref_exact and (canon is None or canon not in ref_canon):
                        unmatched.append(c)
                if unmatched:
                    findings.append(
                        Finding(
                            code="CHR-010",
                            severity="WARNING",
                            title="Contigs absent from reference dictionary",
                            message=(
                                f"{p.path} references contig(s) not present in {ref_p.path}: "
                                f"{', '.join(unmatched[:5])}."
                            ),
                            files=[ref_p.path, p.path],
                            remediation="Filter out unmapped/decoy contigs or confirm the reference FASTA matches.",
                        )
                    )

    # CHR-011: Sex chromosome or alt/decoy scaffold mismatch across reference-backed files
    if len(profiles) >= 2:
        alt_flags = [(p.path, p.has_alt_or_decoy) for p in profiles if p.header_order]
        if len(alt_flags) >= 2 and len({flag for _, flag in alt_flags}) > 1:
            with_alt = [path for path, flag in alt_flags if flag]
            without_alt = [path for path, flag in alt_flags if not flag]
            findings.append(
                Finding(
                    code="CHR-011",
                    severity="WARNING",
                    title="Alt/decoy scaffold presence mismatch",
                    message=(
                        f"{', '.join(with_alt)} include alt/decoy/unlocalized contigs while "
                        f"{', '.join(without_alt)} contain primary contigs only."
                    ),
                    files=with_alt + without_alt,
                    remediation="Ensure downstream tools ignore unplaced scaffolds or filter to primary chromosomes.",
                )
            )

        # Sex chromosome mismatch across reference-backed inputs (e.g. Y missing in one female-only or primary ref)
        ref_profiles = [p for p in profiles if p.header_order and len(p.header_order) >= 20]
        if len(ref_profiles) >= 2:
            sex_sets = [
                (
                    p.path,
                    {
                        classify_contig(c)[0]
                        for c in p.header_order
                        if classify_contig(c)[0] in ("X", "Y") and not is_alt_or_decoy(c)
                    },
                )
                for p in ref_profiles
            ]
            if len({frozenset(s) for _, s in sex_sets}) > 1:
                desc = "; ".join(f"{path}: {sorted(s) or ['none']}" for path, s in sex_sets)
                findings.append(
                    Finding(
                        code="CHR-011",
                        severity="WARNING",
                        title="Sex chromosome representation mismatch",
                        message=f"Reference dictionaries differ on X/Y inclusion ({desc}).",
                        files=[path for path, _ in sex_sets],
                        remediation="Verify whether chrY/Y is intentionally excluded or filter regions accordingly.",
                    )
                )

    if not findings:
        findings.append(
            Finding(
                code="CHR-012",
                severity="INFO",
                title="All inputs are compatible",
                message=f"Verified {len(profiles)} file(s): genome build, contig naming, order, and bounds agree.",
                files=file_paths,
            )
        )

    has_error = any(f.severity == "ERROR" for f in findings)
    has_warning = any(f.severity == "WARNING" for f in findings)
    compatible = not has_error
    if fail_on == "warning":
        exit_code = 1 if (has_error or has_warning) else 0
    else:
        exit_code = 1 if has_error else 0

    return CheckReport(
        files=profiles,
        findings=findings,
        compatible=compatible,
        exit_code=exit_code,
    )
