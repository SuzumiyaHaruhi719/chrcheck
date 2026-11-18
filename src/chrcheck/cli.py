"""Command-line interface and table/JSON/Markdown formatters for chrcheck."""

import argparse
import json
import os
import sys
from typing import List, Optional

from chrcheck import __version__
from chrcheck.checker import check_profiles
from chrcheck.harmonize import (
    UnsafeBuildMappingError,
    build_chr_map,
    format_chr_map,
    rewrite_bed_or_gff,
)
from chrcheck.models import CheckReport, FileProfile
from chrcheck.parsers import parse_file
from chrcheck.sarif import report_to_sarif

_SUBCOMMANDS = {"check", "inspect", "rename-map", "fix-bed"}


def _use_color(no_color: bool) -> bool:
    if no_color or os.environ.get("NO_COLOR") is not None:
        return False
    return sys.stdout.isatty()


def _c(text: str, code: str, enabled: bool) -> str:
    if not enabled:
        return text
    return f"\033[{code}m{text}\033[0m"


def _format_assembly_cell(assemblies: List[str]) -> str:
    if not assemblies:
        return "custom/unknown"
    if len(assemblies) <= 2:
        return ", ".join(assemblies)
    return f"{assemblies[0]} (+{len(assemblies) - 1})"


def _format_summary_table(report: CheckReport, color: bool) -> str:
    lines: List[str] = []
    headers = ["File", "Format", "Style", "Mito", "Assembly", "Contigs", "Records"]
    rows = []
    for p in report.files:
        asm = _format_assembly_cell(p.inferred_assemblies)
        contigs_n = str(len(p.header_order) or len(p.observed_contigs))
        rec_n = str(p.total_records_scanned) + ("+" if p.truncated_scan else "")
        rows.append([
            os.path.basename(p.path),
            p.fmt,
            p.naming_style,
            p.mito_style or "-",
            asm,
            contigs_n,
            rec_n,
        ])

    widths = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]
    sep = "  "
    lines.append(sep.join(h.ljust(widths[i]) for i, h in enumerate(headers)))
    lines.append(sep.join("-" * widths[i] for i in range(len(headers))))
    for r in rows:
        lines.append(sep.join(r[i].ljust(widths[i]) for i in range(len(headers))))

    lines.append("")
    for f in report.findings:
        if f.severity == "ERROR":
            badge = _c(f"[{f.code} ERROR]", "1;31", color)
        elif f.severity == "WARNING":
            badge = _c(f"[{f.code} WARN ]", "1;33", color)
        else:
            badge = _c(f"[{f.code} OK   ]", "1;32", color)
        lines.append(f"{badge} {f.title}")
        lines.append(f"  {f.message}")
        if f.remediation:
            lines.append(f"  Fix: {f.remediation}")

    if report.suggested_commands:
        lines.append("")
        lines.append(_c("Suggested actions:", "1;36", color))
        for cmd in report.suggested_commands:
            lines.append(f"  {cmd}")
    return "\n".join(lines) + "\n"


def _format_markdown(report: CheckReport) -> str:
    lines = ["## chrcheck Preflight Report", ""]
    lines.append("| File | Format | Style | Mito | Assembly | Contigs | Records |")
    lines.append("|---|---|---|---|---|---:|---:|")
    for p in report.files:
        asm = ", ".join(p.inferred_assemblies) if p.inferred_assemblies else "custom/unknown"
        contigs_n = len(p.header_order) or len(p.observed_contigs)
        rec_n = f"{p.total_records_scanned}{'+' if p.truncated_scan else ''}"
        lines.append(
            f"| `{os.path.basename(p.path)}` | `{p.fmt}` | `{p.naming_style}` | "
            f"`{p.mito_style or '-'}` | `{asm}` | {contigs_n} | {rec_n} |"
        )
    lines.append("")
    for f in report.findings:
        lines.append(f"- **{f.code} ({f.severity}) — {f.title}**: {f.message}")
        if f.remediation:
            lines.append(f"  - Fix: `{f.remediation}`")
    return "\n".join(lines) + "\n"


def _report_to_dict(report: CheckReport) -> dict:
    return {
        "compatible": report.compatible,
        "exit_code": report.exit_code,
        "files": [
            {
                "path": p.path,
                "format": p.fmt,
                "naming_style": p.naming_style,
                "mito_style": p.mito_style,
                "inferred_assemblies": p.inferred_assemblies,
                "header_contig_count": len(p.header_order),
                "observed_contig_count": len(p.observed_contigs),
                "records_scanned": p.total_records_scanned,
                "truncated_scan": p.truncated_scan,
                "header_contigs": [
                    {"name": c.name, "length": c.length, "canonical": c.canonical, "style": c.style}
                    for c in p.header_contigs
                ],
                "observed_contigs": p.observed_contigs,
            }
            for p in report.files
        ],
        "findings": [
            {
                "code": f.code,
                "severity": f.severity,
                "title": f.title,
                "message": f.message,
                "files": f.files,
                "remediation": f.remediation,
            }
            for f in report.findings
        ],
        "suggested_commands": report.suggested_commands,
    }


def _resolve_target_arg(target_arg: str, sample_records: Optional[int]):
    lower = target_arg.lower()
    if lower in ("ucsc", "ensembl", "refseq"):
        return None, lower
    return parse_file(target_arg, sample_records=sample_records), None


def _cmd_check(args: argparse.Namespace) -> int:
    sample_records = None if args.full else args.sample_records
    profiles: List[FileProfile] = []
    for path in args.files:
        profiles.append(parse_file(path, sample_records=sample_records))

    report = check_profiles(profiles, fail_on=args.fail_on)

    if args.emit_chr_map and len(profiles) >= 2:
        ref = profiles[0]
        mismatched = [
            p for p in profiles[1:]
            if p.naming_style in ("ucsc", "ensembl", "refseq", "broken_mito", "mixed")
            and p.naming_style != ref.naming_style
        ]
        if mismatched:
            combined_pairs = []
            seen_pairs = set()
            mapped_sources = []
            for src_profile in mismatched:
                try:
                    for pair in build_chr_map(src_profile, target=ref, force=args.force):
                        if pair not in seen_pairs:
                            seen_pairs.add(pair)
                            combined_pairs.append(pair)
                    mapped_sources.append(os.path.basename(src_profile.path))
                except UnsafeBuildMappingError as exc:
                    sys.stderr.write(f"Warning: skipped --emit-chr-map for {src_profile.path}: {exc}\n")
            if combined_pairs:
                with open(args.emit_chr_map, "wt", encoding="utf-8") as fh:
                    fh.write(format_chr_map(combined_pairs))
                report.suggested_commands.append(
                    f"Wrote {len(combined_pairs)} contig mappings ({', '.join(mapped_sources)} -> {os.path.basename(ref.path)}) to {args.emit_chr_map}"
                )

    if not args.quiet:
        if args.format == "json":
            sys.stdout.write(json.dumps(_report_to_dict(report), indent=2) + "\n")
        elif args.format == "markdown":
            sys.stdout.write(_format_markdown(report))
        elif args.format == "sarif":
            sys.stdout.write(json.dumps(report_to_sarif(report), indent=2) + "\n")
        else:
            sys.stdout.write(_format_summary_table(report, _use_color(args.no_color)))
    return report.exit_code


def _cmd_inspect(args: argparse.Namespace) -> int:
    sample_records = None if args.full else args.sample_records
    profiles = [parse_file(path, sample_records=sample_records) for path in args.files]
    if args.format == "json":
        payload = {
            "files": [
                {
                    "path": p.path,
                    "format": p.fmt,
                    "naming_style": p.naming_style,
                    "mito_style": p.mito_style,
                    "inferred_assemblies": p.inferred_assemblies,
                    "header_contigs": [
                        {"name": c.name, "length": c.length, "canonical": c.canonical, "style": c.style}
                        for c in p.header_contigs
                    ],
                    "observed_contigs": p.observed_contigs,
                    "records_scanned": p.total_records_scanned,
                }
                for p in profiles
            ]
        }
        sys.stdout.write(json.dumps(payload, indent=2) + "\n")
        return 0

    for p in profiles:
        asm = ", ".join(p.inferred_assemblies) if p.inferred_assemblies else "custom/unknown"
        sys.stdout.write(f"File: {p.path} ({p.fmt})\n")
        sys.stdout.write(f"  Style: {p.naming_style} | Mito: {p.mito_style or '-'} | Assembly: {asm}\n")
        contigs = p.header_order or p.observed_contigs
        preview = ", ".join(contigs[:10]) + ("..." if len(contigs) > 10 else "")
        sys.stdout.write(f"  Contigs ({len(contigs)}): {preview}\n")
    return 0


def _cmd_rename_map(args: argparse.Namespace) -> int:
    sample_records = None if args.full else args.sample_records
    src = parse_file(args.source, sample_records=sample_records)
    target_profile, target_style = _resolve_target_arg(args.to, sample_records=sample_records)
    pairs = build_chr_map(src, target=target_profile, target_style=target_style, force=args.force)
    content = format_chr_map(pairs)
    if args.output:
        with open(args.output, "wt", encoding="utf-8") as fh:
            fh.write(content)
        sys.stderr.write(f"Wrote {len(pairs)} contig mappings to {args.output}\n")
    else:
        sys.stdout.write(content)
    return 0


def _cmd_fix_bed(args: argparse.Namespace) -> int:
    sample_records = None if args.full else args.sample_records
    src = parse_file(args.source, sample_records=None)
    if src.fmt not in ("bed", "gtf"):
        raise ValueError(f"fix-bed expects a BED/narrowPeak/GTF/GFF3 file, got {src.fmt}")
    target_profile, target_style = _resolve_target_arg(args.to, sample_records=sample_records)
    rewritten, total = rewrite_bed_or_gff(
        args.source,
        args.output,
        src,
        target=target_profile,
        target_style=target_style,
        force=args.force,
    )
    sys.stderr.write(f"Rewrote {rewritten}/{total} record contig names -> {args.output}\n")
    return 0


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chrcheck",
        description="Genome build and chromosome naming preflight CLI for BAM, VCF, BED, GTF, and FASTA.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    sub = parser.add_subparsers(dest="command")

    p_check = sub.add_parser("check", help="Check build, naming, order, and bounds across files.")
    p_check.add_argument("files", nargs="+", help="Genomics input files (reference/BAM first recommended).")
    p_check.add_argument("--format", choices=("table", "json", "markdown", "sarif"), default="table")
    p_check.add_argument("--fail-on", choices=("error", "warning"), default="error")
    p_check.add_argument("--emit-chr-map", metavar="PATH", help="Write bcftools --rename-chrs TSV map.")
    p_check.add_argument("--sample-records", type=int, default=50000)
    p_check.add_argument("--full", action="store_true", help="Scan all records instead of sampling.")
    p_check.add_argument("--force", action="store_true", help="Emit rename map even on build mismatch.")
    p_check.add_argument("--no-color", action="store_true")
    p_check.add_argument("--quiet", action="store_true")

    p_inspect = sub.add_parser("inspect", help="Inspect header and contig profile of one or more files.")
    p_inspect.add_argument("files", nargs="+")
    p_inspect.add_argument("--format", choices=("table", "json"), default="table")
    p_inspect.add_argument("--sample-records", type=int, default=50000)
    p_inspect.add_argument("--full", action="store_true")

    p_map = sub.add_parser("rename-map", help="Generate bcftools annotate --rename-chrs TSV map.")
    p_map.add_argument("source", help="Input file whose contigs need renaming.")
    p_map.add_argument("--to", required=True, help="Target file (.bam/.fai/.vcf) or style (ucsc|ensembl|refseq).")
    p_map.add_argument("-o", "--output", help="Output TSV file (defaults to stdout).")
    p_map.add_argument("--sample-records", type=int, default=50000)
    p_map.add_argument("--full", action="store_true")
    p_map.add_argument("--force", action="store_true")

    p_fix = sub.add_parser("fix-bed", help="Safely remap column 1 of a BED, GTF, or GFF3 file.")
    p_fix.add_argument("source", help="Input BED/narrowPeak/GTF/GFF3 file.")
    p_fix.add_argument("--to", required=True, help="Target file (.bam/.fai/.vcf) or style (ucsc|ensembl|refseq).")
    p_fix.add_argument("-o", "--output", required=True, help="Output file path.")
    p_fix.add_argument("--sample-records", type=int, default=50000)
    p_fix.add_argument("--full", action="store_true")
    p_fix.add_argument("--force", action="store_true")

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    if raw_argv and raw_argv[0] not in _SUBCOMMANDS and raw_argv[0] not in ("-h", "--help", "--version"):
        raw_argv = ["check"] + raw_argv

    parser = _build_parser()
    if not raw_argv:
        parser.print_help()
        return 0

    args = parser.parse_args(raw_argv)
    try:
        if args.command == "check":
            return _cmd_check(args)
        if args.command == "inspect":
            return _cmd_inspect(args)
        if args.command == "rename-map":
            return _cmd_rename_map(args)
        if args.command == "fix-bed":
            return _cmd_fix_bed(args)
    except (FileNotFoundError, ValueError, UnsafeBuildMappingError) as exc:
        sys.stderr.write(f"chrcheck error: {exc}\n")
        return 1 if isinstance(exc, UnsafeBuildMappingError) else 2
    return 0
