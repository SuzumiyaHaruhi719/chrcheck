"""Core data models for contig profiling and compatibility findings."""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class ContigInfo:
    """Metadata for a single contig declared in a header or index."""

    name: str
    length: Optional[int] = None
    canonical: Optional[str] = None
    style: str = "other"
    assembly_hint: Optional[str] = None
    md5: Optional[str] = None


@dataclass
class FileProfile:
    """Parsed summary of one genomics input file."""

    path: str
    fmt: str
    header_contigs: List[ContigInfo] = field(default_factory=list)
    header_order: List[str] = field(default_factory=list)
    header_lengths: Dict[str, int] = field(default_factory=dict)
    observed_contigs: List[str] = field(default_factory=list)
    observed_counts: Dict[str, int] = field(default_factory=dict)
    observed_max_pos: Dict[str, int] = field(default_factory=dict)
    total_records_scanned: int = 0
    truncated_scan: bool = False
    naming_style: str = "unknown"
    observed_styles: List[str] = field(default_factory=list)
    mito_style: Optional[str] = None
    inferred_assemblies: List[str] = field(default_factory=list)
    reference_hints: List[str] = field(default_factory=list)
    malformed_intervals: List[str] = field(default_factory=list)
    zero_length_bed_count: int = 0
    zero_length_examples: List[str] = field(default_factory=list)
    has_alt_or_decoy: bool = False
    parse_warnings: List[str] = field(default_factory=list)

    def all_contigs_ordered(self) -> List[str]:
        """Return header contigs plus any additional observed contigs in encounter order."""
        if not self.header_order:
            return list(self.observed_contigs)
        seen = set(self.header_order)
        extra = [c for c in self.observed_contigs if c not in seen]
        return list(self.header_order) + extra

    def contig_length(self, name: str) -> Optional[int]:
        """Look up length by exact name."""
        return self.header_lengths.get(name)


@dataclass
class Finding:
    """Single diagnostic rule violation or compatibility confirmation."""

    code: str
    severity: str  # ERROR, WARNING, INFO
    title: str
    message: str
    files: List[str] = field(default_factory=list)
    remediation: str = ""
    details: Dict[str, object] = field(default_factory=dict)


@dataclass
class CheckReport:
    """Combined cross-file check result."""

    files: List[FileProfile]
    findings: List[Finding]
    compatible: bool
    exit_code: int
    suggested_commands: List[str] = field(default_factory=list)
