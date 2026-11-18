"""Reference genome assembly fingerprints and contig equivalence rules."""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass(frozen=True)
class AssemblySpec:
    """Canonical contig definitions and lengths for a reference assembly."""

    name: str
    aliases: Tuple[str, ...]
    species: str
    lengths: Dict[str, int]  # canonical name ("1".."22", "X", "Y", "MT") -> length
    ucsc_names: Dict[str, str]
    ensembl_names: Dict[str, str]
    refseq_names: Dict[str, str]
    default_mito: str


def _build_human_maps(
    lengths: Dict[str, int],
    refseq: Dict[str, str],
    default_mito: str = "chrM",
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    ucsc: Dict[str, str] = {}
    ensembl: Dict[str, str] = {}
    for canon in lengths:
        if canon == "MT":
            ucsc[canon] = default_mito
            ensembl[canon] = "MT"
        else:
            ucsc[canon] = f"chr{canon}"
            ensembl[canon] = canon
    return ucsc, ensembl, refseq


def _build_mouse_maps(
    lengths: Dict[str, int],
    refseq: Dict[str, str],
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    ucsc: Dict[str, str] = {}
    ensembl: Dict[str, str] = {}
    for canon in lengths:
        if canon == "MT":
            ucsc[canon] = "chrM"
            ensembl[canon] = "MT"
        else:
            ucsc[canon] = f"chr{canon}"
            ensembl[canon] = canon
    return ucsc, ensembl, refseq


def _build_organism_maps(
    lengths: Dict[str, int],
    refseq: Dict[str, str],
    ensembl_mito: str = "MT",
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """UCSC (``chr`` prefix, ``chrM``) and Ensembl names for non-human, non-mouse organisms."""
    ucsc: Dict[str, str] = {}
    ensembl: Dict[str, str] = {}
    for canon in lengths:
        if canon == "MT":
            ucsc[canon] = "chrM"
            ensembl[canon] = ensembl_mito
        else:
            ucsc[canon] = f"chr{canon}"
            ensembl[canon] = canon
    return ucsc, ensembl, refseq


_GRCH38_LENGTHS: Dict[str, int] = {
    "1": 248956422,
    "2": 242193529,
    "3": 198295559,
    "4": 190214555,
    "5": 181538259,
    "6": 170805979,
    "7": 159345973,
    "8": 145138636,
    "9": 138394717,
    "10": 133797422,
    "11": 135086622,
    "12": 133275309,
    "13": 114364328,
    "14": 107043718,
    "15": 101991189,
    "16": 90338345,
    "17": 83257441,
    "18": 80373285,
    "19": 58617616,
    "20": 64444167,
    "21": 46709983,
    "22": 50818468,
    "X": 156040895,
    "Y": 57227415,
    "MT": 16569,
}

_GRCH38_REFSEQ: Dict[str, str] = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "X": "NC_000023.11",
    "Y": "NC_000024.10",
    "MT": "NC_012920.1",
}

_GRCH37_BASE_LENGTHS: Dict[str, int] = {
    "1": 249250621,
    "2": 243199373,
    "3": 198022430,
    "4": 191154276,
    "5": 180915260,
    "6": 171115067,
    "7": 159138663,
    "8": 146364022,
    "9": 141213431,
    "10": 135534747,
    "11": 135006516,
    "12": 133851895,
    "13": 115169878,
    "14": 107349540,
    "15": 102531392,
    "16": 90354753,
    "17": 81195210,
    "18": 78077248,
    "19": 59128983,
    "20": 63025520,
    "21": 48129895,
    "22": 51304566,
    "X": 155270560,
    "Y": 59373566,
    "MT": 16569,
}

_HG19_LENGTHS: Dict[str, int] = {**_GRCH37_BASE_LENGTHS, "MT": 16571}

_GRCH37_REFSEQ: Dict[str, str] = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "X": "NC_000023.10",
    "Y": "NC_000024.9",
    "MT": "NC_012920.1",
}

_HG19_REFSEQ: Dict[str, str] = {**_GRCH37_REFSEQ, "MT": "NC_001807.4"}

_CHM13_LENGTHS: Dict[str, int] = {
    "1": 248387328,
    "2": 242696752,
    "3": 201105948,
    "4": 193574945,
    "5": 182045439,
    "6": 172126628,
    "7": 160567428,
    "8": 146259331,
    "9": 150617247,
    "10": 134758134,
    "11": 135127769,
    "12": 133324548,
    "13": 113566686,
    "14": 101161492,
    "15": 99753195,
    "16": 96330374,
    "17": 84276897,
    "18": 80542538,
    "19": 61707364,
    "20": 66210255,
    "21": 45090682,
    "22": 51324926,
    "X": 154259566,
    "Y": 62460029,
    "MT": 16569,
}

_CHM13_REFSEQ: Dict[str, str] = {
    "1": "NC_060925.1",
    "2": "NC_060926.1",
    "3": "NC_060927.1",
    "4": "NC_060928.1",
    "5": "NC_060929.1",
    "6": "NC_060930.1",
    "7": "NC_060931.1",
    "8": "NC_060932.1",
    "9": "NC_060933.1",
    "10": "NC_060934.1",
    "11": "NC_060935.1",
    "12": "NC_060936.1",
    "13": "NC_060937.1",
    "14": "NC_060938.1",
    "15": "NC_060939.1",
    "16": "NC_060940.1",
    "17": "NC_060941.1",
    "18": "NC_060942.1",
    "19": "NC_060943.1",
    "20": "NC_060944.1",
    "21": "NC_060945.1",
    "22": "NC_060946.1",
    "X": "NC_060947.1",
    "Y": "NC_060948.1",
    "MT": "NC_012920.1",
}

_GRCM39_LENGTHS: Dict[str, int] = {
    "1": 195154279,
    "2": 181755017,
    "3": 159745316,
    "4": 156860686,
    "5": 151758149,
    "6": 149588044,
    "7": 144995196,
    "8": 130127694,
    "9": 124359700,
    "10": 130530862,
    "11": 121973369,
    "12": 120092757,
    "13": 120883175,
    "14": 125139656,
    "15": 104073951,
    "16": 98008968,
    "17": 95294699,
    "18": 90720763,
    "19": 61420004,
    "X": 169476592,
    "Y": 91455967,
    "MT": 16299,
}

_GRCM39_REFSEQ: Dict[str, str] = {
    "1": "NC_000067.7",
    "2": "NC_000068.8",
    "3": "NC_000069.7",
    "4": "NC_000070.7",
    "5": "NC_000071.7",
    "6": "NC_000072.7",
    "7": "NC_000073.7",
    "8": "NC_000074.7",
    "9": "NC_000075.7",
    "10": "NC_000076.7",
    "11": "NC_000077.7",
    "12": "NC_000078.7",
    "13": "NC_000079.7",
    "14": "NC_000080.7",
    "15": "NC_000081.7",
    "16": "NC_000082.7",
    "17": "NC_000083.7",
    "18": "NC_000084.7",
    "19": "NC_000085.7",
    "X": "NC_000086.8",
    "Y": "NC_000087.8",
    "MT": "NC_005089.1",
}

_GRCM38_LENGTHS: Dict[str, int] = {
    "1": 195471971,
    "2": 182113224,
    "3": 160039680,
    "4": 156508116,
    "5": 151834684,
    "6": 149736546,
    "7": 145441459,
    "8": 129401213,
    "9": 124595110,
    "10": 130694993,
    "11": 122082543,
    "12": 120129022,
    "13": 120421639,
    "14": 124902244,
    "15": 104043685,
    "16": 98207768,
    "17": 94987271,
    "18": 90702639,
    "19": 61431566,
    "X": 171031299,
    "Y": 91744698,
    "MT": 16299,
}

_GRCM38_REFSEQ: Dict[str, str] = {
    "1": "NC_000067.6",
    "2": "NC_000068.7",
    "3": "NC_000069.6",
    "4": "NC_000070.6",
    "5": "NC_000071.6",
    "6": "NC_000072.6",
    "7": "NC_000073.6",
    "8": "NC_000074.6",
    "9": "NC_000075.6",
    "10": "NC_000076.6",
    "11": "NC_000077.6",
    "12": "NC_000078.6",
    "13": "NC_000079.6",
    "14": "NC_000080.6",
    "15": "NC_000081.6",
    "16": "NC_000082.6",
    "17": "NC_000083.6",
    "18": "NC_000084.6",
    "19": "NC_000085.6",
    "X": "NC_000086.7",
    "Y": "NC_000087.7",
    "MT": "NC_005089.1",
}


_GRCZ11_LENGTHS: Dict[str, int] = {
    "1": 59578282,
    "2": 59640629,
    "3": 62628489,
    "4": 78093715,
    "5": 72500376,
    "6": 60270059,
    "7": 74282399,
    "8": 54304671,
    "9": 56459846,
    "10": 45420867,
    "11": 45484837,
    "12": 49182954,
    "13": 52186027,
    "14": 52660232,
    "15": 48040578,
    "16": 55266484,
    "17": 53461100,
    "18": 51023478,
    "19": 48449771,
    "20": 55201332,
    "21": 45934066,
    "22": 39133080,
    "23": 46223584,
    "24": 42172926,
    "25": 37502051,
    "MT": 16596,
}

_GRCZ11_REFSEQ: Dict[str, str] = {
    "1": "NC_007112.7",
    "2": "NC_007113.7",
    "3": "NC_007114.7",
    "4": "NC_007115.7",
    "5": "NC_007116.7",
    "6": "NC_007117.7",
    "7": "NC_007118.7",
    "8": "NC_007119.7",
    "9": "NC_007120.7",
    "10": "NC_007121.7",
    "11": "NC_007122.7",
    "12": "NC_007123.7",
    "13": "NC_007124.7",
    "14": "NC_007125.7",
    "15": "NC_007126.7",
    "16": "NC_007127.7",
    "17": "NC_007128.7",
    "18": "NC_007129.7",
    "19": "NC_007130.7",
    "20": "NC_007131.7",
    "21": "NC_007132.7",
    "22": "NC_007133.7",
    "23": "NC_007134.7",
    "24": "NC_007135.7",
    "25": "NC_007136.7",
    "MT": "NC_002333.2",
}

_BDGP6_LENGTHS: Dict[str, int] = {
    "X": 23542271,
    "2L": 23513712,
    "2R": 25286936,
    "3L": 28110227,
    "3R": 32079331,
    "4": 1348131,
    "Y": 3667352,
    "MT": 19524,
}

_BDGP6_REFSEQ: Dict[str, str] = {
    "X": "NC_004354.4",
    "2L": "NT_033779.5",
    "2R": "NT_033778.4",
    "3L": "NT_037436.4",
    "3R": "NT_033777.3",
    "4": "NC_004353.4",
    "Y": "NC_024512.1",
    "MT": "NC_024511.2",
}

_WBCEL235_LENGTHS: Dict[str, int] = {
    "I": 15072434,
    "II": 15279421,
    "III": 13783801,
    "IV": 17493829,
    "V": 20924180,
    "X": 17718942,
    "MT": 13794,
}

_WBCEL235_REFSEQ: Dict[str, str] = {
    "I": "NC_003279.8",
    "II": "NC_003280.10",
    "III": "NC_003281.10",
    "IV": "NC_003282.8",
    "V": "NC_003283.11",
    "X": "NC_003284.9",
    "MT": "NC_001328.1",
}


def _make_spec(
    name: str,
    aliases: Tuple[str, ...],
    species: str,
    lengths: Dict[str, int],
    refseq: Dict[str, str],
    default_mito: str,
) -> AssemblySpec:
    if species == "human":
        ucsc, ensembl, rs = _build_human_maps(lengths, refseq, default_mito=default_mito)
    elif species == "worm":
        ucsc, ensembl, rs = _build_organism_maps(lengths, refseq, ensembl_mito="MtDNA")
    elif species in ("zebrafish", "fly"):
        ucsc, ensembl, rs = _build_organism_maps(lengths, refseq)
    else:
        ucsc, ensembl, rs = _build_mouse_maps(lengths, refseq)
    return AssemblySpec(
        name=name,
        aliases=aliases,
        species=species,
        lengths=lengths,
        ucsc_names=ucsc,
        ensembl_names=ensembl,
        refseq_names=rs,
        default_mito=default_mito,
    )


ASSEMBLIES: Dict[str, AssemblySpec] = {
    "GRCh38": _make_spec(
        "GRCh38",
        ("hg38", "grch38", "hs38dh", "gca_000001405.15"),
        "human",
        _GRCH38_LENGTHS,
        _GRCH38_REFSEQ,
        "chrM",
    ),
    "GRCh37": _make_spec(
        "GRCh37",
        ("b37", "grch37", "hs37d5", "g1kv37"),
        "human",
        _GRCH37_BASE_LENGTHS,
        _GRCH37_REFSEQ,
        "MT",
    ),
    "hg19": _make_spec(
        "hg19",
        ("hg19",),
        "human",
        _HG19_LENGTHS,
        _HG19_REFSEQ,
        "chrM",
    ),
    "T2T-CHM13v2.0": _make_spec(
        "T2T-CHM13v2.0",
        ("chm13", "chm13v2.0", "t2t-chm13v2.0", "hs1"),
        "human",
        _CHM13_LENGTHS,
        _CHM13_REFSEQ,
        "chrM",
    ),
    "GRCm39": _make_spec(
        "GRCm39",
        ("mm39", "grcm39"),
        "mouse",
        _GRCM39_LENGTHS,
        _GRCM39_REFSEQ,
        "chrM",
    ),
    "GRCm38": _make_spec(
        "GRCm38",
        ("mm10", "grcm38"),
        "mouse",
        _GRCM38_LENGTHS,
        _GRCM38_REFSEQ,
        "chrM",
    ),
    "GRCz11": _make_spec(
        "GRCz11",
        ("danrer11", "grcz11"),
        "zebrafish",
        _GRCZ11_LENGTHS,
        _GRCZ11_REFSEQ,
        "chrM",
    ),
    "BDGP6": _make_spec(
        "BDGP6",
        ("dm6", "bdgp6"),
        "fly",
        _BDGP6_LENGTHS,
        _BDGP6_REFSEQ,
        "chrM",
    ),
    "WBcel235": _make_spec(
        "WBcel235",
        ("ce11", "wbcel235"),
        "worm",
        _WBCEL235_LENGTHS,
        _WBCEL235_REFSEQ,
        "chrM",
    ),
}

# Build reverse lookup maps from RefSeq accession (exact and versionless base)
_REFSEQ_EXACT_TO_CANON_AND_BUILD: Dict[str, Tuple[str, str]] = {}
_REFSEQ_BASE_TO_CANON: Dict[str, str] = {}

for _build_name, _spec in ASSEMBLIES.items():
    for _canon, _acc in _spec.refseq_names.items():
        _REFSEQ_EXACT_TO_CANON_AND_BUILD[_acc] = (_canon, _build_name)
        _base = _acc.split(".", 1)[0]
        _REFSEQ_BASE_TO_CANON[_base] = _canon

_HUMAN_PRIMARY_CANONICAL = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
_MOUSE_PRIMARY_CANONICAL = [str(i) for i in range(1, 20)] + ["X", "Y", "MT"]

# Non-numeric primary chromosome names: X/Y, Drosophila arms (2L, 2R, 3L, 3R) and
# C. elegans Roman numerals (I-V). Digit-only names are matched separately.
_NAMED_CHROM = r"[XY]|2L|2R|3L|3R|III|II|IV|I|V"

_ALT_OR_DECOY_RE = re.compile(
    r"(_alt$|_random$|^chrUn_|^Un_|^HLA-|^EBV$|^hs37d5$|^chrEBV$|_decoy$|^GL\d+|^KI\d+|^JH\d+|^KZ\d+)",
    re.IGNORECASE,
)


def is_alt_or_decoy(contig: str) -> bool:
    """Return True when contig looks like an alt, unlocalized, or decoy scaffold."""
    return bool(_ALT_OR_DECOY_RE.search(contig))


def classify_contig(name: str) -> Tuple[Optional[str], str, Optional[str]]:
    """Return (canonical_name, style, assembly_hint) for a contig name.

    Styles: 'ucsc', 'ensembl', 'refseq', 'broken_mito', 'other'.
    """
    raw = name.strip()
    if not raw:
        return None, "other", None

    if raw in _REFSEQ_EXACT_TO_CANON_AND_BUILD:
        canon, build_hint = _REFSEQ_EXACT_TO_CANON_AND_BUILD[raw]
        return canon, "refseq", build_hint

    base = raw.split(".", 1)[0]
    if raw.startswith("NC_") and base in _REFSEQ_BASE_TO_CANON:
        return _REFSEQ_BASE_TO_CANON[base], "refseq", None

    if raw == "chrMT":
        return "MT", "ucsc", None
    if raw == "chrM":
        return "MT", "ucsc", None
    if raw in ("MT", "MtDNA"):
        return "MT", "ensembl", None
    if raw == "M":
        return "MT", "broken_mito", None

    m_ucsc = re.fullmatch(rf"chr([1-9][0-9]?|{_NAMED_CHROM})", raw)
    if m_ucsc:
        return m_ucsc.group(1), "ucsc", None

    m_ens = re.fullmatch(rf"([1-9][0-9]?|{_NAMED_CHROM})", raw)
    if m_ens:
        return m_ens.group(1), "ensembl", None

    # Handle case-insensitive prefix anomalies such as CHR1 or Chr1
    m_case = re.fullmatch(r"chr([1-9][0-9]?|[XYM])", raw, re.IGNORECASE)
    if m_case:
        sym = m_case.group(1).upper()
        return ("MT" if sym == "M" else sym), "other", None

    return None, "other", None


# Within the autosome group, ranks are scaled by 10 so that Drosophila arms sort
# between "2"/"3"/"4": 2L, 2R, 3L, 3R, 4.
_ARM_RANKS: Dict[str, int] = {"2L": 20, "2R": 21, "3L": 30, "3R": 31}
_ROMAN_RANKS: Dict[str, int] = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5}


def karyotypic_rank(name: str) -> Tuple[int, int, str]:
    """Return deterministic karyotypic sort key for a contig name."""
    canon, _, _ = classify_contig(name)
    if canon is not None:
        if canon.isdigit():
            return (0, int(canon) * 10, name)
        if canon in _ARM_RANKS:
            return (0, _ARM_RANKS[canon], name)
        if canon in _ROMAN_RANKS:
            return (0, _ROMAN_RANKS[canon] * 10, name)
        if canon == "X":
            return (1, 0, name)
        if canon == "Y":
            return (1, 1, name)
        if canon == "MT":
            return (1, 2, name)
    return (2, 0, name)


def infer_assemblies_from_header(
    header_lengths: Dict[str, int],
    reference_hints: Optional[List[str]] = None,
) -> List[str]:
    """Identify matching reference assemblies from contig lengths and metadata hints."""
    if not header_lengths:
        return _infer_from_reference_hints(reference_hints or [])

    canon_lengths: Dict[str, int] = {}
    mito_name: Optional[str] = None
    refseq_build_votes: Dict[str, int] = {}
    for contig, length in header_lengths.items():
        canon, style, build_hint = classify_contig(contig)
        if canon is not None:
            # Keep first canonical length unless MT duplicates occur
            canon_lengths.setdefault(canon, length)
            if canon == "MT":
                mito_name = contig
        if style == "refseq" and build_hint is not None:
            refseq_build_votes[build_hint] = refseq_build_votes.get(build_hint, 0) + 1

    if not canon_lengths:
        return _infer_from_reference_hints(reference_hints or [])

    matches: List[str] = []
    for build_name, spec in ASSEMBLIES.items():
        overlap = [c for c in canon_lengths if c in spec.lengths]
        if not overlap:
            continue
        if all(canon_lengths[c] == spec.lengths[c] for c in overlap):
            # Distinguish hg19 (chrM=16571) from GRCh37 (MT=16569) when MT isn't present
            if build_name == "GRCh37" and "MT" not in canon_lengths and mito_name is None:
                matches.append(build_name)
            elif build_name == "hg19" and "MT" not in canon_lengths:
                matches.append(build_name)
            else:
                matches.append(build_name)

    # Disambiguate GRCh37 vs hg19 when only one autosomal hint or naming style is present
    if "GRCh37" in matches and "hg19" in matches:
        # Look for RefSeq hints or reference URI hints
        if "hg19" in refseq_build_votes and "GRCh37" not in refseq_build_votes:
            matches = [m for m in matches if m != "GRCh37"]
        elif "GRCh37" in refseq_build_votes and "hg19" not in refseq_build_votes:
            matches = [m for m in matches if m != "hg19"]
        else:
            hinted = _infer_from_reference_hints(reference_hints or [])
            if hinted == ["hg19"]:
                matches = ["hg19"]
            elif hinted == ["GRCh37"]:
                matches = ["GRCh37"]

    # If multiple builds matched on lengths but RefSeq accessions uniquely identify one
    if len(matches) > 1 and refseq_build_votes:
        narrowed = [m for m in matches if m in refseq_build_votes]
        if narrowed:
            matches = narrowed

    # Distinguish GRCh37 vs hg19 by mito naming if MT wasn't in header_lengths? No,
    # if MT IS in header_lengths (16569 vs 16571), the length loop already split them!
    return matches


def _infer_from_reference_hints(hints: List[str]) -> List[str]:
    found: List[str] = []
    joined = " ".join(hints).lower()
    if not joined:
        return found
    for build_name, spec in ASSEMBLIES.items():
        if any(alias in joined for alias in spec.aliases):
            found.append(build_name)
    return found


def infer_assemblies_from_coordinates(
    observed_max_pos: Dict[str, int],
    reference_hints: Optional[List[str]] = None,
) -> List[str]:
    """Narrow candidate assemblies for headerless files using observed coordinates."""
    if not observed_max_pos:
        return _infer_from_reference_hints(reference_hints or [])

    canon_max: Dict[str, int] = {}
    saw_chr_above_19 = False
    for contig, pos in observed_max_pos.items():
        canon, _, _ = classify_contig(contig)
        if canon is None:
            continue
        if canon in ("20", "21", "22"):
            saw_chr_above_19 = True
        prev = canon_max.get(canon, 0)
        if pos > prev:
            canon_max[canon] = pos

    if not canon_max:
        return _infer_from_reference_hints(reference_hints or [])

    candidates: List[str] = []
    for build_name, spec in ASSEMBLIES.items():
        # Chromosomes 20-22 exist in human and zebrafish only (not mouse).
        if saw_chr_above_19 and spec.species not in ("human", "zebrafish"):
            continue
        if all(c in spec.lengths and max_pos <= spec.lengths[c] for c, max_pos in canon_max.items()):
            candidates.append(build_name)

    hinted = _infer_from_reference_hints(reference_hints or [])
    if hinted:
        narrowed = [c for c in candidates if c in hinted]
        if narrowed:
            return narrowed
    return candidates


def builds_compatible(builds_a: List[str], builds_b: List[str]) -> bool:
    """Return True if two candidate build sets overlap or one is empty/unknown."""
    if not builds_a or not builds_b:
        return True
    # Note: GRCh37 and hg19 share autosomal coordinates (only MT differs 16569 vs 16571),
    # which is handled specifically by CHR-002. Here, check strict overlap first.
    if set(builds_a) & set(builds_b):
        return True
    if set(builds_a) <= {"GRCh37", "hg19"} and set(builds_b) <= {"GRCh37", "hg19"}:
        return True
    return False


def map_contig_to_style(
    contig: str,
    target_style: str,
    target_build: Optional[str] = None,
    target_mito: Optional[str] = None,
    target_contigs: Optional[Dict[str, str]] = None,
) -> Optional[str]:
    """Map a contig name to the requested target style or target contig dictionary."""
    canon, _, build_hint = classify_contig(contig)
    if canon is None:
        return None

    # Exact target contig dictionary (canonical -> actual target contig name) takes precedence
    if target_contigs and canon in target_contigs:
        return target_contigs[canon]

    if target_style == "ucsc":
        if canon == "MT":
            return target_mito or "chrM"
        return f"chr{canon}"

    if target_style == "ensembl":
        if canon == "MT":
            spec = ASSEMBLIES.get(target_build or build_hint or "")
            return spec.ensembl_names.get("MT", "MT") if spec else "MT"
        return canon

    if target_style == "refseq":
        build = target_build or build_hint or "GRCh38"
        spec = ASSEMBLIES.get(build)
        if spec and canon in spec.refseq_names:
            return spec.refseq_names[canon]
        return None

    return None
