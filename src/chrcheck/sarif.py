"""SARIF 2.1.0 serialisation of chrcheck findings (standard library only)."""

from typing import Dict, List

from chrcheck import __version__
from chrcheck.models import CheckReport

SARIF_SCHEMA = "https://json.schemastore.org/sarif-2.1.0.json"
SARIF_VERSION = "2.1.0"
INFORMATION_URI = "https://github.com/SuzumiyaHaruhi719/chrcheck"

_LEVELS = {"ERROR": "error", "WARNING": "warning", "INFO": "note"}


def _level(severity: str) -> str:
    return _LEVELS.get(severity.upper(), "note")


def _uri(path: str) -> str:
    return path.replace("\\", "/")


def report_to_sarif(report: CheckReport) -> Dict[str, object]:
    """Convert a :class:`CheckReport` into a SARIF 2.1.0 log dictionary.

    One SARIF rule is emitted per distinct ``CHR-`` code (in first-seen order);
    every finding becomes one result that references its rule by id and index.
    """
    rules: List[Dict[str, object]] = []
    rule_index: Dict[str, int] = {}
    results: List[Dict[str, object]] = []

    for finding in report.findings:
        if finding.code not in rule_index:
            rule_index[finding.code] = len(rules)
            rules.append(
                {
                    "id": finding.code,
                    "name": finding.code,
                    "shortDescription": {"text": finding.title},
                    "defaultConfiguration": {"level": _level(finding.severity)},
                }
            )
        text = finding.message
        if finding.remediation:
            text = f"{text} Fix: {finding.remediation}"
        result: Dict[str, object] = {
            "ruleId": finding.code,
            "ruleIndex": rule_index[finding.code],
            "level": _level(finding.severity),
            "message": {"text": text},
        }
        if finding.files:
            result["locations"] = [
                {"physicalLocation": {"artifactLocation": {"uri": _uri(path)}}}
                for path in finding.files
            ]
        results.append(result)

    return {
        "$schema": SARIF_SCHEMA,
        "version": SARIF_VERSION,
        "runs": [
            {
                "tool": {
                    "driver": {
                        "name": "chrcheck",
                        "version": __version__,
                        "informationUri": INFORMATION_URI,
                        "rules": rules,
                    }
                },
                "results": results,
            }
        ],
    }
