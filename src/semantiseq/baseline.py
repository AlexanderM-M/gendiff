"""Compact semantic baselines that do not retain genomic records."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Sequence, Tuple

from semantiseq.alignment import _scan as scan_alignment
from semantiseq.batch import BatchEntry, read_manifest
from semantiseq.compare import _kind
from semantiseq.fingerprint import Fingerprint, digest
from semantiseq.profiles import fields_for
from semantiseq.variant import _normalize as normalize_variant
from semantiseq.variant import _scan as scan_variant


class BaselineError(ValueError):
    """Raised when a baseline cannot be created or checked."""


def _fingerprint(value: Fingerprint) -> Dict[str, Any]:
    return {
        "count": value.count,
        "total": f"{value.total:064x}",
        "squares": f"{value.squares:064x}",
        "xor": f"{value.xor:064x}",
    }


@lru_cache(maxsize=8)
def _file_digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)
            if not block:
                break
            hasher.update(block)
    return hasher.hexdigest()


def _scan_entry(
    entry: BatchEntry,
    path: Path,
    threads: int,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    progress: bool,
) -> Dict[str, Any]:
    if not path.is_file():
        raise BaselineError(f"file not found: {path}")
    kind = _kind(path)
    fields = fields_for(kind, entry.profile)
    reference_digest = _file_digest(entry.reference) if entry.reference else None
    with TemporaryDirectory(prefix="semantiseq-baseline-") as work:
        scan_path = path
        if kind == "variant" and entry.normalize:
            if entry.reference is None:
                raise BaselineError("variant normalization requires a reference")
            scan_path = Path(work) / "normalized.bcf"
            normalize_variant(path, scan_path, entry.reference)
        if kind == "alignment":
            if entry.normalize:
                raise BaselineError("normalization is only valid for variant files")
            scan = scan_alignment(
                scan_path,
                entry.reference,
                threads,
                fields,
                ignore_tags,
                None,
                progress,
                f"{entry.sample} / {entry.stage}",
            )
        else:
            scan = scan_variant(
                scan_path,
                threads,
                fields,
                ignore_info,
                None,
                progress,
                f"{entry.sample} / {entry.stage}",
            )
    return {
        "sample": entry.sample,
        "stage": entry.stage,
        "kind": kind,
        "profile": entry.profile,
        "normalize": entry.normalize,
        "reference_sha256": reference_digest,
        "record_count": scan.count,
        "fields": list(fields),
        "fingerprints": {
            name: _fingerprint(value) for name, value in scan.fingerprints.items()
        },
        "structural_header": f"{digest(scan.structural):064x}",
    }


def create_baseline(
    manifest: Path,
    *,
    threads: int,
    profile: str,
    reference: Optional[Path],
    normalize_variants: bool,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    progress: bool,
) -> Dict[str, Any]:
    entries = read_manifest(manifest, profile, reference, normalize_variants)
    items = [
        _scan_entry(entry, entry.before, threads, ignore_tags, ignore_info, progress)
        for entry in entries
    ]
    return {
        "semantiseq_baseline_format": 1,
        "privacy": "Contains semantic hashes and counts; no genomic records or loci.",
        "settings": {
            "ignore_tags": list(ignore_tags),
            "ignore_info": list(ignore_info),
        },
        "comparisons": items,
    }


@dataclass(frozen=True)
class BaselineCheck:
    sample: str
    stage: str
    passed: bool
    record_count: int
    changed_fields: Tuple[str, ...]
    structure_equal: bool
    reasons: Tuple[str, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample": self.sample,
            "stage": self.stage,
            "passed": self.passed,
            "record_count": self.record_count,
            "changed_fields": list(self.changed_fields),
            "structural_header_equal": self.structure_equal,
            "reasons": list(self.reasons),
        }


def _read_baseline(path: Path) -> Dict[str, Any]:
    if not path.is_file():
        raise BaselineError(f"baseline not found: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise BaselineError(f"invalid baseline JSON: {error}") from error
    if not isinstance(payload, dict) or payload.get("semantiseq_baseline_format") != 1:
        raise BaselineError("unsupported baseline format")
    if not isinstance(payload.get("comparisons"), list):
        raise BaselineError("baseline comparisons must be an array")
    settings = payload.get("settings", {})
    if not isinstance(settings, dict):
        raise BaselineError("baseline settings must be an object")
    for name in ("ignore_tags", "ignore_info"):
        value = settings.get(name, [])
        if not isinstance(value, list) or not all(
            isinstance(item, str) for item in value
        ):
            raise BaselineError(f"baseline {name} must be an array of strings")
    identities = set()
    for index, item in enumerate(payload["comparisons"], 1):
        if not isinstance(item, dict):
            raise BaselineError(f"baseline comparison {index} must be an object")
        identity = (item.get("sample"), item.get("stage"))
        if not all(isinstance(value, str) and value for value in identity):
            raise BaselineError(
                f"baseline comparison {index} needs sample and stage strings"
            )
        if identity in identities:
            raise BaselineError(f"duplicate baseline comparison: {identity}")
        identities.add(identity)
        if not isinstance(item.get("fields"), list) or not isinstance(
            item.get("fingerprints"), dict
        ):
            raise BaselineError(f"baseline comparison {index} has invalid fingerprints")
    return payload


def check_baseline(
    baseline: Path,
    manifest: Path,
    *,
    threads: int,
    profile: str,
    reference: Optional[Path],
    normalize_variants: bool,
    progress: bool,
) -> Tuple[BaselineCheck, ...]:
    payload = _read_baseline(baseline)
    settings = payload.get("settings", {})
    ignore_tags = settings.get("ignore_tags", [])
    ignore_info = settings.get("ignore_info", [])
    entries = read_manifest(manifest, profile, reference, normalize_variants)
    expected = {
        (item.get("sample"), item.get("stage")): item
        for item in payload["comparisons"]
        if isinstance(item, dict)
    }
    actual_keys = {(entry.sample, entry.stage) for entry in entries}
    if set(expected) != actual_keys:
        missing = sorted(set(expected) - actual_keys)
        extra = sorted(actual_keys - set(expected))
        messages = []
        if missing:
            messages.append(f"missing comparisons: {missing}")
        if extra:
            messages.append(f"unexpected comparisons: {extra}")
        raise BaselineError("; ".join(messages))
    checks = []
    for entry in entries:
        previous = expected[(entry.sample, entry.stage)]
        if entry.profile != previous.get("profile"):
            raise BaselineError(
                f"{entry.sample} / {entry.stage}: profile differs from baseline"
            )
        if entry.normalize != previous.get("normalize"):
            raise BaselineError(
                f"{entry.sample} / {entry.stage}: normalization differs from baseline"
            )
        current = _scan_entry(
            entry, entry.after, threads, ignore_tags, ignore_info, progress
        )
        if current["reference_sha256"] != previous.get("reference_sha256"):
            raise BaselineError(
                f"{entry.sample} / {entry.stage}: reference differs from baseline"
            )
        previous_fingerprints = previous.get("fingerprints", {})
        changed = tuple(
            field
            for field in previous.get("fields", [])
            if current["fingerprints"].get(field) != previous_fingerprints.get(field)
        )
        content_equal = current["fingerprints"].get(
            "record"
        ) == previous_fingerprints.get("record")
        structure_equal = current["structural_header"] == previous.get(
            "structural_header"
        )
        reasons = []
        if not content_equal:
            reasons.append("logical records differ")
        if not structure_equal:
            reasons.append("structural header differs")
        checks.append(
            BaselineCheck(
                sample=entry.sample,
                stage=entry.stage,
                passed=not reasons,
                record_count=current["record_count"],
                changed_fields=changed,
                structure_equal=structure_equal,
                reasons=tuple(reasons),
            )
        )
    return tuple(checks)


def baseline_report(checks: Sequence[BaselineCheck]) -> Dict[str, Any]:
    passed = sum(item.passed for item in checks)
    return {
        "semantiseq_baseline_check_format": 1,
        "passed": passed == len(checks),
        "summary": {
            "comparisons": len(checks),
            "passed": passed,
            "failed": len(checks) - passed,
        },
        "comparisons": [item.to_dict() for item in checks],
    }
