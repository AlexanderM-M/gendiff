"""Pipeline regression comparisons driven by a small TSV manifest."""

from __future__ import annotations

import csv
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from gendiff.artifacts import output_workspace, safe_slug, write_manifest
from gendiff.compare import GenDiffError, compare_files
from gendiff.model import ComparisonResult

_REQUIRED_COLUMNS = {"sample", "stage", "before", "after"}
_OPTIONAL_COLUMNS = {"profile", "reference", "normalize"}


class BatchError(ValueError):
    """Raised when a batch manifest or policy is invalid."""


@dataclass(frozen=True)
class BatchEntry:
    sample: str
    stage: str
    before: Path
    after: Path
    profile: str
    reference: Optional[Path]
    normalize: bool
    position: int


@dataclass(frozen=True)
class BatchPolicy:
    max_modified: int = 0
    max_only: int = 0
    min_overlap: float = 0.0
    allow_structural_differences: bool = False
    fail_transitions: Tuple[str, ...] = ()

    def to_dict(self) -> Dict[str, Any]:
        return {
            "max_modified": self.max_modified,
            "max_only": self.max_only,
            "min_overlap": self.min_overlap,
            "allow_structural_differences": self.allow_structural_differences,
            "fail_transitions": list(self.fail_transitions),
        }


@dataclass(frozen=True)
class BatchItem:
    entry: BatchEntry
    result: ComparisonResult
    passed: bool
    reasons: Tuple[str, ...]
    duration_seconds: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample": self.entry.sample,
            "stage": self.entry.stage,
            "passed": self.passed,
            "reasons": list(self.reasons),
            "duration_seconds": round(self.duration_seconds, 6),
            "comparison": self.result.to_dict(),
        }


@dataclass(frozen=True)
class BatchResult:
    manifest: Path
    policy: BatchPolicy
    items: Tuple[BatchItem, ...]
    transitions: Dict[str, int] = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        return all(item.passed for item in self.items)

    @property
    def earliest_divergence(self) -> Dict[str, str]:
        earliest: Dict[str, str] = {}
        for item in self.items:
            if not item.result.equivalent and item.entry.sample not in earliest:
                earliest[item.entry.sample] = item.entry.stage
        return earliest

    @property
    def stage_order(self) -> Tuple[str, ...]:
        return tuple(dict.fromkeys(item.entry.stage for item in self.items))

    @property
    def sample_order(self) -> Tuple[str, ...]:
        return tuple(dict.fromkeys(item.entry.sample for item in self.items))

    def to_dict(self) -> Dict[str, Any]:
        passed = sum(item.passed for item in self.items)
        return {
            "passed": self.passed,
            "manifest": str(self.manifest),
            "summary": {
                "comparisons": len(self.items),
                "passed": passed,
                "failed": len(self.items) - passed,
            },
            "policy": self.policy.to_dict(),
            "stage_order": list(self.stage_order),
            "earliest_divergence": self.earliest_divergence,
            "transitions": self.transitions,
            "comparisons": [item.to_dict() for item in self.items],
        }


def _manifest_path(manifest: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = manifest.parent / path
    return path.absolute()


def _boolean(value: str, line: int) -> bool:
    normalized = value.strip().lower()
    if normalized in {"", "false", "no", "0"}:
        return False
    if normalized in {"true", "yes", "1"}:
        return True
    raise BatchError(f"line {line}: normalize must be true or false")


def read_manifest(
    manifest: Path,
    default_profile: str,
    default_reference: Optional[Path],
    normalize_variants: bool,
) -> Tuple[BatchEntry, ...]:
    if not manifest.is_file():
        raise BatchError(f"manifest not found: {manifest}")
    with manifest.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise BatchError("manifest is empty")
        columns = {column.strip() for column in reader.fieldnames}
        missing = sorted(_REQUIRED_COLUMNS - columns)
        unknown = sorted(columns - _REQUIRED_COLUMNS - _OPTIONAL_COLUMNS)
        if missing:
            raise BatchError(f"manifest is missing columns: {', '.join(missing)}")
        if unknown:
            raise BatchError(f"manifest has unknown columns: {', '.join(unknown)}")

        entries: List[BatchEntry] = []
        seen = set()
        for line, raw in enumerate(reader, 2):
            if None in raw:
                raise BatchError(f"line {line}: too many columns")
            row = {key.strip(): (value or "").strip() for key, value in raw.items()}
            if not any(row.values()):
                continue
            if row["sample"].startswith("#"):
                continue
            for column in _REQUIRED_COLUMNS:
                if not row[column]:
                    raise BatchError(f"line {line}: {column} cannot be empty")
            identity = (row["sample"], row["stage"])
            if identity in seen:
                raise BatchError(
                    f"line {line}: duplicate sample and stage: "
                    f"{row['sample']} / {row['stage']}"
                )
            seen.add(identity)
            reference_text = row.get("reference", "")
            reference = (
                _manifest_path(manifest, reference_text)
                if reference_text
                else default_reference
            )
            entries.append(
                BatchEntry(
                    sample=row["sample"],
                    stage=row["stage"],
                    before=_manifest_path(manifest, row["before"]),
                    after=_manifest_path(manifest, row["after"]),
                    profile=row.get("profile", "") or default_profile,
                    reference=reference,
                    normalize=(
                        _boolean(row["normalize"], line)
                        if row.get("normalize", "")
                        else normalize_variants
                    ),
                    position=len(entries),
                )
            )
    if not entries:
        raise BatchError("manifest contains no comparisons")
    return tuple(entries)


def _evaluate(result: ComparisonResult, policy: BatchPolicy) -> Tuple[str, ...]:
    reasons = []
    details = result.details
    if details is None:
        raise BatchError("batch comparisons require detailed record matching")
    if details.modified > policy.max_modified:
        reasons.append(
            f"modified records {details.modified:,} exceed {policy.max_modified:,}"
        )
    only = details.left_only + details.right_only
    if only > policy.max_only:
        reasons.append(f"sample-only records {only:,} exceed {policy.max_only:,}")
    if result.identity_overlap < policy.min_overlap:
        reasons.append(
            f"identity overlap {result.identity_overlap:.1%} is below "
            f"{policy.min_overlap:.1%}"
        )
    if not result.structure_equal and not policy.allow_structural_differences:
        reasons.append("structural header differs")
    for pattern in policy.fail_transitions:
        count = sum(
            value
            for transition, value in details.transitions.items()
            if pattern.casefold() in transition.casefold()
        )
        if count:
            reasons.append(f"forbidden transition {pattern!r} occurred {count:,} times")
    return tuple(reasons)


def _compare_entry(
    entry: BatchEntry,
    policy: BatchPolicy,
    threads: int,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    max_examples: int,
    progress: bool,
    temp_dir: Optional[Path],
    diff_dir: Optional[Path],
) -> BatchItem:
    started = time.perf_counter()
    try:
        result = compare_files(
            entry.before,
            entry.after,
            entry.reference,
            threads,
            profile=entry.profile,
            ignore_tags=ignore_tags,
            ignore_info=ignore_info,
            normalize_variants=entry.normalize,
            explain=True,
            max_examples=max_examples,
            progress=progress,
            temp_dir=temp_dir,
            left_label=f"{entry.sample} baseline",
            right_label=f"{entry.sample} candidate",
            diff_dir=diff_dir,
        )
    except (GenDiffError, OSError, ValueError) as error:
        raise BatchError(f"{entry.sample} / {entry.stage}: {error}") from error
    reasons = _evaluate(result, policy)
    return BatchItem(
        entry=entry,
        result=result,
        passed=not reasons,
        reasons=reasons,
        duration_seconds=time.perf_counter() - started,
    )


def _run_entries(
    entries: Sequence[BatchEntry],
    policy: BatchPolicy,
    threads: int,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    max_examples: int,
    progress: bool,
    temp_dir: Optional[Path],
    diff_workspace: Optional[Path],
) -> Tuple[BatchItem, ...]:
    workers = min(len(entries), max(1, threads // 2))
    threads_per_comparison = max(1, threads // workers)

    def compare(entry: BatchEntry) -> BatchItem:
        output = None
        if diff_workspace is not None:
            name = (
                f"{entry.position + 1:03d}-{safe_slug(entry.sample)}-"
                f"{safe_slug(entry.stage)}"
            )
            output = diff_workspace / name
        return _compare_entry(
            entry,
            policy,
            threads_per_comparison,
            ignore_tags,
            ignore_info,
            max_examples,
            progress,
            temp_dir,
            output,
        )

    if workers == 1:
        return tuple(compare(entry) for entry in entries)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        return tuple(executor.map(compare, entries))


def _transition_totals(items: Sequence[BatchItem]) -> Dict[str, int]:
    transitions = Counter()
    for item in items:
        if item.result.details is not None:
            transitions.update(item.result.details.transitions)
    return dict(transitions.most_common(100))


def _published_items(
    items: Sequence[BatchItem], workspace: Path, target: Path
) -> Tuple[BatchItem, ...]:
    published = target.absolute()
    updated = []
    for item in items:
        artifacts = {
            key: str(published / Path(value).relative_to(workspace))
            for key, value in item.result.artifacts.items()
        }
        updated.append(replace(item, result=replace(item.result, artifacts=artifacts)))
    return tuple(updated)


def compare_manifest(
    manifest: Path,
    *,
    threads: int,
    profile: str,
    reference: Optional[Path],
    normalize_variants: bool,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    max_examples: int,
    progress: bool,
    temp_dir: Optional[Path],
    diff_dir: Optional[Path],
    force: bool,
    policy: BatchPolicy,
) -> BatchResult:
    entries = read_manifest(manifest, profile, reference, normalize_variants)
    if diff_dir is None:
        items = _run_entries(
            entries,
            policy,
            threads,
            ignore_tags,
            ignore_info,
            max_examples,
            progress,
            temp_dir,
            None,
        )
        return BatchResult(
            manifest=manifest,
            policy=policy,
            items=items,
            transitions=_transition_totals(items),
        )

    with output_workspace(diff_dir, force) as workspace:
        items = _run_entries(
            entries,
            policy,
            threads,
            ignore_tags,
            ignore_info,
            max_examples,
            progress,
            temp_dir,
            workspace,
        )
        items = _published_items(items, workspace, diff_dir)
        write_manifest(
            workspace / "manifest.json",
            {
                "gendiff_diff_format": 1,
                "kind": "batch",
                "comparisons": [
                    {
                        "sample": item.entry.sample,
                        "stage": item.entry.stage,
                        "directory": str(
                            Path(item.result.artifacts["manifest"]).parent.relative_to(
                                diff_dir.absolute()
                            )
                        ),
                    }
                    for item in items
                ],
            },
        )
    return BatchResult(
        manifest=manifest,
        policy=policy,
        items=items,
        transitions=_transition_totals(items),
    )
