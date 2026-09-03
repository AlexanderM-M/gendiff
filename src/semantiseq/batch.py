"""Pipeline regression comparisons driven by a small TSV manifest."""

from __future__ import annotations

import csv
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from semantiseq.artifacts import output_workspace, safe_slug, write_manifest
from semantiseq.compare import SemantiSeqError, compare_files
from semantiseq.model import ComparisonResult
from semantiseq.policy import (
    BatchPolicy,
    PolicyDocument,
    evaluate_policy,
    resolve_policy,
)

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
class BatchItem:
    entry: BatchEntry
    result: ComparisonResult
    passed: bool
    reasons: Tuple[str, ...]
    warnings: Tuple[str, ...]
    policy_trace: Tuple[str, ...]
    matched_rules: Tuple[str, ...]
    effective_policy: BatchPolicy
    duration_seconds: float

    def to_dict(self) -> Dict[str, Any]:
        return {
            "sample": self.entry.sample,
            "stage": self.entry.stage,
            "passed": self.passed,
            "reasons": list(self.reasons),
            "warnings": list(self.warnings),
            "matched_policy_rules": list(self.matched_rules),
            "policy_trace": list(self.policy_trace),
            "effective_policy": self.effective_policy.to_dict(),
            "duration_seconds": round(self.duration_seconds, 6),
            "comparison": self.result.to_dict(),
        }


@dataclass(frozen=True)
class BatchResult:
    manifest: Path
    policy: BatchPolicy
    items: Tuple[BatchItem, ...]
    transitions: Dict[str, int] = field(default_factory=dict)
    policy_document: Optional[PolicyDocument] = None

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
            "semantiseq_batch_format": 1,
            "passed": self.passed,
            "manifest": str(self.manifest),
            "summary": {
                "comparisons": len(self.items),
                "passed": passed,
                "failed": len(self.items) - passed,
            },
            "policy": self.policy.to_dict(),
            "policy_document": (
                self.policy_document.to_dict() if self.policy_document else None
            ),
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
    track_dir: Optional[Path],
    policy_document: Optional[PolicyDocument],
    cache_dir: Optional[Path],
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
            track_dir=track_dir,
            cache_dir=cache_dir,
        )
    except (SemantiSeqError, OSError, ValueError) as error:
        raise BatchError(f"{entry.sample} / {entry.stage}: {error}") from error
    effective, matched = resolve_policy(
        policy,
        policy_document,
        sample=entry.sample,
        stage=entry.stage,
        kind=result.kind,
        profile=result.profile,
    )
    decision = evaluate_policy(result, effective, matched)
    return BatchItem(
        entry=entry,
        result=result,
        passed=not decision.errors,
        reasons=decision.errors,
        warnings=decision.warnings,
        policy_trace=decision.trace,
        matched_rules=decision.matched_rules,
        effective_policy=decision.effective,
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
    track_workspace: Optional[Path],
    policy_document: Optional[PolicyDocument],
    cache_dir: Optional[Path],
) -> Tuple[BatchItem, ...]:
    def output_for(workspace: Optional[Path], entry: BatchEntry) -> Optional[Path]:
        output = None
        if workspace is not None:
            name = (
                f"{entry.position + 1:03d}-{safe_slug(entry.sample)}-"
                f"{safe_slug(entry.stage)}"
            )
            output = workspace / name
        return output

    if len(entries) == 1:
        entry = entries[0]
        return (
            _compare_entry(
                entry,
                policy,
                threads,
                ignore_tags,
                ignore_info,
                max_examples,
                progress,
                temp_dir,
                output_for(diff_workspace, entry),
                output_for(track_workspace, entry),
                policy_document,
                cache_dir,
            ),
        )
    workers = min(len(entries), threads)
    if workers == 1:
        return tuple(
            _compare_entry(
                entry,
                policy,
                1,
                ignore_tags,
                ignore_info,
                max_examples,
                progress,
                temp_dir,
                output_for(diff_workspace, entry),
                output_for(track_workspace, entry),
                policy_document,
                cache_dir,
            )
            for entry in entries
        )
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _compare_entry,
                entry,
                policy,
                1,
                ignore_tags,
                ignore_info,
                max_examples,
                progress,
                temp_dir,
                output_for(diff_workspace, entry),
                output_for(track_workspace, entry),
                policy_document,
                cache_dir,
            )
            for entry in entries
        ]
        return tuple(future.result() for future in futures)


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
        artifacts = dict(item.result.artifacts)
        for key, value in item.result.artifacts.items():
            try:
                relative = Path(value).relative_to(workspace)
            except ValueError:
                continue
            artifacts[key] = str(published / relative)
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
    policy_document: Optional[PolicyDocument] = None,
    track_dir: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
) -> BatchResult:
    entries = read_manifest(manifest, profile, reference, normalize_variants)
    if diff_dir is None and track_dir is None:
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
            None,
            policy_document,
            cache_dir,
        )
        return BatchResult(
            manifest=manifest,
            policy=policy,
            items=items,
            transitions=_transition_totals(items),
            policy_document=policy_document,
        )

    with ExitStack() as stack:
        diff_workspace = (
            stack.enter_context(output_workspace(diff_dir, force)) if diff_dir else None
        )
        track_workspace = (
            stack.enter_context(output_workspace(track_dir, force))
            if track_dir
            else None
        )
        items = _run_entries(
            entries,
            policy,
            threads,
            ignore_tags,
            ignore_info,
            max_examples,
            progress,
            temp_dir,
            diff_workspace,
            track_workspace,
            policy_document,
            cache_dir,
        )
        if diff_workspace is not None and diff_dir is not None:
            items = _published_items(items, diff_workspace, diff_dir)
            write_manifest(
                diff_workspace / "manifest.json",
                {
                    "semantiseq_diff_format": 1,
                    "kind": "batch",
                    "comparisons": [
                        {
                            "sample": item.entry.sample,
                            "stage": item.entry.stage,
                            "directory": str(
                                Path(
                                    item.result.artifacts["manifest"]
                                ).parent.relative_to(diff_dir.absolute())
                            ),
                        }
                        for item in items
                    ],
                },
            )
        if track_workspace is not None and track_dir is not None:
            items = _published_items(items, track_workspace, track_dir)
            write_manifest(
                track_workspace / "manifest.json",
                {
                    "semantiseq_diff_format": 1,
                    "kind": "batch-tracks",
                    "comparisons": [
                        {
                            "sample": item.entry.sample,
                            "stage": item.entry.stage,
                            "directory": str(
                                Path(
                                    item.result.artifacts["track_manifest"]
                                ).parent.relative_to(track_dir.absolute())
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
        policy_document=policy_document,
    )
