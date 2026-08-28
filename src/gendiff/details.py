"""Disk-backed record matching for explainable comparisons."""

from __future__ import annotations

import json
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

from gendiff.fingerprint import normalize
from gendiff.model import DifferenceDetails

_DIGEST_BYTES = 32


class DetailWriter:
    """Write compact comparison records without retaining them in memory."""

    def __init__(self, path: Path) -> None:
        self._connection = sqlite3.connect(str(path))
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute(
            "CREATE TABLE records ("
            "identity TEXT NOT NULL, full BLOB NOT NULL, "
            "fields BLOB NOT NULL, summary TEXT NOT NULL)"
        )
        self._batch: List[Tuple[str, bytes, bytes, str]] = []

    def add(
        self,
        identity: Any,
        full_digest: int,
        field_digests: Sequence[int],
        summary: Dict[str, Any],
    ) -> None:
        self._batch.append(
            (
                json.dumps(normalize(identity), separators=(",", ":")),
                full_digest.to_bytes(_DIGEST_BYTES, "big"),
                b"".join(
                    value.to_bytes(_DIGEST_BYTES, "big") for value in field_digests
                ),
                json.dumps(normalize(summary), separators=(",", ":"), sort_keys=True),
            )
        )
        if len(self._batch) >= 1000:
            self._flush()

    def _flush(self) -> None:
        if not self._batch:
            return
        self._connection.executemany(
            "INSERT INTO records VALUES (?, ?, ?, ?)", self._batch
        )
        self._batch.clear()

    def close(self) -> None:
        self._flush()
        self._connection.execute("CREATE INDEX records_identity ON records(identity)")
        self._connection.commit()
        self._connection.close()


@dataclass(frozen=True)
class _StoredRecord:
    full: bytes
    fields: bytes
    summary: str


def _groups(path: Path) -> Iterator[Tuple[str, List[_StoredRecord]]]:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        cursor = connection.execute(
            "SELECT identity, full, fields, summary "
            "FROM records ORDER BY identity, full"
        )
        current: Optional[str] = None
        records: List[_StoredRecord] = []
        for identity, full, fields, summary in cursor:
            if current is not None and identity != current:
                yield current, records
                records = []
            current = identity
            records.append(_StoredRecord(full, fields, summary))
        if current is not None:
            yield current, records
    finally:
        connection.close()


def _without_exact(
    left: Sequence[_StoredRecord], right: Sequence[_StoredRecord]
) -> Tuple[int, List[_StoredRecord], List[_StoredRecord]]:
    right_by_digest: Dict[bytes, List[int]] = defaultdict(list)
    for index, record in enumerate(right):
        right_by_digest[record.full].append(index)
    used_right = set()
    remaining_left = []
    identical = 0
    for record in left:
        matches = right_by_digest.get(record.full)
        if matches:
            used_right.add(matches.pop())
            identical += 1
        else:
            remaining_left.append(record)
    remaining_right = [
        record for index, record in enumerate(right) if index not in used_right
    ]
    return identical, remaining_left, remaining_right


def _field_parts(record: _StoredRecord) -> Tuple[bytes, ...]:
    return tuple(
        record.fields[offset : offset + _DIGEST_BYTES]
        for offset in range(0, len(record.fields), _DIGEST_BYTES)
    )


def _pair_modified(
    left: List[_StoredRecord], right: List[_StoredRecord]
) -> Tuple[
    List[Tuple[_StoredRecord, _StoredRecord]], List[_StoredRecord], List[_StoredRecord]
]:
    pairs = []
    remaining_left = list(left)
    remaining_right = list(right)
    while remaining_left and remaining_right:
        left_record = remaining_left.pop()
        left_fields = _field_parts(left_record)
        best_index = min(
            range(len(remaining_right)),
            key=lambda index: sum(
                first != second
                for first, second in zip(
                    left_fields, _field_parts(remaining_right[index])
                )
            ),
        )
        pairs.append((left_record, remaining_right.pop(best_index)))
    return pairs, remaining_left, remaining_right


def _location(summary: Dict[str, Any]) -> Optional[Tuple[str, int]]:
    contig = summary.get("reference", summary.get("contig"))
    position = summary.get("position")
    if contig is None or position is None or position < 0:
        return None
    return str(contig), int(position)


def _record_location(
    summary: Dict[str, Any], regions: Counter, loci: List[str], max_loci: int
) -> None:
    location = _location(summary)
    if location is None:
        return
    contig, position = location
    start = ((position - 1) // 1_000_000) * 1_000_000
    regions[f"{contig}:{start + 1}-{start + 1_000_000}"] += 1
    if len(loci) < max_loci:
        loci.append(f"{contig}:{max(1, position - 100)}-{position + 100}")


def _example(
    identity: str,
    status: str,
    left: Optional[_StoredRecord],
    right: Optional[_StoredRecord],
    changed_fields: Sequence[str],
) -> Dict[str, Any]:
    return {
        "identity": json.loads(identity),
        "status": status,
        "changed_fields": list(changed_fields),
        "left": json.loads(left.summary) if left else None,
        "right": json.loads(right.summary) if right else None,
    }


def analyze_details(
    left_path: Path,
    right_path: Path,
    fields: Sequence[str],
    max_examples: int,
) -> DifferenceDetails:
    left_groups = iter(_groups(left_path))
    right_groups = iter(_groups(right_path))
    left_group = next(left_groups, None)
    right_group = next(right_groups, None)
    identical = modified = left_only = right_only = 0
    field_changes = Counter()
    regions = Counter()
    sample_changes = Counter()
    examples: List[Dict[str, Any]] = []
    loci: List[str] = []

    while left_group is not None or right_group is not None:
        if right_group is None or (
            left_group is not None and left_group[0] < right_group[0]
        ):
            identity, left_records = left_group
            for record in left_records:
                left_only += 1
                summary = json.loads(record.summary)
                _record_location(summary, regions, loci, max_examples)
                if len(examples) < max_examples:
                    examples.append(_example(identity, "left_only", record, None, ()))
            left_group = next(left_groups, None)
            continue
        if left_group is None or right_group[0] < left_group[0]:
            identity, right_records = right_group
            for record in right_records:
                right_only += 1
                summary = json.loads(record.summary)
                _record_location(summary, regions, loci, max_examples)
                if len(examples) < max_examples:
                    examples.append(_example(identity, "right_only", None, record, ()))
            right_group = next(right_groups, None)
            continue

        identity = left_group[0]
        matched, remaining_left, remaining_right = _without_exact(
            left_group[1], right_group[1]
        )
        identical += matched
        pairs, remaining_left, remaining_right = _pair_modified(
            remaining_left, remaining_right
        )
        for left_record, right_record in pairs:
            modified += 1
            changed = [
                name
                for name, first, second in zip(
                    fields, _field_parts(left_record), _field_parts(right_record)
                )
                if first != second
            ]
            field_changes.update(changed)
            left_summary = json.loads(left_record.summary)
            right_summary = json.loads(right_record.summary)
            _record_location(right_summary, regions, loci, max_examples)
            left_genotypes = left_summary.get("genotypes", {})
            right_genotypes = right_summary.get("genotypes", {})
            for sample in set(left_genotypes) | set(right_genotypes):
                if left_genotypes.get(sample) != right_genotypes.get(sample):
                    sample_changes[sample] += 1
            if len(examples) < max_examples:
                examples.append(
                    _example(identity, "modified", left_record, right_record, changed)
                )
        for record, side in [
            *((record, "left_only") for record in remaining_left),
            *((record, "right_only") for record in remaining_right),
        ]:
            if side == "left_only":
                left_only += 1
            else:
                right_only += 1
            summary = json.loads(record.summary)
            _record_location(summary, regions, loci, max_examples)
            if len(examples) < max_examples:
                examples.append(
                    _example(
                        identity,
                        side,
                        record if side == "left_only" else None,
                        record if side == "right_only" else None,
                        (),
                    )
                )
        left_group = next(left_groups, None)
        right_group = next(right_groups, None)

    return DifferenceDetails(
        identical=identical,
        modified=modified,
        left_only=left_only,
        right_only=right_only,
        field_changes=dict(field_changes.most_common()),
        top_regions=[
            {"region": region, "changes": changes}
            for region, changes in regions.most_common(20)
        ],
        sample_changes=dict(sample_changes.most_common()),
        examples=examples,
        loci=loci,
    )


def relationship_label(overlap: float, content_equal: bool, empty: bool) -> str:
    if empty:
        return "empty datasets"
    if content_equal:
        return "identical logical records"
    if overlap >= 0.98:
        return "likely the same dataset"
    if overlap >= 0.20:
        return "related datasets"
    return "different datasets"
