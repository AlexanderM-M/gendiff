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


def identity_text(identity: Any) -> str:
    return json.dumps(normalize(identity), separators=(",", ":"))


def digest_bytes(value: int) -> bytes:
    return value.to_bytes(_DIGEST_BYTES, "big")


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
                identity_text(identity),
                digest_bytes(full_digest),
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


class SelectionWriter:
    """Store the exact record multiplicities selected for diff outputs."""

    def __init__(self, path: Path) -> None:
        self._connection = sqlite3.connect(str(path))
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute(
            "CREATE TABLE selections ("
            "side INTEGER NOT NULL, identity TEXT NOT NULL, full BLOB NOT NULL, "
            "status TEXT NOT NULL, count INTEGER NOT NULL, "
            "PRIMARY KEY (side, identity, full, status))"
        )

    def add(
        self,
        side: int,
        identity: str,
        record: "_StoredRecord",
        status: str,
    ) -> None:
        self._connection.execute(
            "INSERT INTO selections VALUES (?, ?, ?, ?, 1) "
            "ON CONFLICT(side, identity, full, status) "
            "DO UPDATE SET count=count+1",
            (side, identity, record.full, status),
        )

    def close(self) -> None:
        self._connection.commit()
        self._connection.close()


class SelectionReader:
    """Consume selected records while rescanning an original input."""

    def __init__(self, path: Path, side: int) -> None:
        self._connection = sqlite3.connect(str(path))
        self._side = side

    def take(self, identity: Any, full_digest: int) -> Optional[str]:
        encoded_identity = identity_text(identity)
        encoded_digest = digest_bytes(full_digest)
        row = self._connection.execute(
            "SELECT status, count FROM selections "
            "WHERE side=? AND identity=? AND full=? AND count>0 "
            "ORDER BY status LIMIT 1",
            (self._side, encoded_identity, encoded_digest),
        ).fetchone()
        if row is None:
            return None
        status, count = row
        self._connection.execute(
            "UPDATE selections SET count=? "
            "WHERE side=? AND identity=? AND full=? AND status=?",
            (
                count - 1,
                self._side,
                encoded_identity,
                encoded_digest,
                status,
            ),
        )
        return str(status)

    def remaining(self) -> int:
        row = self._connection.execute(
            "SELECT COALESCE(SUM(count), 0) FROM selections WHERE side=? AND count>0",
            (self._side,),
        ).fetchone()
        return int(row[0])

    def close(self) -> None:
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
    def visible(record: Optional[_StoredRecord]) -> Optional[Dict[str, Any]]:
        if record is None:
            return None
        summary = json.loads(record.summary)
        return {key: value for key, value in summary.items() if not key.startswith("_")}

    return {
        "identity": json.loads(identity),
        "status": status,
        "changed_fields": list(changed_fields),
        "left": visible(left),
        "right": visible(right),
    }


def _value_text(value: Any) -> str:
    if value is None:
        return "missing"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, list):
        return "/".join("." if item is None else str(item) for item in value)
    return str(value)


def _flag_states(flags: int) -> Dict[str, str]:
    alignment_class = "primary"
    if flags & 0x800:
        alignment_class = "supplementary"
    elif flags & 0x100:
        alignment_class = "secondary"
    return {
        "mapping status": "unmapped" if flags & 0x4 else "mapped",
        "alignment class": alignment_class,
        "duplicate flag": "set" if flags & 0x400 else "clear",
        "QC status": "failed" if flags & 0x200 else "passed",
        "strand": "reverse" if flags & 0x10 else "forward",
    }


def _transition_counts(
    left: Dict[str, Any],
    right: Dict[str, Any],
    changed_fields: Sequence[str],
) -> Counter:
    transitions = Counter()
    changed = set(changed_fields)
    if "flags" in changed:
        left_states = _flag_states(int(left.get("flags", 0)))
        right_states = _flag_states(int(right.get("flags", 0)))
        for name, first in left_states.items():
            second = right_states[name]
            if first != second:
                transitions[f"{name}: {first} → {second}"] += 1

    direct = {
        "mapping_quality": "MAPQ",
        "reference": "reference",
        "cigar": "CIGAR",
        "quality": "QUAL",
        "filters": "filters",
    }
    for field, label in direct.items():
        if field in changed:
            first = _value_text(left.get(field))
            second = _value_text(right.get(field))
            transitions[f"{label}: {first} → {second}"] += 1

    if "position" in changed:
        transitions["position changed"] += 1
    for field in ("sequence", "qualities", "mate", "template_length"):
        if field in changed:
            transitions[f"{field.replace('_', ' ')} changed"] += 1

    if "tags" in changed:
        left_values = left.get("_tag_digests", {})
        right_values = right.get("_tag_digests", {})
        left_tags = set(left_values)
        right_tags = set(right_values)
        for tag in sorted(right_tags - left_tags):
            transitions[f"tag added: {tag}"] += 1
        for tag in sorted(left_tags - right_tags):
            transitions[f"tag removed: {tag}"] += 1
        for tag in sorted(left_tags & right_tags):
            if left_values[tag] != right_values[tag]:
                transitions[f"tag value changed: {tag}"] += 1

    if "info" in changed:
        left_info = left.get("_info", {})
        right_info = right.get("_info", {})
        left_keys = set(left_info)
        right_keys = set(right_info)
        for key in sorted(right_keys - left_keys):
            transitions[f"INFO added: {key}"] += 1
        for key in sorted(left_keys - right_keys):
            transitions[f"INFO removed: {key}"] += 1
        for key in sorted(left_keys & right_keys):
            if left_info[key] != right_info[key]:
                transitions[f"INFO value changed: {key}"] += 1

    if "samples" in changed:
        left_genotypes = left.get("genotypes", {})
        right_genotypes = right.get("genotypes", {})
        genotype_transition = False
        for sample in sorted(set(left_genotypes) | set(right_genotypes)):
            first = left_genotypes.get(sample, {}).get("GT")
            second = right_genotypes.get(sample, {}).get("GT")
            if first != second:
                transitions[
                    f"genotype {sample}: {_value_text(first)} → {_value_text(second)}"
                ] += 1
                genotype_transition = True
        if not genotype_transition:
            transitions["sample FORMAT values changed"] += 1

    handled = set(direct) | {
        "flags",
        "position",
        "sequence",
        "qualities",
        "mate",
        "template_length",
        "tags",
        "info",
        "samples",
    }
    for field in changed - handled:
        transitions[f"{field.replace('_', ' ')} changed"] += 1
    return transitions


def analyze_details(
    left_path: Path,
    right_path: Path,
    fields: Sequence[str],
    max_examples: int,
    selection_path: Optional[Path] = None,
) -> DifferenceDetails:
    left_groups = iter(_groups(left_path))
    right_groups = iter(_groups(right_path))
    left_group = next(left_groups, None)
    right_group = next(right_groups, None)
    identical = modified = left_only = right_only = 0
    field_changes = Counter()
    regions = Counter()
    sample_changes = Counter()
    transitions = Counter()
    examples: List[Dict[str, Any]] = []
    loci: List[str] = []
    selector = SelectionWriter(selection_path) if selection_path else None

    try:
        while left_group is not None or right_group is not None:
            if right_group is None or (
                left_group is not None and left_group[0] < right_group[0]
            ):
                identity, left_records = left_group
                for record in left_records:
                    left_only += 1
                    if selector:
                        selector.add(0, identity, record, "only")
                    summary = json.loads(record.summary)
                    _record_location(summary, regions, loci, max_examples)
                    if len(examples) < max_examples:
                        examples.append(
                            _example(identity, "left_only", record, None, ())
                        )
                left_group = next(left_groups, None)
                continue
            if left_group is None or right_group[0] < left_group[0]:
                identity, right_records = right_group
                for record in right_records:
                    right_only += 1
                    if selector:
                        selector.add(1, identity, record, "only")
                    summary = json.loads(record.summary)
                    _record_location(summary, regions, loci, max_examples)
                    if len(examples) < max_examples:
                        examples.append(
                            _example(identity, "right_only", None, record, ())
                        )
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
                if selector:
                    selector.add(0, identity, left_record, "modified")
                    selector.add(1, identity, right_record, "modified")
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
                transitions.update(
                    _transition_counts(left_summary, right_summary, changed)
                )
                _record_location(right_summary, regions, loci, max_examples)
                left_genotypes = left_summary.get("genotypes", {})
                right_genotypes = right_summary.get("genotypes", {})
                for sample in set(left_genotypes) | set(right_genotypes):
                    if left_genotypes.get(sample) != right_genotypes.get(sample):
                        sample_changes[sample] += 1
                if len(examples) < max_examples:
                    examples.append(
                        _example(
                            identity,
                            "modified",
                            left_record,
                            right_record,
                            changed,
                        )
                    )
            for record, side in [
                *((record, "left_only") for record in remaining_left),
                *((record, "right_only") for record in remaining_right),
            ]:
                if side == "left_only":
                    left_only += 1
                    if selector:
                        selector.add(0, identity, record, "only")
                else:
                    right_only += 1
                    if selector:
                        selector.add(1, identity, record, "only")
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
    finally:
        if selector:
            selector.close()

    return DifferenceDetails(
        identical=identical,
        modified=modified,
        left_only=left_only,
        right_only=right_only,
        field_changes=dict(field_changes.most_common()),
        transitions=dict(transitions.most_common(50)),
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
