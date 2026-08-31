"""Disk-backed record matching for explainable comparisons."""

from __future__ import annotations

import json
import math
import sqlite3
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple

from gendiff.difference_table import DifferenceTableWriter
from gendiff.fingerprint import normalize
from gendiff.metrics import MetricCounts, distribution_shifts
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


def merge_detail_databases(sources: Sequence[Path], target: Path) -> None:
    """Merge independently scanned detail shards into one indexed database."""
    output = sqlite3.connect(str(target))
    output.execute("PRAGMA journal_mode=OFF")
    output.execute("PRAGMA synchronous=OFF")
    output.execute(
        "CREATE TABLE records (identity TEXT NOT NULL, full BLOB NOT NULL, "
        "fields BLOB NOT NULL, summary TEXT NOT NULL)"
    )
    try:
        for source in sources:
            connection = sqlite3.connect(f"file:{source}?mode=ro", uri=True)
            try:
                cursor = connection.execute(
                    "SELECT identity, full, fields, summary FROM records"
                )
                while rows := cursor.fetchmany(2000):
                    output.executemany("INSERT INTO records VALUES (?, ?, ?, ?)", rows)
            finally:
                connection.close()
        output.execute("CREATE INDEX records_identity ON records(identity)")
        output.commit()
    finally:
        output.close()


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


class TrackWriter:
    """Aggregate changed genomic positions in a disk-backed table."""

    def __init__(self, path: Path) -> None:
        self._connection = sqlite3.connect(str(path))
        self._connection.execute("PRAGMA journal_mode=OFF")
        self._connection.execute("PRAGMA synchronous=OFF")
        self._connection.execute(
            "CREATE TABLE loci ("
            "contig TEXT NOT NULL, start INTEGER NOT NULL, status TEXT NOT NULL, "
            "count INTEGER NOT NULL, PRIMARY KEY (contig, start, status))"
        )

    def add(self, summary: Dict[str, Any], status: str) -> None:
        location = _location(summary)
        if location is None:
            return
        contig, position = location
        self._connection.execute(
            "INSERT INTO loci VALUES (?, ?, ?, 1) "
            "ON CONFLICT(contig, start, status) DO UPDATE SET count=count+1",
            (contig, position - 1, status),
        )

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
    regions[(contig, start)] += 1
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
    track_path: Optional[Path] = None,
    table_path: Optional[Path] = None,
    left_contigs: Optional[Dict[str, int]] = None,
    right_contigs: Optional[Dict[str, int]] = None,
    lengths: Optional[Dict[str, int]] = None,
    left_windows: Optional[Dict[Tuple[str, int], int]] = None,
    right_windows: Optional[Dict[Tuple[str, int], int]] = None,
    left_metrics: Optional[MetricCounts] = None,
    right_metrics: Optional[MetricCounts] = None,
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
    tracks = TrackWriter(track_path) if track_path else None
    table = DifferenceTableWriter(table_path) if table_path else None

    def capture(
        identity: str,
        status: str,
        left: Optional[_StoredRecord],
        right: Optional[_StoredRecord],
        changed: Sequence[str],
    ) -> None:
        if table is None and len(examples) >= max_examples:
            return
        item = _example(identity, status, left, right, changed)
        if table:
            table.add(
                {"left_only": "only_in_first", "right_only": "only_in_second"}.get(
                    status, status
                ),
                item["identity"],
                item["changed_fields"],
                item["left"],
                item["right"],
            )
        if len(examples) < max_examples:
            examples.append(item)

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
                    if tracks:
                        tracks.add(summary, "only_in_first")
                    _record_location(summary, regions, loci, max_examples)
                    capture(identity, "left_only", record, None, ())
                left_group = next(left_groups, None)
                continue
            if left_group is None or right_group[0] < left_group[0]:
                identity, right_records = right_group
                for record in right_records:
                    right_only += 1
                    if selector:
                        selector.add(1, identity, record, "only")
                    summary = json.loads(record.summary)
                    if tracks:
                        tracks.add(summary, "only_in_second")
                    _record_location(summary, regions, loci, max_examples)
                    capture(identity, "right_only", None, record, ())
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
                if tracks:
                    tracks.add(right_summary, "modified")
                    if _location(left_summary) != _location(right_summary):
                        tracks.add(left_summary, "modified_previous")
                transitions.update(
                    _transition_counts(left_summary, right_summary, changed)
                )
                _record_location(right_summary, regions, loci, max_examples)
                left_genotypes = left_summary.get("genotypes", {})
                right_genotypes = right_summary.get("genotypes", {})
                for sample in set(left_genotypes) | set(right_genotypes):
                    if left_genotypes.get(sample) != right_genotypes.get(sample):
                        sample_changes[sample] += 1
                capture(identity, "modified", left_record, right_record, changed)
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
                if tracks:
                    tracks.add(
                        summary,
                        f"only_in_{'first' if side == 'left_only' else 'second'}",
                    )
                _record_location(summary, regions, loci, max_examples)
                capture(
                    identity,
                    side,
                    record if side == "left_only" else None,
                    record if side == "right_only" else None,
                    (),
                )
            left_group = next(left_groups, None)
            right_group = next(right_groups, None)
    finally:
        if selector:
            selector.close()
        if tracks:
            tracks.close()
        if table:
            table.close()

    left_counts = left_contigs or {}
    right_counts = right_contigs or {}
    contig_changes = Counter()
    for (contig, _), changes in regions.items():
        contig_changes[contig] += changes
    contig_stats = []
    for contig, changes in contig_changes.most_common():
        first = left_counts.get(contig, 0)
        second = right_counts.get(contig, 0)
        contig_stats.append(
            {
                "contig": contig,
                "changes": changes,
                "first_records": first,
                "second_records": second,
                "difference_fraction": changes / max(1, first + second),
                "coverage_ratio": math.log2((second + 1) / (first + 1)),
                "length": (lengths or {}).get(contig, 0),
            }
        )
    shifts = distribution_shifts(left_metrics or {}, right_metrics or {})
    findings = _findings(
        field_changes,
        modified,
        left_only,
        right_only,
        contig_changes,
        shifts,
    )
    first_windows = left_windows or {}
    second_windows = right_windows or {}
    window_keys = (
        set(regions)
        | {(contig, window * 1_000_000) for contig, window in first_windows}
        | {(contig, window * 1_000_000) for contig, window in second_windows}
    )

    return DifferenceDetails(
        identical=identical,
        modified=modified,
        left_only=left_only,
        right_only=right_only,
        field_changes=dict(field_changes.most_common()),
        transitions=dict(transitions.most_common()),
        top_regions=[
            {
                "region": f"{contig}:{start + 1}-{start + 1_000_000}",
                "changes": changes,
            }
            for (contig, start), changes in regions.most_common(20)
        ],
        sample_changes=dict(sample_changes.most_common()),
        examples=examples,
        loci=loci,
        region_density=[
            {
                "contig": contig,
                "start": start,
                "end": min(
                    start + 1_000_000,
                    (lengths or {}).get(contig, start + 1_000_000),
                ),
                "changes": regions.get((contig, start), 0),
                "first_records": first_windows.get((contig, start // 1_000_000), 0),
                "second_records": second_windows.get((contig, start // 1_000_000), 0),
                "difference_fraction": regions.get((contig, start), 0)
                / max(
                    1,
                    first_windows.get((contig, start // 1_000_000), 0)
                    + second_windows.get((contig, start // 1_000_000), 0),
                ),
                "coverage_ratio": math.log2(
                    (second_windows.get((contig, start // 1_000_000), 0) + 1)
                    / (first_windows.get((contig, start // 1_000_000), 0) + 1)
                ),
            }
            for contig, start in sorted(window_keys)
        ],
        contig_stats=contig_stats,
        distribution_shifts=shifts,
        findings=findings,
    )


def _findings(
    fields: Counter,
    modified: int,
    left_only: int,
    right_only: int,
    contigs: Counter,
    shifts: Sequence[Dict[str, Any]],
) -> List[str]:
    findings = []
    if fields:
        field, count = fields.most_common(1)[0]
        findings.append(
            f"Primary record change: {field.replace('_', ' ')} "
            f"({count / max(1, modified):.0%} of modified records)."
        )
    elif left_only or right_only:
        findings.append(
            f"Record membership changed: {left_only:,} only in the first input and "
            f"{right_only:,} only in the second."
        )
    if shifts:
        shift = shifts[0]
        if shift["kind"] == "rate":
            direction = (
                "increased"
                if shift["second_value"] > shift["first_value"]
                else "decreased"
            )
            findings.append(
                f"{shift['label']} {direction} from {shift['first_value']:.1%} to "
                f"{shift['second_value']:.1%}."
            )
        else:
            if shift["first_median"] == shift["second_median"]:
                findings.append(
                    f"{shift['label']} distribution changed by "
                    f"{shift['distance']:.0%}; median remained "
                    f"{shift['first_median']:g}."
                )
            else:
                findings.append(
                    f"{shift['label']} shifted: median {shift['first_median']:g} → "
                    f"{shift['second_median']:g}."
                )
    if contigs:
        contig, count = contigs.most_common(1)[0]
        total = sum(contigs.values())
        fraction = count / total
        if fraction >= 0.5:
            findings.append(
                f"Differences are concentrated on {contig} "
                f"({fraction:.0%} of located changes)."
            )
        elif fraction <= 0.2:
            findings.append(
                "Differences are genome-wide; no contig contributes more than "
                f"{fraction:.0%}."
            )
        else:
            findings.append(
                f"Differences span {len(contigs):,} contigs; {contig} contributes "
                f"the most ({fraction:.0%})."
            )
    return findings[:3]


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
