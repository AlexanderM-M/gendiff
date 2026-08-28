from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Sequence, Set

import pysam

from gendiff.details import DetailWriter, analyze_details, relationship_label
from gendiff.fingerprint import (
    Fingerprint,
    FingerprintBuilder,
    Sketch,
    SketchBuilder,
    digest,
    digest_parts,
    merge_sketches,
    normalize,
    sketch_containment,
)
from gendiff.model import ComparisonResult
from gendiff.progress import ProgressReporter


@dataclass(frozen=True)
class _ScanResult:
    fingerprints: Dict[str, Fingerprint]
    count: int
    structural: Dict[str, Any]
    metadata: Dict[str, Any]
    sketch: Sketch


def _mode(path: Path) -> str:
    return "rc" if path.name.lower().endswith(".cram") else "rb"


def _ordered_header(header: pysam.AlignmentHeader) -> Dict[str, Any]:
    value = header.to_dict()
    result: Dict[str, Any] = {}
    for key, item in value.items():
        if key == "SQ":
            result[key] = sorted(item, key=lambda entry: entry.get("SN", ""))
        elif key in {"RG", "PG"}:
            result[key] = sorted(item, key=lambda entry: entry.get("ID", ""))
        elif key == "CO":
            result[key] = sorted(item)
        else:
            result[key] = item
    return normalize(result)


def _structural_header(header: pysam.AlignmentHeader) -> Dict[str, Any]:
    ordered = _ordered_header(header)
    return {key: ordered[key] for key in ("SQ", "RG") if key in ordered}


def _identity(record: pysam.AlignedSegment) -> Any:
    role_flags = record.flag & (0x40 | 0x80 | 0x100 | 0x800)
    return record.query_name, role_flags


def _record_values(
    record: pysam.AlignedSegment, ignore_tags: Set[str]
) -> Dict[str, Any]:
    tags = sorted(
        (
            tag
            for tag in record.get_tags(with_value_type=True)
            if tag[0] not in ignore_tags
        ),
        key=lambda tag: tag[0],
    )
    return {
        "read_name": record.query_name,
        "flags": record.flag,
        "reference": record.reference_name,
        "position": record.reference_start,
        "mapping_quality": record.mapping_quality,
        "cigar": record.cigartuples,
        "mate": (record.next_reference_name, record.next_reference_start),
        "template_length": record.template_length,
        "sequence": record.query_sequence,
        "qualities": record.query_qualities,
        "tags": tags,
    }


def _record_summary(
    record: pysam.AlignedSegment, values: Dict[str, Any]
) -> Dict[str, Any]:
    return {
        "read_name": record.query_name,
        "flags": record.flag,
        "reference": record.reference_name,
        "position": record.reference_start + 1 if record.reference_start >= 0 else -1,
        "mapping_quality": record.mapping_quality,
        "cigar": record.cigarstring,
        "query_length": record.query_length,
        "tags": [tag[0] for tag in values["tags"]],
    }


def _scan(
    path: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    ignore_tags: Sequence[str],
    details_path: Optional[Path],
    progress: bool,
    label: str,
    regions: Optional[Sequence[str]] = None,
) -> _ScanResult:
    builders = {name: FingerprintBuilder() for name in ("record",) + tuple(fields)}
    sketch = SketchBuilder()
    writer = DetailWriter(details_path) if details_path else None
    reporter = ProgressReporter(label, progress)
    ignored = set(ignore_tags)
    kwargs = {"reference_filename": str(reference)} if reference else {}
    kwargs["threads"] = threads
    count = 0
    try:
        with pysam.AlignmentFile(str(path), _mode(path), **kwargs) as handle:
            structural = _structural_header(handle.header)
            metadata = _ordered_header(handle.header)
            iterators = (
                (handle.fetch(region) for region in regions)
                if regions is not None
                else (handle.fetch(until_eof=True),)
            )
            for iterator in iterators:
                for record in iterator:
                    identity = _identity(record)
                    values = _record_values(record, ignored)
                    field_digests = []
                    for name in fields:
                        field_digest = digest(values[name])
                        builders[name].add_digest(field_digest)
                        field_digests.append(field_digest)
                    record_digest = digest_parts(field_digests)
                    builders["record"].add_digest(record_digest)
                    sketch.add(identity)
                    if writer:
                        writer.add(
                            identity,
                            record_digest,
                            field_digests,
                            _record_summary(record, values),
                        )
                    count += 1
                    reporter.update(count)
    finally:
        if writer:
            writer.close()
        reporter.finish(count)
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return _ScanResult(fingerprints, count, structural, metadata, sketch.finish())


def _region_chunks(
    path: Path, reference: Optional[Path], workers: int
) -> Optional[list[list[str]]]:
    kwargs = {"reference_filename": str(reference)} if reference else {}
    with pysam.AlignmentFile(str(path), _mode(path), **kwargs) as handle:
        if not handle.has_index():
            return None
        items = sorted(
            zip(handle.references, handle.lengths),
            key=lambda item: item[1],
            reverse=True,
        )
    bins: list[list[str]] = [[] for _ in range(min(workers, max(1, len(items))))]
    sizes = [0] * len(bins)
    for name, length in items:
        index = min(range(len(bins)), key=sizes.__getitem__)
        bins[index].append(name)
        sizes[index] += length
    bins[-1].append("*")
    return bins


def _merge_scans(scans: Sequence[_ScanResult], fields: Sequence[str]) -> _ScanResult:
    builders = {name: FingerprintBuilder() for name in ("record",) + tuple(fields)}
    for scan in scans:
        for name, fingerprint in scan.fingerprints.items():
            builders[name].merge(fingerprint)
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return _ScanResult(
        fingerprints,
        sum(scan.count for scan in scans),
        scans[0].structural,
        scans[0].metadata,
        merge_sketches(scan.sketch for scan in scans),
    )


def _scan_indexed_pair(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    ignore_tags: Sequence[str],
    progress: bool,
) -> Optional[tuple[_ScanResult, _ScanResult]]:
    workers_per_input = max(1, threads // 2)
    left_chunks = _region_chunks(left, reference, workers_per_input)
    right_chunks = _region_chunks(right, reference, workers_per_input)
    if left_chunks is None or right_chunks is None:
        return None
    task_count = len(left_chunks) + len(right_chunks)
    with ProcessPoolExecutor(max_workers=min(threads, task_count)) as executor:
        left_futures = [
            executor.submit(
                _scan,
                left,
                reference,
                1,
                fields,
                ignore_tags,
                None,
                progress,
                f"left:{index + 1}",
                regions,
            )
            for index, regions in enumerate(left_chunks)
        ]
        right_futures = [
            executor.submit(
                _scan,
                right,
                reference,
                1,
                fields,
                ignore_tags,
                None,
                progress,
                f"right:{index + 1}",
                regions,
            )
            for index, regions in enumerate(right_chunks)
        ]
        left_scans = [future.result() for future in left_futures]
        right_scans = [future.result() for future in right_futures]
    return _merge_scans(left_scans, fields), _merge_scans(right_scans, fields)


def _scan_pair(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    ignore_tags: Sequence[str],
    left_details: Optional[Path],
    right_details: Optional[Path],
    progress: bool,
) -> tuple[_ScanResult, _ScanResult]:
    if threads > 2 and left_details is None and right_details is None:
        indexed = _scan_indexed_pair(
            left, right, reference, threads, fields, ignore_tags, progress
        )
        if indexed is not None:
            return indexed
    input_threads = max(1, threads // 2)
    arguments = (reference, input_threads, fields, ignore_tags)
    if threads == 1:
        return (
            _scan(left, *arguments, left_details, progress, "left"),
            _scan(right, *arguments, right_details, progress, "right"),
        )
    with ProcessPoolExecutor(max_workers=2) as executor:
        left_future = executor.submit(
            _scan, left, *arguments, left_details, progress, "left"
        )
        right_future = executor.submit(
            _scan, right, *arguments, right_details, progress, "right"
        )
        return left_future.result(), right_future.result()


def compare_alignments(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    profile: str,
    ignore_tags: Sequence[str],
    explain: bool,
    max_examples: int,
    progress: bool,
    temp_dir: Optional[Path],
) -> ComparisonResult:
    with TemporaryDirectory(prefix="gendiff-", dir=temp_dir) as work:
        left_details = Path(work) / "left.sqlite" if explain else None
        right_details = Path(work) / "right.sqlite" if explain else None
        left_scan, right_scan = _scan_pair(
            left,
            right,
            reference,
            threads,
            fields,
            ignore_tags,
            left_details,
            right_details,
            progress,
        )
        changed = [
            name
            for name in fields
            if left_scan.fingerprints[name] != right_scan.fingerprints[name]
        ]
        content_equal = (
            left_scan.fingerprints["record"] == right_scan.fingerprints["record"]
        )
        overlap = sketch_containment(left_scan.sketch, right_scan.sketch)
        details = (
            analyze_details(left_details, right_details, fields, max_examples)
            if explain and left_details and right_details
            else None
        )
    return ComparisonResult(
        kind="alignment",
        left=left,
        right=right,
        left_records=left_scan.count,
        right_records=right_scan.count,
        content_equal=content_equal,
        structure_equal=digest(left_scan.structural) == digest(right_scan.structural),
        metadata_equal=digest(left_scan.metadata) == digest(right_scan.metadata),
        changed_fields=changed,
        relationship=relationship_label(
            overlap, content_equal, left_scan.count == right_scan.count == 0
        ),
        identity_overlap=overlap,
        profile=profile,
        details=details,
    )
