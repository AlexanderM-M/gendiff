from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Sequence, Set, Tuple

import pysam

from gendiff.artifacts import output_workspace, sample_slugs, write_manifest
from gendiff.details import (
    DetailWriter,
    SelectionReader,
    analyze_details,
    merge_detail_databases,
    relationship_label,
)
from gendiff.difference_table import publish_table
from gendiff.fingerprint import (
    Fingerprint,
    FingerprintBuilder,
    Sketch,
    SketchBuilder,
    digest,
    digest_native,
    digest_parts,
    merge_sketches,
    normalize,
    sketch_containment,
)
from gendiff.model import ComparisonResult
from gendiff.progress import ProgressReporter
from gendiff.regions import RegionFilter
from gendiff.tracks import contig_lengths, write_tracks


@dataclass(frozen=True)
class _ScanResult:
    fingerprints: Dict[str, Fingerprint]
    count: int
    structural: Dict[str, Any]
    metadata: Dict[str, Any]
    sketch: Sketch
    contig_counts: Dict[str, int]


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
        "_tag_digests": {
            tag: digest_native((value, value_type))
            for tag, value, value_type in values["tags"]
        },
    }


def _record_digests(
    values: Dict[str, Any], fields: Sequence[str]
) -> Tuple[list[int], int]:
    field_digests = [digest_native(values[name]) for name in fields]
    return field_digests, digest_parts(field_digests)


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
    region_filter: Optional[RegionFilter] = None,
) -> _ScanResult:
    builders = {name: FingerprintBuilder() for name in ("record",) + tuple(fields)}
    sketch = SketchBuilder()
    writer = DetailWriter(details_path) if details_path else None
    reporter = ProgressReporter(label, progress)
    ignored = set(ignore_tags)
    kwargs = {"reference_filename": str(reference)} if reference else {}
    kwargs["threads"] = threads
    count = 0
    contig_counts = Counter()
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
                    start = record.reference_start
                    end = record.reference_end or (start + 1 if start >= 0 else start)
                    if region_filter and not region_filter.allows(
                        record.reference_name, start, end
                    ):
                        continue
                    identity = _identity(record)
                    values = _record_values(record, ignored)
                    field_digests, record_digest = _record_digests(values, fields)
                    for name, field_digest in zip(fields, field_digests):
                        builders[name].add_digest(field_digest)
                    builders["record"].add_digest(record_digest)
                    sketch.add_digest(digest_native(identity))
                    if writer:
                        writer.add(
                            identity,
                            record_digest,
                            field_digests,
                            _record_summary(record, values),
                        )
                    count += 1
                    if record.reference_name is not None:
                        contig_counts[record.reference_name] += 1
                    reporter.update(count)
    finally:
        if writer:
            writer.close()
        reporter.finish(count)
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return _ScanResult(
        fingerprints, count, structural, metadata, sketch.finish(), dict(contig_counts)
    )


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
        dict(sum((Counter(scan.contig_counts) for scan in scans), Counter())),
    )


def _scan_indexed_pair(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    ignore_tags: Sequence[str],
    progress: bool,
    left_label: str,
    right_label: str,
    left_details: Optional[Path],
    right_details: Optional[Path],
    region_filter: Optional[RegionFilter],
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
                left_details.parent / f"left-{index}.sqlite" if left_details else None,
                progress,
                f"{left_label}:{index + 1}",
                regions,
                region_filter,
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
                right_details.parent / f"right-{index}.sqlite"
                if right_details
                else None,
                progress,
                f"{right_label}:{index + 1}",
                regions,
                region_filter,
            )
            for index, regions in enumerate(right_chunks)
        ]
        left_scans = [future.result() for future in left_futures]
        right_scans = [future.result() for future in right_futures]
    if left_details:
        merge_detail_databases(
            [
                left_details.parent / f"left-{index}.sqlite"
                for index in range(len(left_chunks))
            ],
            left_details,
        )
    if right_details:
        merge_detail_databases(
            [
                right_details.parent / f"right-{index}.sqlite"
                for index in range(len(right_chunks))
            ],
            right_details,
        )
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
    left_label: str,
    right_label: str,
    region_filter: Optional[RegionFilter],
) -> tuple[_ScanResult, _ScanResult]:
    if threads > 2:
        indexed = _scan_indexed_pair(
            left,
            right,
            reference,
            threads,
            fields,
            ignore_tags,
            progress,
            left_label,
            right_label,
            left_details,
            right_details,
            region_filter,
        )
        if indexed is not None:
            return indexed
    input_threads = max(1, threads // 2)
    arguments = (reference, input_threads, fields, ignore_tags)
    if threads == 1:
        return (
            _scan(
                left,
                *arguments,
                left_details,
                progress,
                left_label,
                None,
                region_filter,
            ),
            _scan(
                right,
                *arguments,
                right_details,
                progress,
                right_label,
                None,
                region_filter,
            ),
        )
    with ProcessPoolExecutor(max_workers=2) as executor:
        left_future = executor.submit(
            _scan,
            left,
            *arguments,
            left_details,
            progress,
            left_label,
            None,
            region_filter,
        )
        right_future = executor.submit(
            _scan,
            right,
            *arguments,
            right_details,
            progress,
            right_label,
            None,
            region_filter,
        )
        return left_future.result(), right_future.result()


def _write_diff_side(
    source: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    ignore_tags: Sequence[str],
    selection_path: Path,
    side: int,
    only_path: Path,
    modified_path: Path,
) -> Counter:
    kwargs = {"reference_filename": str(reference)} if reference else {}
    kwargs["threads"] = max(1, threads)
    reader = SelectionReader(selection_path, side)
    counts = Counter()
    try:
        with pysam.AlignmentFile(str(source), _mode(source), **kwargs) as handle:
            with (
                pysam.AlignmentFile(
                    str(only_path), "wb", header=handle.header, threads=max(1, threads)
                ) as only_output,
                pysam.AlignmentFile(
                    str(modified_path),
                    "wb",
                    header=handle.header,
                    threads=max(1, threads),
                ) as modified_output,
            ):
                ignored = set(ignore_tags)
                for record in handle.fetch(until_eof=True):
                    values = _record_values(record, ignored)
                    _, record_digest = _record_digests(values, fields)
                    status = reader.take(_identity(record), record_digest)
                    if status == "only":
                        only_output.write(record)
                        counts["only"] += 1
                    elif status == "modified":
                        modified_output.write(record)
                        counts["modified"] += 1
        remaining = reader.remaining()
        if remaining:
            raise ValueError(
                f"could not recover {remaining:,} selected records from {source}"
            )
    finally:
        reader.close()
    return counts


def _write_diff_outputs(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    profile: str,
    ignore_tags: Sequence[str],
    selection_path: Path,
    target: Path,
    force: bool,
    left_label: str,
    right_label: str,
) -> Dict[str, str]:
    left_slug, right_slug = sample_slugs(left_label, right_label)
    names = {
        "only_in_first": f"only-in-{left_slug}.bam",
        "only_in_second": f"only-in-{right_slug}.bam",
        "modified_first": f"modified-{left_slug}.bam",
        "modified_second": f"modified-{right_slug}.bam",
    }
    with output_workspace(target, force) as workspace:
        per_input_threads = max(1, threads // 2)
        left_counts = _write_diff_side(
            left,
            reference,
            per_input_threads,
            fields,
            ignore_tags,
            selection_path,
            0,
            workspace / names["only_in_first"],
            workspace / names["modified_first"],
        )
        right_counts = _write_diff_side(
            right,
            reference,
            per_input_threads,
            fields,
            ignore_tags,
            selection_path,
            1,
            workspace / names["only_in_second"],
            workspace / names["modified_second"],
        )
        write_manifest(
            workspace / "manifest.json",
            {
                "gendiff_diff_format": 1,
                "format": "BAM",
                "profile": profile,
                "compared_fields": fields,
                "inputs": [
                    {
                        "label": left_label,
                        "only": left_counts["only"],
                        "modified": left_counts["modified"],
                    },
                    {
                        "label": right_label,
                        "only": right_counts["only"],
                        "modified": right_counts["modified"],
                    },
                ],
                "files": names,
            },
        )
    published = target.absolute()
    return {key: str(published / name) for key, name in names.items()} | {
        "manifest": str(published / "manifest.json")
    }


def compare_alignments(
    left: Path,
    right: Path,
    reference: Optional[Path],
    threads: int,
    fields: Sequence[str],
    profile: str,
    left_label: str,
    right_label: str,
    ignore_tags: Sequence[str],
    explain: bool,
    max_examples: int,
    progress: bool,
    temp_dir: Optional[Path],
    diff_dir: Optional[Path],
    track_dir: Optional[Path],
    difference_table: Optional[Path],
    region_filter: Optional[RegionFilter],
    force: bool,
) -> ComparisonResult:
    artifacts: Dict[str, str] = {}
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
            left_label,
            right_label,
            region_filter,
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
        selection_path = Path(work) / "selections.sqlite" if diff_dir else None
        track_path = Path(work) / "tracks.sqlite" if track_dir else None
        table_path = (
            Path(work)
            / (
                "differences.tsv.gz"
                if difference_table and difference_table.suffix == ".gz"
                else "differences.tsv"
            )
            if difference_table
            else None
        )
        details = (
            analyze_details(
                left_details,
                right_details,
                fields,
                max_examples,
                selection_path,
                track_path,
                table_path,
                left_scan.contig_counts,
                right_scan.contig_counts,
                contig_lengths(left_scan.structural, right_scan.structural),
            )
            if explain and left_details and right_details
            else None
        )
        if diff_dir and selection_path:
            artifacts.update(
                _write_diff_outputs(
                    left,
                    right,
                    reference,
                    threads,
                    fields,
                    profile,
                    ignore_tags,
                    selection_path,
                    diff_dir,
                    force,
                    left_label,
                    right_label,
                )
            )
        if track_dir and track_path:
            artifacts.update(
                write_tracks(
                    track_path,
                    track_dir,
                    force,
                    contig_lengths(left_scan.structural, right_scan.structural),
                )
            )
        if difference_table and table_path:
            artifacts["difference_table"] = publish_table(
                table_path, difference_table, force
            )
    return ComparisonResult(
        kind="alignment",
        left=left,
        right=right,
        left_label=left_label,
        right_label=right_label,
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
        artifacts=artifacts,
    )
