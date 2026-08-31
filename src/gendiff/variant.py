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


def _structural_header(header: pysam.VariantHeader) -> Dict[str, Any]:
    contigs = {
        name: {"length": header.contigs[name].length} for name in sorted(header.contigs)
    }
    return {"contigs": contigs, "samples": sorted(header.samples)}


def _metadata_header(header: pysam.VariantHeader) -> Dict[str, Any]:
    records = sorted(str(record).strip() for record in header.records)
    return {"records": records, "samples": sorted(header.samples)}


def _sample_values(record: pysam.VariantRecord) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for name in sorted(record.samples):
        sample = record.samples[name]
        result[name] = {
            "values": {key: sample[key] for key in sorted(sample.keys())},
            "phased": sample.phased,
        }
    return result


def _record_values(
    record: pysam.VariantRecord, ignore_info: Set[str]
) -> Dict[str, Any]:
    return {
        "locus": (record.contig, record.pos, record.stop),
        "identifiers": record.id,
        "alleles": record.alleles,
        "quality": record.qual,
        "filters": sorted(record.filter.keys()),
        "info": {
            key: record.info[key]
            for key in sorted(record.info.keys())
            if key not in ignore_info
        },
        "samples": _sample_values(record),
    }


def _identity(record: pysam.VariantRecord) -> Any:
    return record.contig, record.pos, record.alleles


def _record_summary(
    record: pysam.VariantRecord, values: Dict[str, Any]
) -> Dict[str, Any]:
    return {
        "contig": record.contig,
        "position": record.pos,
        "identifiers": record.id,
        "alleles": record.alleles,
        "quality": record.qual,
        "filters": sorted(record.filter.keys()),
        "_info": values["info"],
        "genotypes": {
            sample: {
                "GT": record.samples[sample].get("GT"),
                "phased": record.samples[sample].phased,
            }
            for sample in sorted(record.samples)
        },
    }


def _record_digests(
    values: Dict[str, Any], fields: Sequence[str]
) -> Tuple[list[int], int]:
    field_digests = [digest_native(values[name]) for name in fields]
    return field_digests, digest_parts(field_digests)


def _scan(
    path: Path,
    threads: int,
    fields: Sequence[str],
    ignore_info: Sequence[str],
    details_path: Optional[Path],
    progress: bool,
    label: str,
    region_filter: Optional[RegionFilter] = None,
) -> _ScanResult:
    builders = {name: FingerprintBuilder() for name in ("record",) + tuple(fields)}
    sketch = SketchBuilder()
    writer = DetailWriter(details_path) if details_path else None
    reporter = ProgressReporter(label, progress)
    ignored = set(ignore_info)
    count = 0
    contig_counts = Counter()
    try:
        with pysam.VariantFile(str(path), threads=threads) as handle:
            structural = _structural_header(handle.header)
            metadata = _metadata_header(handle.header)
            for record in handle:
                if region_filter and not region_filter.allows(
                    record.contig, record.start, record.stop
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
                contig_counts[record.contig] += 1
                reporter.update(count)
    finally:
        if writer:
            writer.close()
        reporter.finish(count)
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return _ScanResult(
        fingerprints, count, structural, metadata, sketch.finish(), dict(contig_counts)
    )


def _scan_pair(
    left: Path,
    right: Path,
    threads: int,
    fields: Sequence[str],
    ignore_info: Sequence[str],
    left_details: Optional[Path],
    right_details: Optional[Path],
    progress: bool,
    left_label: str,
    right_label: str,
    region_filter: Optional[RegionFilter],
) -> tuple[_ScanResult, _ScanResult]:
    input_threads = max(1, threads // 2)
    arguments = (input_threads, fields, ignore_info)
    if threads == 1:
        return (
            _scan(left, *arguments, left_details, progress, left_label, region_filter),
            _scan(
                right, *arguments, right_details, progress, right_label, region_filter
            ),
        )
    with ProcessPoolExecutor(max_workers=2) as executor:
        left_future = executor.submit(
            _scan, left, *arguments, left_details, progress, left_label, region_filter
        )
        right_future = executor.submit(
            _scan,
            right,
            *arguments,
            right_details,
            progress,
            right_label,
            region_filter,
        )
        return left_future.result(), right_future.result()


def _write_diff_side(
    source: Path,
    threads: int,
    fields: Sequence[str],
    ignore_info: Sequence[str],
    selection_path: Path,
    side: int,
    only_path: Path,
    modified_path: Path,
) -> Counter:
    reader = SelectionReader(selection_path, side)
    counts = Counter()
    try:
        with pysam.VariantFile(str(source), threads=max(1, threads)) as handle:
            with (
                pysam.VariantFile(
                    str(only_path), "w", header=handle.header
                ) as only_output,
                pysam.VariantFile(
                    str(modified_path), "w", header=handle.header
                ) as modified_output,
            ):
                ignored = set(ignore_info)
                for record in handle:
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
    threads: int,
    fields: Sequence[str],
    profile: str,
    ignore_info: Sequence[str],
    selection_path: Path,
    target: Path,
    force: bool,
    left_label: str,
    right_label: str,
    normalized: bool,
) -> Dict[str, str]:
    left_slug, right_slug = sample_slugs(left_label, right_label)
    names = {
        "only_in_first": f"only-in-{left_slug}.vcf",
        "only_in_second": f"only-in-{right_slug}.vcf",
        "modified_first": f"modified-{left_slug}.vcf",
        "modified_second": f"modified-{right_slug}.vcf",
    }
    with output_workspace(target, force) as workspace:
        per_input_threads = max(1, threads // 2)
        left_counts = _write_diff_side(
            left,
            per_input_threads,
            fields,
            ignore_info,
            selection_path,
            0,
            workspace / names["only_in_first"],
            workspace / names["modified_first"],
        )
        right_counts = _write_diff_side(
            right,
            per_input_threads,
            fields,
            ignore_info,
            selection_path,
            1,
            workspace / names["only_in_second"],
            workspace / names["modified_second"],
        )
        write_manifest(
            workspace / "manifest.json",
            {
                "gendiff_diff_format": 1,
                "format": "VCF",
                "normalized": normalized,
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


def _normalize(path: Path, output: Path, reference: Path) -> None:
    from pysam import bcftools

    bcftools.norm(
        "-f",
        str(reference),
        "-m",
        "-any",
        "-Ob",
        "-o",
        str(output),
        str(path),
        catch_stdout=False,
    )


def compare_variants(
    left: Path,
    right: Path,
    reference: Optional[Path],
    normalize_variants: bool,
    threads: int,
    fields: Sequence[str],
    profile: str,
    left_label: str,
    right_label: str,
    ignore_info: Sequence[str],
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
        workspace = Path(work)
        scan_left, scan_right = left, right
        if normalize_variants:
            if reference is None:
                raise ValueError("--normalize requires --reference")
            scan_left = workspace / "left.normalized.bcf"
            scan_right = workspace / "right.normalized.bcf"
            _normalize(left, scan_left, reference)
            _normalize(right, scan_right, reference)
        left_details = workspace / "left.sqlite" if explain else None
        right_details = workspace / "right.sqlite" if explain else None
        left_scan, right_scan = _scan_pair(
            scan_left,
            scan_right,
            threads,
            fields,
            ignore_info,
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
        selection_path = workspace / "selections.sqlite" if diff_dir else None
        track_path = workspace / "tracks.sqlite" if track_dir else None
        table_path = (
            workspace
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
                    scan_left,
                    scan_right,
                    threads,
                    fields,
                    profile,
                    ignore_info,
                    selection_path,
                    diff_dir,
                    force,
                    left_label,
                    right_label,
                    normalize_variants,
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
        kind="variant",
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
