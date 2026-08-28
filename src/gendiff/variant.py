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


def _record_summary(record: pysam.VariantRecord) -> Dict[str, Any]:
    return {
        "contig": record.contig,
        "position": record.pos,
        "identifiers": record.id,
        "alleles": record.alleles,
        "quality": record.qual,
        "filters": sorted(record.filter.keys()),
        "genotypes": {
            sample: {
                "GT": record.samples[sample].get("GT"),
                "phased": record.samples[sample].phased,
            }
            for sample in sorted(record.samples)
        },
    }


def _scan(
    path: Path,
    threads: int,
    fields: Sequence[str],
    ignore_info: Sequence[str],
    details_path: Optional[Path],
    progress: bool,
    label: str,
) -> _ScanResult:
    builders = {name: FingerprintBuilder() for name in ("record",) + tuple(fields)}
    sketch = SketchBuilder()
    writer = DetailWriter(details_path) if details_path else None
    reporter = ProgressReporter(label, progress)
    ignored = set(ignore_info)
    count = 0
    try:
        with pysam.VariantFile(str(path), threads=threads) as handle:
            structural = _structural_header(handle.header)
            metadata = _metadata_header(handle.header)
            for record in handle:
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
                        _record_summary(record),
                    )
                count += 1
                reporter.update(count)
    finally:
        if writer:
            writer.close()
        reporter.finish(count)
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return _ScanResult(fingerprints, count, structural, metadata, sketch.finish())


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
) -> tuple[_ScanResult, _ScanResult]:
    input_threads = max(1, threads // 2)
    arguments = (input_threads, fields, ignore_info)
    if threads == 1:
        return (
            _scan(left, *arguments, left_details, progress, left_label),
            _scan(right, *arguments, right_details, progress, right_label),
        )
    with ProcessPoolExecutor(max_workers=2) as executor:
        left_future = executor.submit(
            _scan, left, *arguments, left_details, progress, left_label
        )
        right_future = executor.submit(
            _scan, right, *arguments, right_details, progress, right_label
        )
        return left_future.result(), right_future.result()


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
) -> ComparisonResult:
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
    )
