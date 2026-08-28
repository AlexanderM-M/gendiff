from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import pysam

from gendiff.fingerprint import (
    Fingerprint,
    FingerprintBuilder,
    digest,
    digest_parts,
    normalize,
)
from gendiff.model import ComparisonResult


_FIELDS = (
    "read_name",
    "flags",
    "reference",
    "position",
    "mapping_quality",
    "cigar",
    "mate",
    "template_length",
    "sequence",
    "qualities",
    "tags",
)


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


def _record_values(record: pysam.AlignedSegment) -> Dict[str, Any]:
    values: Dict[str, Any] = {
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
        "tags": sorted(record.get_tags(with_value_type=True), key=lambda tag: tag[0]),
    }
    return values


def _scan(
    path: Path, reference: Optional[Path], threads: int
) -> Tuple[Dict[str, Fingerprint], int, Dict[str, Any], Dict[str, Any]]:
    builders = {name: FingerprintBuilder() for name in ("record",) + _FIELDS}
    kwargs = {"reference_filename": str(reference)} if reference else {}
    kwargs["threads"] = threads
    with pysam.AlignmentFile(str(path), _mode(path), **kwargs) as handle:
        structural = _structural_header(handle.header)
        metadata = _ordered_header(handle.header)
        for record in handle.fetch(until_eof=True):
            values = _record_values(record)
            parts = []
            for name in _FIELDS:
                field_digest = digest(values[name])
                builders[name].add_digest(field_digest)
                parts.append(field_digest)
            builders["record"].add_digest(digest_parts(parts))
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return fingerprints, fingerprints["record"].count, structural, metadata


def compare_alignments(
    left: Path, right: Path, reference: Optional[Path], threads: int
) -> ComparisonResult:
    if threads == 1:
        left_scan = _scan(left, reference, 1)
        right_scan = _scan(right, reference, 1)
    else:
        input_threads = max(1, threads // 2)
        with ProcessPoolExecutor(max_workers=2) as executor:
            left_future = executor.submit(_scan, left, reference, input_threads)
            right_future = executor.submit(_scan, right, reference, input_threads)
            left_scan = left_future.result()
            right_scan = right_future.result()

    left_fp, left_count, left_structure, left_metadata = left_scan
    right_fp, right_count, right_structure, right_metadata = right_scan
    changed = [name for name in _FIELDS if left_fp[name] != right_fp[name]]
    return ComparisonResult(
        kind="alignment",
        left=left,
        right=right,
        left_records=left_count,
        right_records=right_count,
        content_equal=left_fp["record"] == right_fp["record"],
        structure_equal=digest(left_structure) == digest(right_structure),
        metadata_equal=digest(left_metadata) == digest(right_metadata),
        changed_fields=changed,
    )
