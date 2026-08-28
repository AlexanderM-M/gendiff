from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

import pysam

from gendiff.fingerprint import Fingerprint, FingerprintBuilder, digest, normalize
from gendiff.model import ComparisonResult


_FIELDS = ("locus", "identifiers", "alleles", "quality", "filters", "info", "samples")


def _structural_header(header: pysam.VariantHeader) -> Dict[str, Any]:
    contigs = {
        name: {"length": header.contigs[name].length}
        for name in sorted(header.contigs)
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


def _record_values(record: pysam.VariantRecord) -> Dict[str, Any]:
    return {
        "locus": (record.contig, record.pos, record.stop),
        "identifiers": record.id,
        "alleles": record.alleles,
        "quality": record.qual,
        "filters": sorted(record.filter.keys()),
        "info": {key: record.info[key] for key in sorted(record.info.keys())},
        "samples": _sample_values(record),
    }


def _scan(
    path: Path,
) -> Tuple[Dict[str, Fingerprint], int, Dict[str, Any], Dict[str, Any]]:
    builders = {name: FingerprintBuilder() for name in ("record",) + _FIELDS}
    with pysam.VariantFile(str(path)) as handle:
        structural = _structural_header(handle.header)
        metadata = _metadata_header(handle.header)
        for record in handle:
            values = normalize(_record_values(record))
            builders["record"].add(values)
            for name in _FIELDS:
                builders[name].add(values[name])
    fingerprints = {name: builder.finish() for name, builder in builders.items()}
    return fingerprints, fingerprints["record"].count, structural, metadata


def compare_variants(left: Path, right: Path) -> ComparisonResult:
    left_fp, left_count, left_structure, left_metadata = _scan(left)
    right_fp, right_count, right_structure, right_metadata = _scan(right)
    changed = [name for name in _FIELDS if left_fp[name] != right_fp[name]]
    return ComparisonResult(
        kind="variant",
        left=left,
        right=right,
        left_records=left_count,
        right_records=right_count,
        content_equal=left_fp["record"] == right_fp["record"],
        structure_equal=digest(left_structure) == digest(right_structure),
        metadata_equal=digest(left_metadata) == digest(right_metadata),
        changed_fields=changed,
    )
