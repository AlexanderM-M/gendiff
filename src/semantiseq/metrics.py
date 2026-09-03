from __future__ import annotations

import math
from collections import Counter
from typing import Any, Dict, Iterable, Mapping, MutableMapping, Sequence

MetricCounts = Dict[str, Dict[float, int]]

_META = {
    "mapping_quality": ("Mapping quality", "distribution"),
    "read_length": ("Read length", "distribution"),
    "insert_size": ("Insert size", "distribution"),
    "unmapped_rate": ("Unmapped rate", "rate"),
    "duplicate_rate": ("Duplicate rate", "rate"),
    "secondary_rate": ("Secondary alignment rate", "rate"),
    "supplementary_rate": ("Supplementary alignment rate", "rate"),
    "variant_quality": ("Variant quality", "distribution"),
    "depth": ("Depth", "distribution"),
    "allele_frequency": ("Allele frequency", "distribution"),
    "filtered_rate": ("Filtered variant rate", "rate"),
    "heterozygous_rate": ("Heterozygous genotype rate", "rate"),
}


def _add(metrics: MutableMapping[str, Counter], name: str, value: float) -> None:
    metrics.setdefault(name, Counter())[value] += 1


def _bucket(value: float, step: float, maximum: float) -> float:
    return min(maximum, math.floor(value / step) * step)


def add_alignment_metrics(
    metrics: MutableMapping[str, Counter], record: Any, fields: Sequence[str]
) -> None:
    selected = set(fields)
    if "mapping_quality" in selected:
        _add(metrics, "mapping_quality", _bucket(record.mapping_quality, 5, 255))
    if "sequence" in selected and record.query_length is not None:
        _add(metrics, "read_length", _bucket(record.query_length, 250, 100_000))
    if "template_length" in selected and record.template_length:
        _add(metrics, "insert_size", _bucket(abs(record.template_length), 100, 10_000))
    if "flags" in selected:
        _add(metrics, "unmapped_rate", float(record.is_unmapped))
        _add(metrics, "duplicate_rate", float(record.is_duplicate))
        _add(metrics, "secondary_rate", float(record.is_secondary))
        _add(metrics, "supplementary_rate", float(record.is_supplementary))


def _number(value: Any) -> Iterable[float]:
    values = value if isinstance(value, tuple) else (value,)
    for item in values:
        if isinstance(item, (int, float)) and math.isfinite(float(item)):
            yield float(item)


def add_variant_metrics(
    metrics: MutableMapping[str, Counter],
    record: Any,
    fields: Sequence[str],
    ignored_info: Sequence[str],
) -> None:
    selected = set(fields)
    ignored = set(ignored_info)
    if "quality" in selected and record.qual is not None:
        _add(metrics, "variant_quality", _bucket(record.qual, 10, 1_000))
    if "info" in selected:
        if "DP" not in ignored and "DP" in record.info:
            for value in _number(record.info["DP"]):
                _add(metrics, "depth", _bucket(value, 10, 10_000))
        if "AF" not in ignored and "AF" in record.info:
            for value in _number(record.info["AF"]):
                _add(metrics, "allele_frequency", round(value * 20) / 20)
    if "filters" in selected:
        filters = set(record.filter.keys())
        _add(metrics, "filtered_rate", float(bool(filters - {"PASS"})))
    if "samples" in selected:
        for sample in record.samples.values():
            genotype = tuple(
                allele for allele in sample.get("GT", ()) if allele is not None
            )
            if genotype:
                _add(metrics, "heterozygous_rate", float(len(set(genotype)) > 1))


def merge_metric_counts(values: Iterable[MetricCounts]) -> MetricCounts:
    merged: Dict[str, Counter] = {}
    for metrics in values:
        for name, counts in metrics.items():
            merged.setdefault(name, Counter()).update(counts)
    return {name: dict(counts) for name, counts in merged.items()}


def _median(counts: Mapping[float, int]) -> float:
    total = sum(counts.values())
    targets = ((total - 1) // 2, total // 2)
    answers = []
    seen = 0
    for value, count in sorted(counts.items()):
        while len(answers) < 2 and targets[len(answers)] < seen + count:
            answers.append(value)
        seen += count
        if len(answers) == 2:
            break
    return sum(answers) / len(answers) if answers else 0.0


def _bins(
    first: Mapping[float, int], second: Mapping[float, int], limit: int = 32
) -> list[Dict[str, Any]]:
    values = sorted(set(first) | set(second))
    size = max(1, math.ceil(len(values) / limit))
    first_total = sum(first.values()) or 1
    second_total = sum(second.values()) or 1
    result = []
    for offset in range(0, len(values), size):
        group = values[offset : offset + size]
        label = f"{group[0]:g}" if len(group) == 1 else f"{group[0]:g}–{group[-1]:g}"
        result.append(
            {
                "label": label,
                "first": sum(first.get(value, 0) for value in group) / first_total,
                "second": sum(second.get(value, 0) for value in group) / second_total,
            }
        )
    return result


def distribution_shifts(
    first: MetricCounts, second: MetricCounts
) -> list[Dict[str, Any]]:
    shifts = []
    for name in sorted(set(first) | set(second)):
        before, after = first.get(name, {}), second.get(name, {})
        first_total, second_total = sum(before.values()), sum(after.values())
        if not first_total or not second_total:
            continue
        values = set(before) | set(after)
        distance = (
            sum(
                abs(
                    before.get(value, 0) / first_total
                    - after.get(value, 0) / second_total
                )
                for value in values
            )
            / 2
        )
        if distance < 0.01:
            continue
        label, kind = _META[name]
        item: Dict[str, Any] = {
            "metric": name,
            "label": label,
            "kind": kind,
            "distance": distance,
        }
        if kind == "rate":
            item["first_value"] = before.get(1.0, 0) / first_total
            item["second_value"] = after.get(1.0, 0) / second_total
        else:
            item["first_median"] = _median(before)
            item["second_median"] = _median(after)
            item["bins"] = _bins(before, after)
        shifts.append(item)
    return sorted(shifts, key=lambda item: item["distance"], reverse=True)[:6]
