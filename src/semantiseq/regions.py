from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Tuple


@dataclass(frozen=True)
class _Intervals:
    starts: Tuple[int, ...]
    ends: Tuple[int, ...]

    def overlaps(self, start: int, end: int) -> bool:
        index = bisect_left(self.starts, end) - 1
        return index >= 0 and self.ends[index] > start


@dataclass(frozen=True)
class RegionFilter:
    include: Dict[str, _Intervals]
    exclude: Dict[str, _Intervals]

    def allows(self, contig: Optional[str], start: int, end: int) -> bool:
        if contig is None or start < 0:
            return not self.include
        included = self.include.get(contig)
        if self.include and (included is None or not included.overlaps(start, end)):
            return False
        excluded = self.exclude.get(contig)
        return excluded is None or not excluded.overlaps(start, end)


def _merge(values: Iterable[Tuple[str, int, int]]) -> Dict[str, _Intervals]:
    grouped: Dict[str, list[Tuple[int, int]]] = {}
    for contig, start, end in values:
        grouped.setdefault(contig, []).append((start, end))
    result = {}
    for contig, intervals in grouped.items():
        merged: list[list[int]] = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        result[contig] = _Intervals(
            tuple(item[0] for item in merged), tuple(item[1] for item in merged)
        )
    return result


def _parse_region(value: str) -> Tuple[str, int, int]:
    try:
        contig, coordinates = value.rsplit(":", 1)
        start_text, end_text = coordinates.replace(",", "").split("-", 1)
        start, end = int(start_text), int(end_text)
    except (TypeError, ValueError) as error:
        raise ValueError(
            f"invalid region: {value}; expected CONTIG:START-END"
        ) from error
    if not contig or start < 1 or end < start:
        raise ValueError(f"invalid region: {value}; expected CONTIG:START-END")
    return contig, start - 1, end


def _read_bed(path: Path) -> list[Tuple[str, int, int]]:
    intervals = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            columns = stripped.split()
            try:
                contig, start, end = columns[0], int(columns[1]), int(columns[2])
            except (IndexError, ValueError) as error:
                raise ValueError(f"invalid BED row {line_number} in {path}") from error
            if start < 0 or end <= start:
                raise ValueError(f"invalid BED interval on row {line_number} in {path}")
            intervals.append((contig, start, end))
    return intervals


def load_region_filter(
    regions: Sequence[str] = (),
    include_bed: Optional[Path] = None,
    exclude_bed: Optional[Path] = None,
) -> Optional[RegionFilter]:
    include = [_parse_region(value) for value in regions]
    if include_bed:
        include.extend(_read_bed(include_bed))
    exclude = _read_bed(exclude_bed) if exclude_bed else []
    if not include and not exclude:
        return None
    return RegionFilter(_merge(include), _merge(exclude))
