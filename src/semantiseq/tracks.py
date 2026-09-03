"""Write standard genomic tracks from disk-backed changed-locus counts."""

from __future__ import annotations

import math
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

from semantiseq.artifacts import output_workspace, write_manifest

_WINDOW = 1_000_000


def _rows(path: Path) -> Iterator[Tuple[str, int, str, int]]:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        yield from connection.execute(
            "SELECT contig, start, status, count FROM loci "
            "ORDER BY contig, start, status"
        )
    finally:
        connection.close()


def contig_lengths(*headers: Mapping[str, Any]) -> Dict[str, int]:
    """Extract the largest known length for each contig from either header family."""
    lengths: Dict[str, int] = {}
    for header in headers:
        for item in header.get("SQ", []):
            name, length = item.get("SN"), item.get("LN")
            if name is not None and isinstance(length, int):
                lengths[str(name)] = max(lengths.get(str(name), 0), length)
        for name, item in header.get("contigs", {}).items():
            length = item.get("length")
            if isinstance(length, int):
                lengths[str(name)] = max(lengths.get(str(name), 0), length)
    return lengths


def write_tracks(
    database: Path,
    target: Path,
    force: bool,
    lengths: Optional[Mapping[str, int]] = None,
    first_windows: Optional[Mapping[Tuple[str, int], int]] = None,
    second_windows: Optional[Mapping[Tuple[str, int], int]] = None,
) -> Dict[str, str]:
    names = {
        "bed": "changes.bed",
        "bedgraph": "difference-density.bedgraph",
        "rate": "difference-rate.bedgraph",
        "coverage": "coverage-ratio.bedgraph",
        "contigs": "contig-summary.tsv",
    }
    with output_workspace(target, force) as workspace:
        contigs = Counter()
        windows = Counter()
        with (workspace / names["bed"]).open("w", encoding="utf-8") as bed:
            for contig, start, status, count in _rows(database):
                bed.write(f"{contig}\t{start}\t{start + 1}\t{status};x{count}\n")
                contigs[contig] += count
                windows[(contig, start // _WINDOW)] += count
        with (workspace / names["bedgraph"]).open("w", encoding="utf-8") as graph:
            for (contig, window), count in sorted(windows.items()):
                start = window * _WINDOW
                end = start + _WINDOW
                if lengths and contig in lengths:
                    end = min(end, lengths[contig])
                graph.write(f"{contig}\t{start}\t{end}\t{count}\n")
        first = first_windows or {}
        second = second_windows or {}
        window_keys = set(windows) | set(first) | set(second)
        with (
            (workspace / names["rate"]).open("w", encoding="utf-8") as rates,
            (workspace / names["coverage"]).open("w", encoding="utf-8") as coverage,
        ):
            for contig, window in sorted(window_keys):
                start = window * _WINDOW
                end = min(start + _WINDOW, (lengths or {}).get(contig, start + _WINDOW))
                before = first.get((contig, window), 0)
                after = second.get((contig, window), 0)
                rate = min(
                    1.0, windows.get((contig, window), 0) / max(1, before + after)
                )
                ratio = math.log2((after + 1) / (before + 1))
                rates.write(f"{contig}\t{start}\t{end}\t{rate:.6g}\n")
                coverage.write(f"{contig}\t{start}\t{end}\t{ratio:.6g}\n")
        with (workspace / names["contigs"]).open("w", encoding="utf-8") as summary:
            summary.write("contig\tchanged_records\n")
            for contig, count in contigs.most_common():
                summary.write(f"{contig}\t{count}\n")
        write_manifest(
            workspace / "manifest.json",
            {
                "semantiseq_diff_format": 1,
                "kind": "genomic-tracks",
                "coordinate_systems": {
                    "bed": "zero-based, half-open",
                    "bedgraph": "zero-based, half-open; one-megabase windows",
                    "rate": "changed records / compared records by window",
                    "coverage": "log2(second / first) record count by window",
                },
                "files": names,
            },
        )
    published = target.absolute()
    return {f"track_{key}": str(published / name) for key, name in names.items()} | {
        "track_manifest": str(published / "manifest.json")
    }
