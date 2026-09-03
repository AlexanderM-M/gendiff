from __future__ import annotations

import statistics
from dataclasses import dataclass
from html import escape
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Sequence, Tuple

from semantiseq.compare import compare_files


@dataclass(frozen=True)
class MatrixResult:
    files: Tuple[Path, ...]
    labels: Tuple[str, ...]
    values: Tuple[Tuple[float, ...], ...]
    order: Tuple[int, ...]
    outliers: Tuple[int, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "semantiseq_matrix_format": 1,
            "files": [
                {"label": label, "path": str(path)}
                for label, path in zip(self.labels, self.files)
            ],
            "similarity": [list(row) for row in self.values],
            "cluster_order": list(self.order),
            "outliers": [self.labels[index] for index in self.outliers],
        }


def _cluster(values: Sequence[Sequence[float]]) -> Tuple[int, ...]:
    clusters = [[index] for index in range(len(values))]
    while len(clusters) > 1:
        first, second = max(
            (
                (left, right)
                for left in range(len(clusters))
                for right in range(left + 1, len(clusters))
            ),
            key=lambda pair: (
                sum(values[a][b] for a in clusters[pair[0]] for b in clusters[pair[1]])
                / (len(clusters[pair[0]]) * len(clusters[pair[1]]))
            ),
        )
        left, right = clusters[first], clusters[second]
        if values[left[-1]][right[0]] < values[right[-1]][left[0]]:
            left, right = right, left
        merged = left + right
        clusters.pop(second)
        clusters.pop(first)
        clusters.append(merged)
    return tuple(clusters[0])


def compare_matrix(
    files: Sequence[Path],
    *,
    threads: int,
    profile: str,
    reference: Optional[Path],
    normalize_variants: bool,
    ignore_tags: Sequence[str],
    ignore_info: Sequence[str],
    cache_dir: Optional[Path],
) -> MatrixResult:
    if len(files) < 2:
        raise ValueError("matrix comparison needs at least two files")
    resolved = tuple(path.absolute() for path in files)
    if len(set(resolved)) != len(resolved):
        raise ValueError("matrix inputs must be unique")

    def run(active_cache: Path) -> MatrixResult:
        size = len(resolved)
        values = [
            [1.0 if row == column else 0.0 for column in range(size)]
            for row in range(size)
        ]
        for first in range(size):
            for second in range(first + 1, size):
                result = compare_files(
                    resolved[first],
                    resolved[second],
                    reference,
                    threads,
                    profile=profile,
                    normalize_variants=normalize_variants,
                    ignore_tags=ignore_tags,
                    ignore_info=ignore_info,
                    cache_dir=active_cache,
                )
                values[first][second] = values[second][first] = (
                    result.content_similarity
                )
        averages = [
            sum(value for column, value in enumerate(row) if column != index)
            / (size - 1)
            for index, row in enumerate(values)
        ]
        median = statistics.median(averages)
        outliers = tuple(
            index for index, value in enumerate(averages) if value < median - 0.15
        )
        labels = tuple(path.name for path in resolved)
        return MatrixResult(
            resolved,
            labels,
            tuple(tuple(row) for row in values),
            _cluster(values),
            outliers,
        )

    if cache_dir:
        return run(cache_dir)
    with TemporaryDirectory(prefix="semantiseq-matrix-") as temporary:
        return run(Path(temporary))


def render_matrix_text(result: MatrixResult) -> str:
    lines = [f"Files: {len(result.files)}", "Similarity order:"]
    lines.extend(f"  {result.labels[index]}" for index in result.order)
    if result.outliers:
        lines.append(
            "Outliers: " + ", ".join(result.labels[i] for i in result.outliers)
        )
    return "\n".join(lines)


def render_matrix_html(result: MatrixResult) -> str:
    order = result.order
    headers = "".join(
        f"<th><span>{escape(result.labels[index])}</span></th>" for index in order
    )
    rows = []
    for first in order:
        cells = []
        for second in order:
            value = result.values[first][second]
            hue = value * 120
            cells.append(
                f"<td style='background:hsl({hue:.1f} 65% 82%)' "
                f"title='{escape(result.labels[first])} and "
                f"{escape(result.labels[second])}: {value:.1%} semantic similarity'>"
                f"{value:.0%}</td>"
            )
        rows.append(f"<tr><th>{escape(result.labels[first])}</th>{''.join(cells)}</tr>")
    outliers = (
        "<p class='warning'>Outliers: "
        + ", ".join(escape(result.labels[index]) for index in result.outliers)
        + "</p>"
        if result.outliers
        else ""
    )
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SemantiSeq similarity matrix</title><style>
body{{font:14px/1.4 system-ui,sans-serif;color:#1f2328;margin:32px}}
main{{max-width:1200px;margin:auto}}h1{{margin-bottom:4px}}p{{color:#59636e}}
.matrix{{overflow:auto;border:1px solid #d0d7de;border-radius:8px;padding:12px}}
table{{border-collapse:separate;border-spacing:3px}}th{{text-align:right;max-width:180px;
overflow:hidden;text-overflow:ellipsis;white-space:nowrap}}thead th span{{display:block;
writing-mode:vertical-rl;transform:rotate(180deg);height:150px}}td{{min-width:58px;
height:42px;text-align:center;border-radius:4px;font-weight:600}}.warning{{color:#cf222e}}
</style></head><body><main><h1>SemantiSeq similarity matrix</h1>
<p>{len(result.files)} files, grouped by semantic content similarity.</p>{outliers}
<div class="matrix"><table><thead><tr><th></th>{headers}</tr></thead>
<tbody>{"".join(rows)}</tbody></table></div></main></body></html>\n"""
