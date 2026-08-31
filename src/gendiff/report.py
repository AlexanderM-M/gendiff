from __future__ import annotations

import json
import os
import re
from collections import Counter
from html import escape
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Dict, Iterable, List

from gendiff.model import ComparisonResult


def _bar_rows(values: Dict[str, int]) -> str:
    if not values:
        return "<p>No changed field groups.</p>"
    maximum = max(values.values()) or 1
    rows = []
    for name, value in values.items():
        width = max(1.0, value * 100.0 / maximum)
        rows.append(
            "<div class='bar-row'>"
            f"<span>{escape(name.replace('_', ' '))}</span>"
            f"<div class='bar'><i style='width:{width:.2f}%'></i></div>"
            f"<strong>{value:,}</strong></div>"
        )
    return "".join(rows)


def _table(headers: Iterable[str], rows: Iterable[Iterable[str]]) -> str:
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(value)}</td>" for value in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _example_status(status: str, result: ComparisonResult) -> str:
    return {
        "left_only": f"Only in {result.left_label}",
        "right_only": f"Only in {result.right_label}",
        "modified": "Modified",
    }.get(status, status.replace("_", " ").title())


def _composition(result: ComparisonResult) -> str:
    details = result.details
    if details is None:
        return ""
    values = (
        ("Identical", details.identical, "#1a7f37"),
        ("Modified", details.modified, "#bf8700"),
        (f"Only in {result.left_label}", details.left_only, "#0969da"),
        (f"Only in {result.right_label}", details.right_only, "#8250df"),
    )
    total = sum(value for _, value, _ in values) or 1
    segments = "".join(
        f"<i style='width:{value * 100 / total:.3f}%;background:{color}' "
        f"title='{escape(name)}: {value:,}'></i>"
        for name, value, color in values
        if value
    )
    legend = "".join(
        f"<span><b style='background:{color}'></b>{escape(name)} "
        f"<strong>{value:,}</strong></span>"
        for name, value, color in values
    )
    return (
        f"<div class='composition'>{segments}</div><div class='legend'>{legend}</div>"
    )


def _genome_map(result: ComparisonResult) -> str:
    details = result.details
    if details is None or not details.contig_stats:
        return ""
    density = {}
    for item in details.region_density:
        density.setdefault(item["contig"], []).append(item)
    maximum = max((item["changes"] for item in details.region_density), default=1)
    rows = []
    for stat in details.contig_stats[:12]:
        contig = str(stat["contig"])
        length = int(
            stat["length"]
            or max((x["end"] for x in density.get(contig, [])), default=1)
        )
        cells = "".join(
            f"<i style='left:{100 * item['start'] / length:.3f}%;"
            f"width:{max(0.35, 100 * (item['end'] - item['start']) / length):.3f}%;"
            f"opacity:{0.18 + 0.82 * item['changes'] / maximum:.3f}' "
            f"title='{escape(contig)}:{item['start'] + 1:,} — "
            f"{item['changes']:,} changes'></i>"
            for item in density.get(contig, [])
        )
        rows.append(
            f"<div class='genome-row'><strong>{escape(contig)}</strong>"
            f"<div class='genome-track'>{cells}</div>"
            f"<span>{stat['changes']:,}</span></div>"
        )
    return (
        "<section><h2>Differences across the genome</h2>" + "".join(rows) + "</section>"
    )


def _transition_matrix(result: ComparisonResult) -> str:
    details = result.details
    if details is None:
        return ""
    patterns = (
        ("Genotype transitions", re.compile(r"^genotype .*?: (.+) → (.+)$")),
        ("MAPQ transitions", re.compile(r"^MAPQ: (.+) → (.+)$")),
    )
    for title, pattern in patterns:
        pairs = Counter()
        for transition, count in details.transitions.items():
            match = pattern.match(transition)
            if match:
                pairs[match.groups()] += count
        if not pairs:
            continue
        state_counts = Counter()
        for (first, second), count in pairs.items():
            state_counts[first] += count
            state_counts[second] += count
        states = [value for value, _ in state_counts.most_common(8)]
        maximum = max(pairs.values())
        header = "".join(f"<th>{escape(value)}</th>" for value in states)
        rows = []
        for first in states:
            cells = "".join(
                f"<td class='heat' style='--heat:"
                f"{pairs[(first, second)] / maximum:.3f}'>"
                f"{pairs[(first, second)] or ''}</td>"
                for second in states
            )
            rows.append(f"<tr><th>{escape(first)}</th>{cells}</tr>")
        return (
            f"<section><h2>{title}</h2><table class='matrix'><tr>"
            f"<th>From ↓ / to →</th>{header}</tr>{''.join(rows)}</table></section>"
        )
    return ""


def render_html(result: ComparisonResult) -> str:
    details = result.details
    status = "Equivalent" if result.equivalent else "Different"
    status_class = "ok" if result.equivalent else "different"
    detail_sections: List[str] = []
    if details is not None:
        detail_sections.append(
            "<section><h2>Record comparison</h2>" + _composition(result) + "</section>"
        )
        detail_sections.append(_genome_map(result))
        detail_sections.append(_transition_matrix(result))
        detail_sections.append(
            "<section><h2>Changed fields</h2>"
            f"{_bar_rows(details.field_changes)}</section>"
        )
        if details.transitions:
            detail_sections.append(
                "<details><summary>Change transitions</summary><section>"
                + _table(
                    ("Transition", "Records"),
                    (
                        (transition, f"{count:,}")
                        for transition, count in list(details.transitions.items())[:50]
                    ),
                )
                + "</section></details>"
            )
        if details.top_regions:
            detail_sections.append(
                "<details><summary>Most affected regions</summary><section>"
                + _table(
                    ("Region", "Changes"),
                    (
                        (str(item["region"]), f"{item['changes']:,}")
                        for item in details.top_regions
                    ),
                )
                + "</section></details>"
            )
        if details.sample_changes:
            detail_sections.append(
                "<section><h2>Genotype changes by sample</h2>"
                + _table(
                    ("Sample", "Changes"),
                    (
                        (sample, f"{count:,}")
                        for sample, count in details.sample_changes.items()
                    ),
                )
                + "</section>"
            )
        if details.examples:
            detail_sections.append(
                "<details><summary>Record examples</summary><section>"
                + _table(
                    (
                        "Status",
                        "Identity",
                        "Changed fields",
                        result.left_label,
                        result.right_label,
                    ),
                    (
                        (
                            _example_status(str(item["status"]), result),
                            json.dumps(item["identity"], separators=(",", ":")),
                            ", ".join(item["changed_fields"]),
                            json.dumps(
                                item["left"], separators=(",", ":"), sort_keys=True
                            ),
                            json.dumps(
                                item["right"], separators=(",", ":"), sort_keys=True
                            ),
                        )
                        for item in details.examples
                    ),
                )
                + "</section></details>"
            )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GenDiff report</title>
<style>
:root{{--bg:#f6f8fa;--panel:#fff;--text:#1f2328;--muted:#59636e;--line:#d0d7de;
--accent:#0969da;--ok:#1a7f37;--bad:#cf222e}}*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--text);
font:14px/1.5 system-ui,sans-serif}}
main{{max-width:1100px;margin:0 auto;padding:32px 20px 60px}}
h1{{margin:0 0 4px;font-size:30px}}h2{{font-size:18px;margin:0 0 16px}}
section{{background:var(--panel);border:1px solid var(--line);border-radius:8px;
padding:20px;margin-top:16px;overflow:auto}}.status{{font-size:20px;font-weight:700}}
.status.ok{{color:var(--ok)}}.status.different{{color:var(--bad)}}.muted{{color:var(--muted)}}
.summary{{display:grid;grid-template-columns:180px 1fr;gap:7px 16px;margin-top:18px}}
.summary dt{{color:var(--muted)}}.summary dd{{margin:0;overflow-wrap:anywhere}}
.cards{{display:grid;grid-template-columns:repeat(4,minmax(120px,1fr));gap:12px}}
.cards article{{border:1px solid var(--line);border-radius:6px;padding:14px}}
.cards strong{{display:block;font-size:22px}}.cards span{{color:var(--muted)}}
.composition{{display:flex;height:28px;border-radius:7px;overflow:hidden;
background:#eaeef2}}.composition i{{display:block;min-width:2px}}
.legend{{display:flex;flex-wrap:wrap;gap:14px;margin-top:12px}}
.legend span{{color:var(--muted)}}.legend b{{display:inline-block;width:10px;
height:10px;border-radius:2px;margin-right:6px}}
.legend strong{{color:var(--text);margin-left:3px}}
.genome-row{{display:grid;grid-template-columns:110px 1fr 70px;gap:12px;
align-items:center;margin:8px 0}}.genome-row span{{text-align:right}}
.genome-track{{height:14px;background:#eaeef2;border-radius:4px;
position:relative;overflow:hidden}}.genome-track i{{position:absolute;height:100%;
background:var(--bad)}}details{{margin-top:16px}}
summary{{cursor:pointer;font-weight:700;font-size:16px}}
.heat{{text-align:center;background:rgba(207,34,46,var(--heat))}}
.matrix th{{text-align:center}}
.bar-row{{display:grid;grid-template-columns:150px 1fr 80px;gap:12px;
align-items:center;margin:9px 0}}
.bar{{height:10px;background:#eaeef2;border-radius:5px;overflow:hidden}}
.bar i{{display:block;height:100%;background:var(--accent)}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border-bottom:1px solid var(--line);padding:8px;text-align:left;
vertical-align:top;max-width:360px;overflow-wrap:anywhere}}
th{{background:#f6f8fa}}
@media(max-width:700px){{.cards{{grid-template-columns:1fr 1fr}}
.summary{{grid-template-columns:1fr}}}}
</style>
</head>
<body><main>
<h1>GenDiff</h1><div class="status {status_class}">{status}</div>
<section><h2>Summary</h2><dl class="summary">
<dt>Identity overlap</dt><dd>{result.identity_overlap:.1%}</dd>
<dt>{escape(result.left_label)}</dt>
<dd>{result.left_records:,} records</dd>
<dt>{escape(result.right_label)}</dt>
<dd>{result.right_records:,} records</dd>
<dt>Logical records</dt><dd>{"same" if result.content_equal else "different"}</dd>
<dt>Structural header</dt><dd>{"same" if result.structure_equal else "different"}</dd>
<dt>Metadata header</dt>
<dd>{"same" if result.metadata_equal else "different"} (informational)</dd>
</dl></section>
{"".join(detail_sections)}
</main></body></html>
"""


def _short(value: str, limit: int = 52) -> str:
    return value if len(value) <= limit else value[: limit - 1] + "…"


def render_svg(result: ComparisonResult) -> str:
    """Render a compact, standalone comparison figure."""
    details = result.details
    record_values = []
    if details is not None:
        record_values = [
            ("Identical", details.identical, "#1a7f37"),
            ("Modified", details.modified, "#bf8700"),
            (f"Only in {result.left_label}", details.left_only, "#0969da"),
            (f"Only in {result.right_label}", details.right_only, "#8250df"),
        ]
    transition_values = list((details.transitions if details else {}).items())[:8]
    field_values = (
        transition_values
        or list((details.field_changes if details else {}).items())[:8]
    )
    height = 420 if not field_values else 444 + len(field_values) * 36
    status = "Equivalent" if result.equivalent else "Different"
    status_color = "#1a7f37" if result.equivalent else "#cf222e"
    overlap_width = 920 * max(0.0, min(1.0, result.identity_overlap))
    field_maximum = max((value for _, value in field_values), default=1) or 1

    record_rows = []
    total_records = sum(value for _, value, _ in record_values) or 1
    offset = 48.0
    for index, (name, value, color) in enumerate(record_values):
        width = 920 * value / total_records
        if value:
            record_rows.append(
                f'<rect x="{offset:.1f}" y="286" width="{width:.1f}" '
                f'height="24" fill="{color}"/>'
            )
        offset += width
        x = 54 + (index % 2) * 510
        y = 342 + (index // 2) * 28
        record_rows.append(
            f'<rect x="{x}" y="{y - 11}" width="10" height="10" rx="2" '
            f'fill="{color}"/><text x="{x + 18}" y="{y}" class="label">'
            f"{escape(_short(name))}: {value:,}</text>"
        )

    field_rows = []
    field_start = 430
    for index, (name, value) in enumerate(field_values):
        y = field_start + index * 36
        width = 550 * value / field_maximum
        field_rows.append(
            f'<text x="54" y="{y}" class="label">'
            f"{escape(name.replace('_', ' '))}</text>"
            f'<rect x="310" y="{y - 15}" width="550" height="14" rx="4" '
            'fill="#eaeef2"/>'
            f'<rect x="310" y="{y - 15}" width="{width:.1f}" height="14" '
            'rx="4" fill="#0969da"/>'
            f'<text x="1040" y="{y}" class="value">{value:,}</text>'
        )
    field_section = ""
    if field_rows:
        section_title = "Top transitions" if transition_values else "Changed fields"
        field_section = (
            f'<text x="40" y="402" class="section">{section_title}</text>'
            + "".join(field_rows)
        )
    title_text = escape(
        f"GenDiff comparison: {result.left_label} and {result.right_label}"
    )
    left_name = escape(_short(result.left_label))
    right_name = escape(_short(result.right_label))
    record_section = "".join(record_rows)
    if not record_section:
        record_section = (
            '<text x="54" y="332" class="muted">'
            "Run with --explain to match individual records.</text>"
        )

    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="{height}"
 viewBox="0 0 1120 {height}" role="img" aria-labelledby="title description">
<title id="title">{title_text}</title>
<desc id="description">{status}; {result.identity_overlap:.1%} identity overlap.</desc>
<style>
text{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;fill:#1f2328}}
.title{{font-size:28px;font-weight:700}}.status{{font-size:17px;font-weight:700}}
.sample{{font-size:19px;font-weight:700}}.path{{font-size:12px;fill:#59636e}}
.count{{font-size:24px;font-weight:700}}.muted{{font-size:13px;fill:#59636e}}
.section{{font-size:17px;font-weight:700}}.label{{font-size:13px}}
.value{{font-size:13px;font-weight:700;text-anchor:end}}
</style>
<rect width="1120" height="{height}" fill="#f6f8fa"/>
<rect x="24" y="24" width="1072" height="{height - 48}" rx="12" fill="#fff"
 stroke="#d0d7de"/>
<text x="48" y="70" class="title">GenDiff</text>
<text x="1048" y="66" class="status" text-anchor="end"
 style="fill:{status_color}">{status}</text>
<rect x="48" y="88" width="488" height="86" rx="8" fill="#f6f8fa"
 stroke="#d0d7de"/>
<text x="68" y="118" class="sample">{left_name}</text>
<text x="68" y="148" class="count">{result.left_records:,}</text>
<text x="68" y="165" class="muted">records</text>
<rect x="560" y="88" width="488" height="86" rx="8" fill="#f6f8fa"
 stroke="#d0d7de"/>
<text x="580" y="118" class="sample">{right_name}</text>
<text x="580" y="148" class="count">{result.right_records:,}</text>
<text x="580" y="165" class="muted">records</text>
<text x="48" y="206" class="section">Identity overlap</text>
<rect x="48" y="220" width="920" height="18" rx="6" fill="#eaeef2"/>
<rect x="48" y="220" width="{overlap_width:.1f}" height="18" rx="6" fill="#0969da"/>
<text x="1048" y="235" class="value">{result.identity_overlap:.1%}</text>
<text x="48" y="270" class="section">Record comparison</text>
<rect x="48" y="286" width="920" height="24" rx="6" fill="#eaeef2"/>
{record_section}
{field_section}
</svg>
'''


def render_igv_batch(result: ComparisonResult, reference: Path) -> str:
    if result.details is None or not result.details.loci:
        raise ValueError("no changed genomic loci are available for IGV")
    lines = [
        "new",
        f"genome {reference.resolve().as_uri()}",
        f"load {result.left.resolve().as_uri()}",
        f"load {result.right.resolve().as_uri()}",
    ]
    for locus in result.details.loci:
        lines.append(f"goto {locus}")
    lines.append("collapse")
    return "\n".join(lines) + "\n"


def write_text(path: Path, content: str, force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"output exists: {path}; use --force to replace it")
    if not path.parent.is_dir():
        raise FileNotFoundError(f"output directory not found: {path.parent}")
    with NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(content)
    try:
        os.replace(temporary, path)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
