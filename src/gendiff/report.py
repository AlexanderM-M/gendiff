from __future__ import annotations

import json
import os
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


def render_html(result: ComparisonResult) -> str:
    details = result.details
    status = "Equivalent" if result.equivalent else "Different"
    status_class = "ok" if result.equivalent else "different"
    detail_sections: List[str] = []
    if details is not None:
        detail_sections.append(
            "<section><h2>Record comparison</h2><div class='cards'>"
            f"<article><strong>{details.identical:,}</strong>"
            "<span>Identical</span></article>"
            f"<article><strong>{details.modified:,}</strong>"
            "<span>Modified</span></article>"
            f"<article><strong>{details.left_only:,}</strong>"
            "<span>Only in left</span></article>"
            f"<article><strong>{details.right_only:,}</strong>"
            "<span>Only in right</span></article>"
            "</div></section>"
        )
        detail_sections.append(
            "<section><h2>Changed fields</h2>"
            f"{_bar_rows(details.field_changes)}</section>"
        )
        if details.top_regions:
            detail_sections.append(
                "<section><h2>Most affected regions</h2>"
                + _table(
                    ("Region", "Changes"),
                    (
                        (str(item["region"]), f"{item['changes']:,}")
                        for item in details.top_regions
                    ),
                )
                + "</section>"
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
                "<section><h2>Examples</h2>"
                + _table(
                    ("Status", "Identity", "Changed fields", "Left", "Right"),
                    (
                        (
                            str(item["status"]),
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
                + "</section>"
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
<dt>Relationship</dt><dd>{escape(result.relationship)}</dd>
<dt>Identity overlap</dt><dd>{result.identity_overlap:.1%}</dd>
<dt>Profile</dt><dd>{escape(result.profile)}</dd>
<dt>Left</dt><dd>{escape(str(result.left))}</dd>
<dt>Right</dt><dd>{escape(str(result.right))}</dd>
<dt>Records</dt><dd>{result.left_records:,} / {result.right_records:,}</dd>
<dt>Logical records</dt><dd>{"same" if result.content_equal else "different"}</dd>
<dt>Structural header</dt><dd>{"same" if result.structure_equal else "different"}</dd>
<dt>Metadata header</dt>
<dd>{"same" if result.metadata_equal else "different"} (informational)</dd>
</dl></section>
{"".join(detail_sections)}
</main></body></html>
"""


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
