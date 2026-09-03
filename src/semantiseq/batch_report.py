"""Terminal, HTML, JUnit, and MultiQC output for batch comparisons."""

from __future__ import annotations

import json
import os
import xml.etree.ElementTree as ET
from html import escape
from pathlib import Path
from typing import Iterable, Optional

from semantiseq.batch import BatchItem, BatchResult


def render_batch_text(result: BatchResult) -> str:
    passed = sum(item.passed for item in result.items)
    lines = [
        f"Passed: {'yes' if result.passed else 'no'}",
        f"Comparisons: {len(result.items)} ({passed} passed, "
        f"{len(result.items) - passed} failed)",
        "Results:",
    ]
    for item in result.items:
        details = item.result.details
        modified = details.modified if details else 0
        only = details.left_only + details.right_only if details else 0
        lines.append(
            f"  {'PASS' if item.passed else 'FAIL'}  {item.entry.sample} / "
            f"{item.entry.stage}: overlap {item.result.identity_overlap:.1%}, "
            f"modified {modified:,}, only {only:,}"
        )
        for reason in item.reasons:
            lines.append(f"        {reason}")
        for warning in item.warnings:
            lines.append(f"        WARNING: {warning}")
        if item.matched_rules:
            lines.append(f"        Policy: {', '.join(item.matched_rules)}")
    if result.earliest_divergence:
        lines.append("Earliest divergence:")
        for sample, stage in result.earliest_divergence.items():
            lines.append(f"  {sample}: {stage}")
    return "\n".join(lines)


def _table(headers: Iterable[str], rows: Iterable[Iterable[str]]) -> str:
    head = "".join(f"<th>{escape(value)}</th>" for value in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{value}</td>" for value in row) + "</tr>" for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _artifact_link(item: BatchItem, report_path: Optional[Path]) -> str:
    manifests = [
        ("diff files", item.result.artifacts.get("manifest")),
        ("tracks", item.result.artifacts.get("track_manifest")),
    ]
    manifests = [(label, value) for label, value in manifests if value is not None]
    if not manifests:
        return "—"
    links = []
    for label, manifest in manifests:
        directory = Path(manifest).parent
        if report_path is None:
            links.append(escape(label))
            continue
        relative = os.path.relpath(directory, report_path.absolute().parent)
        href = escape(Path(relative).as_posix(), quote=True)
        links.append(f'<a href="{href}/manifest.json">{escape(label)}</a>')
    return " · ".join(links)


def render_batch_html(result: BatchResult, report_path: Optional[Path] = None) -> str:
    passed = sum(item.passed for item in result.items)
    earliest = result.earliest_divergence
    comparison_rows = []
    item_lookup = {(item.entry.sample, item.entry.stage): item for item in result.items}
    for item in result.items:
        details = item.result.details
        modified = details.modified if details else 0
        only = details.left_only + details.right_only if details else 0
        reasons = "; ".join(item.reasons) or "—"
        warnings = "; ".join(item.warnings) or "—"
        policy_rules = ", ".join(item.matched_rules) or "defaults"
        comparison_rows.append(
            (
                escape(item.entry.sample),
                escape(item.entry.stage),
                (
                    '<strong class="pass">PASS</strong>'
                    if item.passed
                    else '<strong class="fail">FAIL</strong>'
                ),
                f"{item.result.identity_overlap:.1%}",
                f"{modified:,}",
                f"{only:,}",
                escape(reasons),
                escape(warnings),
                escape(policy_rules),
                _artifact_link(item, report_path),
            )
        )

    matrix_rows = []
    for sample in result.sample_order:
        cells = [escape(sample)]
        for stage in result.stage_order:
            item = item_lookup.get((sample, stage))
            if item is None:
                cells.append("—")
                continue
            status = "PASS" if item.passed else "FAIL"
            status_class = "pass" if item.passed else "fail"
            first = earliest.get(sample) == stage
            suffix = " · first" if first else ""
            cells.append(f'<strong class="{status_class}">{status}{suffix}</strong>')
        matrix_rows.append(tuple(cells))
    matrix = _table(("Sample", *result.stage_order), matrix_rows)

    earliest_section = ""
    if earliest:
        earliest_section = (
            "<section><h2>Earliest divergence</h2>"
            + _table(
                ("Sample", "Stage"),
                ((escape(sample), escape(stage)) for sample, stage in earliest.items()),
            )
            + "</section>"
        )
    transition_section = ""
    if result.transitions:
        transition_section = (
            "<section><h2>Top transitions</h2>"
            + _table(
                ("Transition", "Records"),
                (
                    (escape(name), f"{count:,}")
                    for name, count in list(result.transitions.items())[:20]
                ),
            )
            + "</section>"
        )
    comparison_table = _table(
        (
            "Sample",
            "Stage",
            "Status",
            "Overlap",
            "Modified",
            "Only",
            "Reason",
            "Warning",
            "Policy",
            "Artifacts",
        ),
        comparison_rows,
    )
    status = "Passed" if result.passed else "Failed"
    status_class = "pass" if result.passed else "fail"
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SemantiSeq pipeline regression</title>
<style>
:root{{--bg:#f6f8fa;--panel:#fff;--text:#1f2328;--muted:#59636e;
--line:#d0d7de;--ok:#1a7f37;--bad:#cf222e}}*{{box-sizing:border-box}}
body{{margin:0;background:var(--bg);color:var(--text);
font:14px/1.5 system-ui,sans-serif}}
main{{max-width:1200px;margin:auto;padding:32px 20px}}
h1{{margin:0 0 4px;font-size:28px}}h2{{font-size:18px;margin:0 0 14px}}
.muted{{color:var(--muted)}}.pass{{color:var(--ok)}}.fail{{color:var(--bad)}}
.cards{{display:grid;grid-template-columns:repeat(3,minmax(120px,1fr));gap:12px}}
.cards article,section{{background:var(--panel);border:1px solid var(--line);
border-radius:8px;padding:18px}}section{{margin-top:16px;overflow:auto}}
.cards strong{{display:block;font-size:24px}}.cards span{{color:var(--muted)}}
table{{border-collapse:collapse;width:100%;font-size:13px}}th,td{{padding:8px;
border-bottom:1px solid var(--line);text-align:left;vertical-align:top}}
th{{background:#f6f8fa;white-space:nowrap}}a{{color:#0969da}}
@media(max-width:700px){{.cards{{grid-template-columns:1fr}}}}
</style></head><body><main>
<h1>SemantiSeq pipeline regression</h1>
<div class="{status_class}"><strong>{status}</strong></div>
<p class="muted">{escape(result.manifest.name)}</p>
<div class="cards">
<article><strong>{len(result.items)}</strong><span>Comparisons</span></article>
<article><strong class="pass">{passed}</strong><span>Passed</span></article>
<article><strong class="fail">{len(result.items) - passed}</strong>
<span>Failed</span></article>
</div>
<section><h2>Stage overview</h2>{matrix}</section>
<section><h2>Comparisons</h2>
{comparison_table}
</section>{earliest_section}{transition_section}
</main></body></html>
"""


def render_junit(result: BatchResult) -> str:
    failures = sum(not item.passed for item in result.items)
    duration = sum(item.duration_seconds for item in result.items)
    suite = ET.Element(
        "testsuite",
        {
            "name": "semantiseq pipeline regression",
            "tests": str(len(result.items)),
            "failures": str(failures),
            "errors": "0",
            "time": f"{duration:.6f}",
        },
    )
    properties = ET.SubElement(suite, "properties")
    for name, value in result.policy.to_dict().items():
        ET.SubElement(
            properties,
            "property",
            {"name": name, "value": json.dumps(value, separators=(",", ":"))},
        )
    for item in result.items:
        case = ET.SubElement(
            suite,
            "testcase",
            {
                "classname": f"semantiseq.batch.{item.entry.sample}",
                "name": item.entry.stage,
                "time": f"{item.duration_seconds:.6f}",
            },
        )
        if not item.passed:
            message = "; ".join(item.reasons)
            failure = ET.SubElement(
                case,
                "failure",
                {"message": message, "type": "SemantiSeqRegression"},
            )
            failure.text = message
        decision_text = "\n".join(
            [
                *(f"rule: {name}" for name in item.matched_rules),
                *(f"warning: {warning}" for warning in item.warnings),
                *item.policy_trace,
            ]
        )
        if decision_text:
            ET.SubElement(case, "system-out").text = decision_text
    ET.indent(suite, space="  ")
    return ET.tostring(suite, encoding="utf-8", xml_declaration=True).decode("utf-8")


def render_multiqc(result: BatchResult) -> str:
    earliest = result.earliest_divergence
    data = {}
    for item in result.items:
        details = item.result.details
        key = f"{item.entry.sample} — {item.entry.stage}"
        data[key] = {
            "stage": item.entry.stage,
            "status": "PASS" if item.passed else "FAIL",
            "identity_overlap": round(item.result.identity_overlap * 100.0, 4),
            "modified": details.modified if details else 0,
            "only": details.left_only + details.right_only if details else 0,
            "first_divergence": earliest.get(item.entry.sample) == item.entry.stage,
        }
    payload = {
        "id": "semantiseq_pipeline_regression",
        "section_name": "SemantiSeq pipeline regression",
        "description": "Semantic differences between baseline and candidate outputs.",
        "plot_type": "table",
        "pconfig": {"id": "semantiseq_regression", "title": "SemantiSeq comparisons"},
        "headers": {
            "stage": {"title": "Stage"},
            "status": {"title": "Status"},
            "identity_overlap": {
                "title": "Identity overlap",
                "suffix": "%",
                "min": 0,
                "max": 100,
            },
            "modified": {"title": "Modified"},
            "only": {"title": "Only in one input"},
            "first_divergence": {"title": "First divergence"},
        },
        "data": data,
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"
