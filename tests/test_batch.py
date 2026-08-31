import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pysam
import pytest

from gendiff.batch import BatchError, BatchPolicy, compare_manifest, read_manifest
from gendiff.cli import main


def _write_vcf(path: Path, genotype=(0, 1)) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.contigs.add("chr1", length=1000)
    header.formats.add("GT", 1, "String", "Genotype")
    header.add_sample("sample")
    with pysam.VariantFile(path, "w", header=header) as output:
        record = output.new_record(contig="chr1", start=9, alleles=("A", "C"), qual=60)
        record.filter.add("PASS")
        record.samples["sample"]["GT"] = genotype
        output.write(record)


def _manifest(path: Path, rows) -> None:
    lines = ["sample\tstage\tbefore\tafter"]
    lines.extend("\t".join(row) for row in rows)
    path.write_text("\n".join(lines) + "\n")


def test_batch_detects_first_divergent_stage_and_applies_policy(
    tmp_path: Path,
) -> None:
    baseline = tmp_path / "baseline.vcf"
    same = tmp_path / "same.vcf"
    changed = tmp_path / "changed.vcf"
    _write_vcf(baseline)
    _write_vcf(same)
    _write_vcf(changed, genotype=(1, 1))
    manifest = tmp_path / "comparisons.tsv"
    _manifest(
        manifest,
        [
            ("Sample A", "alignment", baseline.name, same.name),
            ("Sample A", "variants", baseline.name, changed.name),
            ("Sample B", "variants", baseline.name, same.name),
        ],
    )

    strict = compare_manifest(
        manifest,
        threads=4,
        profile="strict",
        reference=None,
        normalize_variants=False,
        ignore_tags=(),
        ignore_info=(),
        max_examples=0,
        progress=False,
        temp_dir=None,
        diff_dir=None,
        force=False,
        policy=BatchPolicy(),
    )

    assert not strict.passed
    assert [item.passed for item in strict.items] == [True, False, True]
    assert strict.earliest_divergence == {"Sample A": "variants"}
    assert strict.transitions == {"genotype sample: 0/1 → 1/1": 1}

    allowed = compare_manifest(
        manifest,
        threads=2,
        profile="strict",
        reference=None,
        normalize_variants=False,
        ignore_tags=(),
        ignore_info=(),
        max_examples=0,
        progress=False,
        temp_dir=None,
        diff_dir=None,
        force=False,
        policy=BatchPolicy(max_modified=1),
    )
    assert allowed.passed
    assert allowed.earliest_divergence == {"Sample A": "variants"}

    forbidden = compare_manifest(
        manifest,
        threads=2,
        profile="strict",
        reference=None,
        normalize_variants=False,
        ignore_tags=(),
        ignore_info=(),
        max_examples=0,
        progress=False,
        temp_dir=None,
        diff_dir=None,
        force=False,
        policy=BatchPolicy(max_modified=1, fail_transitions=("genotype",)),
    )
    assert not forbidden.passed
    assert "forbidden transition" in forbidden.items[1].reasons[0]


def test_batch_cli_writes_all_outputs_and_real_diffs(capsys, tmp_path: Path) -> None:
    baseline = tmp_path / "baseline.vcf"
    candidate = tmp_path / "candidate.vcf"
    _write_vcf(baseline)
    _write_vcf(candidate, genotype=(1, 1))
    manifest = tmp_path / "comparisons.tsv"
    _manifest(
        manifest,
        [("Sample A", "variants", baseline.name, candidate.name)],
    )
    html = tmp_path / "regression.html"
    json_report = tmp_path / "regression.json"
    junit = tmp_path / "regression.xml"
    multiqc = tmp_path / "gendiff_mqc.json"
    diffs = tmp_path / "diffs"
    tracks = tmp_path / "tracks"

    status = main(
        [
            "batch",
            str(manifest),
            "--max-modified",
            "1",
            "--html",
            str(html),
            "--json",
            str(json_report),
            "--junit",
            str(junit),
            "--multiqc",
            str(multiqc),
            "--write-diff",
            str(diffs),
            "--tracks",
            str(tracks),
        ]
    )

    assert status == 0
    assert "Comparisons: 1 (1 passed, 0 failed)" in capsys.readouterr().out
    assert "GenDiff pipeline regression" in html.read_text()
    assert "Stage overview" in html.read_text()
    assert "diff files" in html.read_text()
    aggregate = json.loads(json_report.read_text())
    assert aggregate["passed"] is True
    assert aggregate["earliest_divergence"] == {"Sample A": "variants"}
    suite = ET.parse(junit).getroot()
    assert suite.attrib["tests"] == "1"
    assert suite.attrib["failures"] == "0"
    multiqc_data = json.loads(multiqc.read_text())
    assert multiqc_data["id"] == "gendiff_pipeline_regression"

    batch_manifest = json.loads((diffs / "manifest.json").read_text())
    assert batch_manifest["kind"] == "batch"
    child = diffs / batch_manifest["comparisons"][0]["directory"]
    child_manifest = json.loads((child / "manifest.json").read_text())
    modified = child / child_manifest["files"]["modified_second"]
    with pysam.VariantFile(modified) as handle:
        assert [record.pos for record in handle] == [10]

    track_manifest = json.loads((tracks / "manifest.json").read_text())
    assert track_manifest["kind"] == "batch-tracks"
    track_child = tracks / track_manifest["comparisons"][0]["directory"]
    assert (track_child / "changes.bed").read_text().startswith("chr1\t9\t10")

    assert main(["batch", str(manifest), "--html", str(html)]) == 2


def test_manifest_validation_is_explicit(tmp_path: Path) -> None:
    manifest = tmp_path / "bad.tsv"
    manifest.write_text("sample\tstage\tbefore\nA\tvariants\ta.vcf\n")
    with pytest.raises(BatchError, match="missing columns: after"):
        read_manifest(manifest, "strict", None, False)

    manifest.write_text(
        "sample\tstage\tbefore\tafter\n"
        "A\tvariants\ta.vcf\tb.vcf\n"
        "A\tvariants\ta.vcf\tb.vcf\n"
    )
    with pytest.raises(BatchError, match="duplicate sample and stage"):
        read_manifest(manifest, "strict", None, False)
