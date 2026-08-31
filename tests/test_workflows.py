import json
from pathlib import Path

import pysam

from gendiff.baseline import baseline_report, check_baseline, create_baseline
from gendiff.batch import BatchPolicy, compare_manifest
from gendiff.cli import main
from gendiff.policy import load_policy


def _write_vcf(path: Path, genotype=(0, 1)) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.contigs.add("chr1", length=1000)
    header.formats.add("GT", 1, "String", "Genotype")
    header.add_sample("sample")
    with pysam.VariantFile(path, "w", header=header) as output:
        record = output.new_record(contig="chr1", start=9, alleles=("A", "C"))
        record.filter.add("PASS")
        record.samples["sample"]["GT"] = genotype
        output.write(record)


def _manifest(path: Path, before: Path, after: Path, stage: str = "variants") -> None:
    path.write_text(
        f"sample\tstage\tbefore\tafter\nExample\t{stage}\t{before.name}\t{after.name}\n"
    )


def test_stage_policy_records_warnings_and_trace(tmp_path: Path) -> None:
    before = tmp_path / "before.vcf"
    after = tmp_path / "after.vcf"
    manifest = tmp_path / "comparisons.tsv"
    policy_path = tmp_path / "policy.json"
    _write_vcf(before)
    _write_vcf(after, genotype=(1, 1))
    _manifest(manifest, before, after)
    policy_path.write_text(
        json.dumps(
            {
                "gendiff_policy": 1,
                "defaults": {
                    "max_modified": None,
                    "transition_default": "allow",
                },
                "rules": [
                    {
                        "name": "variant genotype review",
                        "match": {"stage": "variants", "kind": "variant"},
                        "set": {"deny_fields": ["samples"], "severity": "warning"},
                    }
                ],
            }
        )
    )
    result = compare_manifest(
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
        policy=BatchPolicy(),
        policy_document=load_policy(policy_path),
    )

    assert result.passed
    assert result.items[0].matched_rules == ("variant genotype review",)
    assert result.items[0].warnings == ("forbidden changed fields: samples",)
    assert "matched rule: variant genotype review" in result.items[0].policy_trace


def test_manifest_generation_pairs_relative_paths(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline" / "alignment"
    candidate = tmp_path / "candidate" / "alignment"
    baseline.mkdir(parents=True)
    candidate.mkdir(parents=True)
    _write_vcf(baseline / "example.vcf")
    _write_vcf(candidate / "example.vcf")
    output = tmp_path / "comparisons.tsv"

    assert (
        main(
            [
                "manifest",
                str(baseline.parent),
                str(candidate.parent),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    content = output.read_text()
    assert "example\talignment" in content
    assert "baseline/alignment/example.vcf" in content
    assert "candidate/alignment/example.vcf" in content


def test_compact_baseline_detects_candidate_change(tmp_path: Path) -> None:
    before = tmp_path / "before.vcf"
    same = tmp_path / "same.vcf"
    changed = tmp_path / "changed.vcf"
    manifest = tmp_path / "comparisons.tsv"
    baseline_path = tmp_path / "baseline.gendiff.json"
    _write_vcf(before)
    _write_vcf(same)
    _write_vcf(changed, genotype=(1, 1))
    _manifest(manifest, before, same)
    payload = create_baseline(
        manifest,
        threads=2,
        profile="strict",
        reference=None,
        normalize_variants=False,
        ignore_tags=(),
        ignore_info=(),
        progress=False,
    )
    baseline_path.write_text(json.dumps(payload))
    assert "read_name" not in baseline_path.read_text()
    assert baseline_report(
        check_baseline(
            baseline_path,
            manifest,
            threads=2,
            profile="strict",
            reference=None,
            normalize_variants=False,
            progress=False,
        )
    )["passed"]

    _manifest(manifest, before, changed)
    report = baseline_report(
        check_baseline(
            baseline_path,
            manifest,
            threads=2,
            profile="strict",
            reference=None,
            normalize_variants=False,
            progress=False,
        )
    )
    assert not report["passed"]
    assert report["comparisons"][0]["changed_fields"] == ["samples"]


def test_track_outputs_are_standard_text_files(tmp_path: Path) -> None:
    before = tmp_path / "before.vcf"
    after = tmp_path / "after.vcf"
    tracks = tmp_path / "tracks"
    _write_vcf(before)
    _write_vcf(after, genotype=(1, 1))

    assert main([str(before), str(after), "--tracks", str(tracks)]) == 1
    assert (tracks / "changes.bed").read_text() == "chr1\t9\t10\tmodified;x1\n"
    assert (tracks / "difference-density.bedgraph").read_text() == (
        "chr1\t0\t1000\t1\n"
    )
    assert "chr1\t1" in (tracks / "contig-summary.tsv").read_text()
