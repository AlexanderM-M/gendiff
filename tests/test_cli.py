import json
from pathlib import Path

import pytest

from gendiff.cli import main


def test_missing_file_returns_error(capsys, tmp_path: Path) -> None:
    status = main([str(tmp_path / "one.bam"), str(tmp_path / "two.bam")])

    assert status == 2
    assert "file not found" in capsys.readouterr().err


def test_json_is_machine_readable(capsys, tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    status = main([str(vcf), str(vcf), "--json"])
    output = json.loads(capsys.readouterr().out)

    assert status == 0
    assert output["equivalent"] is True
    assert output["type"] == "variant"


def test_threads_must_be_positive(tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    with pytest.raises(SystemExit) as error:
        main([str(vcf), str(vcf), "--threads", "0"])

    assert error.value.code == 2


def test_html_report_is_written_and_protected(capsys, tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    report = tmp_path / "report.html"
    vcf.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    assert main([str(vcf), str(vcf), "--html", str(report)]) == 0
    assert "<!doctype html>" in report.read_text()
    assert main([str(vcf), str(vcf), "--html", str(report)]) == 2
    assert "output exists" in capsys.readouterr().err


def test_igv_batch_contains_inputs_and_changed_locus(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    reference = tmp_path / "reference.fa"
    batch = tmp_path / "differences.igv.batch"
    header = (
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=chr1,length=20>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    left.write_text(header + "chr1\t2\t.\tA\tC\t.\tPASS\t.\n")
    right.write_text(header + "chr1\t4\t.\tA\tC\t.\tPASS\t.\n")
    reference.write_text(">chr1\nAAAAAAAAAAAAAAAAAAAA\n")

    status = main(
        [
            str(left),
            str(right),
            "--reference",
            str(reference),
            "--igv-batch",
            str(batch),
        ]
    )

    content = batch.read_text()
    assert status == 1
    assert f"genome {reference.resolve().as_uri()}" in content
    assert f"load {left.resolve().as_uri()}" in content
    assert "goto chr1:" in content


def test_igv_validation_happens_before_writing_reports(tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    report = tmp_path / "report.html"
    batch = tmp_path / "differences.igv.batch"
    vcf.write_text(
        "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    status = main(
        [str(vcf), str(vcf), "--html", str(report), "--igv-batch", str(batch)]
    )

    assert status == 2
    assert not report.exists()
    assert not batch.exists()
