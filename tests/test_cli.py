import json
from pathlib import Path

from gendiff.cli import main


def test_missing_file_returns_error(capsys, tmp_path: Path) -> None:
    status = main([str(tmp_path / "one.bam"), str(tmp_path / "two.bam")])

    assert status == 2
    assert "file not found" in capsys.readouterr().err


def test_json_is_machine_readable(capsys, tmp_path: Path) -> None:
    vcf = tmp_path / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )

    status = main([str(vcf), str(vcf), "--json"])
    output = json.loads(capsys.readouterr().out)

    assert status == 0
    assert output["equivalent"] is True
    assert output["type"] == "variant"
