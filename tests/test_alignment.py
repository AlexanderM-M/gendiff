from pathlib import Path

import pysam

from gendiff.compare import compare_files


def _write_bam(path: Path, records, program: str) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "PG": [{"ID": program, "PN": program}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for index, mapq in records:
            record = pysam.AlignedSegment(output.header)
            record.query_name = f"read-{index}"
            record.query_sequence = "ACGT"
            record.flag = 0
            record.reference_id = 0
            record.reference_start = 10 + index
            record.mapping_quality = mapq
            record.cigarstring = "4M"
            record.query_qualities = pysam.qualitystring_to_array("IIII")
            output.write(record)


def test_alignment_ignores_order_and_provenance(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20), (1, 30)], "old-tool")
    _write_bam(right, [(1, 30), (0, 20)], "new-tool")

    result = compare_files(left, right)

    assert result.equivalent
    assert result.content_equal
    assert result.metadata_equal is False


def test_alignment_detects_changed_mapping_quality(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20)], "tool")
    _write_bam(right, [(0, 40)], "tool")

    result = compare_files(left, right)

    assert not result.equivalent
    assert result.changed_fields == ["mapping_quality"]
