from pathlib import Path

import pysam

from gendiff.compare import compare_files


def _write_bam(
    path: Path,
    records,
    program: str,
    *,
    tag_offset: int = 0,
    sort_order: str = "unsorted",
    include_unmapped: bool = False,
) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": sort_order},
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
            record.set_tag("NM", index + tag_offset)
            output.write(record)
        if include_unmapped:
            record = pysam.AlignedSegment(output.header)
            record.query_name = "unmapped-read"
            record.query_sequence = "ACGT"
            record.flag = 4
            record.reference_id = -1
            record.reference_start = -1
            record.mapping_quality = 0
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


def test_single_and_parallel_scans_agree(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20), (1, 30)], "tool")
    _write_bam(right, [(0, 20), (1, 40)], "tool")

    sequential = compare_files(left, right, threads=1)
    parallel = compare_files(left, right, threads=2)

    assert sequential == parallel


def test_explain_matches_modified_records(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20), (1, 30)], "tool")
    _write_bam(right, [(0, 20), (1, 40)], "tool")

    result = compare_files(left, right, explain=True)

    assert result.relationship == "likely the same dataset"
    assert result.identity_overlap == 1.0
    assert result.details is not None
    assert result.details.identical == 1
    assert result.details.modified == 1
    assert result.details.field_changes == {"mapping_quality": 1}


def test_core_profile_ignores_tags(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20)], "tool", tag_offset=0)
    _write_bam(right, [(0, 20)], "tool", tag_offset=1)

    assert not compare_files(left, right).equivalent
    assert compare_files(left, right, profile="core").equivalent
    assert compare_files(left, right, ignore_tags=("NM",)).equivalent


def test_indexed_sharding_preserves_results(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    records = [(0, 20), (1, 30)]
    _write_bam(left, records, "tool", sort_order="coordinate", include_unmapped=True)
    _write_bam(right, records, "tool", sort_order="coordinate", include_unmapped=True)
    pysam.index(str(left))
    pysam.index(str(right))

    sequential = compare_files(left, right, threads=1)
    sharded = compare_files(left, right, threads=4)

    assert sharded.left_records == 3
    assert sequential == sharded
