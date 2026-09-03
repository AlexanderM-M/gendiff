import gzip
from pathlib import Path

import pysam

from semantiseq.compare import compare_files
from semantiseq.report import render_html, render_svg


def _write_bam(
    path: Path,
    records,
    program: str,
    *,
    tag_offset: int = 0,
    sort_order: str = "unsorted",
    include_unmapped: bool = False,
    read_group: str = "",
) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": sort_order},
        "SQ": [{"SN": "chr1", "LN": 1000}],
        "PG": [{"ID": program, "PN": program}],
    }
    if read_group:
        header["RG"] = [{"ID": read_group, "SM": "sample"}]
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
            if read_group:
                record.set_tag("RG", read_group)
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
    assert result.details.transitions == {"MAPQ: 30 → 40": 1}
    assert result.details.distribution_shifts[0]["label"] == "Mapping quality"
    assert result.details.region_density[0]["difference_fraction"] == 0.25
    assert len(result.details.findings) == 3
    assert "Change transitions" in render_html(result)
    assert "Distribution shifts" in render_svg(result)
    assert "Mapping quality" in render_svg(result)


def test_writes_record_level_bam_diffs(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    output = tmp_path / "diff"
    _write_bam(left, [(0, 20), (1, 30), (2, 50)], "tool")
    _write_bam(right, [(0, 20), (1, 40), (3, 60)], "tool")

    result = compare_files(
        left,
        right,
        diff_dir=output,
        left_label="Before",
        right_label="After",
    )

    assert result.details is not None
    assert result.details.modified == 1
    assert result.details.left_only == 1
    assert result.details.right_only == 1
    expected = {
        "only_in_first": ["read-2"],
        "only_in_second": ["read-3"],
        "modified_first": ["read-1"],
        "modified_second": ["read-1"],
    }
    for key, names in expected.items():
        with pysam.AlignmentFile(result.artifacts[key], "rb") as handle:
            assert [record.query_name for record in handle] == names
    assert (output / "manifest.json").is_file()


def test_core_profile_ignores_tags(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20)], "tool", tag_offset=0)
    _write_bam(right, [(0, 20)], "tool", tag_offset=1)

    strict = compare_files(left, right, explain=True)
    assert not strict.equivalent
    assert strict.details is not None
    assert strict.details.transitions == {"tag value changed: NM": 1}
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


def test_region_filter_and_compressed_difference_table(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    table = tmp_path / "differences.tsv.gz"
    _write_bam(left, [(0, 20), (1, 30)], "tool")
    _write_bam(right, [(0, 20), (1, 40)], "tool")

    result = compare_files(
        left,
        right,
        regions=("chr1:12-12",),
        difference_table=table,
    )

    assert result.details is not None
    assert result.details.modified == 1
    assert result.details.contig_stats[0]["contig"] == "chr1"
    with gzip.open(table, "rt", encoding="utf-8") as handle:
        content = handle.read()
    assert "mapping_quality" in content
    assert result.artifacts["difference_table"] == str(table.absolute())


def test_indexed_detailed_matching_uses_shards(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    _write_bam(left, [(0, 20), (1, 30)], "tool", sort_order="coordinate")
    _write_bam(right, [(0, 20), (1, 40)], "tool", sort_order="coordinate")
    pysam.index(str(left))
    pysam.index(str(right))

    result = compare_files(left, right, threads=4, explain=True)

    assert result.details is not None
    assert result.details.modified == 1


def test_cache_reuses_detailed_scans_and_attributes_read_groups(tmp_path: Path) -> None:
    left = tmp_path / "left.bam"
    right = tmp_path / "right.bam"
    cache = tmp_path / "cache"
    _write_bam(left, [(0, 20)], "tool", read_group="lane-a")
    _write_bam(right, [(0, 40)], "tool", read_group="lane-a")

    first = compare_files(left, right, explain=True, cache_dir=cache)
    second = compare_files(left, right, explain=True, cache_dir=cache)

    assert first.cache == {"first": "miss", "second": "miss"}
    assert second.cache == {"first": "hit", "second": "hit"}
    assert second.details is not None
    assert second.details.group_changes == {"lane-a": 1}
