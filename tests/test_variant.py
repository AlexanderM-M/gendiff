from pathlib import Path

import pysam

from gendiff.compare import compare_files


def _write_vcf(
    path: Path,
    positions,
    genotype=(0, 1),
    source="caller",
    quality=60,
    depth=20,
) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value=source)
    header.contigs.add("chr1", length=1000)
    header.formats.add("GT", 1, "String", "Genotype")
    header.info.add("DP", 1, "Integer", "Depth")
    header.add_sample("sample")
    with pysam.VariantFile(path, "w", header=header) as output:
        for position in positions:
            record = output.new_record(
                contig="chr1", start=position - 1, alleles=("A", "C"), qual=quality
            )
            record.filter.add("PASS")
            record.info["DP"] = depth
            record.samples["sample"]["GT"] = genotype
            output.write(record)


def test_variant_ignores_record_order_and_provenance(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    _write_vcf(left, [10, 20], source="old-caller")
    _write_vcf(right, [20, 10], source="new-caller")

    result = compare_files(left, right)

    assert result.equivalent
    assert result.content_equal
    assert not result.metadata_equal


def test_variant_detects_changed_genotype(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    _write_vcf(left, [10])
    _write_vcf(right, [10], genotype=(1, 1))

    result = compare_files(left, right)

    assert not result.equivalent
    assert result.changed_fields == ["samples"]


def test_variant_explain_reports_sample_changes(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    _write_vcf(left, [10])
    _write_vcf(right, [10], genotype=(1, 1))

    result = compare_files(left, right, explain=True)

    assert result.details is not None
    assert result.details.modified == 1
    assert result.details.sample_changes == {"sample": 1}


def test_calls_profile_ignores_quality_and_info(tmp_path: Path) -> None:
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    _write_vcf(left, [10], quality=60, depth=20)
    _write_vcf(right, [10], quality=10, depth=5)

    assert not compare_files(left, right).equivalent
    assert compare_files(left, right, profile="calls").equivalent

    same_quality = tmp_path / "same-quality.vcf"
    _write_vcf(same_quality, [10], quality=60, depth=5)
    assert compare_files(left, same_quality, ignore_info=("DP",)).equivalent


def _write_indel(path: Path, position: int) -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.contigs.add("chr1", length=20)
    with pysam.VariantFile(path, "w", header=header) as output:
        output.write(
            output.new_record(contig="chr1", start=position - 1, alleles=("AA", "A"))
        )


def test_reference_normalization_equates_indel_representations(tmp_path: Path) -> None:
    reference = tmp_path / "reference.fa"
    reference.write_text(">chr1\nAAAAAAAAAAAAAAAAAAAA\n")
    pysam.faidx(str(reference))
    left = tmp_path / "left.vcf"
    right = tmp_path / "right.vcf"
    _write_indel(left, 2)
    _write_indel(right, 1)

    assert not compare_files(left, right).equivalent
    normalized = compare_files(
        left,
        right,
        reference=reference,
        normalize_variants=True,
    )

    assert normalized.equivalent
