from pathlib import Path

import pysam

from gendiff.compare import compare_files


def _write_vcf(path: Path, positions, genotype=(0, 1), source="caller") -> None:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value=source)
    header.contigs.add("chr1", length=1000)
    header.formats.add("GT", 1, "String", "Genotype")
    header.add_sample("sample")
    with pysam.VariantFile(path, "w", header=header) as output:
        for position in positions:
            record = output.new_record(
                contig="chr1", start=position - 1, alleles=("A", "C"), qual=60
            )
            record.filter.add("PASS")
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
