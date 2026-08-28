from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

from gendiff.alignment import compare_alignments
from gendiff.model import ComparisonResult
from gendiff.profiles import fields_for
from gendiff.variant import compare_variants


class GenDiffError(Exception):
    """Raised when inputs cannot be compared."""


def _kind(path: Path) -> str:
    name = path.name.lower()
    if name.endswith((".bam", ".cram")):
        return "alignment"
    if name.endswith((".vcf", ".vcf.gz", ".bcf")):
        return "variant"
    raise GenDiffError(f"unsupported file type: {path}")


def compare_files(
    left: Path,
    right: Path,
    reference: Optional[Path] = None,
    threads: int = 2,
    *,
    profile: str = "strict",
    ignore_tags: Sequence[str] = (),
    ignore_info: Sequence[str] = (),
    normalize_variants: bool = False,
    explain: bool = False,
    max_examples: int = 10,
    progress: bool = False,
    temp_dir: Optional[Path] = None,
) -> ComparisonResult:
    if threads < 1:
        raise GenDiffError("threads must be at least 1")
    for path in (left, right):
        if not path.is_file():
            raise GenDiffError(f"file not found: {path}")
    if reference is not None and not reference.is_file():
        raise GenDiffError(f"reference not found: {reference}")
    if max_examples < 0:
        raise GenDiffError("max examples must be at least 0")
    if temp_dir is not None and not temp_dir.is_dir():
        raise GenDiffError(f"temporary directory not found: {temp_dir}")

    left_kind = _kind(left)
    right_kind = _kind(right)
    if left_kind != right_kind:
        raise GenDiffError("inputs must be the same file type family")
    try:
        fields = fields_for(left_kind, profile)
    except ValueError as error:
        raise GenDiffError(str(error)) from error

    if left_kind == "alignment":
        if normalize_variants:
            raise GenDiffError("--normalize is only valid for VCF/BCF comparisons")
        if ignore_info:
            raise GenDiffError("--ignore-info is only valid for VCF/BCF comparisons")
        return compare_alignments(
            left,
            right,
            reference,
            threads,
            fields,
            profile,
            ignore_tags,
            explain,
            max_examples,
            progress,
            temp_dir,
        )
    if ignore_tags:
        raise GenDiffError("--ignore-tag is only valid for BAM/CRAM comparisons")
    if normalize_variants:
        if reference is None:
            raise GenDiffError("--normalize requires --reference")
        reference_index = Path(f"{reference}.fai")
        if not reference_index.is_file():
            raise GenDiffError(
                f"reference index not found: {reference_index} "
                "(create it with samtools faidx)"
            )
    return compare_variants(
        left,
        right,
        reference,
        normalize_variants,
        threads,
        fields,
        profile,
        ignore_info,
        explain,
        max_examples,
        progress,
        temp_dir,
    )
