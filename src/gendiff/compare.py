from __future__ import annotations

from pathlib import Path
from typing import Optional

from gendiff.alignment import compare_alignments
from gendiff.model import ComparisonResult
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
    left: Path, right: Path, reference: Optional[Path] = None
) -> ComparisonResult:
    for path in (left, right):
        if not path.is_file():
            raise GenDiffError(f"file not found: {path}")
    if reference is not None and not reference.is_file():
        raise GenDiffError(f"reference not found: {reference}")

    left_kind = _kind(left)
    right_kind = _kind(right)
    if left_kind != right_kind:
        raise GenDiffError("inputs must be the same file type family")
    if reference is not None and left_kind != "alignment":
        raise GenDiffError("--reference is only valid for BAM/CRAM comparisons")

    if left_kind == "alignment":
        return compare_alignments(left, right, reference)
    return compare_variants(left, right)
