from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

from gendiff.alignment import compare_alignments
from gendiff.model import ComparisonResult
from gendiff.profiles import fields_for
from gendiff.variant import compare_variants


class GenDiffError(Exception):
    """Raised when inputs cannot be compared."""


_FORMAT_SUFFIXES = (".vcf.gz", ".cram", ".bam", ".bcf", ".vcf")
_PROCESSING_SUFFIXES = (
    ".merged",
    ".sorted",
    ".deduplicated",
    ".dedup",
    ".markdup",
    ".aligned",
)


def _inferred_label(path: Path) -> str:
    label = path.name
    lowered = label.lower()
    for suffix in _FORMAT_SUFFIXES:
        if lowered.endswith(suffix):
            label = label[: -len(suffix)]
            break
    lowered = label.lower()
    for suffix in _PROCESSING_SUFFIXES:
        if lowered.endswith(suffix):
            label = label[: -len(suffix)]
            break
    return label or path.name


def _input_labels(
    left: Path,
    right: Path,
    left_label: Optional[str],
    right_label: Optional[str],
) -> tuple[str, str]:
    for value, option in (
        (left_label, "--name-a"),
        (right_label, "--name-b"),
    ):
        if value is not None and not value.strip():
            raise GenDiffError(f"{option} cannot be empty")
    first = left_label.strip() if left_label else _inferred_label(left)
    second = right_label.strip() if right_label else _inferred_label(right)
    if first == second:
        if left.resolve() != right.resolve():
            first = f"{left.parent.name}/{first}"
            second = f"{right.parent.name}/{second}"
        else:
            first = f"{first} (A)"
            second = f"{second} (B)"
    return first, second


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
    left_label: Optional[str] = None,
    right_label: Optional[str] = None,
    diff_dir: Optional[Path] = None,
    force: bool = False,
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
    if diff_dir is not None:
        explain = True
        resolved_diff = diff_dir.absolute()
        if resolved_diff in {left.absolute(), right.absolute()}:
            raise GenDiffError("diff output directory cannot be an input file")
    labels = _input_labels(left, right, left_label, right_label)

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
            *labels,
            ignore_tags,
            explain,
            max_examples,
            progress,
            temp_dir,
            diff_dir,
            force,
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
        *labels,
        ignore_info,
        explain,
        max_examples,
        progress,
        temp_dir,
        diff_dir,
        force,
    )
