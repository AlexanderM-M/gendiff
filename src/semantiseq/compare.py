from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

from semantiseq.alignment import compare_alignments
from semantiseq.model import ComparisonResult
from semantiseq.profiles import fields_for
from semantiseq.regions import load_region_filter
from semantiseq.variant import compare_variants


class SemantiSeqError(Exception):
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
            raise SemantiSeqError(f"{option} cannot be empty")
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
    raise SemantiSeqError(f"unsupported file type: {path}")


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
    track_dir: Optional[Path] = None,
    difference_table: Optional[Path] = None,
    regions: Sequence[str] = (),
    regions_file: Optional[Path] = None,
    exclude_regions: Optional[Path] = None,
    cache_dir: Optional[Path] = None,
    force: bool = False,
) -> ComparisonResult:
    if threads < 1:
        raise SemantiSeqError("threads must be at least 1")
    for path in (left, right):
        if not path.is_file():
            raise SemantiSeqError(f"file not found: {path}")
    if reference is not None and not reference.is_file():
        raise SemantiSeqError(f"reference not found: {reference}")
    if max_examples < 0:
        raise SemantiSeqError("max examples must be at least 0")
    if temp_dir is not None and not temp_dir.is_dir():
        raise SemantiSeqError(f"temporary directory not found: {temp_dir}")
    if cache_dir is not None and cache_dir.exists() and not cache_dir.is_dir():
        raise SemantiSeqError(f"cache path is not a directory: {cache_dir}")
    for path, label in (
        (regions_file, "regions file"),
        (exclude_regions, "excluded regions file"),
    ):
        if path is not None and not path.is_file():
            raise SemantiSeqError(f"{label} not found: {path}")
    try:
        region_filter = load_region_filter(regions, regions_file, exclude_regions)
    except (OSError, ValueError) as error:
        raise SemantiSeqError(str(error)) from error
    if diff_dir is not None:
        explain = True
        resolved_diff = diff_dir.absolute()
        if resolved_diff in {left.absolute(), right.absolute()}:
            raise SemantiSeqError("diff output directory cannot be an input file")
    if track_dir is not None:
        explain = True
        resolved_tracks = track_dir.absolute()
        if resolved_tracks in {left.absolute(), right.absolute()}:
            raise SemantiSeqError("track output directory cannot be an input file")
        if diff_dir is not None and resolved_tracks == diff_dir.absolute():
            raise SemantiSeqError("diff and track output directories must differ")
    if difference_table is not None:
        explain = True
        if difference_table.absolute() in {left.absolute(), right.absolute()}:
            raise SemantiSeqError("difference table cannot be an input file")
    labels = _input_labels(left, right, left_label, right_label)

    left_kind = _kind(left)
    right_kind = _kind(right)
    if left_kind != right_kind:
        raise SemantiSeqError("inputs must be the same file type family")
    try:
        fields = fields_for(left_kind, profile)
    except ValueError as error:
        raise SemantiSeqError(str(error)) from error

    if left_kind == "alignment":
        if normalize_variants:
            raise SemantiSeqError("--normalize is only valid for VCF/BCF comparisons")
        if ignore_info:
            raise SemantiSeqError("--ignore-info is only valid for VCF/BCF comparisons")
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
            track_dir,
            difference_table,
            region_filter,
            cache_dir,
            force,
        )
    if ignore_tags:
        raise SemantiSeqError("--ignore-tag is only valid for BAM/CRAM comparisons")
    if normalize_variants:
        if reference is None:
            raise SemantiSeqError("--normalize requires --reference")
        reference_index = Path(f"{reference}.fai")
        if not reference_index.is_file():
            raise SemantiSeqError(
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
        track_dir,
        difference_table,
        region_filter,
        cache_dir,
        force,
    )
