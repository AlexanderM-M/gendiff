from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from gendiff import __version__
from gendiff.batch import BatchError, BatchPolicy, compare_manifest
from gendiff.batch_report import (
    render_batch_html,
    render_batch_text,
    render_junit,
    render_multiqc,
)
from gendiff.compare import GenDiffError, compare_files
from gendiff.model import ComparisonResult
from gendiff.profiles import ALL_PROFILES
from gendiff.report import render_html, render_igv_batch, render_svg, write_text


def _positive_int(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return number


def _nonnegative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("must be at least 0")
    return number


def _proportion(value: str) -> float:
    number = float(value)
    if not 0.0 <= number <= 1.0:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return number


def _single_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gendiff",
        description="Compare the logical contents of two genomics files.",
        epilog="Use 'gendiff batch --help' for pipeline regression comparisons.",
    )
    parser.add_argument("sample_a", type=Path, help="first BAM, CRAM, VCF, or BCF")
    parser.add_argument("sample_b", type=Path, help="second file")
    parser.add_argument("--name-a", help="display name for the first sample")
    parser.add_argument("--name-b", help="display name for the second sample")
    parser.add_argument(
        "--reference", type=Path, help="reference FASTA for CRAM, normalization, or IGV"
    )
    parser.add_argument(
        "--threads",
        type=_positive_int,
        default=2,
        help="total worker threads (default: 2)",
    )
    parser.add_argument(
        "--profile",
        choices=ALL_PROFILES,
        default="strict",
        help="comparison policy (default: strict)",
    )
    parser.add_argument(
        "--ignore-tag",
        action="append",
        default=[],
        metavar="TAG",
        help="ignore a BAM tag; repeat as needed",
    )
    parser.add_argument(
        "--ignore-info",
        action="append",
        default=[],
        metavar="FIELD",
        help="ignore a VCF INFO field; repeat as needed",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="normalize VCF/BCF alleles using --reference",
    )
    parser.add_argument(
        "--explain", action="store_true", help="match records and explain differences"
    )
    parser.add_argument(
        "--max-examples",
        type=_nonnegative_int,
        default=10,
        help="maximum detailed examples (default: 10)",
    )
    parser.add_argument(
        "--progress", action="store_true", help="report progress to standard error"
    )
    parser.add_argument(
        "--temp-dir",
        type=Path,
        help="temporary storage for detailed record matching",
    )
    parser.add_argument("--html", type=Path, help="write a self-contained HTML report")
    parser.add_argument("--svg", type=Path, help="write a standalone SVG summary")
    parser.add_argument("--igv-batch", type=Path, help="write an IGV batch file")
    parser.add_argument(
        "--write-diff",
        type=Path,
        metavar="DIR",
        help="write only-in and modified records as BAM or VCF files",
    )
    parser.add_argument("--force", action="store_true", help="replace outputs")
    parser.add_argument("--json", action="store_true", help="write JSON output")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return parser


def _batch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gendiff batch",
        description="Compare baseline and candidate pipeline outputs from a TSV.",
    )
    parser.add_argument(
        "manifest",
        type=Path,
        help="TSV with sample, stage, before, and after columns",
    )
    parser.add_argument(
        "--threads",
        type=_positive_int,
        default=2,
        help="total worker threads across comparisons (default: 2)",
    )
    parser.add_argument(
        "--profile",
        choices=ALL_PROFILES,
        default="strict",
        help="default comparison profile (default: strict)",
    )
    parser.add_argument("--reference", type=Path, help="default reference FASTA")
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="normalize variant inputs unless overridden by the manifest",
    )
    parser.add_argument(
        "--ignore-tag",
        action="append",
        default=[],
        metavar="TAG",
        help="ignore a BAM tag; repeat as needed",
    )
    parser.add_argument(
        "--ignore-info",
        action="append",
        default=[],
        metavar="FIELD",
        help="ignore a VCF INFO field; repeat as needed",
    )
    parser.add_argument(
        "--max-modified",
        type=_nonnegative_int,
        default=0,
        help="maximum modified records per comparison (default: 0)",
    )
    parser.add_argument(
        "--max-only",
        type=_nonnegative_int,
        default=0,
        help="maximum combined sample-only records per comparison (default: 0)",
    )
    parser.add_argument(
        "--min-overlap",
        type=_proportion,
        default=0.0,
        help="minimum identity overlap from 0 to 1 (default: 0)",
    )
    parser.add_argument(
        "--allow-structural-differences",
        action="store_true",
        help="do not fail when structural headers differ",
    )
    parser.add_argument(
        "--fail-transition",
        action="append",
        default=[],
        metavar="TEXT",
        help="fail if a transition contains this text; repeat as needed",
    )
    parser.add_argument(
        "--max-examples",
        type=_nonnegative_int,
        default=0,
        help="examples retained per comparison (default: 0)",
    )
    parser.add_argument(
        "--progress", action="store_true", help="report progress to standard error"
    )
    parser.add_argument("--temp-dir", type=Path, help="temporary matching storage")
    parser.add_argument("--html", type=Path, help="write aggregate HTML")
    parser.add_argument("--json", type=Path, help="write aggregate JSON")
    parser.add_argument("--junit", type=Path, help="write JUnit XML")
    parser.add_argument(
        "--multiqc",
        type=Path,
        metavar="FILE",
        help="write MultiQC custom-content JSON",
    )
    parser.add_argument(
        "--write-diff",
        type=Path,
        metavar="DIR",
        help="write record-level diff files for every comparison",
    )
    parser.add_argument("--force", action="store_true", help="replace outputs")
    return parser


def _status(value: bool) -> str:
    return "same" if value else "different"


def _render(result: ComparisonResult) -> str:
    lines = [
        f"Equivalent: {'yes' if result.equivalent else 'no'}",
        f"Type: {result.kind}",
        f"Relationship: {result.relationship}",
        f"Identity overlap: {result.identity_overlap:.1%}",
        f"Profile: {result.profile}",
        "Inputs:",
        f"  {result.left_label}: {result.left} ({result.left_records:,} records)",
        f"  {result.right_label}: {result.right} ({result.right_records:,} records)",
        f"Logical records: {_status(result.content_equal)}",
        f"Structural header: {_status(result.structure_equal)}",
        f"Metadata header: {_status(result.metadata_equal)} (informational)",
    ]
    if result.changed_fields:
        fields = ", ".join(field.replace("_", " ") for field in result.changed_fields)
        lines.append(f"Changed fields: {fields}")
    if result.details is not None:
        details = result.details
        lines.extend(
            [
                "Record differences:",
                f"  Identical: {details.identical:,}",
                f"  Modified: {details.modified:,}",
                f"  Only in {result.left_label}: {details.left_only:,}",
                f"  Only in {result.right_label}: {details.right_only:,}",
            ]
        )
        if details.field_changes:
            changes = ", ".join(
                f"{name.replace('_', ' ')}={count:,}"
                for name, count in details.field_changes.items()
            )
            lines.append(f"Field change counts: {changes}")
        if details.transitions:
            transitions = ", ".join(
                f"{name}={count:,}"
                for name, count in list(details.transitions.items())[:10]
            )
            lines.append(f"Top transitions: {transitions}")
        if details.top_regions:
            regions = ", ".join(
                f"{item['region']} ({item['changes']:,})"
                for item in details.top_regions[:5]
            )
            lines.append(f"Top regions: {regions}")
    if result.artifacts:
        lines.append(f"Diff files: {Path(result.artifacts['manifest']).parent}")
    return "\n".join(lines)


def _single_main(argv: List[str]) -> int:
    args = _single_parser().parse_args(argv)
    try:
        outputs = [output for output in (args.html, args.svg, args.igv_batch) if output]
        for output in outputs:
            if output and output.exists() and not args.force:
                raise FileExistsError(
                    f"output exists: {output}; use --force to replace it"
                )
        if len(set(outputs)) != len(outputs):
            raise ValueError("report outputs must use different paths")
        if args.write_diff and args.write_diff.absolute() in {
            output.absolute() for output in outputs
        }:
            raise ValueError(
                "diff directory and report outputs must use different paths"
            )
        if args.igv_batch and args.reference is None:
            raise ValueError("--igv-batch requires --reference")
        explain = any(
            (
                args.explain,
                args.html is not None,
                args.svg is not None,
                args.igv_batch,
                args.write_diff is not None,
            )
        )
        result = compare_files(
            args.sample_a,
            args.sample_b,
            args.reference,
            args.threads,
            profile=args.profile,
            ignore_tags=args.ignore_tag,
            ignore_info=args.ignore_info,
            normalize_variants=args.normalize,
            explain=explain,
            max_examples=args.max_examples,
            progress=args.progress,
            temp_dir=args.temp_dir,
            left_label=args.name_a,
            right_label=args.name_b,
            diff_dir=args.write_diff,
            force=args.force,
        )
        html_report = render_html(result) if args.html else None
        svg_report = render_svg(result) if args.svg else None
        igv_report = (
            render_igv_batch(result, args.reference)
            if args.igv_batch and args.reference is not None
            else None
        )
        if args.html:
            write_text(args.html, html_report or "", args.force)
            print(f"Wrote {args.html}", file=sys.stderr)
        if args.svg:
            write_text(args.svg, svg_report or "", args.force)
            print(f"Wrote {args.svg}", file=sys.stderr)
        if args.igv_batch:
            write_text(args.igv_batch, igv_report or "", args.force)
            print(f"Wrote {args.igv_batch}", file=sys.stderr)
        if args.write_diff:
            print(f"Wrote {args.write_diff}", file=sys.stderr)
    except (GenDiffError, OSError, ValueError) as error:
        print(f"gendiff: error: {error}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
    else:
        print(_render(result))
    return 0 if result.equivalent else 1


def _batch_main(argv: List[str]) -> int:
    args = _batch_parser().parse_args(argv)
    try:
        outputs = [
            output
            for output in (args.html, args.json, args.junit, args.multiqc)
            if output
        ]
        if len(set(output.absolute() for output in outputs)) != len(outputs):
            raise ValueError("batch outputs must use different paths")
        for output in outputs:
            if output.exists() and not args.force:
                raise FileExistsError(
                    f"output exists: {output}; use --force to replace it"
                )
        if args.write_diff and args.write_diff.absolute() in {
            output.absolute() for output in outputs
        }:
            raise ValueError(
                "diff directory and report outputs must use different paths"
            )
        if args.multiqc and not args.multiqc.name.endswith("_mqc.json"):
            raise ValueError("--multiqc filename must end with _mqc.json")
        if not all(pattern.strip() for pattern in args.fail_transition):
            raise ValueError("--fail-transition cannot be empty")
        policy = BatchPolicy(
            max_modified=args.max_modified,
            max_only=args.max_only,
            min_overlap=args.min_overlap,
            allow_structural_differences=args.allow_structural_differences,
            fail_transitions=tuple(args.fail_transition),
        )
        result = compare_manifest(
            args.manifest,
            threads=args.threads,
            profile=args.profile,
            reference=args.reference,
            normalize_variants=args.normalize,
            ignore_tags=args.ignore_tag,
            ignore_info=args.ignore_info,
            max_examples=args.max_examples,
            progress=args.progress,
            temp_dir=args.temp_dir,
            diff_dir=args.write_diff,
            force=args.force,
            policy=policy,
        )
        reports = {
            args.html: render_batch_html(result, args.html) if args.html else None,
            args.json: (
                json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n"
                if args.json
                else None
            ),
            args.junit: render_junit(result) if args.junit else None,
            args.multiqc: render_multiqc(result) if args.multiqc else None,
        }
        for output, content in reports.items():
            if output is not None and content is not None:
                write_text(output, content, args.force)
                print(f"Wrote {output}", file=sys.stderr)
        if args.write_diff:
            print(f"Wrote {args.write_diff}", file=sys.stderr)
    except (BatchError, GenDiffError, OSError, ValueError) as error:
        print(f"gendiff: error: {error}", file=sys.stderr)
        return 2
    print(render_batch_text(result))
    return 0 if result.passed else 1


def main(argv: Optional[List[str]] = None) -> int:
    values = list(sys.argv[1:] if argv is None else argv)
    if values and values[0] == "batch":
        return _batch_main(values[1:])
    return _single_main(values)
