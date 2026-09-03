from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List, Optional

from semantiseq import __version__
from semantiseq.baseline import (
    BaselineError,
    baseline_report,
    check_baseline,
    create_baseline,
)
from semantiseq.batch import BatchError, BatchPolicy, compare_manifest
from semantiseq.batch_report import (
    render_batch_html,
    render_batch_text,
    render_junit,
    render_multiqc,
)
from semantiseq.compare import SemantiSeqError, compare_files
from semantiseq.identity import compare_identity, render_identity
from semantiseq.manifest import ManifestError, pair_directories, render_manifest
from semantiseq.matrix import compare_matrix, render_matrix_html, render_matrix_text
from semantiseq.model import ComparisonResult
from semantiseq.policy import PolicyError, load_policy
from semantiseq.profiles import ALL_PROFILES
from semantiseq.provenance import reproducibility_manifest
from semantiseq.report import render_html, render_igv_batch, render_svg, write_text


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
        prog="semantiseq",
        description="Compare the logical contents of two genomics files.",
        epilog=(
            "Subcommands: batch, baseline, manifest, matrix, and identity. "
            "Run 'semantiseq COMMAND --help' for details."
        ),
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
    parser.add_argument(
        "--cache",
        type=Path,
        metavar="DIR",
        help="reuse validated scan and matching data",
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
    parser.add_argument(
        "--tracks",
        type=Path,
        metavar="DIR",
        help="write BED, bedGraph, and per-contig difference tracks",
    )
    parser.add_argument(
        "--difference-table",
        type=Path,
        metavar="FILE",
        help="write every difference as TSV or TSV.GZ",
    )
    parser.add_argument(
        "--region",
        action="append",
        default=[],
        metavar="CONTIG:START-END",
        help="compare only this region; repeat as needed",
    )
    parser.add_argument(
        "--regions", type=Path, metavar="BED", help="compare only BED regions"
    )
    parser.add_argument(
        "--exclude-regions", type=Path, metavar="BED", help="exclude BED regions"
    )
    parser.add_argument("--force", action="store_true", help="replace outputs")
    parser.add_argument("--json", action="store_true", help="write JSON output")
    parser.add_argument(
        "--reproducibility",
        type=Path,
        metavar="FILE",
        help="write checksums and settings for exact reproduction",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return parser


def _batch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="semantiseq batch",
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
        "--max-modified-fraction",
        type=_proportion,
        help="maximum modified fraction per comparison",
    )
    parser.add_argument(
        "--max-only",
        type=_nonnegative_int,
        default=0,
        help="maximum combined sample-only records per comparison (default: 0)",
    )
    parser.add_argument(
        "--max-only-fraction",
        type=_proportion,
        help="maximum sample-only fraction per comparison",
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
        "--policy",
        type=Path,
        help="stage-aware JSON policy file",
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
    parser.add_argument("--cache", type=Path, metavar="DIR", help="shared scan cache")
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
    parser.add_argument(
        "--tracks",
        type=Path,
        metavar="DIR",
        help="write BED, bedGraph, and per-contig tracks for every comparison",
    )
    parser.add_argument("--force", action="store_true", help="replace outputs")
    return parser


def _manifest_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="semantiseq manifest",
        description="Pair matching genomics files in two directory trees.",
    )
    parser.add_argument("baseline", type=Path, help="baseline output directory")
    parser.add_argument("candidate", type=Path, help="candidate output directory")
    parser.add_argument("--output", type=Path, required=True, help="output TSV")
    parser.add_argument(
        "--dry-run", action="store_true", help="print the manifest without writing it"
    )
    parser.add_argument("--force", action="store_true", help="replace the output")
    return parser


def _matrix_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="semantiseq matrix",
        description="Compare many genomics files in a similarity matrix.",
    )
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument("--threads", type=_positive_int, default=2)
    parser.add_argument("--profile", choices=ALL_PROFILES, default="strict")
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--normalize", action="store_true")
    parser.add_argument("--ignore-tag", action="append", default=[])
    parser.add_argument("--ignore-info", action="append", default=[])
    parser.add_argument("--cache", type=Path, metavar="DIR")
    parser.add_argument("--html", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def _identity_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="semantiseq identity",
        description="Check whether two files represent the same biological sample.",
    )
    parser.add_argument("sample_a", type=Path)
    parser.add_argument("sample_b", type=Path)
    parser.add_argument("--sites", type=Path, help="known SNPs for BAM/CRAM inputs")
    parser.add_argument("--reference", type=Path, help="reference FASTA for CRAM")
    parser.add_argument("--min-depth", type=_positive_int, default=5)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--force", action="store_true")
    return parser


def _baseline_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="semantiseq baseline",
        description="Create or check compact semantic baselines.",
    )
    commands = parser.add_subparsers(dest="baseline_command", required=True)
    create = commands.add_parser("create", help="create a baseline from before files")
    check = commands.add_parser("check", help="check after files against a baseline")
    create.add_argument("manifest", type=Path)
    check.add_argument("baseline", type=Path)
    check.add_argument("manifest", type=Path)
    for command in (create, check):
        command.add_argument("--threads", type=_positive_int, default=2)
        command.add_argument("--profile", choices=ALL_PROFILES, default="strict")
        command.add_argument("--reference", type=Path)
        command.add_argument("--normalize", action="store_true")
        command.add_argument("--progress", action="store_true")
    create.add_argument("--output", type=Path, required=True)
    create.add_argument("--ignore-tag", action="append", default=[])
    create.add_argument("--ignore-info", action="append", default=[])
    create.add_argument("--force", action="store_true")
    check.add_argument("--json", type=Path, help="write machine-readable results")
    check.add_argument("--force", action="store_true")
    return parser


def _status(value: bool) -> str:
    return "same" if value else "different"


def _render(result: ComparisonResult) -> str:
    lines = [
        f"Equivalent: {'yes' if result.equivalent else 'no'}",
        f"Type: {result.kind}",
        f"Relationship: {result.relationship}",
        f"Identity overlap: {result.identity_overlap:.1%}",
        f"Content similarity: {result.content_similarity:.1%}",
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
    if result.header_differences:
        lines.append("Header differences:")
        lines.extend(f"  - {value}" for value in result.header_differences[:10])
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
        if details.findings:
            lines.append("What changed:")
            lines.extend(f"  - {finding}" for finding in details.findings)
        if details.group_changes:
            groups = ", ".join(
                f"{name}={count:,}"
                for name, count in list(details.group_changes.items())[:5]
            )
            lines.append(f"Read-group attribution: {groups}")
    if "manifest" in result.artifacts:
        lines.append(f"Diff files: {Path(result.artifacts['manifest']).parent}")
    if "track_manifest" in result.artifacts:
        lines.append(
            f"Genomic tracks: {Path(result.artifacts['track_manifest']).parent}"
        )
    if result.cache:
        lines.append(
            f"Cache: {result.cache.get('first')} / {result.cache.get('second')}"
        )
    return "\n".join(lines)


def _single_main(argv: List[str]) -> int:
    started = time.perf_counter()
    args = _single_parser().parse_args(argv)
    try:
        outputs = [
            output
            for output in (
                args.html,
                args.svg,
                args.igv_batch,
                args.difference_table,
                args.reproducibility,
            )
            if output
        ]
        for output in outputs:
            if output and output.exists() and not args.force:
                raise FileExistsError(
                    f"output exists: {output}; use --force to replace it"
                )
        if len(set(outputs)) != len(outputs):
            raise ValueError("report outputs must use different paths")
        if any(
            output.absolute() in {args.sample_a.absolute(), args.sample_b.absolute()}
            for output in outputs
        ):
            raise ValueError("report output cannot replace an input file")
        directories = [path for path in (args.write_diff, args.tracks) if path]
        if len({path.absolute() for path in directories}) != len(directories):
            raise ValueError("diff and track output directories must differ")
        if args.cache and any(
            args.cache.absolute() == directory.absolute() for directory in directories
        ):
            raise ValueError("cache and artifact directories must differ")
        if any(
            directory.absolute() in {output.absolute() for output in outputs}
            for directory in directories
        ):
            raise ValueError("artifact directories and report outputs must differ")
        if args.igv_batch and args.reference is None:
            raise ValueError("--igv-batch requires --reference")
        explain = any(
            (
                args.explain,
                args.html is not None,
                args.svg is not None,
                args.igv_batch,
                args.write_diff is not None,
                args.tracks is not None,
                args.difference_table is not None,
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
            track_dir=args.tracks,
            difference_table=args.difference_table,
            regions=args.region,
            regions_file=args.regions,
            exclude_regions=args.exclude_regions,
            cache_dir=args.cache,
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
        if args.tracks:
            print(f"Wrote {args.tracks}", file=sys.stderr)
        if args.difference_table:
            print(f"Wrote {args.difference_table}", file=sys.stderr)
        if args.reproducibility:
            manifest = reproducibility_manifest(
                result,
                argv,
                args.reference,
                time.perf_counter() - started,
            )
            write_text(
                args.reproducibility,
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                args.force,
            )
            print(f"Wrote {args.reproducibility}", file=sys.stderr)
    except (SemantiSeqError, OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
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
        directories = [path for path in (args.write_diff, args.tracks) if path]
        if len({path.absolute() for path in directories}) != len(directories):
            raise ValueError("diff and track output directories must differ")
        if any(
            directory.absolute() in {output.absolute() for output in outputs}
            for directory in directories
        ):
            raise ValueError("artifact directories and report outputs must differ")
        if args.multiqc and not args.multiqc.name.endswith("_mqc.json"):
            raise ValueError("--multiqc filename must end with _mqc.json")
        if not all(pattern.strip() for pattern in args.fail_transition):
            raise ValueError("--fail-transition cannot be empty")
        policy = BatchPolicy(
            max_modified=args.max_modified,
            max_only=args.max_only,
            max_modified_fraction=args.max_modified_fraction,
            max_only_fraction=args.max_only_fraction,
            min_overlap=args.min_overlap,
            allow_structural_differences=args.allow_structural_differences,
            fail_transitions=tuple(args.fail_transition),
        )
        policy_document = load_policy(args.policy) if args.policy else None
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
            policy_document=policy_document,
            track_dir=args.tracks,
            cache_dir=args.cache,
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
        if args.tracks:
            print(f"Wrote {args.tracks}", file=sys.stderr)
    except (BatchError, SemantiSeqError, PolicyError, OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
        return 2
    print(render_batch_text(result))
    return 0 if result.passed else 1


def _manifest_main(argv: List[str]) -> int:
    args = _manifest_parser().parse_args(argv)
    try:
        pairs = pair_directories(args.baseline, args.candidate)
        content = render_manifest(pairs, args.output)
        if args.dry_run:
            print(content, end="")
        else:
            write_text(args.output, content, args.force)
            print(f"Wrote {args.output} ({len(pairs)} comparisons)")
    except (ManifestError, OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
        return 2
    return 0


def _matrix_main(argv: List[str]) -> int:
    args = _matrix_parser().parse_args(argv)
    try:
        outputs = [value for value in (args.html, args.json) if value]
        if len({value.absolute() for value in outputs}) != len(outputs):
            raise ValueError("matrix outputs must use different paths")
        for output in outputs:
            if output.exists() and not args.force:
                raise FileExistsError(
                    f"output exists: {output}; use --force to replace it"
                )
            if output.absolute() in {path.absolute() for path in args.files}:
                raise ValueError("matrix output cannot replace an input file")
        result = compare_matrix(
            args.files,
            threads=args.threads,
            profile=args.profile,
            reference=args.reference,
            normalize_variants=args.normalize,
            ignore_tags=args.ignore_tag,
            ignore_info=args.ignore_info,
            cache_dir=args.cache,
        )
        if args.html:
            write_text(args.html, render_matrix_html(result), args.force)
        if args.json:
            write_text(
                args.json,
                json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
                args.force,
            )
        print(render_matrix_text(result))
        return 0
    except (SemantiSeqError, OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
        return 2


def _identity_main(argv: List[str]) -> int:
    args = _identity_parser().parse_args(argv)
    try:
        for path in (args.sample_a, args.sample_b):
            if not path.is_file():
                raise ValueError(f"file not found: {path}")
        if args.sites and not args.sites.is_file():
            raise ValueError(f"known-sites file not found: {args.sites}")
        if args.reference and not args.reference.is_file():
            raise ValueError(f"reference not found: {args.reference}")
        if args.json and args.json.exists() and not args.force:
            raise FileExistsError(
                f"output exists: {args.json}; use --force to replace it"
            )
        if args.json and args.json.absolute() in {
            args.sample_a.absolute(),
            args.sample_b.absolute(),
        }:
            raise ValueError("identity output cannot replace an input file")
        result = compare_identity(
            args.sample_a,
            args.sample_b,
            sites=args.sites,
            reference=args.reference,
            min_depth=args.min_depth,
        )
        if args.json:
            write_text(
                args.json,
                json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n",
                args.force,
            )
        print(render_identity(result))
        return 0 if result.passed else 1
    except (OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
        return 2


def _baseline_main(argv: List[str]) -> int:
    args = _baseline_parser().parse_args(argv)
    try:
        if args.baseline_command == "create":
            payload = create_baseline(
                args.manifest,
                threads=args.threads,
                profile=args.profile,
                reference=args.reference,
                normalize_variants=args.normalize,
                ignore_tags=args.ignore_tag,
                ignore_info=args.ignore_info,
                progress=args.progress,
            )
            write_text(
                args.output,
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                args.force,
            )
            print(f"Wrote {args.output} ({len(payload['comparisons'])} comparisons)")
            return 0
        checks = check_baseline(
            args.baseline,
            args.manifest,
            threads=args.threads,
            profile=args.profile,
            reference=args.reference,
            normalize_variants=args.normalize,
            progress=args.progress,
        )
        report = baseline_report(checks)
        if args.json:
            write_text(
                args.json,
                json.dumps(report, indent=2, sort_keys=True) + "\n",
                args.force,
            )
        print(
            f"Passed: {'yes' if report['passed'] else 'no'}\n"
            f"Comparisons: {report['summary']['comparisons']} "
            f"({report['summary']['passed']} passed, "
            f"{report['summary']['failed']} failed)"
        )
        for item in checks:
            fields = (
                f"; fields: {', '.join(item.changed_fields)}"
                if item.changed_fields
                else ""
            )
            print(
                f"  {'PASS' if item.passed else 'FAIL'}  {item.sample} / "
                f"{item.stage}{fields}"
            )
        return 0 if report["passed"] else 1
    except (BaselineError, BatchError, SemantiSeqError, OSError, ValueError) as error:
        print(f"semantiseq: error: {error}", file=sys.stderr)
        return 2


def main(argv: Optional[List[str]] = None) -> int:
    values = list(sys.argv[1:] if argv is None else argv)
    if values and values[0] == "batch":
        return _batch_main(values[1:])
    if values and values[0] == "manifest":
        return _manifest_main(values[1:])
    if values and values[0] == "matrix":
        return _matrix_main(values[1:])
    if values and values[0] == "identity":
        return _identity_main(values[1:])
    if values and values[0] == "baseline":
        return _baseline_main(values[1:])
    return _single_main(values)
