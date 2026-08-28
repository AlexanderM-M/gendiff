from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import List, Optional

from gendiff import __version__
from gendiff.compare import GenDiffError, compare_files
from gendiff.model import ComparisonResult


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="gendiff",
        description="Compare the logical contents of two genomics files.",
    )
    parser.add_argument("left", type=Path, help="first BAM, CRAM, VCF, or BCF")
    parser.add_argument("right", type=Path, help="second file")
    parser.add_argument(
        "--reference", type=Path, help="reference FASTA used to decode CRAM"
    )
    parser.add_argument("--json", action="store_true", help="write JSON output")
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    return parser


def _status(value: bool) -> str:
    return "same" if value else "different"


def _render(result: ComparisonResult) -> str:
    lines = [
        f"Equivalent: {'yes' if result.equivalent else 'no'}",
        f"Type: {result.kind}",
        f"Records: {result.left_records} / {result.right_records}",
        f"Logical records: {_status(result.content_equal)}",
        f"Structural header: {_status(result.structure_equal)}",
        f"Metadata header: {_status(result.metadata_equal)} (informational)",
    ]
    if result.changed_fields:
        fields = ", ".join(field.replace("_", " ") for field in result.changed_fields)
        lines.append(f"Changed fields: {fields}")
    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        result = compare_files(args.left, args.right, args.reference)
    except (GenDiffError, OSError, ValueError) as error:
        print(f"gendiff: error: {error}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(result.to_dict(), indent=2, sort_keys=True))
    else:
        print(_render(result))
    return 0 if result.equivalent else 1
