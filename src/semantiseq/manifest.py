"""Create deterministic batch manifests from two matching directory trees."""

from __future__ import annotations

import csv
import io
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple


class ManifestError(ValueError):
    """Raised when directory pairing is incomplete or ambiguous."""


_SUFFIXES = (".vcf.gz", ".cram", ".bam", ".bcf", ".vcf")
_PROCESSING_SUFFIXES = (
    ".merged",
    ".sorted",
    ".deduplicated",
    ".dedup",
    ".markdup",
    ".aligned",
)


@dataclass(frozen=True)
class ManifestPair:
    sample: str
    stage: str
    before: Path
    after: Path
    relative: Path


def _strip_suffix(path: Path) -> str:
    value = path.name
    lowered = value.lower()
    for suffix in _SUFFIXES:
        if lowered.endswith(suffix):
            value = value[: -len(suffix)]
            break
    lowered = value.lower()
    for suffix in _PROCESSING_SUFFIXES:
        if lowered.endswith(suffix):
            value = value[: -len(suffix)]
            break
    return value or path.name


def _files(root: Path) -> Dict[str, Path]:
    if not root.is_dir():
        raise ManifestError(f"directory not found: {root}")
    found = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file() or not path.name.lower().endswith(_SUFFIXES):
            continue
        relative = path.relative_to(root).as_posix()
        found[relative] = path.absolute()
    if not found:
        raise ManifestError(f"no BAM, CRAM, VCF, or BCF files found in {root}")
    return found


def pair_directories(before: Path, after: Path) -> Tuple[ManifestPair, ...]:
    baseline = _files(before)
    candidate = _files(after)
    missing_after = sorted(set(baseline) - set(candidate))
    missing_before = sorted(set(candidate) - set(baseline))
    if missing_after or missing_before:
        parts = []
        if missing_after:
            parts.append("missing from candidate: " + ", ".join(missing_after[:10]))
        if missing_before:
            parts.append("missing from baseline: " + ", ".join(missing_before[:10]))
        raise ManifestError("directory trees do not match; " + "; ".join(parts))
    pairs = []
    identities = set()
    for relative_text in sorted(baseline):
        relative = Path(relative_text)
        sample = _strip_suffix(relative)
        stage = relative.parent.as_posix() if relative.parent != Path(".") else "output"
        identity = (sample, stage)
        if identity in identities:
            raise ManifestError(
                f"ambiguous sample and stage inferred for {relative}; edit the layout "
                "or create the manifest manually"
            )
        identities.add(identity)
        pairs.append(
            ManifestPair(
                sample=sample,
                stage=stage,
                before=baseline[relative_text],
                after=candidate[relative_text],
                relative=relative,
            )
        )
    return tuple(pairs)


def render_manifest(pairs: Tuple[ManifestPair, ...], output: Path) -> str:
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(("sample", "stage", "before", "after"))
    for pair in pairs:
        writer.writerow(
            (
                pair.sample,
                pair.stage,
                os.path.relpath(pair.before, output.absolute().parent),
                os.path.relpath(pair.after, output.absolute().parent),
            )
        )
    return stream.getvalue()
