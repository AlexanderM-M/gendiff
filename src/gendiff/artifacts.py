"""Safe publication helpers for generated comparison artifacts."""

from __future__ import annotations

import json
import os
import re
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Iterator, Tuple


def sample_slugs(left_label: str, right_label: str) -> Tuple[str, str]:
    """Return short, filesystem-safe, distinct sample names."""

    def slug(value: str) -> str:
        cleaned = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
        return cleaned or "sample"

    left = slug(left_label)
    right = slug(right_label)
    if left == right:
        return f"{left}-a", f"{right}-b"
    return left, right


@contextmanager
def output_workspace(target: Path, force: bool) -> Iterator[Path]:
    """Build an output directory privately, then publish it as one operation."""
    target = target.absolute()
    if target == Path(target.anchor) or target == Path.home():
        raise ValueError(f"unsafe diff output directory: {target}")
    if not target.parent.is_dir():
        raise FileNotFoundError(f"output directory not found: {target.parent}")
    if target.exists() and not force:
        raise FileExistsError(f"output exists: {target}; use --force to replace it")
    if target.exists() and force:
        manifest = target / "manifest.json"
        try:
            previous = json.loads(manifest.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, TypeError):
            raise ValueError(
                f"refusing to replace unmanaged directory: {target}"
            ) from None
        if previous.get("gendiff_diff_format") != 1:
            raise ValueError(f"refusing to replace unmanaged directory: {target}")

    with TemporaryDirectory(prefix=f".{target.name}.", dir=target.parent) as root:
        root_path = Path(root)
        workspace = root_path / "new"
        workspace.mkdir()
        yield workspace

        previous = root_path / "previous"
        if target.exists():
            os.replace(target, previous)
        try:
            os.replace(workspace, target)
        except Exception:
            if previous.exists():
                os.replace(previous, target)
            raise


def write_manifest(path: Path, payload: Dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
