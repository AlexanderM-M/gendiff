from __future__ import annotations

import csv
import gzip
import json
import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, Optional, Sequence, TextIO


class DifferenceTableWriter:
    """Stream record differences to a compact TSV ledger."""

    def __init__(self, path: Path) -> None:
        self._handle: TextIO = (
            gzip.open(path, "wt", encoding="utf-8", newline="")
            if path.suffix == ".gz"
            else path.open("w", encoding="utf-8", newline="")
        )
        self._writer = csv.writer(self._handle, delimiter="\t", lineterminator="\n")
        self._writer.writerow(
            (
                "status",
                "identity",
                "contig",
                "position",
                "changed_fields",
                "first",
                "second",
            )
        )

    def add(
        self,
        status: str,
        identity: Any,
        changed_fields: Sequence[str],
        first: Optional[Dict[str, Any]],
        second: Optional[Dict[str, Any]],
    ) -> None:
        record = second or first or {}
        self._writer.writerow(
            (
                status,
                json.dumps(identity, separators=(",", ":")),
                record.get("reference", record.get("contig", "")),
                record.get("position", ""),
                ",".join(changed_fields),
                json.dumps(first, separators=(",", ":"), sort_keys=True)
                if first
                else "",
                json.dumps(second, separators=(",", ":"), sort_keys=True)
                if second
                else "",
            )
        )

    def close(self) -> None:
        self._handle.close()


def publish_table(source: Path, target: Path, force: bool) -> str:
    if target.exists() and not force:
        raise FileExistsError(f"output exists: {target}; use --force to replace it")
    if not target.parent.is_dir():
        raise FileNotFoundError(f"output directory not found: {target.parent}")
    with NamedTemporaryFile(
        dir=target.parent, prefix=f".{target.name}.", delete=False
    ) as handle:
        temporary = Path(handle.name)
        with source.open("rb") as source_handle:
            shutil.copyfileobj(source_handle, handle)
    try:
        os.replace(temporary, target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return str(target.absolute())
