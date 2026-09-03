from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any, Dict, Optional

from semantiseq import __version__
from semantiseq.fingerprint import Fingerprint, Sketch


def cache_key(path: Path, settings: Dict[str, Any]) -> str:
    stat = path.stat()
    payload = {
        "version": __version__,
        "path": str(path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "settings": settings,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encode_scan(scan: Any) -> Dict[str, Any]:
    return {
        "format": 1,
        "fingerprints": {
            name: [value.count, value.total, value.squares, value.xor]
            for name, value in scan.fingerprints.items()
        },
        "count": scan.count,
        "structural": scan.structural,
        "metadata": scan.metadata,
        "sketch": {"values": list(scan.sketch.values), "size": scan.sketch.size},
        "content_sketch": {
            "values": list(scan.content_sketch.values),
            "size": scan.content_sketch.size,
        },
        "contig_counts": scan.contig_counts,
        "window_counts": [
            [contig, window, count]
            for (contig, window), count in scan.window_counts.items()
        ],
        "metric_counts": {
            name: [[value, count] for value, count in counts.items()]
            for name, counts in scan.metric_counts.items()
        },
    }


def decode_scan(payload: Dict[str, Any]) -> Dict[str, Any]:
    if payload.get("format") != 1:
        raise ValueError("unsupported cache format")
    return {
        "fingerprints": {
            name: Fingerprint(*values)
            for name, values in payload["fingerprints"].items()
        },
        "count": payload["count"],
        "structural": payload["structural"],
        "metadata": payload["metadata"],
        "sketch": Sketch(tuple(payload["sketch"]["values"]), payload["sketch"]["size"]),
        "content_sketch": Sketch(
            tuple(payload["content_sketch"]["values"]),
            payload["content_sketch"]["size"],
        ),
        "contig_counts": payload["contig_counts"],
        "window_counts": {
            (contig, int(window)): count
            for contig, window, count in payload["window_counts"]
        },
        "metric_counts": {
            name: {float(value): count for value, count in values}
            for name, values in payload["metric_counts"].items()
        },
    }


def load_entry(
    cache_dir: Path, key: str, details_target: Optional[Path]
) -> Optional[Dict[str, Any]]:
    summary = cache_dir / f"{key}.json"
    details = cache_dir / f"{key}.sqlite"
    if not summary.is_file() or (details_target is not None and not details.is_file()):
        return None
    try:
        payload = json.loads(summary.read_text(encoding="utf-8"))
        decoded = decode_scan(payload)
        if details_target is not None:
            shutil.copyfile(details, details_target)
        return decoded
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def store_entry(
    cache_dir: Path,
    key: str,
    scan: Any,
    details_source: Optional[Path],
) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    summary = cache_dir / f"{key}.json"
    with NamedTemporaryFile(
        "w", encoding="utf-8", dir=cache_dir, prefix=".semantiseq-", delete=False
    ) as handle:
        temporary = Path(handle.name)
        json.dump(encode_scan(scan), handle, separators=(",", ":"), sort_keys=True)
    os.replace(temporary, summary)
    if details_source is not None:
        target = cache_dir / f"{key}.sqlite"
        with NamedTemporaryFile(
            dir=cache_dir, prefix=".semantiseq-", delete=False
        ) as handle:
            temporary_details = Path(handle.name)
        try:
            shutil.copyfile(details_source, temporary_details)
            os.replace(temporary_details, target)
        except Exception:
            temporary_details.unlink(missing_ok=True)
            raise
