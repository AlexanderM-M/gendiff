from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional, Sequence

from semantiseq import __version__
from semantiseq.model import ComparisonResult


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def reproducibility_manifest(
    result: ComparisonResult,
    argv: Sequence[str],
    reference: Optional[Path],
    duration_seconds: float,
) -> Dict[str, Any]:
    inputs = [
        {"path": str(path.resolve()), "sha256": _sha256(path)}
        for path in (result.left, result.right)
    ]
    return {
        "semantiseq_reproducibility_format": 1,
        "semantiseq_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": ["semantiseq", *argv],
        "duration_seconds": round(duration_seconds, 6),
        "inputs": inputs,
        "reference": (
            {"path": str(reference.resolve()), "sha256": _sha256(reference)}
            if reference
            else None
        ),
        "comparison": result.to_dict(),
    }
