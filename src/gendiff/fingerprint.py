"""Order-independent fingerprints for streams of logical records."""

from __future__ import annotations

from array import array
from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any


_BITS = 256
_MODULUS = 1 << _BITS


def normalize(value: Any) -> Any:
    """Convert pysam and Python values to a stable JSON-compatible form."""
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return {"float": "nan"}
        if math.isinf(value):
            return {"float": "inf" if value > 0 else "-inf"}
        return value
    if isinstance(value, bytes):
        return {"bytes": value.hex()}
    if isinstance(value, (list, tuple, array)):
        return [normalize(item) for item in value]
    if isinstance(value, dict):
        return {
            str(key): normalize(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    return str(value)


def digest(value: Any) -> int:
    payload = json.dumps(
        normalize(value), ensure_ascii=True, separators=(",", ":"), sort_keys=True
    ).encode("utf-8")
    return int.from_bytes(
        hashlib.blake2b(payload, digest_size=32, person=b"gendiff-v1").digest(),
        "big",
    )


def digest_parts(parts: list[int]) -> int:
    """Hash already-computed field digests into one logical record digest."""
    hasher = hashlib.blake2b(digest_size=32, person=b"gendiff-rec-v1")
    for part in parts:
        hasher.update(part.to_bytes(32, "big"))
    return int.from_bytes(hasher.digest(), "big")


@dataclass(frozen=True)
class Fingerprint:
    count: int
    total: int
    squares: int
    xor: int


class FingerprintBuilder:
    """Build a bounded-memory multiset fingerprint."""

    def __init__(self) -> None:
        self._count = 0
        self._total = 0
        self._squares = 0
        self._xor = 0

    def add(self, value: Any) -> None:
        self.add_digest(digest(value))

    def add_digest(self, item: int) -> None:
        self._count += 1
        self._total = (self._total + item) % _MODULUS
        self._squares = (self._squares + item * item) % _MODULUS
        self._xor ^= item

    def finish(self) -> Fingerprint:
        return Fingerprint(self._count, self._total, self._squares, self._xor)
