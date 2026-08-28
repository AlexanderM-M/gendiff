"""Order-independent fingerprints for streams of logical records."""

from __future__ import annotations

import hashlib
import heapq
import math
import pickle
from array import array
from dataclasses import dataclass
from typing import Any, Iterable

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
    """Hash a canonical value; suitable for headers and public helpers."""
    payload = pickle.dumps(normalize(value), protocol=4)
    return int.from_bytes(
        hashlib.blake2b(payload, digest_size=32, person=b"gendiff-v2").digest(),
        "big",
    )


def digest_native(value: Any) -> int:
    """Hash an internally constructed value without redundant normalization."""
    payload = pickle.dumps(value, protocol=5)
    return int.from_bytes(
        hashlib.blake2b(payload, digest_size=32, person=b"gendiff-v2").digest(),
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

    def merge(self, fingerprint: Fingerprint) -> None:
        self._count += fingerprint.count
        self._total = (self._total + fingerprint.total) % _MODULUS
        self._squares = (self._squares + fingerprint.squares) % _MODULUS
        self._xor ^= fingerprint.xor

    def finish(self) -> Fingerprint:
        return Fingerprint(self._count, self._total, self._squares, self._xor)


@dataclass(frozen=True)
class Sketch:
    values: tuple[int, ...]
    size: int = 512

    @property
    def cardinality(self) -> float:
        if len(self.values) < self.size:
            return float(len(self.values))
        return (self.size - 1) * _MODULUS / self.values[-1]


class SketchBuilder:
    """Build a fixed-size KMV sketch of unique identities."""

    def __init__(self, size: int = 512) -> None:
        self._size = size
        self._heap: list[int] = []
        self._values: set[int] = set()

    def add(self, value: Any) -> None:
        self.add_digest(digest(value))

    def add_digest(self, value: int) -> None:
        if value in self._values:
            return
        if len(self._heap) < self._size:
            heapq.heappush(self._heap, -value)
            self._values.add(value)
            return
        largest = -self._heap[0]
        if value >= largest:
            return
        removed = -heapq.heapreplace(self._heap, -value)
        self._values.remove(removed)
        self._values.add(value)

    def merge(self, sketch: Sketch) -> None:
        for value in sketch.values:
            self.add_digest(value)

    def finish(self) -> Sketch:
        return Sketch(tuple(sorted(self._values)), self._size)


def sketch_containment(left: Sketch, right: Sketch) -> float:
    """Estimate how much of the smaller identity set occurs in the larger."""
    if not left.values and not right.values:
        return 1.0
    if not left.values or not right.values:
        return 0.0
    sample = sorted(set(left.values) | set(right.values))[: left.size]
    shared = sum(value in left.values and value in right.values for value in sample)
    jaccard = shared / len(sample)
    if jaccard == 0:
        return 0.0
    intersection = jaccard * (left.cardinality + right.cardinality) / (1.0 + jaccard)
    return min(1.0, intersection / min(left.cardinality, right.cardinality))


def merge_sketches(sketches: Iterable[Sketch]) -> Sketch:
    builder = SketchBuilder()
    for sketch in sketches:
        builder.merge(sketch)
    return builder.finish()
