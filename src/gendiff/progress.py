from __future__ import annotations

import sys
import time


class ProgressReporter:
    def __init__(self, label: str, enabled: bool, interval: float = 2.0) -> None:
        self._label = label
        self._enabled = enabled
        self._interval = interval
        self._started = time.monotonic()
        self._last = self._started

    def update(self, records: int) -> None:
        if not self._enabled:
            return
        now = time.monotonic()
        if now - self._last < self._interval:
            return
        elapsed = max(now - self._started, 0.001)
        rate = records / elapsed
        print(
            f"[{self._label}] {records:,} records ({rate:,.0f} records/s)",
            file=sys.stderr,
            flush=True,
        )
        self._last = now

    def finish(self, records: int) -> None:
        if not self._enabled:
            return
        elapsed = max(time.monotonic() - self._started, 0.001)
        print(
            f"[{self._label}] {records:,} records in {elapsed:.1f}s",
            file=sys.stderr,
            flush=True,
        )
