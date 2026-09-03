"""Configuration hook for the optional MultiQC module."""

from __future__ import annotations


def configure() -> None:
    from multiqc import config

    config.sp["semantiseq"] = {
        "fn": "*semantiseq-batch.json",
        "max_filesize": 10_000_000,
    }
    names = [
        item if isinstance(item, str) else next(iter(item), "")
        for item in config.module_order
    ]
    if "semantiseq" not in names:
        config.module_order.append("semantiseq")
