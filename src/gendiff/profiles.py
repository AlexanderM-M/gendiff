"""Named comparison policies."""

from __future__ import annotations

from typing import Dict, Tuple

ALIGNMENT_FIELDS = (
    "read_name",
    "flags",
    "reference",
    "position",
    "mapping_quality",
    "cigar",
    "mate",
    "template_length",
    "sequence",
    "qualities",
    "tags",
)

VARIANT_FIELDS = (
    "locus",
    "identifiers",
    "alleles",
    "quality",
    "filters",
    "info",
    "samples",
)

ALIGNMENT_PROFILES: Dict[str, Tuple[str, ...]] = {
    "strict": ALIGNMENT_FIELDS,
    "core": tuple(field for field in ALIGNMENT_FIELDS if field != "tags"),
    "mapping": tuple(
        field
        for field in ALIGNMENT_FIELDS
        if field not in {"sequence", "qualities", "tags"}
    ),
}

VARIANT_PROFILES: Dict[str, Tuple[str, ...]] = {
    "strict": VARIANT_FIELDS,
    "calls": ("locus", "alleles", "filters", "samples"),
    "genotypes": ("locus", "alleles", "samples"),
}

ALL_PROFILES = tuple(sorted(set(ALIGNMENT_PROFILES) | set(VARIANT_PROFILES)))


def fields_for(kind: str, profile: str) -> Tuple[str, ...]:
    profiles = ALIGNMENT_PROFILES if kind == "alignment" else VARIANT_PROFILES
    if profile not in profiles:
        valid = ", ".join(profiles)
        raise ValueError(
            f"profile {profile!r} is not valid for {kind} files; use {valid}"
        )
    return profiles[profile]
