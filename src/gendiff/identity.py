from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, Optional, Tuple

import pysam


@dataclass(frozen=True)
class IdentityPair:
    first: str
    second: str
    shared_calls: int
    concordance: float

    @property
    def status(self) -> str:
        if self.shared_calls < 20:
            return "insufficient calls"
        if self.concordance >= 0.98:
            return "match"
        if self.concordance >= 0.90:
            return "possible match"
        return "mismatch"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "first": self.first,
            "second": self.second,
            "shared_calls": self.shared_calls,
            "concordance": self.concordance,
            "status": self.status,
        }


@dataclass(frozen=True)
class IdentityResult:
    pairs: Tuple[IdentityPair, ...]
    best_matches: Tuple[IdentityPair, ...]

    @property
    def passed(self) -> bool:
        return bool(self.best_matches) and all(
            pair.status in {"match", "possible match"} for pair in self.best_matches
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "gendiff_identity_format": 1,
            "passed": self.passed,
            "best_matches": [pair.to_dict() for pair in self.best_matches],
            "pairs": [pair.to_dict() for pair in self.pairs],
        }


def _variant_database(path: Path, database: Path) -> Tuple[str, ...]:
    connection = sqlite3.connect(str(database))
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("CREATE TABLE calls (sample TEXT, locus TEXT, genotype TEXT)")
    batch = []
    with pysam.VariantFile(str(path)) as handle:
        samples = tuple(handle.header.samples)
        if not samples:
            raise ValueError(f"VCF has no samples: {path}")
        for record in handle:
            locus = json.dumps(
                (record.contig, record.pos, tuple(record.alleles or ())),
                separators=(",", ":"),
            )
            for sample in samples:
                genotype = record.samples[sample].get("GT")
                if genotype and all(allele is not None for allele in genotype):
                    batch.append(
                        (
                            sample,
                            locus,
                            json.dumps(tuple(sorted(genotype)), separators=(",", ":")),
                        )
                    )
                    if len(batch) >= 2000:
                        connection.executemany(
                            "INSERT INTO calls VALUES (?, ?, ?)", batch
                        )
                        batch.clear()
    if batch:
        connection.executemany("INSERT INTO calls VALUES (?, ?, ?)", batch)
    connection.execute("CREATE INDEX calls_locus ON calls(locus)")
    connection.commit()
    connection.close()
    return samples


def _variant_identity(first: Path, second: Path) -> IdentityResult:
    with TemporaryDirectory(prefix="gendiff-identity-") as work:
        workspace = Path(work)
        first_db, second_db = workspace / "first.sqlite", workspace / "second.sqlite"
        first_samples = _variant_database(first, first_db)
        second_samples = _variant_database(second, second_db)
        connection = sqlite3.connect(str(first_db))
        connection.execute("ATTACH DATABASE ? AS second", (str(second_db),))
        rows = {
            (first_name, second_name): (shared, matches)
            for first_name, second_name, shared, matches in connection.execute(
                "SELECT a.sample, b.sample, COUNT(*), "
                "SUM(CASE WHEN a.genotype=b.genotype THEN 1 ELSE 0 END) "
                "FROM calls a JOIN second.calls b ON a.locus=b.locus "
                "GROUP BY a.sample, b.sample"
            )
        }
        connection.close()
    pairs = []
    for first_name in first_samples:
        for second_name in second_samples:
            shared, matches = rows.get((first_name, second_name), (0, 0))
            pairs.append(
                IdentityPair(
                    first_name,
                    second_name,
                    shared,
                    matches / shared if shared else 0.0,
                )
            )
    return _identity_result(tuple(pairs), first_samples)


def _sites(path: Path) -> list[Tuple[str, int, str, str]]:
    result = []
    with pysam.VariantFile(str(path)) as handle:
        for record in handle:
            if (
                len(record.ref) == 1
                and record.alts
                and len(record.alts) == 1
                and len(record.alts[0]) == 1
            ):
                result.append((record.contig, record.pos, record.ref, record.alts[0]))
    if not result:
        raise ValueError(f"known-sites file contains no biallelic SNPs: {path}")
    return result


def _alignment_calls(
    path: Path,
    sites: list[Tuple[str, int, str, str]],
    reference: Optional[Path],
    min_depth: int,
) -> Dict[str, Dict[Tuple[Any, ...], Tuple[int, ...]]]:
    mode = "rc" if path.name.lower().endswith(".cram") else "rb"
    kwargs = {"reference_filename": str(reference)} if reference else {}
    calls = {}
    with pysam.AlignmentFile(str(path), mode, **kwargs) as handle:
        if not handle.has_index():
            raise ValueError(f"identity checking requires an index: {path}")
        for contig, position, reference_base, alternate in sites:
            counts = {reference_base.upper(): 0, alternate.upper(): 0}
            try:
                columns = handle.pileup(
                    contig,
                    position - 1,
                    position,
                    truncate=True,
                    min_base_quality=13,
                    min_mapping_quality=20,
                )
                for column in columns:
                    if column.reference_pos != position - 1:
                        continue
                    for read in column.pileups:
                        if read.query_position is None:
                            continue
                        base = read.alignment.query_sequence[
                            read.query_position
                        ].upper()
                        if base in counts:
                            counts[base] += 1
            except ValueError:
                continue
            depth = sum(counts.values())
            if depth < min_depth:
                continue
            fraction = counts[alternate.upper()] / depth
            genotype = (
                (0, 0) if fraction < 0.2 else (1, 1) if fraction > 0.8 else (0, 1)
            )
            calls[(contig, position, reference_base, alternate)] = genotype
    return {path.name: calls}


def _pairs(
    first: Dict[str, Dict[Tuple[Any, ...], Tuple[int, ...]]],
    second: Dict[str, Dict[Tuple[Any, ...], Tuple[int, ...]]],
) -> IdentityResult:
    pairs = []
    for first_name, first_calls in first.items():
        for second_name, second_calls in second.items():
            shared = set(first_calls) & set(second_calls)
            matches = sum(first_calls[key] == second_calls[key] for key in shared)
            pairs.append(
                IdentityPair(
                    first_name,
                    second_name,
                    len(shared),
                    matches / len(shared) if shared else 0.0,
                )
            )
    return _identity_result(tuple(pairs), tuple(first))


def _identity_result(
    pairs: Tuple[IdentityPair, ...], first_names: Tuple[str, ...]
) -> IdentityResult:
    ordered = tuple(
        sorted(
            pairs, key=lambda pair: (-pair.concordance, -pair.shared_calls, pair.second)
        )
    )
    best = []
    used_first, used_second = set(), set()
    for pair in ordered:
        if pair.first not in used_first and pair.second not in used_second:
            best.append(pair)
            used_first.add(pair.first)
            used_second.add(pair.second)
    for name in first_names:
        if name not in used_first:
            best.append(
                max(
                    (pair for pair in ordered if pair.first == name),
                    key=lambda pair: (pair.concordance, pair.shared_calls),
                )
            )
    return IdentityResult(ordered, tuple(best))


def compare_identity(
    first: Path,
    second: Path,
    *,
    sites: Optional[Path],
    reference: Optional[Path],
    min_depth: int,
) -> IdentityResult:
    names = (first.name.lower(), second.name.lower())
    if all(name.endswith((".vcf", ".vcf.gz", ".bcf")) for name in names):
        return _variant_identity(first, second)
    if all(name.endswith((".bam", ".cram")) for name in names):
        if sites is None:
            raise ValueError("BAM/CRAM identity checking requires --sites")
        variants = _sites(sites)
        return _pairs(
            _alignment_calls(first, variants, reference, min_depth),
            _alignment_calls(second, variants, reference, min_depth),
        )
    raise ValueError("identity inputs must both be variant files or both alignments")


def render_identity(result: IdentityResult) -> str:
    lines = [f"Identity check: {'pass' if result.passed else 'warning'}"]
    for pair in result.best_matches:
        lines.append(
            f"  {pair.first} → {pair.second}: {pair.concordance:.1%} "
            f"across {pair.shared_calls:,} calls ({pair.status})"
        )
    return "\n".join(lines)
