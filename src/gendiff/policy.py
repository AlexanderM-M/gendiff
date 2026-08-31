"""Stage-aware, auditable policy evaluation for batch comparisons."""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from fnmatch import fnmatchcase
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

from gendiff.model import ComparisonResult


class PolicyError(ValueError):
    """Raised when a policy file is invalid."""


@dataclass(frozen=True)
class BatchPolicy:
    max_modified: Optional[int] = 0
    max_only: Optional[int] = 0
    max_modified_fraction: Optional[float] = None
    max_only_fraction: Optional[float] = None
    min_overlap: float = 0.0
    allow_structural_differences: bool = False
    fail_transitions: Tuple[str, ...] = ()
    deny_transitions: Tuple[str, ...] = ()
    allow_fields: Optional[Tuple[str, ...]] = None
    deny_fields: Tuple[str, ...] = ()
    allow_transitions: Tuple[str, ...] = ()
    transition_default: str = "allow"
    severity: str = "error"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "max_modified": self.max_modified,
            "max_only": self.max_only,
            "max_modified_fraction": self.max_modified_fraction,
            "max_only_fraction": self.max_only_fraction,
            "min_overlap": self.min_overlap,
            "allow_structural_differences": self.allow_structural_differences,
            "fail_transitions": list(self.fail_transitions),
            "deny_transitions": list(self.deny_transitions),
            "allow_fields": (
                list(self.allow_fields) if self.allow_fields is not None else None
            ),
            "deny_fields": list(self.deny_fields),
            "allow_transitions": list(self.allow_transitions),
            "transition_default": self.transition_default,
            "severity": self.severity,
        }


@dataclass(frozen=True)
class PolicyRule:
    name: str
    match: Mapping[str, str]
    settings: Mapping[str, Any]


@dataclass(frozen=True)
class PolicyDocument:
    source: Path
    defaults: Mapping[str, Any]
    rules: Tuple[PolicyRule, ...]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "source": str(self.source),
            "defaults": dict(self.defaults),
            "rules": [
                {
                    "name": rule.name,
                    "match": dict(rule.match),
                    "set": dict(rule.settings),
                }
                for rule in self.rules
            ],
        }


@dataclass(frozen=True)
class PolicyDecision:
    errors: Tuple[str, ...]
    warnings: Tuple[str, ...]
    trace: Tuple[str, ...]
    matched_rules: Tuple[str, ...]
    effective: BatchPolicy


_SETTING_KEYS = set(BatchPolicy.__dataclass_fields__)
_MATCH_KEYS = {"sample", "stage", "kind", "profile"}
_SEVERITIES = {"error", "warning", "info"}
_TRANSITION_DEFAULTS = {"allow", "fail"}


def _mapping(value: Any, location: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise PolicyError(f"{location} must be an object")
    return value


def _validated_settings(raw: Mapping[str, Any], location: str) -> Dict[str, Any]:
    unknown = sorted(set(raw) - _SETTING_KEYS)
    if unknown:
        raise PolicyError(f"{location} has unknown settings: {', '.join(unknown)}")
    values = dict(raw)
    for name in ("max_modified", "max_only"):
        value = values.get(name)
        if value is not None and (
            isinstance(value, bool) or not isinstance(value, int) or value < 0
        ):
            raise PolicyError(
                f"{location}.{name} must be null or a non-negative integer"
            )
    for name in ("max_modified_fraction", "max_only_fraction"):
        value = values.get(name)
        if value is not None and (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not 0 <= value <= 1
        ):
            raise PolicyError(f"{location}.{name} must be null or between 0 and 1")
    overlap = values.get("min_overlap")
    if overlap is not None and (
        isinstance(overlap, bool)
        or not isinstance(overlap, (int, float))
        or not 0 <= overlap <= 1
    ):
        raise PolicyError(f"{location}.min_overlap must be between 0 and 1")
    if "allow_structural_differences" in values and not isinstance(
        values["allow_structural_differences"], bool
    ):
        raise PolicyError(f"{location}.allow_structural_differences must be boolean")
    for name in (
        "fail_transitions",
        "deny_transitions",
        "allow_fields",
        "deny_fields",
        "allow_transitions",
    ):
        if name not in values or values[name] is None:
            continue
        value = values[name]
        if not isinstance(value, list) or not all(
            isinstance(item, str) and item.strip() for item in value
        ):
            raise PolicyError(
                f"{location}.{name} must be an array of non-empty strings"
            )
        values[name] = tuple(value)
    if values.get("transition_default", "allow") not in _TRANSITION_DEFAULTS:
        raise PolicyError(f"{location}.transition_default must be 'allow' or 'fail'")
    if values.get("severity", "error") not in _SEVERITIES:
        raise PolicyError(f"{location}.severity must be error, warning, or info")
    return values


def load_policy(path: Path) -> PolicyDocument:
    if not path.is_file():
        raise PolicyError(f"policy not found: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise PolicyError(f"invalid policy JSON: {error}") from error
    root = _mapping(payload, "policy")
    unknown = sorted(set(root) - {"gendiff_policy", "defaults", "rules"})
    if unknown:
        raise PolicyError(f"policy has unknown keys: {', '.join(unknown)}")
    if root.get("gendiff_policy") != 1:
        raise PolicyError("gendiff_policy must be 1")
    defaults = {"transition_default": "fail"}
    defaults.update(
        _validated_settings(_mapping(root.get("defaults", {}), "defaults"), "defaults")
    )
    raw_rules = root.get("rules", [])
    if not isinstance(raw_rules, list):
        raise PolicyError("rules must be an array")
    rules = []
    names = set()
    for index, item in enumerate(raw_rules, 1):
        raw = _mapping(item, f"rule {index}")
        unknown = sorted(set(raw) - {"name", "match", "set"})
        if unknown:
            raise PolicyError(f"rule {index} has unknown keys: {', '.join(unknown)}")
        name = raw.get("name")
        if not isinstance(name, str) or not name.strip():
            raise PolicyError(f"rule {index}.name must be a non-empty string")
        if name in names:
            raise PolicyError(f"duplicate rule name: {name}")
        names.add(name)
        match = _mapping(raw.get("match", {}), f"rule {index}.match")
        unknown_match = sorted(set(match) - _MATCH_KEYS)
        if unknown_match:
            raise PolicyError(
                f"rule {index}.match has unknown keys: {', '.join(unknown_match)}"
            )
        if not match:
            raise PolicyError(f"rule {index}.match cannot be empty")
        if not all(
            isinstance(value, str) and value.strip() for value in match.values()
        ):
            raise PolicyError(f"rule {index}.match values must be non-empty strings")
        settings = _validated_settings(
            _mapping(raw.get("set", {}), f"rule {index}.set"), f"rule {index}.set"
        )
        if not settings:
            raise PolicyError(f"rule {index}.set cannot be empty")
        rules.append(PolicyRule(name.strip(), dict(match), settings))
    return PolicyDocument(path.absolute(), defaults, tuple(rules))


def _matches(pattern: str, value: str) -> bool:
    return fnmatchcase(value.casefold(), pattern.casefold())


def resolve_policy(
    base: BatchPolicy,
    document: Optional[PolicyDocument],
    *,
    sample: str,
    stage: str,
    kind: str,
    profile: str,
) -> Tuple[BatchPolicy, Tuple[str, ...]]:
    effective = base
    matched = []
    if document is None:
        return effective, ()
    effective = replace(effective, **document.defaults)
    context = {"sample": sample, "stage": stage, "kind": kind, "profile": profile}
    for rule in document.rules:
        if all(_matches(pattern, context[key]) for key, pattern in rule.match.items()):
            effective = replace(effective, **rule.settings)
            matched.append(rule.name)
    return effective, tuple(matched)


def _pattern_matches(pattern: str, value: str) -> bool:
    folded_pattern = pattern.casefold()
    folded_value = value.casefold()
    return fnmatchcase(folded_value, folded_pattern) or folded_pattern in folded_value


def evaluate_policy(
    result: ComparisonResult,
    policy: BatchPolicy,
    matched_rules: Sequence[str] = (),
) -> PolicyDecision:
    details = result.details
    if details is None:
        raise PolicyError("batch comparisons require detailed record matching")
    violations = []
    trace = [f"matched rule: {name}" for name in matched_rules]
    modified_denominator = max(1, min(result.left_records, result.right_records))
    only_denominator = max(1, result.left_records + result.right_records)
    modified_fraction = details.modified / modified_denominator
    only = details.left_only + details.right_only
    only_fraction = only / only_denominator

    def limit(
        actual: float, maximum: Optional[float], label: str, formatted: str
    ) -> None:
        if maximum is None:
            trace.append(f"{label}: no limit")
        elif actual > maximum:
            violations.append(f"{label} {formatted} exceed {maximum:.1%}")
        else:
            trace.append(f"{label}: {formatted} within {maximum:.1%}")

    if policy.max_modified is not None:
        if details.modified > policy.max_modified:
            violations.append(
                f"modified records {details.modified:,} exceed {policy.max_modified:,}"
            )
        else:
            trace.append(
                f"modified records: {details.modified:,} within {policy.max_modified:,}"
            )
    limit(
        modified_fraction,
        policy.max_modified_fraction,
        "modified fraction",
        f"{modified_fraction:.1%}",
    )
    if policy.max_only is not None:
        if only > policy.max_only:
            violations.append(
                f"sample-only records {only:,} exceed {policy.max_only:,}"
            )
        else:
            trace.append(f"sample-only records: {only:,} within {policy.max_only:,}")
    limit(
        only_fraction,
        policy.max_only_fraction,
        "sample-only fraction",
        f"{only_fraction:.1%}",
    )
    if result.identity_overlap < policy.min_overlap:
        violations.append(
            f"identity overlap {result.identity_overlap:.1%} is below "
            f"{policy.min_overlap:.1%}"
        )
    else:
        trace.append(
            f"identity overlap: {result.identity_overlap:.1%} meets "
            f"{policy.min_overlap:.1%}"
        )
    if not result.structure_equal and not policy.allow_structural_differences:
        violations.append("structural header differs")
    elif not result.structure_equal:
        trace.append("structural header difference allowed")

    changed_fields = set(details.field_changes)
    denied = changed_fields & set(policy.deny_fields)
    if denied:
        violations.append(f"forbidden changed fields: {', '.join(sorted(denied))}")
    if policy.allow_fields is not None:
        unexpected = changed_fields - set(policy.allow_fields)
        if unexpected:
            violations.append(
                f"unexpected changed fields: {', '.join(sorted(unexpected))}"
            )

    for transition, count in details.transitions.items():
        denied_transition = next(
            (
                pattern
                for pattern in tuple(policy.fail_transitions)
                + tuple(policy.deny_transitions)
                if _pattern_matches(pattern, transition)
            ),
            None,
        )
        if denied_transition is not None:
            violations.append(
                f"forbidden transition {denied_transition!r} occurred {count:,} times"
            )
            continue
        allowed = any(
            _pattern_matches(pattern, transition)
            for pattern in policy.allow_transitions
        )
        if policy.transition_default == "fail" and not allowed:
            violations.append(
                f"unapproved transition {transition!r} occurred {count:,} times"
            )
        elif allowed:
            trace.append(f"allowed transition: {transition} ({count:,})")

    if policy.fail_transitions or policy.deny_transitions:
        trace.append("deny rules take precedence over allow rules")
    unique_violations = tuple(dict.fromkeys(violations))
    if policy.severity == "error":
        errors, warnings = unique_violations, ()
    elif policy.severity == "warning":
        errors, warnings = (), unique_violations
    else:
        errors, warnings = (), ()
        trace.extend(f"informational: {item}" for item in unique_violations)
    return PolicyDecision(errors, warnings, tuple(trace), tuple(matched_rules), policy)
