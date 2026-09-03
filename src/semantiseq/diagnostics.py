from __future__ import annotations

from typing import Any, Dict, Mapping, Sequence


def _contigs(header: Mapping[str, Any]) -> Dict[str, Any]:
    if "contigs" in header:
        return {
            str(name): item.get("length") for name, item in header["contigs"].items()
        }
    return {str(item.get("SN")): item.get("LN") for item in header.get("SQ", [])}


def _set_changes(label: str, first: Sequence[str], second: Sequence[str]) -> list[str]:
    messages = []
    added = sorted(set(second) - set(first))
    removed = sorted(set(first) - set(second))
    if added:
        messages.append(f"{label} added: {', '.join(added[:10])}")
    if removed:
        messages.append(f"{label} removed: {', '.join(removed[:10])}")
    return messages


def header_differences(
    kind: str, first: Mapping[str, Any], second: Mapping[str, Any]
) -> list[str]:
    messages = []
    first_contigs, second_contigs = _contigs(first), _contigs(second)
    messages.extend(_set_changes("Contigs", first_contigs, second_contigs))
    changed_lengths = [
        name
        for name in sorted(set(first_contigs) & set(second_contigs))
        if first_contigs[name] != second_contigs[name]
    ]
    if changed_lengths:
        messages.append(f"Contig lengths differ: {', '.join(changed_lengths[:10])}")

    if kind == "alignment":
        first_groups = {item.get("ID"): item for item in first.get("RG", [])}
        second_groups = {item.get("ID"): item for item in second.get("RG", [])}
        messages.extend(_set_changes("Read groups", first_groups, second_groups))
        for name in sorted(set(first_groups) & set(second_groups)):
            changed = [
                key
                for key in ("SM", "LB", "PL", "PU")
                if first_groups[name].get(key) != second_groups[name].get(key)
            ]
            if changed:
                messages.append(f"Read group {name} changed: {', '.join(changed)}")
    else:
        messages.extend(
            _set_changes("Samples", first.get("samples", []), second.get("samples", []))
        )
        for section, label in (
            ("info", "INFO definitions"),
            ("formats", "FORMAT definitions"),
            ("filters", "FILTER definitions"),
        ):
            first_values = first.get(section, {})
            second_values = second.get(section, {})
            messages.extend(_set_changes(label, first_values, second_values))
            changed = [
                name
                for name in sorted(set(first_values) & set(second_values))
                if first_values[name] != second_values[name]
            ]
            if changed:
                messages.append(f"{label} changed: {', '.join(changed[:10])}")

    shared = set(first_contigs) & set(second_contigs)
    union = set(first_contigs) | set(second_contigs)
    if union and (len(shared) / len(union) < 0.8 or changed_lengths):
        messages.insert(
            0,
            "Reference compatibility warning: contig names or lengths indicate "
            "different reference definitions.",
        )
    return messages
