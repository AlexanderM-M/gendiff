from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass(frozen=True)
class DifferenceDetails:
    identical: int
    modified: int
    left_only: int
    right_only: int
    field_changes: Dict[str, int]
    top_regions: List[Dict[str, Any]] = field(default_factory=list)
    sample_changes: Dict[str, int] = field(default_factory=dict)
    examples: List[Dict[str, Any]] = field(default_factory=list)
    loci: List[str] = field(default_factory=list)

    def to_dict(self, left_label: str, right_label: str) -> Dict[str, Any]:
        statuses = {
            "left_only": "only_in_first",
            "right_only": "only_in_second",
        }
        examples = [
            {
                "identity": example["identity"],
                "status": statuses.get(example["status"], example["status"]),
                "changed_fields": example["changed_fields"],
                "first": {
                    "label": left_label,
                    "record": example["left"],
                },
                "second": {
                    "label": right_label,
                    "record": example["right"],
                },
            }
            for example in self.examples
        ]
        return {
            "records": {
                "identical": self.identical,
                "modified": self.modified,
                "only_in_first": {
                    "label": left_label,
                    "count": self.left_only,
                },
                "only_in_second": {
                    "label": right_label,
                    "count": self.right_only,
                },
            },
            "field_changes": self.field_changes,
            "top_regions": self.top_regions,
            "sample_changes": self.sample_changes,
            "examples": examples,
        }


@dataclass(frozen=True)
class ComparisonResult:
    kind: str
    left: Path
    right: Path
    left_label: str
    right_label: str
    left_records: int
    right_records: int
    content_equal: bool
    structure_equal: bool
    metadata_equal: bool
    changed_fields: List[str]
    relationship: str = "unknown"
    identity_overlap: float = 0.0
    profile: str = "strict"
    details: Optional[DifferenceDetails] = None

    @property
    def equivalent(self) -> bool:
        return self.content_equal and self.structure_equal

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "equivalent": self.equivalent,
            "type": self.kind,
            "inputs": [
                {
                    "label": self.left_label,
                    "path": str(self.left),
                    "records": self.left_records,
                },
                {
                    "label": self.right_label,
                    "path": str(self.right),
                    "records": self.right_records,
                },
            ],
            "logical_records_equal": self.content_equal,
            "structural_header_equal": self.structure_equal,
            "metadata_header_equal": self.metadata_equal,
            "changed_fields": self.changed_fields,
            "relationship": self.relationship,
            "identity_overlap": self.identity_overlap,
            "profile": self.profile,
        }
        if self.details is not None:
            result["details"] = self.details.to_dict(self.left_label, self.right_label)
        return result
