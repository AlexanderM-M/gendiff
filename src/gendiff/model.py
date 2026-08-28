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

    def to_dict(self) -> Dict[str, Any]:
        return {
            "records": {
                "identical": self.identical,
                "modified": self.modified,
                "left_only": self.left_only,
                "right_only": self.right_only,
            },
            "field_changes": self.field_changes,
            "top_regions": self.top_regions,
            "sample_changes": self.sample_changes,
            "examples": self.examples,
        }


@dataclass(frozen=True)
class ComparisonResult:
    kind: str
    left: Path
    right: Path
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
            "left": str(self.left),
            "right": str(self.right),
            "records": {"left": self.left_records, "right": self.right_records},
            "logical_records_equal": self.content_equal,
            "structural_header_equal": self.structure_equal,
            "metadata_header_equal": self.metadata_equal,
            "changed_fields": self.changed_fields,
            "relationship": self.relationship,
            "identity_overlap": self.identity_overlap,
            "profile": self.profile,
        }
        if self.details is not None:
            result["details"] = self.details.to_dict()
        return result
