from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List


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

    @property
    def equivalent(self) -> bool:
        return self.content_equal and self.structure_equal

    def to_dict(self) -> Dict[str, Any]:
        return {
            "equivalent": self.equivalent,
            "type": self.kind,
            "left": str(self.left),
            "right": str(self.right),
            "records": {"left": self.left_records, "right": self.right_records},
            "logical_records_equal": self.content_equal,
            "structural_header_equal": self.structure_equal,
            "metadata_header_equal": self.metadata_equal,
            "changed_fields": self.changed_fields,
        }
