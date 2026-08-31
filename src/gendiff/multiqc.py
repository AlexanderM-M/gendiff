"""Native MultiQC module for GenDiff batch JSON reports."""

from __future__ import annotations

import json

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table


class MultiqcModule(BaseMultiqcModule):
    def __init__(self) -> None:
        super().__init__(
            name="GenDiff",
            anchor="gendiff",
            href="https://github.com/AlexanderM-M/gendiff",
            info="Semantic regression checks for genomics pipeline outputs.",
        )
        data = {}
        for source in self.find_log_files("gendiff"):
            try:
                payload = json.loads(source["f"])
            except (json.JSONDecodeError, TypeError):
                continue
            if payload.get("gendiff_batch_format") != 1:
                continue
            self.add_data_source(source)
            for item in payload.get("comparisons", []):
                comparison = item.get("comparison", {})
                details = comparison.get("details", {}).get("records", {})
                sample = str(item.get("sample", "sample"))
                stage = str(item.get("stage", "stage"))
                name = self.clean_s_name(f"{sample} — {stage}", source)
                data[name] = {
                    "passed": 1 if item.get("passed") else 0,
                    "overlap": float(comparison.get("identity_overlap", 0)) * 100,
                    "modified": int(details.get("modified", 0)),
                    "only": int(details.get("only_in_first", {}).get("count", 0))
                    + int(details.get("only_in_second", {}).get("count", 0)),
                }
        data = self.ignore_samples(data)
        if not data:
            raise ModuleNoSamplesFound
        self.general_stats_addcols(
            data,
            {
                "passed": {
                    "title": "GenDiff",
                    "description": "GenDiff policy result (1 pass, 0 fail)",
                    "min": 0,
                    "max": 1,
                    "scale": "RdYlGn",
                },
                "overlap": {
                    "title": "Overlap",
                    "description": "GenDiff identity overlap",
                    "suffix": "%",
                    "min": 0,
                    "max": 100,
                    "scale": "RdYlGn",
                    "hidden": True,
                },
            },
        )
        self.add_section(
            name="Pipeline comparisons",
            anchor="gendiff-comparisons",
            description="The most important GenDiff regression metrics by stage.",
            plot=table.plot(
                data,
                {
                    "passed": {"title": "Pass", "min": 0, "max": 1},
                    "overlap": {
                        "title": "Overlap",
                        "suffix": "%",
                        "min": 0,
                        "max": 100,
                    },
                    "modified": {"title": "Modified"},
                    "only": {"title": "Only in one input"},
                },
                {"id": "gendiff_comparisons", "title": "GenDiff comparisons"},
            ),
        )
