from pathlib import Path

import pytest

from gendiff.artifacts import output_workspace


def test_force_does_not_replace_an_unmanaged_directory(tmp_path: Path) -> None:
    target = tmp_path / "existing"
    target.mkdir()
    marker = target / "keep.txt"
    marker.write_text("keep")

    with pytest.raises(ValueError, match="unmanaged"):
        with output_workspace(target, force=True):
            pass

    assert marker.read_text() == "keep"
