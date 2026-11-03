"""Check fixed constants and portable repository paths."""

import importlib
import json
from pathlib import Path

import pytest

from src import config


def test_verified_constants() -> None:
    """Keep the user-supplied scientific constants unchanged."""
    assert config.ELEMENT_LENGTH == 200
    assert config.MOTIF_SPAN == (91, 111)
    assert config.M_VALUE == 20.86
    assert config.CELL_TYPES == ["GM12878", "HepG2", "K562", "SK-N-SH"]
    assert config.JASPAR_ID == "MA0138.2"
    assert config.GENOME_BUILD == "hg38"
    assert config.TARGET_SPECIES == {"rheMac10", "mm39", "canFam6", "bosTau9", "monDom5"}


def test_paths_independent_of_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Resolve all configured paths from the repository, not the current directory."""
    expected = Path(__file__).resolve().parents[1]
    monkeypatch.chdir(tmp_path)
    reloaded = importlib.reload(config)
    assert reloaded.ROOT == expected
    paths = [value for value in vars(reloaded).values() if isinstance(value, Path)]
    paths += list(reloaded.CHAIN_PATHS.values()) + list(reloaded.GENOME_PATHS.values())
    assert paths
    assert all(path.is_absolute() and path.is_relative_to(expected) for path in paths)


def test_notebook_has_no_outputs() -> None:
    """Ensure the starter notebook contains neither results nor execution counts."""
    notebook = json.loads((config.NOTEBOOKS_DIR / "01_explore.ipynb").read_text())
    for cell in notebook["cells"]:
        if cell["cell_type"] == "code":
            assert cell["outputs"] == []
            assert cell["execution_count"] is None
