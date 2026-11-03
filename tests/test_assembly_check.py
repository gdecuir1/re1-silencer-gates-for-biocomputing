"""Data-free tests for exact assembly comparisons and bounded evidence downloads."""

import json
import subprocess
from pathlib import Path
from unittest.mock import Mock

import pandas as pd
import pytest

from src import assembly_check, downloads


@pytest.mark.parametrize(
    "expected,reference,orientation",
    [
        ("ACGA", "acga", "+"),
        ("TCGT", "acga", "-"),
        ("ACGT", "ACGT", "both"),
        ("ACGA", "ACGC", "none"),
    ],
)
def test_exact_comparison(expected: str, reference: str, orientation: str) -> None:
    """Require exact identity while permitting case differences and reverse complement."""
    assert assembly_check.exact_orientation(expected, reference) == orientation


def test_sampling_independent_of_order() -> None:
    """Choose the first available interval per group without using reference outcomes."""
    table = pd.DataFrame(
        {
            "silencer": ["synthetic_late", "synthetic_early", "synthetic_unmatched"],
            "project": ["motif"] * 3,
            "chromosome": ["1"] * 3,
            "start": [20, 10, 1],
            "sequence": ["ACGA", "TCGT", None],
        }
    )
    result = assembly_check.select_samples(table)
    assert result.silencer.tolist() == ["synthetic_early"]
    pd.testing.assert_frame_equal(result, assembly_check.select_samples(table.iloc[::-1]))


def test_shift_diagnostics(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Detect a shifted match without changing the reported original coordinates."""
    sequence = "A" * 100 + "C" * 100
    payload = {"genome": "hg38", "chrom": "chr1", "start": 9, "end": 211, "dna": sequence + "GT"}
    path = tmp_path / "hg38_chr1_9_211.json"
    path.write_text(json.dumps(payload))
    monkeypatch.setattr(assembly_check, "EVIDENCE_DIR", tmp_path.parent)

    def local_cache(url: str, destination: Path, **kwargs: object) -> Path:
        """Supply a synthetic response without permitting network access."""
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(payload))
        return destination

    monkeypatch.setattr(assembly_check, "download_small_file", local_cache)
    row = pd.Series(
        {
            "silencer": "synthetic",
            "chromosome": "1",
            "start": 10,
            "end": 210,
            "project": "motif",
            "sequence": sequence,
        }
    )
    result = assembly_check.compare_sample(row, "hg38", fetch=False)
    assert result["exact_orientation"] == "none"
    assert result["minus_one_orientation"] == "+"
    assert result["plus_one_orientation"] == "none"
    assert result["start"] == 10
    assert result["end"] == 210
    payload["genome"] = "hg19"
    with pytest.raises(ValueError, match="response mismatch: genome"):
        assembly_check.compare_sample(row, "hg38", fetch=False)


def test_download_provenance_and_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Append provenance inside the downloader and reject tampered cached bytes."""
    data = tmp_path / "data"
    data.mkdir()
    ledger = data / "README.md"
    ledger.write_text("Existing provenance must remain intact.\n")
    monkeypatch.setattr(downloads.config, "ROOT", tmp_path)
    monkeypatch.setattr(downloads.config, "DATA_DIR", data)
    destination = data / "synthetic.json"
    url = "https://example.invalid/synthetic"

    def fake_curl(args: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        """Simulate a tiny completed HTTP response; no network is used."""
        Path(args[args.index("--output") + 1]).write_text('{"synthetic": true}')
        return subprocess.CompletedProcess(args, 0, "", "")

    runner = Mock(side_effect=fake_curl)
    monkeypatch.setattr(downloads.subprocess, "run", runner)
    assert downloads.download_small_file(url, destination, fetch=True) == destination
    assert ledger.read_text().startswith("Existing provenance must remain intact.\n")
    assert url in ledger.read_text()
    assert downloads.local_md5(destination) in ledger.read_text()
    initial = ledger.read_text()
    downloads.download_small_file(url, destination)
    assert ledger.read_text() == initial
    assert runner.call_count == 1
    destination.write_text("changed")
    with pytest.raises(ValueError, match="changed"):
        downloads.download_small_file(url, destination)


def test_missing_cache_does_not_fetch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Default offline mode must never trigger a download or silently invent data."""
    monkeypatch.setattr(downloads.config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(downloads.config, "ROOT", tmp_path)
    runner = Mock()
    monkeypatch.setattr(downloads.subprocess, "run", runner)
    with pytest.raises(FileNotFoundError, match="--fetch"):
        downloads.download_small_file("https://example.invalid/synthetic", tmp_path / "absent.json")
    runner.assert_not_called()


def test_download_failure_not_recorded_as_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed transfer leaves neither a successful provenance row nor a partial file."""
    monkeypatch.setattr(downloads.config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(downloads.config, "ROOT", tmp_path)
    ledger = tmp_path / "README.md"
    ledger.write_text("Existing history.\n")
    monkeypatch.setattr(
        downloads.subprocess,
        "run",
        Mock(return_value=subprocess.CompletedProcess([], 22, "", "HTTP 403")),
    )
    with pytest.raises(RuntimeError, match="HTTP 403"):
        downloads.download_small_file(
            "https://example.invalid/synthetic", tmp_path / "bad", fetch=True
        )
    assert ledger.read_text() == "Existing history.\n"
    assert not (tmp_path / "bad").exists()
