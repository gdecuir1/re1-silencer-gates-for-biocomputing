"""Synthetic tests for the mouse-scoring downloader. No network or genomes."""

import subprocess
from pathlib import Path

import pytest

from src import config, fetch_mouse_scoring_inputs, pilot_mouse


def test_chain_url_matches_the_mouse_lift() -> None:
    assert fetch_mouse_scoring_inputs.CHAIN_URL == pilot_mouse.CHAIN_URL


def test_matrix_requires_jaspar_id(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    destination = tmp_path / "MA0138.2.meme"
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(config, "MOTIF_PATH", destination)

    def fake_download(url: str, path: Path, **kwargs: object) -> Path:
        assert url == fetch_mouse_scoring_inputs.MATRIX_URL
        assert kwargs["fetch"] is True
        path.write_text("MEME version 4\n\nMOTIF MA0138.2 REST\n")
        return path

    monkeypatch.setattr(fetch_mouse_scoring_inputs, "download_small_file", fake_download)
    assert fetch_mouse_scoring_inputs.acquire_matrix() == destination
    assert config.JASPAR_ID in fetch_mouse_scoring_inputs.motif_ids(destination)


def test_html_matrix_is_moved_aside(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    destination = tmp_path / "MA0138.2.meme"
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(config, "MOTIF_PATH", destination)

    def fake_download(url: str, path: Path, **kwargs: object) -> Path:
        path.write_text("<html>not a matrix</html>")
        return path

    monkeypatch.setattr(fetch_mouse_scoring_inputs, "download_small_file", fake_download)
    with pytest.raises(ValueError, match="MA0138.2"):
        fetch_mouse_scoring_inputs.acquire_matrix()
    assert not destination.exists()
    assert destination.with_name("MA0138.2.meme.rejected").is_file()


def test_chain_uses_resumable_fetch(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = {}

    monkeypatch.setattr(fetch_mouse_scoring_inputs, "remote_size", lambda url: 73017320)
    monkeypatch.setattr(
        fetch_mouse_scoring_inputs,
        "fetch",
        lambda url, path, size: calls.update(url=url, path=path, size=size) or path,
    )
    fetch_mouse_scoring_inputs.acquire_chain()
    assert calls == {
        "url": fetch_mouse_scoring_inputs.CHAIN_URL,
        "path": config.CHAIN_PATHS["mm39"],
        "size": 73017320,
    }


def test_complete_mm39_does_not_call_the_genome_script(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    genome = tmp_path / "mm39.2bit"
    genome.write_bytes(b"1234")
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "MM39_GENOME", genome)
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "remote_size", lambda url: 4)

    def fail_if_called(*args: object, **kwargs: object) -> subprocess.CompletedProcess:
        raise AssertionError("complete genome must not be downloaded again")

    monkeypatch.setattr(fetch_mouse_scoring_inputs.subprocess, "run", fail_if_called)
    assert fetch_mouse_scoring_inputs.acquire_genome() == genome


def test_missing_mm39_requests_only_that_assembly(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    genome = tmp_path / "data" / "raw" / "genomes" / "mm39.2bit"
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "MM39_GENOME", genome)
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "remote_size", lambda url: 4)
    seen = {}

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess:
        seen["command"] = command
        seen["cwd"] = kwargs["cwd"]
        genome.parent.mkdir(parents=True)
        genome.write_bytes(b"1234")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(fetch_mouse_scoring_inputs.subprocess, "run", fake_run)
    assert fetch_mouse_scoring_inputs.acquire_genome() == genome
    assert seen["command"] == [
        "bash",
        str(tmp_path / "scripts" / "fetch_genomes.sh"),
        "--yes",
        "mm39",
    ]
    assert seen["cwd"] == tmp_path


def test_metadata_failure_stops(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        fetch_mouse_scoring_inputs,
        "fetch_metadata",
        lambda assemblies: ["mm39/md5sum.txt: missing"],
    )
    with pytest.raises(RuntimeError, match="mm39/md5sum.txt"):
        fetch_mouse_scoring_inputs.acquire_metadata()


def test_main_points_at_the_offline_ortholog_step(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "acquire_matrix", lambda: Path("matrix"))
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "acquire_chain", lambda: Path("chain"))
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "acquire_metadata", lambda: None)
    monkeypatch.setattr(fetch_mouse_scoring_inputs, "acquire_genome", lambda: Path("genome"))
    fetch_mouse_scoring_inputs.main([])
    output = capsys.readouterr().out
    assert "does not write orthologs.mm39.csv" in output
    assert "make pilot-orthologs" in output
