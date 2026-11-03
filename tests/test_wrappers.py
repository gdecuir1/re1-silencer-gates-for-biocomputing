"""Exercise external-program contracts with explicit synthetic output fixtures."""

from pathlib import Path

import pandas as pd
import pytest

from src import motifs, orthologs


def test_fimo_best_score_and_unscorable(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Select maximum scores and retain sequences without hits as unknown."""
    fasta, matrix = tmp_path / "synthetic.fa", tmp_path / "synthetic.meme"
    fasta.write_text(">synthetic_a\n" + "A" * 200 + "\n>synthetic_n\n" + "N" * 200 + "\n")
    matrix.write_text("MEME version 4\nMOTIF SYNTHETIC_TEST_ONLY\n")

    def fake_fimo(args: list[str], *, stdout_path: Path) -> str:
        """Write explicitly synthetic FIMO TSV to test best-hit selection."""
        assert args[args.index("--thresh") + 1] == "1"
        assert "--text" in args
        stdout_path.write_text(
            "motif_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\n"
            "SYNTHETIC_TEST_ONLY\tseq_0\t9\t29\t+\t12\t0.2\n"
            "SYNTHETIC_TEST_ONLY\tseq_0\t2\t22\t-\t30\t0.01\n"
            "SYNTHETIC_TEST_ONLY\tseq_0\t3\t23\t+\t30\t0.01\n"
        )
        return ""

    monkeypatch.setattr(motifs, "run_command", fake_fimo)
    result = motifs.score_sequences(fasta, matrix, motif_id="SYNTHETIC_TEST_ONLY")
    assert result["element_id"].tolist() == ["synthetic_a", "synthetic_n"]
    assert result.loc[0, "score"] == 30
    assert result.loc[0, "start"] == 2
    assert result.loc[0, "retains_motif"]
    assert not result.loc[1, "has_hit"]
    assert pd.isna(result.loc[1, "score"])
    assert pd.isna(result.loc[1, "retains_motif"])


def test_lift_returns_failed_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, intervals: pd.DataFrame
) -> None:
    """Send 200 bp intervals and preserve the original failed record and reason."""
    chain = tmp_path / "synthetic.chain"
    chain.touch()

    def fake_liftover(args: list[str]) -> str:
        """Return one synthetic lift and one explicit failed lift."""
        lines = Path(args[2]).read_text().splitlines()
        assert all(int(line.split("\t")[2]) - int(line.split("\t")[1]) == 200 for line in lines)
        Path(args[-2]).write_text(lines[0] + "\n")
        Path(args[-1]).write_text("#Deleted in new\n" + lines[1] + "\n")
        return ""

    monkeypatch.setattr(orthologs, "run_command", fake_liftover)
    result = orthologs.lift_elements(intervals, chain, min_match=0.95)
    assert result.mapped["element_id"].tolist() == ["synthetic_plus"]
    assert result.unmapped["element_id"].tolist() == ["synthetic_minus"]
    assert result.unmapped.loc[0, "reason"] == "Deleted in new"
    assert result.unmapped.loc[0, "start"] == 400
    assert "synthetic_minus" in result.unmapped.loc[0, "raw_record"]


def test_all_failed_lifts_are_valid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, intervals: pd.DataFrame
) -> None:
    """An empty mapped BED is valid when every record is accounted for as failed."""
    chain = tmp_path / "synthetic.chain"
    chain.touch()

    def fake_liftover(args: list[str]) -> str:
        """Preserve every source record in synthetic unmapped output."""
        Path(args[-2]).write_text("")
        Path(args[-1]).write_text(Path(args[2]).read_text())
        return ""

    monkeypatch.setattr(orthologs, "run_command", fake_liftover)
    result = orthologs.lift_elements(intervals, chain, min_match=0.95)
    assert result.mapped.empty
    assert len(result.unmapped) == 2


def test_bedtools_skipped_interval_is_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, intervals: pd.DataFrame
) -> None:
    """Catch bedtools' successful exit with missing chromosome records."""
    genome = tmp_path / "synthetic.fa"
    genome.touch()

    def fake_bedtools(args: list[str], *, stdout_path: Path) -> str:
        """Simulate a silently skipped interval."""
        assert "-s" in args
        stdout_path.write_text("chrSynthetic\t100\t300\tsynthetic_plus\t0\t+\t" + "A" * 200 + "\n")
        return ""

    monkeypatch.setattr(orthologs, "run_command", fake_bedtools)
    with pytest.raises(ValueError, match="every element"):
        orthologs.extract_sequences(intervals, genome)
