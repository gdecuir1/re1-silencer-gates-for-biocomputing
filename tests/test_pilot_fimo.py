"""Synthetic checks for mouse-window FIMO calls. No genome or network is used."""

import pandas as pd
import pytest

from src import config, pilot, pilot_fimo


def _audit() -> pd.DataFrame:
    """One element of each mapping class, plus a second scored element."""
    return pd.DataFrame(
        [
            {
                "element_id": "synthetic_present",
                "call": "MAPPED",
                "call_reason": "both_lifted",
                "mm39_motif_chrom": "chrSynthetic",
                "mm39_motif_start": 100,
                "mm39_motif_end": 121,
                "mm39_motif_strand": "+",
            },
            {
                "element_id": "synthetic_absent",
                "call": "MAPPED",
                "call_reason": "both_lifted",
                "mm39_motif_chrom": "chrSynthetic",
                "mm39_motif_start": 10,
                "mm39_motif_end": 31,
                "mm39_motif_strand": "-",
            },
            {
                "element_id": "synthetic_undetermined",
                "call": "UNDETERMINED",
                "call_reason": "element_only_motif_unmapped",
                "mm39_motif_chrom": pd.NA,
                "mm39_motif_start": pd.NA,
                "mm39_motif_end": pd.NA,
                "mm39_motif_strand": pd.NA,
            },
            {
                "element_id": "synthetic_alignment_gap",
                "call": "ABSENT",
                "call_reason": "neither_lifted",
                "mm39_motif_chrom": pd.NA,
                "mm39_motif_start": pd.NA,
                "mm39_motif_end": pd.NA,
                "mm39_motif_strand": pd.NA,
            },
        ]
    )


def _hits() -> pd.DataFrame:
    """A threshold hit on the plus window and a weak hit on the minus window."""
    return pd.DataFrame(
        [
            {
                "element_id": "synthetic_present",
                "score": config.M_VALUE,
                "p_value": 0.01,
                "has_hit": True,
                "retains_motif": True,
                "start": 41,
                "stop": 61,
                "strand": "+",
            },
            {
                "element_id": "synthetic_absent",
                "score": config.M_VALUE - 1,
                "p_value": 0.2,
                "has_hit": True,
                "retains_motif": False,
                "start": 2,
                "stop": 4,
                "strand": "-",
            },
        ]
    )


def test_flank_is_genomic_and_clips_at_chromosome_ends() -> None:
    windows = pilot_fimo.scoring_windows(_audit(), {"chrSynthetic": 130}, flank=40)
    present = windows.set_index("element_id").loc["synthetic_present"]
    absent = windows.set_index("element_id").loc["synthetic_absent"]
    assert (present.window_start, present.window_end, present.flank_left, present.flank_right) == (
        60,
        130,
        40,
        9,
    )
    assert bool(present.flank_clipped)
    assert (absent.window_start, absent.window_end) == (0, 71)
    assert int(absent.flank_left) == 10
    assert windows.set_index("element_id").loc["synthetic_alignment_gap", "window_chrom"] == ""


def test_threshold_and_unscored_rows_stay_distinct() -> None:
    windows = pilot_fimo.scoring_windows(_audit(), {"chrSynthetic": 500})
    frame = pilot_fimo.assemble(_audit(), windows, _hits())
    status = frame.set_index("element_id")["motif_status"]
    assert status.to_dict() == {
        "synthetic_absent": "absent",
        "synthetic_alignment_gap": "unscored",
        "synthetic_present": "present",
        "synthetic_undetermined": "unscored",
    }
    present = frame.set_index("element_id").loc["synthetic_present"]
    # Plus-strand FIMO positions 41-61 on a window starting at 60.
    assert (int(present.hit_start), int(present.hit_end)) == (100, 121)
    absent = frame.set_index("element_id").loc["synthetic_absent"]
    # Minus-strand positions 2-4 on window [0, 71).
    assert (int(absent.hit_start), int(absent.hit_end)) == (67, 70)
    assert pd.isna(frame.set_index("element_id").loc["synthetic_undetermined", "score"])


def test_minus_strand_hit_uses_the_oriented_sequence() -> None:
    assert pilot_fimo.hit_interval(100, 110, "-", 2, 4) == (106, 109)
    assert pilot_fimo.hit_interval(100, 110, "+", 2, 4) == (101, 104)


def test_spotcheck_is_stable_and_the_report_leaves_observations_blank() -> None:
    windows = pilot_fimo.scoring_windows(_audit(), {"chrSynthetic": 500})
    frame = pilot_fimo.assemble(_audit(), windows, _hits())
    first = pilot_fimo.choose_spotcheck(["synthetic_absent", "synthetic_present"], count=2)
    second = pilot_fimo.choose_spotcheck(["synthetic_present", "synthetic_absent"], count=2)
    assert first == second == ["synthetic_absent", "synthetic_present"]
    text = pilot_fimo.render(frame, first, {"ortholog_sidecar": "abc"})
    assert f">= {config.M_VALUE}" in text
    assert f"{pilot_fimo.FLANK_BP} bp" in text
    assert "Browser observation" in text
    assert "| synthetic_present |" in text
    assert text.split("| synthetic_present |")[1].split("\n")[0].endswith("| |")
    assert "does not prove REST" in text
    assert str(pilot.SAMPLE_SEED) in text


def test_fimo_must_cover_every_mapped_element() -> None:
    windows = pilot_fimo.scoring_windows(_audit(), {"chrSynthetic": 500})
    with pytest.raises(ValueError, match="every MAPPED"):
        pilot_fimo.assemble(_audit(), windows, _hits().iloc[[0]])


def test_missing_sidecar_does_not_publish(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(pilot_fimo.pilot_orthologs, "SIDECAR", tmp_path / "missing.csv")
    monkeypatch.setattr(pilot_fimo, "REPORT", tmp_path / "report.md")
    monkeypatch.setattr(pilot_fimo, "SCORES", tmp_path / "scores.csv")
    with pytest.raises(FileNotFoundError, match="pilot-orthologs"):
        pilot_fimo.main()
    assert not (tmp_path / "report.md").exists()
    assert not (tmp_path / "scores.csv").exists()
