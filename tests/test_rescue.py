"""Flank geometry, missing evidence, shared scoring, and downstream loss safeguards."""

import shutil
from pathlib import Path

import pandas as pd
import pytest

from src import config, rescue
from src.classify import OrthologState, classify_element


def anchor(start=100, end=200, strand="+", chrom="chrTarget"):
    return dict(start=start, end=end, strand=strand, chrom=chrom)


@pytest.mark.parametrize(
    "left,right,reason",
    [
        (None, None, "both_flanks_unmapped"),
        (None, anchor(), "left_flank_unmapped"),
        (anchor(), None, "right_flank_unmapped"),
        (anchor(), anchor(300, 400, chrom="chrOther"), "different_target_chromosomes"),
        (anchor(), anchor(300, 400, strand="-"), "different_target_orientations"),
        (anchor(300, 400), anchor(), "order_not_preserved"),
        (anchor(strand="-"), anchor(300, 400, strand="-"), "order_not_preserved"),
        (anchor(), anchor(150, 250), "negative_gap"),
        (anchor(), anchor(10201, 10301), "gap_too_large"),
        (anchor(-1, 100), anchor(200, 300), "target_out_of_bounds"),
    ],
)
def test_first_failure_is_recorded(left, right, reason):
    pair = rescue.anchor_gap(left, right, {"chrTarget": 20000})
    assert not pair["valid_pair"]
    assert pair["rejection_reason"] == reason


def test_plus_minus_zero_and_maximum_gaps():
    sizes = {"chrTarget": 20000}
    forward = rescue.anchor_gap(anchor(), anchor(400, 600), sizes)
    reverse = rescue.anchor_gap(anchor(400, 600, "-"), anchor(100, 200, "-"), sizes)
    assert forward["valid_pair"] and reverse["valid_pair"]
    assert (forward["gap_start"], forward["gap_end"]) == (200, 400)
    assert (reverse["gap_start"], reverse["gap_end"], reverse["gap_strand"]) == (200, 400, "-")
    assert rescue.anchor_gap(anchor(), anchor(200, 300), sizes)["gap_width"] == 0
    assert rescue.anchor_gap(anchor(), anchor(10200, 10300), sizes)["valid_pair"]


def elements():
    return pd.DataFrame(
        [
            dict(element_id="reverse_oligo", chrom="chrSource", start=1000, end=1200, strand="-"),
            dict(element_id="boundary", chrom="chrSource", start=100, end=300, strand="+"),
        ]
    )


def test_flanks_exclude_full_element_and_track_genomic_orientation():
    left, right, excluded = rescue.build_flanks(elements(), {"chrSource": 5000}, 500)
    assert excluded == {"boundary"}
    assert left[["start", "end", "strand"]].values.tolist() == [[500, 1000, "+"]]
    assert right[["start", "end", "strand"]].values.tolist() == [[1200, 1700, "+"]]
    malformed = elements().assign(end=lambda x: x.start + 150)
    with pytest.raises(ValueError, match="200"):
        rescue.build_flanks(malformed, {"chrSource": 5000}, 500)


def test_source_sidecar_disagreement_is_fatal():
    source = elements().assign(working_build=config.GENOME_BUILD)
    sidecar = source.rename(
        columns={
            k: f"{config.GENOME_BUILD}_element_{k}" for k in ("chrom", "start", "end", "strand")
        }
    )
    sidecar = sidecar.assign(
        call="ABSENT", source_build=config.GENOME_BUILD, target_build="syntheticTarget"
    )
    assert len(rescue.absent_elements(sidecar, source, "syntheticTarget")) == 2
    sidecar.loc[0, f"{config.GENOME_BUILD}_element_start"] += 1
    with pytest.raises(ValueError, match="coordinate mismatch"):
        rescue.absent_elements(sidecar, source, "syntheticTarget")
    with pytest.raises(ValueError, match="target_build"):
        rescue.absent_elements(sidecar, source, "differentTarget")


def geometry():
    return pd.DataFrame(
        [
            dict(
                element_id="wide",
                valid_pair=True,
                gap_chrom="chrTarget",
                gap_start=200,
                gap_end=400,
                gap_width=200,
            ),
            dict(
                element_id="short",
                valid_pair=True,
                gap_chrom="chrTarget",
                gap_start=500,
                gap_end=566,
                gap_width=66,
            ),
            dict(
                element_id="unknown",
                valid_pair=False,
                gap_chrom="",
                gap_start="",
                gap_end="",
                gap_width="",
            ),
        ]
    )


def test_exact_midpoint_windows_and_extension_are_explicit():
    frame = geometry()
    center = rescue.scan_intervals(frame, {"chrTarget": 1000}, "center")
    assert center.scan_bp.tolist() == [101, 101]
    assert center.start.tolist() == [249, 482]
    assert center.end.tolist() == [350, 583]
    assert center.extends_outside_gap.tolist() == [False, True]
    whole = rescue.scan_intervals(frame, {"chrTarget": 1000}, "whole")
    assert whole.scan_bp.tolist() == [200, 66]
    frame.loc[0, ["gap_start", "gap_end"]] = [0, 20]
    with pytest.raises(ValueError, match="unclipped"):
        rescue.scan_intervals(frame, {"chrTarget": 1000}, "center")


def test_three_states_and_uncertainty_survive_classification():
    frame = geometry()
    whole = pd.DataFrame(dict(element_id=["wide", "short"], status=["present", "absent"]))
    center = whole.assign(status=["absent", "absent"])
    resolved = rescue.resolve(frame, {"whole": whole, "center": center}, expected=3).set_index(
        "element_id"
    )
    assert resolved.loc["wide", "whole_resolution"] == "REGION_MOTIF_PRESENT"
    assert resolved.loc["short", "whole_resolution"] == "REGION_BELOW_THRESHOLD"
    assert resolved.loc["wide", "whole_ortholog_state"] == OrthologState.UNSCORABLE
    assert resolved.loc["short", "whole_ortholog_state"] == OrthologState.UNSCORABLE
    assert resolved.loc["unknown", "whole_resolution"] == "NO_SYNTENY"
    assert resolved.loc["wide", "search_space_sensitive"]
    states = {species: OrthologState.NO_MOTIF for species in config.TARGET_SPECIES}
    target = sorted(states)[0]
    states[target] = resolved.loc["unknown", "whole_ortholog_state"]
    classification = classify_element(states)
    assert classification.evolutionary_class == "unresolved"
    assert target in classification.failed_lift_species
    assert target not in classification.no_motif_species
    with pytest.raises(AssertionError):
        rescue.resolve(frame, {"whole": whole, "center": center}, expected=158)
    with pytest.raises(ValueError, match="exactly"):
        rescue.resolve(frame, {"whole": whole.iloc[:1], "center": center}, expected=3)
    with pytest.raises(ValueError, match="unscorable"):
        rescue.resolve(
            frame, {"whole": whole.assign(status="no_hit"), "center": center}, expected=3
        )


@pytest.mark.parametrize(
    "sequence,has_hit,score,expected",
    [
        ("A" * 101, True, 20.86, "present"),
        ("A" * 101, True, 20.85, "absent"),
        ("N" * 101, False, None, "error"),
        ("N" + "A" * 100, True, 1.0, "error"),
        ("A" * 6, False, None, "absent"),
    ],
)
def test_shared_scorer_threshold_and_nohit_safety(
    tmp_path, monkeypatch, sequence, has_hit, score, expected
):
    intervals = pd.DataFrame(
        [
            dict(
                element_id="synthetic",
                chrom="chrTarget",
                start=100,
                end=100 + len(sequence),
                strand="+",
                scan_bp=len(sequence),
                extends_outside_gap=False,
            )
        ]
    )
    matrix = tmp_path / "matrix.meme"
    monkeypatch.setattr(rescue.shutil, "which", lambda _: "twoBitToFa")

    def extract(argv, _log):
        Path(argv[-1]).write_text(f">synthetic\n{sequence}\n")

    monkeypatch.setattr(rescue, "run_tool", extract)

    def scorer(fasta, supplied_matrix):
        assert supplied_matrix == matrix
        assert rescue.read_fasta(fasta) == {"synthetic": sequence}
        return pd.DataFrame(
            [
                dict(
                    element_id="synthetic",
                    has_hit=has_hit,
                    score=score,
                    retains_motif=score >= config.M_VALUE if has_hit else pd.NA,
                    p_value=0.01 if has_hit else None,
                    start=1,
                    stop=21,
                    strand="-",
                )
            ]
        )

    monkeypatch.setattr(rescue.motifs, "score_sequences", scorer)
    if expected == "error":
        with pytest.raises(ValueError, match="Unscorable"):
            rescue.score_intervals(intervals, tmp_path / "genome.2bit", matrix, tmp_path / "work")
    else:
        result = rescue.score_intervals(
            intervals, tmp_path / "genome.2bit", matrix, tmp_path / "work"
        )
        assert result.loc[0, "status"] == expected
        if has_hit:
            assert result.loc[0, "score"] == score
            assert (result.loc[0, "hit_start"], result.loc[0, "hit_end"]) == (100, 121)
        else:
            assert pd.isna(result.loc[0, "score"])
            assert result.loc[0, "score_basis"] == "shorter_than_motif"


def test_zero_gap_is_not_a_fabricated_score(tmp_path):
    intervals = pd.DataFrame(
        [
            dict(
                element_id="deleted",
                chrom="chrTarget",
                start=100,
                end=100,
                strand="+",
                scan_bp=0,
                extends_outside_gap=False,
            )
        ]
    )
    scores = rescue.score_intervals(intervals, tmp_path / "unused", tmp_path / "unused", tmp_path)
    assert scores.loc[0, "status"] == "absent"
    assert pd.isna(scores.loc[0, "score"])
    assert not scores.loc[0, "has_hit"]


@pytest.mark.integration
def test_real_lift_with_other_assembly_and_shared_minmatch(tmp_path, monkeypatch):
    if shutil.which("liftOver") is None:
        pytest.skip("liftOver not installed")
    chain = tmp_path / "synthetic.chain"
    chain.write_text("chain 1000 chrSource 5000 + 0 5000 chrTarget 5000 + 0 5000 1\n5000\n\n")
    monkeypatch.setattr(config, "MINMATCH_CROSS_SPECIES", 0.7)
    frame = rescue.run_anchors(
        elements(),
        {"chrSource": 5000},
        {"chrTarget": 5000},
        "otherAssembly",
        chain,
        tmp_path / "run",
        500,
    )
    mapped = frame.set_index("element_id").loc["reverse_oligo"]
    assert mapped.valid_pair
    assert (mapped.gap_start, mapped.gap_end, mapped.gap_strand) == (1000, 1200, "+")
    assert (
        frame.set_index("element_id").loc["boundary", "rejection_reason"]
        == "source_flank_out_of_bounds"
    )
    assert "-minMatch=0.7" in (tmp_path / "run/left/liftOver.log").read_text()
    assert "otherAssembly" in (tmp_path / "run/right/liftOver.log").read_text()
