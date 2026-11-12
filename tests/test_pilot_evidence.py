"""Published coordinate evidence and explicit exclusions from a human-to-mouse contrast."""

import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq

from src import config
from src import pilot_evidence as ev
from src.classify import OrthologState


def test_published_coordinates_and_bases_independent_of_implementation():
    reference = json.loads(ev.REFERENCE.read_text())
    for row in reference["records"]:
        assert row["fimo_start"] == 91 and row["fimo_stop"] == 111
        a, b = config.MOTIF_SPAN
        core = row["sequence"][a : b + 1]
        matched = row["matched_sequence"]
        if row["fimo_strand"] == "-":
            matched = str(Seq(matched).reverse_complement())
        assert core == matched
        assert row["sequence"][91:112] != matched


@pytest.mark.integration
def test_published_fimo_scores_and_off_by_one_regression(tmp_path, monkeypatch):
    if shutil.which("fimo") is None:
        pytest.skip("FIMO not on PATH")
    monkeypatch.setattr(config, "MOTIF_PATH", Path(__file__).parent / "fixtures/MA0138.2.meme")
    result = ev.check_reference(tmp_path / "correct")
    assert len(result) == 5
    assert result.score.sub(result.score_published).abs().max() == 0
    monkeypatch.setattr(config, "MOTIF_SPAN", (91, 111))
    with pytest.raises(ValueError, match="published FIMO coordinates"):
        ev.check_reference(tmp_path / "wrong")


def case(**kwargs):
    return (
        dict(
            human_core_status="present",
            alignment_status="INTACT_MOTIF",
            mouse_core_status="below_threshold",
            mouse_window_status="below_threshold",
        )
        | kwargs
    )


@pytest.mark.parametrize(
    "overrides,status,state",
    [
        ({}, "NONRETENTION_CANDIDATE", OrthologState.NO_MOTIF),
        (
            {"mouse_core_status": "present", "mouse_window_status": "present"},
            "RETAINED_CANDIDATE",
            OrthologState.RETAINED,
        ),
        ({"human_core_status": "below_threshold"}, "HUMAN_WEAK", OrthologState.UNSCORABLE),
        ({"human_core_status": "unscorable"}, "SOURCE_UNSCORABLE", OrthologState.UNSCORABLE),
        ({"alignment_status": "DISRUPTED_MOTIF"}, "UNRESOLVED", OrthologState.UNSCORABLE),
        ({"alignment_status": "DISCORDANT_CHAINS"}, "UNRESOLVED", OrthologState.UNSCORABLE),
        ({"alignment_status": "AMBIGUOUS_CHAIN"}, "UNRESOLVED", OrthologState.UNSCORABLE),
        ({"alignment_status": "REGIONAL_ANCHORS_ONLY"}, "UNRESOLVED", OrthologState.UNSCORABLE),
        (
            {"alignment_status": "REGIONAL_ANCHORS_ONLY", "mouse_core_status": "present"},
            "UNRESOLVED",
            OrthologState.UNSCORABLE,
        ),
        ({"alignment_status": "NO_MOTIF_CORRESPONDENCE"}, "UNRESOLVED", OrthologState.FAILED_LIFT),
        ({"mouse_core_status": "unscorable"}, "UNRESOLVED", OrthologState.UNSCORABLE),
        ({"mouse_window_status": "present"}, "UNRESOLVED", OrthologState.UNSCORABLE),
    ],
)
def test_negatives_cannot_bypass_source_or_alignment_gates(overrides, status, state):
    result, reason, adapter = ev.comparison(case(**overrides))
    assert result == status and adapter == state and reason


def test_unknown_source_status_is_not_silently_weak():
    with pytest.raises(ValueError, match="Unknown human"):
        ev.comparison(case(human_core_status="missing"))


def test_missing_hits_and_ambiguous_negatives_preserve_uncertainty(tmp_path, monkeypatch):
    records = {"hit": "A" * 21, "weak": "C" * 21, "unknown": "N" * 21, "partial": "A" * 20 + "N"}
    hits = pd.DataFrame(
        dict(
            element_id=list(records),
            has_hit=[True, True, False, True],
            retains_motif=[True, False, pd.NA, False],
            score=[25.0, 10.0, float("nan"), 9.0],
        )
    )
    monkeypatch.setattr(ev.motifs, "score_sequences", lambda *_: hits)
    result = ev.score_records(records, tmp_path).set_index("element_id")
    assert result.status.to_dict() == dict(
        hit="present", weak="below_threshold", unknown="unscorable", partial="unscorable"
    )
    assert pd.isna(result.loc["unknown", "score"])
    assert result.loc["partial", "score"] == 9.0


def test_human_reference_mismatch_stops_publication(tmp_path, monkeypatch):
    source = pd.DataFrame(
        [
            dict(
                element_id="x",
                chrom="chrS",
                start=0,
                end=200,
                strand="+",
                sequence="A" * 200,
                working_build=config.GENOME_BUILD,
            )
        ]
    )
    sidecar = pd.DataFrame(dict(element_id=["x"]))
    monkeypatch.setattr(ev, "extract_intervals", lambda *_: {"x": "C" * 200})
    with pytest.raises(ValueError, match="reference differs"):
        ev.human_scores(source, sidecar, tmp_path / "genome", tmp_path)


def test_all_ids_survive_and_only_supported_negatives_enter_contrast():
    ids = ["weak", "gapped", "candidate", "regional"]
    sidecar = pd.DataFrame(
        dict(
            element_id=ids,
            call=["MAPPED"] * 3 + ["ABSENT"],
            call_reason=["both_lifted"] * 3 + ["neither_lifted"],
        )
    )
    baseline = pd.DataFrame(
        dict(
            element_id=ids,
            motif_status=["absent"] * 3 + ["unscored"],
            score=[1.0, 2.0, 3.0, float("nan")],
        )
    )
    human = pd.DataFrame(
        dict(element_id=ids, human_core_status=["below_threshold"] + ["present"] * 3)
    )
    quality = pd.DataFrame(
        dict(
            element_id=ids,
            alignment_status=[
                "INTACT_MOTIF",
                "DISRUPTED_MOTIF",
                "INTACT_MOTIF",
                "REGIONAL_ANCHORS_ONLY",
            ],
        )
    )
    core = pd.DataFrame(dict(element_id=ids[:3], status=["below_threshold"] * 3))
    regional = pd.DataFrame(
        dict(
            element_id=["regional"], valid_pair=[True], whole_resolution=["REGION_BELOW_THRESHOLD"]
        )
    )
    result = ev.assemble(sidecar, baseline, human, core, quality, regional).set_index("element_id")
    assert set(result.index) == set(ids)
    assert result.comparison_status.to_dict() == dict(
        weak="HUMAN_WEAK",
        gapped="UNRESOLVED",
        candidate="NONRETENTION_CANDIDATE",
        regional="UNRESOLVED",
    )
    assert result.candidate_nonretention.sum() == 1
    assert result.loc[["weak", "gapped", "regional"], "exclude_mouse_downstream"].all()
    assert not result.loc["candidate", "exclude_mouse_downstream"]
    assert result.loc["regional", "ortholog_state"] == OrthologState.UNSCORABLE
    with pytest.raises(ValueError, match="exactly the MAPPED"):
        ev.assemble(sidecar, baseline, human, core.iloc[:2], quality, regional)


def test_manifest_check_detects_changed_inputs_and_outputs(tmp_path, monkeypatch):
    input_path, output_path = tmp_path / "input.csv", tmp_path / "calls.csv"
    input_path.write_text("source\n")
    output_path.write_text("calls\n")
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = dict(
        inputs={"input.csv": digest(input_path)},
        outputs={"calls.csv": digest(output_path)},
        tools={},
    )
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(ev, "OUTPUT", tmp_path)
    monkeypatch.setattr(sys, "argv", ["pilot_evidence", "--check"])
    ev.main()
    input_path.write_text("different source\n")
    with pytest.raises(ValueError, match="stale or modified: input.csv"):
        ev.main()
    input_path.write_text("source\n")
    output_path.write_text("different calls\n")
    with pytest.raises(ValueError, match="stale or modified: calls.csv"):
        ev.main()
