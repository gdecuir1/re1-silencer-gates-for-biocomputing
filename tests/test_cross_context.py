"""Worst-context selection, coverage, rank denominators and separate references."""

import numpy as np
import pandas as pd
import pytest

from src.cross_context import CONTEXTS, analyze, reference_rankings


def measurements(values, scores=None, references=None):
    """Construct paired observations with independently specified context effects."""
    rows = []
    for element, ranges in values.items():
        for (cell, enhancer), dr in zip(CONTEXTS, ranges, strict=True):
            rows.append(
                dict(
                    element_id=element,
                    cell_type=cell,
                    enhancer=enhancer,
                    reference_type=(references or {}).get(element, "scrambled"),
                    fimo_score=(scores or {}).get(element, 25.0),
                    on_activity=dr + 2.0,
                    off_activity=2.0,
                    dynamic_range=dr,
                    separation=dr / 2.0,
                    leakiness=2.0 if pd.notna(dr) else np.nan,
                )
            )
    return pd.DataFrame(rows)


def test_worst_context_precedes_higher_median_and_order_is_deterministic():
    frame = measurements({"consistent": [4] * 5, "fragile": [5, 5, 5, 5, -1], "middle": [2] * 5})
    candidates, _, _, _, _ = analyze(frame)
    assert candidates.element_id.tolist() == ["consistent", "fragile", "middle"]
    a, b = candidates.iloc[0], candidates.iloc[1]
    assert a.minimum_rank_percentile > b.minimum_rank_percentile
    assert a.median_rank_percentile < b.median_rank_percentile
    pd.testing.assert_frame_equal(candidates, analyze(frame.sample(frac=1, random_state=12))[0])


def test_coverage_exclusions_remain_in_context_denominator_and_missing_is_not_worst():
    frame = measurements(
        {
            "complete": [1] * 5,
            "four": [2, 2, 2, 2, np.nan],
            "three": [3, 3, 3, np.nan, np.nan],
            "none": [np.nan] * 5,
        }
    )
    candidates, audit, coverage, denominators, summary = analyze(frame)
    assert candidates.element_id.tolist() == ["four", "complete"]
    four = candidates.iloc[0]
    assert four.n_contexts_measured == 4
    assert four.minimum_rank_percentile == pytest.approx(200 / 3)
    assert pd.isna(four["percentile__SK-N-SH__En19"])
    assert four.missing_contexts == "SK-N-SH/En19"
    assert denominators.n_measured_scrambled.tolist() == [3, 3, 3, 2, 1]
    assert coverage.elements.tolist() == [1, 0, 0, 1, 1, 1]
    assert summary["n_qualified"] == summary["n_excluded_coverage"] == 2
    none = audit.set_index("element_id").loc["none"]
    assert pd.isna(none.minimum_rank_percentile) and pd.isna(none.cross_context_rank)
    assert pd.isna(none.dynamic_range_spread)


def test_percentiles_use_average_ties_and_exclude_other_reference_and_weak_motifs():
    frame = measurements(
        {"a": [1] * 5, "b": [1] * 5, "c": [3] * 5, "fallback": [9] * 5, "weak": [9] * 5},
        scores={"a": 20.86, "weak": 20.85},
        references={"fallback": "NegCtrl"},
    )
    candidates, audit, _, denominators, summary = analyze(frame)
    assert candidates.element_id.tolist() == ["c", "a", "b"]
    assert candidates.minimum_rank_percentile.tolist() == [100, 50, 50]
    assert set(audit.element_id) == {"a", "b", "c"}
    assert denominators.n_measured_scrambled.eq(3).all()
    assert summary["n_above_threshold"] == 4 and summary["n_scrambled_eligible"] == 3


def test_raw_spread_off_levels_and_motif_correlation_are_separate_from_ranks():
    frame = measurements(
        {"a": [-2, -1, 0, 1, 2], "b": [0, 1, 2, 3, 4], "c": [2, 3, 4, 5, 6]},
        scores={"a": 21, "b": 22, "c": 23},
    )
    candidates, _, _, _, summary = analyze(frame)
    assert summary["spearman_min_percentile_fimo"] == pytest.approx(1)
    assert candidates.dynamic_range_spread.eq(4).all()
    assert candidates.dynamic_range_sd.tolist() == pytest.approx([np.sqrt(2.5)] * 3)
    assert candidates.max_off_activity.eq(2).all()
    a = candidates.set_index("element_id").loc["a"]
    assert a["dynamic_range__HepG2__En19"] == -2


def test_missing_raw_pair_cannot_be_counted_from_stale_dynamic_range():
    frame = measurements({"a": [1] * 5})
    frame.loc[0, "on_activity"] = np.nan
    frame.loc[1, "off_activity"] = np.nan
    candidates, audit, _, _, summary = analyze(frame)
    assert candidates.empty and summary["n_excluded_coverage"] == 1
    assert audit.n_contexts_measured.tolist() == [3]
    assert audit["dynamic_range__HepG2__En19"].isna().all()


@pytest.mark.parametrize("fault", ["duplicate", "score", "infinity", "range"])
def test_invalid_measurements_fail(fault):
    frame = measurements({"a": [1] * 5})
    if fault == "duplicate":
        frame = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    elif fault == "score":
        frame.loc[0, "fimo_score"] = 30
    elif fault == "infinity":
        frame.loc[0, "on_activity"] = np.inf
    else:
        frame.loc[0, "dynamic_range"] = 99
    with pytest.raises(ValueError):
        analyze(frame)


def test_reference_specific_ranks_do_not_compete_with_fallback():
    frame = measurements(
        {"a": [1] * 5, "b": [2] * 5, "fallback": [100] * 5},
        references={"fallback": "NegCtrl"},
    )
    frame.loc[frame.element_id.eq("b"), "separation"] = np.nan
    primary = reference_rankings(frame, "scrambled")
    fallback = reference_rankings(frame, "NegCtrl")
    assert primary[primary.element_id.eq("b")].rank_dynamic_range.eq(1).all()
    assert primary[primary.element_id.eq("b")].rank_separation.isna().all()
    assert primary[primary.element_id.eq("a")].rank_dynamic_range.eq(2).all()
    assert primary[primary.element_id.eq("a")].rank_separation.eq(1).all()
    assert fallback.rank_dynamic_range.eq(1).all()
