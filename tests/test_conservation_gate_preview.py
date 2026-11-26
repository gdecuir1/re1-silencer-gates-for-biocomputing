"""Reference exclusions, unchanged robustness ranks and same-case gate correlations."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from src.conservation_gate_preview import analyze, correlate
from src.conservation_pipeline import partial_spearman
from src.cross_context import CONTEXTS


def fixture():
    rows = []
    for i in range(10):
        for j, (cell, enhancer) in enumerate(CONTEXTS):
            off = float((i * 3 + j) % 5)
            dr = float((i * 7 + j) % 11)
            rows.append(
                dict(
                    element_id=f"e{i}",
                    cell_type=cell,
                    enhancer=enhancer,
                    population="motif",
                    reference_type="scrambled",
                    fimo_score=18.0 + i,
                    on_activity=off + dr,
                    off_activity=off,
                    dynamic_range=dr,
                )
            )
    frame = pd.DataFrame(rows)
    features = pd.DataFrame(
        dict(
            element_id=[f"e{i}" for i in range(10)],
            phyloP_motif_mean=[float((i * 3) % 7) for i in range(10)],
            phyloP_motif_coverage_status=["complete"] * 10,
        )
    )
    candidates = pd.DataFrame(
        dict(
            element_id=[f"e{i}" for i in range(3, 10)],
            fimo_score=[18.0 + i for i in range(3, 10)],
            n_contexts_measured=[4] * 7,
            minimum_rank_percentile=[10, 30, 80, 25, 95, 45, 50],
        )
    )
    return features, frame, candidates


def test_reference_and_score_cohorts_are_explicit_and_robustness_is_once_per_element():
    features, frame, candidates = fixture()
    fallback = frame[frame.element_id.eq("e0")].copy()
    fallback.element_id = "fallback"
    fallback.reference_type = "NegCtrl"
    fallback.dynamic_range = 1000
    fallback.on_activity = fallback.off_activity + 1000
    table, cases, _ = analyze(features, pd.concat([frame, fallback]), candidates, partial_spearman)
    assert table.query('cohort == "above_m_value"').n_source.eq(7).all()
    assert table.query('cohort == "all_canonical_scores"').n_source.eq(10).all()
    robust = cases[cases.outcome.eq("minimum_rank_percentile")]
    assert len(robust) == 7 and robust.element_id.is_unique
    np.testing.assert_array_equal(
        robust.minimum_rank_percentile, candidates.minimum_rank_percentile
    )
    assert "fallback" not in set(cases.element_id)


def test_missing_conservation_excludes_cases_without_reranking_and_missing_on_is_not_fallback():
    features, frame, candidates = fixture()
    features.loc[features.element_id.eq("e3"), "phyloP_motif_coverage_status"] = "partial"
    frame.loc[(frame.element_id == "e4") & (frame.cell_type == "HepG2"), "on_activity"] = np.nan
    table, cases, _ = analyze(features, frame, candidates, partial_spearman)
    robust = table.query('outcome == "minimum_rank_percentile"').iloc[0]
    assert robust.n_complete_cases == 6 and robust.n_excluded == 1
    dr = table.query('cohort == "above_m_value" and cell_type == "HepG2"').iloc[0]
    assert dr.n_complete_cases == 5 and dr.missing_pair == 1
    saved = cases.query('outcome == "minimum_rank_percentile"').set_index("element_id")
    assert saved.loc["e3", "minimum_rank_percentile"] == 10 and not saved.loc["e3", "complete_case"]


def test_same_case_rank_adjustment_matches_independent_correlation_formula():
    group = pd.DataFrame(
        dict(
            phyloP_motif_mean=[2, 4, 1, 1, 5, 8, np.nan, 10],
            minimum_rank_percentile=[10, 40, 20, 50, 30, 80, 90, np.nan],
            fimo_score=[21, 24, 22, 26, 23, 25, 28, 27],
            phyloP_motif_coverage_status=["complete"] * 8,
        )
    )
    result, cases = correlate(group, "minimum_rank_percentile", partial_spearman)
    x = group.iloc[:6][["phyloP_motif_mean", "minimum_rank_percentile", "fimo_score"]]
    r = spearmanr(x).statistic
    expected = (r[0, 1] - r[0, 2] * r[1, 2]) / np.sqrt((1 - r[0, 2] ** 2) * (1 - r[1, 2] ** 2))
    assert result["n_complete_cases"] == 6 and cases.complete_case.sum() == 6
    assert result["spearman_phyloP_outcome"] == pytest.approx(r[0, 1])
    assert result["partial_spearman_given_fimo"] == pytest.approx(expected)
    assert result["incremental_rank_r_squared"] == pytest.approx((1 - r[1, 2] ** 2) * expected**2)


def test_duplicate_robustness_rows_and_source_metric_disagreement_fail():
    features, frame, candidates = fixture()
    with pytest.raises(ValueError, match="one value per element"):
        analyze(features, frame, pd.concat([candidates, candidates.iloc[[0]]]), partial_spearman)
    frame.loc[frame.element_id.eq("e3"), "dynamic_range"] = 100
    with pytest.raises(ValueError, match="Dynamic range disagrees"):
        analyze(features, frame, candidates, partial_spearman)
