"""Training-only ranks, paired comparisons and chromosome/cluster leakage guards."""

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata, spearmanr

from src.heldout_validation import empirical_rank, evaluate, fit_predict, summarize


def test_training_midrank_and_unseen_test_tails():
    train = np.array([3, 1, 3, 6], dtype=float)
    np.testing.assert_allclose(empirical_rank(train, train), (rankdata(train) - 0.5) / len(train))
    np.testing.assert_allclose(
        empirical_rank(train, [-1, 1, 2, 3, 5, 6, 10]), [0, 0.125, 0.25, 0.5, 0.75, 0.875, 1]
    )
    with pytest.raises(ValueError, match="impute"):
        empirical_rank(train, [np.nan])


def sample():
    rng = np.random.default_rng(31)
    n = 80
    cases = pd.DataFrame(
        dict(
            element_id=[f"e{i:03d}" for i in range(n)],
            fimo_score=rng.uniform(21, 30, n),
            phyloP_motif_mean=rng.normal(size=n),
            dynamic_range=rng.normal(size=n),
            outcome="dynamic_range",
            cell_type="HepG2",
            enhancer="En19",
        )
    )
    assignments = pd.DataFrame(
        dict(
            element_id=cases.element_id,
            fold=np.arange(n) % 5,
            chrom=[f"chr{i % 5}" for i in range(n)],
            cluster_id=[f"c{i}" for i in range(n)],
            component_id=[f"chr{i % 5}" for i in range(n)],
        )
    )
    return cases, assignments


def test_test_labels_and_other_test_observations_do_not_change_fit_or_prediction():
    cases, _ = sample()
    train, test = cases.iloc[:60], cases.iloc[60:]
    pred, coef = fit_predict(train, test, "dynamic_range", True)
    changed = test.copy()
    changed.dynamic_range = 999
    changed.loc[changed.index[1:], "fimo_score"] = 999
    changed.loc[changed.index[1:], "phyloP_motif_mean"] = 999
    new_pred, new_coef = fit_predict(train, changed, "dynamic_range", True)
    assert coef == new_coef
    assert pred[0] == new_pred[0]
    baseline, baseline_coef = fit_predict(train, test, "dynamic_range", False)
    assert set(baseline_coef) == {"intercept", "fimo_score"}
    assert set(coef) == {"intercept", "fimo_score", "phyloP_motif_mean"}
    assert len(baseline) == len(test)


def test_exactly_once_out_of_fold_and_saved_scores_match_predictions():
    cases, assignments = sample()
    scores, predictions, coefficients = evaluate(cases, assignments)
    assert predictions.element_id.is_unique and len(predictions) == len(cases)
    assert len(scores) == 5 and len(coefficients) == 10
    for row in scores.itertuples():
        p = predictions[predictions.fold == row.fold]
        assert row.n_test == len(p) and row.n_train == len(cases) - len(p)
        for model in ("baseline", "augmented"):
            assert getattr(row, model + "_spearman") == pytest.approx(
                spearmanr(p.observed_outcome, p[model + "_predicted_rank"]).statistic
            )
        assert row.delta == pytest.approx(row.augmented_spearman - row.baseline_spearman)
    new = evaluate(cases.sample(frac=1, random_state=2), assignments)
    for left, right in zip((scores, predictions, coefficients), new, strict=True):
        pd.testing.assert_frame_equal(left, right)


@pytest.mark.parametrize("column", ["chrom", "cluster_id", "component_id"])
def test_leakage_rejected(column):
    cases, assignments = sample()
    assignments.loc[1, column] = assignments.loc[0, column]
    with pytest.raises(ValueError, match="Leakage"):
        evaluate(cases, assignments)


def test_missing_outcomes_not_imputed_and_duplicate_elements_rejected():
    cases, assignments = sample()
    with pytest.raises(ValueError, match="Duplicate element"):
        evaluate(pd.concat([cases, cases.iloc[[0]]]), assignments)
    cases.loc[0, "dynamic_range"] = np.nan
    with pytest.raises(ValueError, match="complete cases"):
        evaluate(cases, assignments)


def test_null_and_inconsistent_results_are_reported_without_selection():
    rows = pd.DataFrame(
        dict(
            outcome="dynamic_range",
            cell_type="HepG2",
            enhancer="En19",
            baseline_spearman=[0.5] * 5,
            augmented_spearman=[0.4, 0.6, 0.6, 0.6, 0.6],
        )
    )
    rows["delta"] = rows.augmented_spearman - rows.baseline_spearman
    rows["improved"] = rows.delta > 0
    summary = summarize(rows).iloc[0]
    assert summary.folds_improved == 4
    assert summary.verdict == "inconsistent improvement across folds"
    assert summary.delta_min == pytest.approx(-0.1)
    rows["delta"], rows["improved"] = 0.0, False
    assert summarize(rows).iloc[0].verdict == "no fold improved"
