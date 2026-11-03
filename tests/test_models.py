"""Test mechanical evaluation and ensure all scientific stubs fail explicitly."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import metrics, models


def test_evaluation_known_ordering() -> None:
    """Perfectly ordered synthetic predictions have known rank and ROC metrics."""
    table = pd.DataFrame(
        {
            "element_id": [f"synthetic_{i}" for i in range(6)],
            "cell_type": ["K562"] * 6,
            "fold": [0, 0, 0, 1, 1, 1],
            "observed": [0, 1, 2, 0, 1, 2],
            "prediction": [0, 1, 2, 2, 1, 0],
            "positive": [0, 0, 1, 0, 0, 1],
        }
    )
    result = models.evaluate_folds(table)
    assert result["spearman"].tolist() == [1.0, -1.0]
    assert result["auroc"].tolist() == [1.0, 0.0]


def test_undefined_metrics_are_missing() -> None:
    """Do not invent scores for constant values or single-class folds."""
    table = pd.DataFrame(
        {
            "element_id": ["synthetic_a", "synthetic_b"],
            "cell_type": ["K562"] * 2,
            "fold": [0, 0],
            "observed": [1, 1],
            "prediction": [1, 1],
            "positive": [0, 0],
        }
    )
    result = models.evaluate_folds(table)
    assert np.isnan(result.loc[0, "spearman"])
    assert np.isnan(result.loc[0, "auroc"])
    assert result.loc[0, "auroc_status"] == "single_class"


def test_stubs_raise() -> None:
    """Every scientific placeholder raises rather than producing plausible output."""
    series, frame = pd.Series(dtype=float), pd.DataFrame()
    for function, args in [
        (models.baseline, (series, series, series)),
        (models.motif_conservation, (frame, series, series)),
        (models.motif_cofactor, (frame, series, series)),
        (models.gapped_kmer_svm, (Path("synthetic_unused.fa"), frame, frame)),
        (metrics.dynamic_range, (series, series)),
        (metrics.off_state_leakiness, (series, series)),
        (metrics.separation, (frame, frame)),
        (metrics.variant_robustness, (series, frame)),
    ]:
        with pytest.raises(NotImplementedError):
            function(*args)
