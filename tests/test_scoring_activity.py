"""Guard missingness, sign direction, and element-level diagnostic denominators."""

import numpy as np
import pandas as pd
import pytest

from src.scoring_validation import CELLS, aggregate, bins, sign_check, stats


def observations():
    return pd.DataFrame(
        {
            "element_id": ["a", "a", "b", "c", "d"],
            "project": ["motif"] * 5,
            "log2FoldChange": [-2.0, 2.0, -1.0, np.nan, 3.0],
            "padj": [0.01, 0.001, np.nan, np.nan, 0.001],
        }
    )


def test_calls_require_negative_activity_and_fdr_and_count_elements_once():
    result = aggregate(observations()).set_index("element_id")
    assert result.loc["a", "activity"] == 0
    assert result.loc["a", "diagnostic_call"]
    assert not result.loc["d", "diagnostic_call"]
    assert pd.isna(result.loc["b", "diagnostic_call"])
    assert pd.isna(result.loc["c", "diagnostic_call"])
    summary = stats(result)
    assert summary["n_elements"] == 4
    assert summary["n_measured"] == 3
    assert summary["n_testable"] == 2
    assert summary["diagnostic_fraction"] == 0.5
    assert np.isnan(summary["published_repressing_fraction"])


def test_bins_include_edges_and_keep_empty_bins_missing():
    frame = aggregate(observations())
    frame["score"] = [0.0, 1.0, 3.0, 3.5]
    table = bins(frame, np.arange(5))
    assert table.n_elements.tolist() == [1, 1, 0, 2]
    assert np.isnan(table.loc[2, "diagnostic_fraction"])
    assert np.isnan(table.loc[2, "median_activity"])


def test_bins_refuse_lost_maximum():
    frame = aggregate(observations()).assign(score=4.0)
    with pytest.raises(ValueError, match="lost or duplicated"):
        bins(frame, np.arange(5))


@pytest.mark.parametrize("bad", [1.0, np.nan])
def test_sign_gate_fails_on_inversion_or_missingness(bad):
    rows = [
        {"cell_type": cell, "population": pop, "median": value}
        for cell in CELLS.values()
        for pop, value in [("motif", -0.2), ("noMotif", 0.0)]
    ]
    sign = pd.DataFrame(rows)
    sign_check(sign)
    sign.loc[0, "median"] = bad
    with pytest.raises(ValueError, match="SIGN GATE FAILED"):
        sign_check(sign)


@pytest.mark.parametrize("field,bad", [("padj", -0.1), ("padj", 1.1), ("log2FoldChange", np.inf)])
def test_invalid_measurements_fail(field, bad):
    frame = observations()
    frame.loc[0, field] = bad
    with pytest.raises(ValueError):
        aggregate(frame)
