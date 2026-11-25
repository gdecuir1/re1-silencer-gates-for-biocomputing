"""Context isolation, sign, missingness, and complete source-key preservation."""

import numpy as np
import pandas as pd
import pytest

from src.dataset import assert_source_rows, attach_metadata, characterize, matched_repression


def source_rows():
    return pd.DataFrame(
        [
            ("c1", "NegCtrl", "A", "E1", 2.0, 0.1),
            ("c2", "NegCtrl", "A", "E1", 4.0, 0.1),
            ("m", "motif", "A", "E1", 4.0, 0.9),
            ("missing", "motif", "A", "E1", np.nan, np.nan),
            ("c1", "NegCtrl", "A", "E2", 10.0, 0.1),
            ("c2", "NegCtrl", "A", "E2", 12.0, 0.1),
            ("m", "motif", "A", "E2", 4.0, 0.9),
            ("c1", "NegCtrl", "B", "E1", 100.0, 0.1),
            ("m", "motif", "B", "E1", 4.0, 0.9),
        ],
        columns=["element_id", "population", "cell_type", "enhancer", "raw_activity", "padj"],
    )


def test_baselines_stay_in_cell_and_enhancer_and_positive_raw_can_repress():
    source = source_rows()
    result = matched_repression(source)
    m = result[result.element_id.eq("m")]
    assert m.negctrl_baseline.tolist() == [3.0, 11.0, 100.0]
    assert m.repression.tolist() == [-1.0, 7.0, 96.0]
    assert m.padj.tolist() == [0.9, 0.9, 0.9]  # No significance filtering.
    missing = result[result.element_id.eq("missing")].iloc[0]
    assert pd.isna(missing.repression)
    assert missing.missing_raw_activity and missing.missing_repression
    assert_source_rows(source, result)
    assert len(result) == 9  # Missing E2/B observations are not synthesized.


def test_all_missing_controls_do_not_become_zero_baseline():
    source = source_rows()
    source.loc[source.population.eq("NegCtrl"), "raw_activity"] = np.nan
    result = matched_repression(source)
    assert result.negctrl_baseline.isna().all()
    assert result.repression.isna().all()


def test_context_without_control_is_retained_with_missing_baseline():
    source = source_rows()
    source = source[~(source.cell_type.eq("B") & source.population.eq("NegCtrl"))]
    result = matched_repression(source)
    assert_source_rows(source, result)
    assert result[result.cell_type.eq("B")].missing_negctrl_baseline.all()


def test_duplicate_keys_fail_instead_of_averaging():
    source = source_rows()
    with pytest.raises(ValueError, match="duplicate"):
        matched_repression(pd.concat([source, source.iloc[[0]]]))


def test_equal_counts_cannot_hide_changed_keys():
    source = source_rows()
    changed = source.copy()
    changed.loc[0, "element_id"] = "wrong"
    with pytest.raises(AssertionError, match="keys differ"):
        assert_source_rows(source, changed)


def test_dropped_rows_fail_count_assertion():
    with pytest.raises(AssertionError, match="Source rows"):
        assert_source_rows(source_rows(), source_rows().iloc[:-1])


def test_metadata_keeps_unverified_and_control_rows_and_nullable_threshold():
    source = source_rows()
    coordinates = pd.DataFrame(
        [("m", "hg38", "chr1", 100, 300, "+")],
        columns=["element_id", "working_build", "chrom", "start", "end", "strand"],
    )
    audit = pd.DataFrame(
        [("m", "lifted+verified"), ("missing", "lifted+mismatch")], columns=["element_id", "bucket"]
    )
    scores = pd.DataFrame([("m", 20.86), ("missing", 10.0)], columns=["element_id", "score"])
    result = attach_metadata(matched_repression(source), coordinates, audit, scores)
    assert_source_rows(source, result)
    assert result[result.element_id.eq("m")].above_m_value.all()
    assert result[result.element_id.eq("missing")].missing_hg38_coordinates.all()
    assert not result[result.element_id.eq("missing")].eligible_for_conservation.any()
    assert result[result.population.eq("NegCtrl")].above_m_value.isna().all()
    assert result[result.population.eq("NegCtrl")].fimo_score.isna().all()


def test_correlations_are_context_specific_with_no_significance_filter():
    rows = []
    for enhancer, activities in [("E1", [3.0, 2.0, 1.0]), ("E2", [1.0, 2.0, 3.0])]:
        for i, activity in enumerate(activities):
            rows.append((str(i), "motif", "A", enhancer, activity, 0.9, float(i)))
        rows.extend(
            [
                ("c", "NegCtrl", "A", enhancer, 4.0, 0.1, np.nan),
                ("n", "noMotif", "A", enhancer, 5.0, 0.1, np.nan),
            ]
        )
    source = pd.DataFrame(
        rows,
        columns=[
            "element_id",
            "population",
            "cell_type",
            "enhancer",
            "raw_activity",
            "padj",
            "fimo_score",
        ],
    )
    table = characterize(matched_repression(source)).set_index("enhancer")
    assert table.loc["E1", "spearman_rho"] == 1.0
    assert table.loc["E2", "spearman_rho"] == -1.0
    assert table.loc["E1", "recommendation"] == "primary_within_cell"
    assert table.loc["E2", "recommendation"] == "do_not_prioritize"
    assert table.n_correlation_pairs.tolist() == [3, 3]
    assert table.motif_vs_noMotif_gap.tolist() == [3.0, 3.0]
