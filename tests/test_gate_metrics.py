"""Matched states, uncertainty, missingness and operating-region selection."""

import numpy as np
import pandas as pd
import pytest

from src import metrics
from src.gate_metrics import (
    assemble,
    control_baselines,
    operating_regions,
    pair_catalog,
    rank_pairs,
)


def observations():
    rows = []
    for enhancer, control, off, on in [("E1", 4.0, 1.0, 3.0), ("E2", 9.0, 2.0, np.nan)]:
        for element, pop, activity, se in [
            ("x_Ref", "motif", off, 0.3),
            ("y_Ref", "motif", 2.0, 0.4),
            ("x_Scr", "scrambled,SNP", on, 0.4),
            ("control", "NegCtrl", control, 0.5),
        ]:
            rows.append(
                (element, "cell", enhancer, pop, activity, se, 0.9, enhancer + "^" + element)
            )
    return pd.DataFrame(
        rows,
        columns=[
            "element_id",
            "cell_type",
            "enhancer",
            "population",
            "raw_activity",
            "lfcSE",
            "padj",
            "source_construct_id",
        ],
    )


def catalog():
    return pd.DataFrame(
        {
            "element_id": ["x_Ref", "y_Ref"],
            "scrambled_id": ["x_Scr", None],
            "has_designed_partner": [True, False],
        }
    )


def built():
    source = observations()
    scores = pd.DataFrame({"element_id": ["x_Ref", "y_Ref"], "score": [20.86, 20.85]})
    return assemble(source, catalog(), control_baselines(source), scores)


def test_matched_reference_and_missing_partner_do_not_silently_switch():
    f = built().set_index(["element_id", "enhancer"])
    assert len(f) == 4
    assert f.loc[("x_Ref", "E1"), "dynamic_range"] == 2.0
    assert f.loc[("x_Ref", "E1"), "leakiness"] == 1.0
    assert f.loc[("x_Ref", "E1"), "reference_type"] == "scrambled"
    assert f.loc[("x_Ref", "E2"), "reference_reason"] == "partner_activity_missing"
    for column in ("on_activity", "dynamic_range", "leakiness", "separation"):
        assert pd.isna(f.loc[("x_Ref", "E2"), column])
    assert f.loc[("y_Ref", "E1"), "reference_type"] == "NegCtrl"
    assert f.loc[("y_Ref", "E1"), "dynamic_range"] == 2.0
    assert f.loc[("y_Ref", "E2"), "dynamic_range"] == 7.0
    assert f.loc[("x_Ref", "E1"), "above_m_value"]
    assert not f.loc[("y_Ref", "E1"), "above_m_value"]
    assert f.off_padj.eq(0.9).all()  # Supplied significance never filters metrics.


def test_missing_context_retains_canonical_row_with_explicit_reason():
    source = observations()
    source = source[~(source.element_id.eq("x_Scr") & source.enhancer.eq("E2"))]
    result = assemble(
        source,
        catalog(),
        control_baselines(source),
        pd.DataFrame({"element_id": ["x_Ref", "y_Ref"], "score": [22.0, 23.0]}),
    )
    row = result[result.element_id.eq("x_Ref") & result.enhancer.eq("E2")].iloc[0]
    assert row.reference_reason == "partner_context_absent"
    assert row.reference_type == "scrambled" and pd.isna(row.dynamic_range)


def test_se_separation_is_rms_uncertainty_and_zero_spread_stays_missing():
    on = pd.DataFrame(
        {"activity": [3.0, 1.0, np.nan, 4.0], "standard_error": [0.4, 0.0, 0.4, np.nan]}
    )
    off = pd.DataFrame({"activity": [1.0, 1.0, -1.0, 2.0], "standard_error": [0.3, 0.0, 0.3, 0.4]})
    result = metrics.separation(on, off)
    assert result.iloc[0] == pytest.approx(2 / np.sqrt((0.4**2 + 0.3**2) / 2))
    assert result.iloc[1:].isna().all()
    assert metrics.dynamic_range(on.activity, off.activity).iloc[1] == 0


def test_negative_off_and_negative_range_are_not_clipped():
    on, off = pd.Series([-2.0]), pd.Series([-1.0])
    assert metrics.dynamic_range(on, off).iloc[0] == -1
    assert metrics.off_state_leakiness(off, on).iloc[0] == -1


def test_metric_inputs_do_not_silently_reindex_or_accept_negative_se():
    with pytest.raises(ValueError, match="identical unique"):
        metrics.dynamic_range(pd.Series([1.0], index=["a"]), pd.Series([1.0], index=["b"]))
    with pytest.raises(ValueError, match="finite or missing"):
        metrics.dynamic_range(pd.Series([np.inf]), pd.Series([1.0]))
    on = pd.DataFrame({"activity": [2.0], "standard_error": [-0.1]})
    with pytest.raises(ValueError, match="negative"):
        metrics.separation(on, on)


def test_duplicate_observations_fail():
    source = observations()
    with pytest.raises(ValueError, match="Duplicate source"):
        assemble(
            pd.concat([source, source.iloc[[0]]]),
            catalog(),
            control_baselines(source),
            pd.DataFrame(),
        )


def test_sequence_pairing_checks_flanks_and_allows_measured_adjacent_edit():
    a = "A" * 200
    b = a[:90] + "C" * 22 + a[112:]
    oligos = pd.DataFrame(
        {"oligo_name": ["x_Ref", "x_Scr"], "sequence": [a, b], "library": ["WG_REST"] * 2}
    )
    pairs, cores = pair_catalog(observations(), oligos)
    x = pairs.set_index("element_id").loc["x_Ref"]
    assert x.has_designed_partner and x.edited_adjacent_base
    assert x.changed_bases == 22 and cores["x_Scr"] == "C" * 21
    oligos.loc[1, "sequence"] = "G" + b[1:]
    with pytest.raises(ValueError, match="outside"):
        pair_catalog(observations(), oligos)


def test_bootstrap_is_context_local_and_order_independent():
    s = observations()
    extra = s[s.population.eq("NegCtrl")].copy()
    extra.element_id = "control2"
    extra.raw_activity += 2
    s = pd.concat([s, extra], ignore_index=True)
    a, b = control_baselines(s), control_baselines(s.iloc[::-1])
    pd.testing.assert_frame_equal(a, b)
    assert a.negctrl_baseline.tolist() == [5.0, 10.0]
    assert a.negctrl_median_se.gt(0).all()


def test_operating_medians_use_complete_high_score_pairs_and_separate_contexts():
    f = built()
    table, _ = operating_regions(f)
    e1, e2 = table.set_index("enhancer").loc["E1"], table.set_index("enhancer").loc["E2"]
    assert e1.n_pairs == 1 and e1.median_on == 3.0 and e1.median_dynamic_range == 2.0
    assert e1.usable
    assert e2.n_pairs == 0 and e2.verdict == "missing_pairs"
    assert pd.isna(e2.median_on)


def test_rank_agreement_is_observed_and_missing_se_is_not_ranked_zero():
    rows = []
    for i in range(30):
        rows.append(
            dict(
                element_id=f"e{i:02}",
                cell_type="A",
                enhancer="E1",
                above_m_value=True,
                dynamic_range=float(i),
                leakiness=float(i),
                separation=float(30 - i) if i != 29 else np.nan,
                reference_type="scrambled",
            )
        )
    f = pd.DataFrame(rows)
    region = pd.DataFrame({"cell_type": ["A"], "enhancer": ["E1"], "usable": [True]})
    ranks, agreement = rank_pairs(f, region)
    assert ranks.iloc[0].element_id == "e29" and ranks.iloc[0].rank_dynamic_range == 1
    assert pd.isna(ranks.iloc[0].rank_separation)
    a = agreement[agreement.reference_type.eq("all_explicit_references")].iloc[0]
    assert a.overlap_dr_separation == 20 and a.common_eligible_overlap == 21
    assert a.spearman_dr_separation == pytest.approx(-1)
    empty = agreement[agreement.reference_type.eq("NegCtrl")].iloc[0]
    assert empty.agreement == "unavailable" and pd.isna(empty.overlap_fraction)


def test_no_usable_context_retains_empty_ranking_schema():
    f = built()
    regions, _ = operating_regions(f)
    regions["usable"] = False
    ranking, agreement = rank_pairs(f, regions)
    assert ranking.empty and agreement.empty
    assert "rank_dynamic_range" in ranking
