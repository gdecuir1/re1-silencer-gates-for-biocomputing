"""Corrected gate population, offline-only replay and secondary feature contract."""

import json
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest
from test_conservation_pipeline import example_features

from src import config, conservation_io
from src.conservation_pipeline import join_features, sanity_checks
from src.conservation_strata import corrected_checks
from src.dataset import digest


def metadata():
    return pd.DataFrame(
        dict(
            element_id=["a", "b", "c"],
            population=["motif", "motif", "noMotif"],
            fimo_score=[20.86, 20.85, np.nan],
        )
    )


def test_corrected_gate_stratifies_exact_threshold_and_ignores_phastcons_failure():
    f = example_features()
    f.loc[f.element_id.eq("b"), "phyloP_motif_mean"] = -20
    f.loc[f.population.eq("motif"), "phastCons_motif_mean"] = 0
    paired, _, old = sanity_checks(f)
    # Use an independently supplied passing original paired check, held fixed.
    paired.loc[
        paired.track.eq("phyloP") & paired.population.eq("motif"),
        ["median_difference", "fraction_motif_greater"],
    ] = [0.1117, 0.571]
    checks, observed, summary = corrected_checks(f, metadata(), paired)
    assert not old.passed.all() and checks.passed.all()
    assert len(checks) == 2 and checks.track.eq("phyloP").all()
    assert observed.set_index("element_id").loc["a", "stratum"] == "canonical_above_m_value"
    assert observed.set_index("element_id").loc["b", "stratum"] == "canonical_below_m_value"
    assert summary.query('scope == "complete"').n.tolist() == [1, 1, 1]
    assert summary.query('stratum == "canonical_below_m_value"')["median"].tolist() == [-20, -20]


def test_corrected_population_failure_still_blocks_join():
    f = example_features()
    f.loc[f.element_id.eq("a"), "phyloP_motif_mean"] = -5
    checks, _, _ = corrected_checks(f, metadata(), sanity_checks(f)[0])
    assert not checks.iloc[1].passed
    with pytest.raises(ValueError, match="join prohibited"):
        join_features(pd.DataFrame(), f, checks)


def test_missing_coordinates_retained_but_missing_fimo_fails():
    f = example_features()
    m = pd.concat(
        [
            metadata(),
            pd.DataFrame(dict(element_id=["unverified"], population=["motif"], fimo_score=[30.0])),
        ],
        ignore_index=True,
    )
    checks, observed, summary = corrected_checks(f, m, sanity_checks(f)[0])
    row = observed.set_index("element_id").loc["unverified"]
    assert row.missing_verified_coordinates and pd.isna(row.phyloP_motif_mean)
    high = summary.query('scope == "complete" and stratum == "canonical_above_m_value"').iloc[0]
    assert high.source_elements == 2 and high.n == 1 and high.missing_coordinates == 1
    m.loc[m.element_id.eq("a"), "fimo_score"] = np.nan
    with pytest.raises(ValueError, match="FIMO score missing"):
        corrected_checks(f, m, sanity_checks(f)[0])


def test_join_preserves_contexts_and_only_includes_full_element_phastcons():
    f = example_features()
    checks = corrected_checks(f, metadata(), sanity_checks(f)[0])[0]
    dataset = metadata().merge(
        pd.DataFrame({"cell_type": ["C", "C"], "enhancer": ["E1", "E2"]}), how="cross"
    )
    dataset["repression"] = [1, np.nan, 3, 4, 5, 6]
    joined = join_features(dataset, f, checks)
    pd.testing.assert_frame_equal(joined[dataset.columns], dataset)
    assert "phastCons_element_mean" in joined and "phyloP_motif_mean" in joined
    assert "phastCons_motif_mean" not in joined and "phastCons_flank_mean" not in joined


def test_offline_loader_refuses_missing_cache_without_network(tmp_path, monkeypatch):
    elements = pd.DataFrame(
        dict(element_id=["a"], chrom=["chr1"], start=[0], end=[200], strand=["+"])
    )
    coordinates = tmp_path / "elements.csv"
    elements.to_csv(coordinates, index=False)
    monkeypatch.setattr(config, "ROOT", tmp_path)
    monkeypatch.setattr(config, "HG38_ELEMENTS", coordinates)
    urls = {"phyloP": "http://example.invalid/p", "phastCons": "http://example.invalid/c"}
    monkeypatch.setattr(config, "BIGWIG_URLS", urls)
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            dict(
                tracks={k: dict(url=v) for k, v in urls.items()},
                inputs={
                    "elements.csv": digest(coordinates),
                    "data/processed/conservation_windows/missing.npy": "abc",
                },
            )
        )
    )
    network = Mock(side_effect=AssertionError("Network must not be called"))
    monkeypatch.setattr(conservation_io, "urlopen", network)
    monkeypatch.setattr(conservation_io.pyBigWig, "open", network)
    with pytest.raises(ValueError, match="remote reads prohibited"):
        conservation_io.load_cached(elements, manifest)
    network.assert_not_called()
