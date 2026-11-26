"""Conservation coverage, orientation, gate enforcement and rank-adjusted preview."""

from unittest.mock import MagicMock, Mock

import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

from src import conservation_io
from src.conservation_pipeline import (
    coverage_table,
    join_features,
    partial_spearman,
    sanity_checks,
    summarize,
    window_summary,
)
from src.conservation_strata import corrected_checks


def elements():
    return pd.DataFrame(
        dict(
            element_id=["a", "b", "c"],
            project=["motif", "motif", "noMotif"],
            chrom=["chrSynthetic"] * 3,
            start=[100, 400, 700],
            end=[300, 600, 900],
            strand=["+", "-", "+"],
        )
    )


def example_features():
    # a and b have the same oriented motif; reverse-strand b maps to offsets 89..109.
    phylo = np.full((3, 200), -1.0)
    phylo[0, 90:111] = 2.0
    phylo[1, 89:110] = 2.0
    phast = np.full((3, 200), 0.1)
    phast[0, 90:111] = 0.9
    phast[1, 89:110] = 0.9
    arrays = dict(phyloP=phylo, phastCons=phast)
    absent = {track: np.zeros(3, dtype=bool) for track in arrays}
    return summarize(elements(), arrays, absent)


def test_partial_missing_and_negative_scores_have_explicit_denominators():
    values = np.array([[-2, 0, 2, np.nan], [np.nan] * 4, [-4, -3, -2, -1]])
    stats = window_summary(values)
    np.testing.assert_allclose(stats["mean"], [0, np.nan, -2.5], equal_nan=True)
    np.testing.assert_allclose(stats["max"], [2, np.nan, -1], equal_nan=True)
    np.testing.assert_allclose(
        stats["positive_fraction_observed"], [1 / 3, np.nan, 0], equal_nan=True
    )
    assert stats["coverage_status"].tolist() == ["partial", "no_data", "complete"]
    assert stats["covered_bases"].tolist() == [3, 0, 4]
    with pytest.raises(ValueError, match="finite-or-missing"):
        window_summary(np.array([[np.inf]]))


def test_both_strands_and_complementary_flanks_are_exact():
    f = example_features().set_index("element_id")
    assert f.phyloP_motif_mean.tolist() == [2, 2, -1]
    assert f.phyloP_flank_mean.eq(-1).all()
    assert f.phyloP_motif_requested_bases.eq(21).all()
    assert f.phyloP_flank_requested_bases.eq(179).all()
    assert f.loc["b", "motif_start"] == 489 and f.loc["b", "motif_end"] == 510
    assert f.loc["a", "phyloP_element_mean"] == pytest.approx((42 - 179) / 200)
    assert f.loc["a", "phyloP_element_positive_fraction_observed"] == pytest.approx(21 / 200)


def test_sanity_checks_require_both_tracks_and_both_expectations():
    f = example_features()
    paired, distributions, checks = sanity_checks(f)
    assert len(checks) == 4 and checks.passed.all()
    assert paired.query('population == "motif"').fraction_motif_greater.eq(1).all()
    assert len(distributions) == 8
    f.loc[f.population.eq("motif"), "phastCons_motif_mean"] = 0
    _, _, failed = sanity_checks(f)
    assert failed[failed.track.eq("phyloP")].passed.all()
    assert not failed[failed.track.eq("phastCons")].passed.any()


def test_partial_coverage_reported_but_excluded_from_primary_checks():
    f = example_features()
    f.loc[0, "phyloP_motif_coverage_status"] = "partial"
    f.loc[0, "phyloP_element_coverage_status"] = "partial"
    coverage = coverage_table(f)
    row = coverage.query('track == "phyloP" and window == "motif" and population == "all"').iloc[0]
    assert row.complete == 2 and row.partial == 1 and row.no_data == 0
    paired, distributions, checks = sanity_checks(f)
    canonical = paired.query('track == "phyloP" and population == "motif"').iloc[0]
    assert canonical.n == 1 and canonical.excluded_partial_or_missing == 1
    assert checks.passed.all()
    assert distributions.query('track == "phyloP" and population == "motif"').n.tolist() == [1, 2]


def test_join_is_blocked_on_failure_and_retains_every_context_and_missing_activity():
    f = example_features()
    checks = corrected_checks(
        f,
        pd.DataFrame(
            {
                "element_id": ["a", "b", "c"],
                "population": ["motif", "motif", "noMotif"],
                "fimo_score": [25.0, 20.0, np.nan],
            }
        ),
        sanity_checks(f)[0],
    )[0]
    dataset = pd.DataFrame(
        dict(
            element_id=["a", "a", "control"],
            cell_type=["C"] * 3,
            enhancer=["E1", "E2", "E1"],
            repression=[2, np.nan, 0.1],
        )
    )
    joined = join_features(dataset, f, checks)
    assert len(joined) == 3
    pd.testing.assert_frame_equal(joined[dataset.columns], dataset)
    assert joined.missing_conservation_coordinates.tolist() == [False, False, True]
    assert pd.isna(joined.loc[2, "phyloP_motif_mean"])
    checks.loc[0, "passed"] = False
    with pytest.raises(ValueError, match="join prohibited"):
        join_features(dataset, f, checks)


def test_duplicate_feature_join_fails_instead_of_expanding_rows():
    f = example_features()
    checks = corrected_checks(
        f,
        pd.DataFrame(
            {
                "element_id": ["a", "b", "c"],
                "population": ["motif", "motif", "noMotif"],
                "fimo_score": [25.0, 20.0, np.nan],
            }
        ),
        sanity_checks(f)[0],
    )[0]
    dataset = pd.DataFrame(dict(element_id=["a"], cell_type=["C"], enhancer=["E"], repression=[1]))
    with pytest.raises(ValueError, match="Duplicate"):
        join_features(dataset, pd.concat([f, f.iloc[[0]]]), checks)


def test_partial_spearman_matches_rank_correlation_identity_with_ties():
    x = [1, 1, 2, 3, 4, 5, 8, 9]
    y = [2, 5, 3, 6, 9, 4, 7, 8]
    z = [5, 2, 3, 1, 6, 7, 8, 9]
    r = spearmanr(np.column_stack([x, y, z])).statistic
    expected = (r[0, 1] - r[0, 2] * r[1, 2]) / np.sqrt((1 - r[0, 2] ** 2) * (1 - r[1, 2] ** 2))
    assert partial_spearman(x, y, z) == pytest.approx(expected)
    assert np.isnan(partial_spearman(z, y, z))
    assert np.isnan(partial_spearman([1, np.nan, 2, 3], [1, 2, 3, 4], [2, 3, 4, 1]))


def test_production_range_probe_refuses_full_response_without_consuming_body(monkeypatch):
    response = MagicMock()
    response.__enter__.return_value = response
    response.status = 200
    response.headers = {}
    monkeypatch.setattr(conservation_io, "urlopen", Mock(return_value=response))
    with pytest.raises(ValueError, match="refusing full"):
        conservation_io.track_identity("https://example.invalid/track.bw")
    response.read.assert_not_called()


def test_cached_windows_are_hash_checked_and_never_requeried(tmp_path, monkeypatch):
    handle = MagicMock()
    handle.__enter__.return_value = handle
    handle.isBigWig.return_value = True
    handle.chroms.return_value = {"chrSynthetic": 1000}
    handle.values.return_value = [-1.0] * 200
    opened = Mock(return_value=handle)
    monkeypatch.setattr(conservation_io.pyBigWig, "open", opened)
    monkeypatch.setattr(conservation_io.pyBigWig, "remote", 1)
    group = elements().iloc[:1]
    identity = dict(url="https://example.invalid/p.bw", etag="fixed")
    first = conservation_io.read_chunk("phyloP", identity, group, tmp_path)
    second = conservation_io.read_chunk("phyloP", identity, group, tmp_path)
    np.testing.assert_array_equal(first[1], second[1])
    assert opened.call_count == 1
    assert handle.values.call_args.args == ("chrSynthetic", 100, 300)
    first[3].write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="identity mismatch"):
        conservation_io.read_chunk("phyloP", identity, group, tmp_path)


def test_track_absent_chromosome_is_flagged_unknown_not_zero(tmp_path, monkeypatch):
    handle = MagicMock()
    handle.__enter__.return_value = handle
    handle.isBigWig.return_value = True
    handle.chroms.return_value = {"chrOther": 1000}
    monkeypatch.setattr(conservation_io.pyBigWig, "open", Mock(return_value=handle))
    monkeypatch.setattr(conservation_io.pyBigWig, "remote", 1)
    result = conservation_io.read_chunk(
        "phyloP", dict(url="https://example.invalid/p.bw"), elements().iloc[:1], tmp_path
    )
    assert result[2] and np.isnan(result[1]).all()
    handle.values.assert_not_called()


def test_transient_read_error_retries_same_window_without_neutral_imputation(tmp_path, monkeypatch):
    handles = [MagicMock(), MagicMock()]
    for handle in handles:
        handle.__enter__.return_value = handle
        handle.isBigWig.return_value = True
        handle.chroms.return_value = {"chrSynthetic": 1000}
    handles[0].values.side_effect = [[-1.0] * 200, RuntimeError("transient network read")]
    handles[1].values.side_effect = [[2.0] * 200, [3.0] * 200]
    monkeypatch.setattr(conservation_io.pyBigWig, "open", Mock(side_effect=handles))
    monkeypatch.setattr(conservation_io.pyBigWig, "remote", 1)
    monkeypatch.setattr(conservation_io.time, "sleep", Mock())
    result = conservation_io.read_chunk(
        "phyloP", dict(url="https://example.invalid/p.bw"), elements(), tmp_path
    )
    assert result[1][:, 0].tolist() == [-1, 2, 3]
    assert handles[1].values.call_args_list[0].args == ("chrSynthetic", 400, 600)
    assert handles[1].values.call_args_list[1].args == ("chrSynthetic", 700, 900)


def test_transport_reuse_requires_identical_file_identity_and_preserves_receipt(
    tmp_path, monkeypatch
):
    import json

    from src import config
    from src.dataset import digest

    monkeypatch.setattr(config, "ROOT", tmp_path)
    cache = tmp_path / "cache"
    group = elements()
    handle = MagicMock()
    handle.__enter__.return_value = handle
    handle.isBigWig.return_value = True
    handle.chroms.return_value = {"chrSynthetic": 1000}
    handle.values.return_value = [2.0] * 200
    monkeypatch.setattr(conservation_io.pyBigWig, "remote", 1)
    monkeypatch.setattr(conservation_io.pyBigWig, "open", Mock(return_value=handle))
    old = dict(
        url="https://example.invalid/p.bw",
        resolved_url="https://example.invalid/p.bw",
        etag="file1",
        last_modified="fixed",
        total_bytes=1000,
        range_status=206,
    )
    original = conservation_io.read_chunk("phyloP", old, group, cache)
    new = dict(old, url="http://example.invalid/p.bw", resolved_url="http://example.invalid/p.bw")
    monkeypatch.setattr(config, "BIGWIG_URLS", {"phyloP": new["url"]})
    monkeypatch.setattr(conservation_io, "track_identity", Mock(return_value=new))
    conservation_io.reuse_transport_cache(group, cache)
    loaded = conservation_io.read_chunk("phyloP", new, group, cache)
    assert digest(loaded[3]) == digest(original[3])
    receipt = json.loads(loaded[4].read_text())
    assert receipt["original_receipt_sha256"] == digest(original[4])
    assert conservation_io.pyBigWig.open.call_count == 1
    monkeypatch.setattr(
        conservation_io, "track_identity", Mock(return_value=dict(new, etag="changed"))
    )
    with pytest.raises(ValueError, match="identities disagree"):
        conservation_io.reuse_transport_cache(group, cache)
