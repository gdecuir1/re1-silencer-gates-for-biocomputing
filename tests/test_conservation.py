"""Test exact windows, remote range safeguards, and missing coverage without HTTP."""

from unittest.mock import MagicMock, Mock

import numpy as np
import pandas as pd
import pytest

from src import conservation

URLS = {
    "phyloP": "https://example.invalid/synthetic-phyloP.bw",
    "phastCons": "https://example.invalid/synthetic-phastCons.bw",
}


def test_motif_coordinates_both_strands(intervals: pd.DataFrame) -> None:
    """Inclusive offsets 90..110 become exactly 21 queried bases on either strand."""
    result = conservation.motif_windows(intervals)
    assert result["start"].tolist() == [190, 489]
    assert result["end"].tolist() == [211, 510]
    assert (result["end"] - result["start"]).eq(21).all()


def test_remote_exact_values_and_missingness(
    intervals: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Query only motif bases and preserve missing coverage instead of filling zero."""
    handles = []
    for _ in range(2):
        handle = Mock()
        handle.chroms.return_value = {"chrSynthetic": 1000}
        handle.isBigWig.return_value = True
        handle.values.side_effect = [[2.0] * 20 + [np.nan], [np.nan] * 21]
        handles.append(handle)
    check = Mock()
    opened = Mock(side_effect=handles)
    monkeypatch.setattr(conservation, "_check_byte_ranges", check)
    monkeypatch.setattr(conservation.pyBigWig, "remote", 1)
    monkeypatch.setattr(conservation.pyBigWig, "open", opened)
    result = conservation.query_conservation(intervals, urls=URLS)
    assert result["phyloP_covered_bases"].tolist() == [20, 0]
    assert result.loc[0, "phyloP_mean"] == 2.0
    assert np.isnan(result.loc[1, "phyloP_mean"])
    for handle in handles:
        assert handle.values.call_args_list[0].args == ("chrSynthetic", 190, 211)
        assert handle.values.call_args_list[1].args == ("chrSynthetic", 489, 510)
        handle.close.assert_called_once()
    assert check.call_count == 2
    assert [call.args[0] for call in opened.call_args_list] == list(URLS.values())


def test_range_server_rejection(monkeypatch: pytest.MonkeyPatch) -> None:
    """Close a server response that ignores Range without reading its body."""
    response = MagicMock()
    response.__enter__.return_value = response
    response.status = 200
    opened = Mock(return_value=response)
    monkeypatch.setattr(conservation, "urlopen", opened)
    with pytest.raises(ValueError, match="refusing full download"):
        conservation._check_byte_ranges(URLS["phyloP"])
    response.read.assert_not_called()
    response.__exit__.assert_called_once()
    assert opened.call_args.args[0].get_header("Range") == "bytes=0-0"


def test_range_server_acceptance(monkeypatch: pytest.MonkeyPatch) -> None:
    """Read just one byte from a valid partial-content probe."""
    response = MagicMock()
    response.__enter__.return_value = response
    response.status = 206
    response.headers = {"Content-Range": "bytes 0-0/100"}
    monkeypatch.setattr(conservation, "urlopen", Mock(return_value=response))
    conservation._check_byte_ranges(URLS["phyloP"])
    response.read.assert_called_once_with(1)


def test_missing_chromosome_closes_track(
    intervals: pd.DataFrame, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Treat mismatched assemblies as errors while releasing the open handle."""
    handle = Mock()
    handle.chroms.return_value = {"otherSyntheticChrom": 1000}
    monkeypatch.setattr(conservation, "_check_byte_ranges", Mock())
    monkeypatch.setattr(conservation.pyBigWig, "remote", 1)
    monkeypatch.setattr(conservation.pyBigWig, "open", Mock(return_value=handle))
    with pytest.raises(ValueError, match="invalid chromosome"):
        conservation.query_conservation(intervals, urls=URLS)
    handle.close.assert_called_once()


def test_unset_urls(intervals: pd.DataFrame) -> None:
    """Unverified track names never become automatic downloads."""
    with pytest.raises(ValueError, match="Fill BIGWIG_URLS"):
        conservation.query_conservation(intervals, urls={"phyloP": None, "phastCons": None})


def test_no_remote_support(intervals: pd.DataFrame, monkeypatch: pytest.MonkeyPatch) -> None:
    """Fail clearly if the installed pyBigWig lacks libcurl support."""
    monkeypatch.setattr(conservation.pyBigWig, "remote", 0)
    with pytest.raises(RuntimeError, match="without remote"):
        conservation.query_conservation(intervals, urls=URLS)
