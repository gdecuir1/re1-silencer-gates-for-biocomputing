"""Read only motif windows from remote bigWigs through pyBigWig/libBigWig."""

import re
from collections.abc import Mapping
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import pyBigWig

from src import config
from src._utils import validate_intervals


def motif_windows(elements: pd.DataFrame) -> pd.DataFrame:
    """Convert oriented 200 bp BED intervals to 21 bp BED motif windows.

    Args:
        elements: element_id, chrom, start, end, strand. Strand is the orientation
            of the sequence in which MOTIF_SPAN is defined, not the FIMO hit strand.

    Returns:
        element_id, chrom, start, end with zero-based, half-open motif coordinates.
    """
    validate_intervals(elements)
    left, right = config.MOTIF_SPAN
    windows = elements[["element_id", "chrom", "start", "end"]].copy()
    reverse = elements["strand"].eq("-")
    windows["start"] = np.where(reverse, elements["end"] - right - 1, elements["start"] + left)
    windows["end"] = np.where(reverse, elements["end"] - left, elements["start"] + right + 1)
    return windows.astype({"start": "int64", "end": "int64"})


def _check_byte_ranges(url: str) -> None:
    """Verify byte-range support without consuming a full response body.

    Args:
        url: HTTP(S) bigWig URL.

    Raises:
        ValueError: URL is not HTTP(S), or server ignores the range request.
    """
    parsed = urlparse(url)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ValueError("bigWig URLs must use HTTP or HTTPS")
    request = Request(url, headers={"Range": "bytes=0-0", "Accept-Encoding": "identity"})
    with urlopen(request, timeout=30) as response:
        if response.status != 206 or not re.fullmatch(
            r"bytes 0-0/\d+", response.headers.get("Content-Range", "")
        ):
            raise ValueError("bigWig server does not honor byte ranges; refusing full download")
        response.read(1)


def query_conservation(
    elements: pd.DataFrame, *, urls: Mapping[str, str | None] | None = None
) -> pd.DataFrame:
    """Return exact per-base means and coverage for phyloP and phastCons.

    Args:
        elements: Oriented hg38 200 bp BED intervals with unique element_id.
        urls: Explicit phyloP and phastCons HTTP(S) URLs, otherwise config settings.

    Returns:
        One row per element with motif coordinates, track means, and covered-base
        counts. Uncovered bases stay missing; all-missing windows have a NaN mean.

    Raises:
        ValueError: Unset URLs, unsupported byte ranges, or invalid coordinates.
        RuntimeError: pyBigWig lacks remote support or cannot open a track.
    """
    windows = motif_windows(elements)
    tracks = dict(config.BIGWIG_URLS if urls is None else urls)
    if set(tracks) != {"phyloP", "phastCons"}:
        raise ValueError("Provide exactly phyloP and phastCons URLs")
    if any(not isinstance(url, str) or not url.strip() for url in tracks.values()):
        raise ValueError("Fill BIGWIG_URLS with verified hg38 phyloP and phastCons URLs")
    if not pyBigWig.remote:
        raise RuntimeError("pyBigWig was built without remote HTTP support (libcurl required)")
    result = windows.copy()
    for track, url in tracks.items():
        _check_byte_ranges(url)
        handle = pyBigWig.open(url)
        if handle is None:
            raise RuntimeError(f"Cannot open remote {track} bigWig")
        try:
            if not handle.isBigWig():
                raise ValueError(f"{track} URL is not a bigWig")
            chromosomes = handle.chroms()
            means, counts = [], []
            for row in windows.itertuples(index=False):
                if row.chrom not in chromosomes or row.end > chromosomes[row.chrom]:
                    raise ValueError(f"{track}: invalid chromosome/window for {row.element_id}")
                values = np.asarray(handle.values(row.chrom, row.start, row.end), dtype=float)
                # pyBigWig can return [] for a bigWig with no entries anywhere.
                if values.size not in (0, row.end - row.start):
                    raise ValueError(f"{track}: incomplete response for {row.element_id}")
                if np.isinf(values).any():
                    raise ValueError(f"{track}: infinite values for {row.element_id}")
                finite = values[np.isfinite(values)]
                means.append(float(finite.mean()) if finite.size else np.nan)
                counts.append(int(finite.size))
            result[f"{track}_mean"] = means
            result[f"{track}_covered_bases"] = counts
        finally:
            handle.close()
    return result
