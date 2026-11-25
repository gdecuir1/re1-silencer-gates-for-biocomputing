"""Continuous gate metrics on matched log2 activities; variant analysis is open.

Positive dynamic range is repression. Leakiness is the OFF level itself.
Separation uses RMS standard error, not biological replicate SD or a test statistic.
"""

import numpy as np
import pandas as pd


def _aligned(*series: pd.Series) -> None:
    """Reject silent reindexing, duplicate keys and infinite measurements."""
    for s in series:
        if not s.index.is_unique or not s.index.equals(series[0].index):
            raise ValueError("Metric inputs must have identical unique indexes")
        if not np.isfinite(s.dropna()).all():
            raise ValueError("Metric inputs must be finite or missing")


def dynamic_range(on_state: pd.Series, off_state: pd.Series) -> pd.Series:
    """ON minus OFF, in log2 activity units; either missing state gives NA."""
    _aligned(on_state, off_state)
    return on_state - off_state


def off_state_leakiness(off_state: pd.Series, reference: pd.Series) -> pd.Series:
    """OFF log2 activity, with no floor/normalization; missing ON or OFF gives NA.

    reference is the matched ON state and supplies the missingness mask only.
    """
    _aligned(off_state, reference)
    return off_state.where(reference.notna())


def separation(on_state: pd.DataFrame, off_state: pd.DataFrame) -> pd.Series:
    """Uncertainty-standardized separation = (ON-OFF)/sqrt((SE_ON²+SE_OFF²)/2).

    Both frames supply activity and standard_error on identical unique indexes.
    SEs describe fitted log2FoldChange uncertainty (or bootstrap uncertainty of a
    control median); they are NOT replicate standard deviations. This dimensionless
    precision ranking is neither Cohen's d nor a paired-test statistic. No covariance
    is assumed or inferred. Missing SE or zero pooled SE yields NA, never infinity.
    """
    on, off = on_state.activity, off_state.activity
    a, b = on_state.standard_error, off_state.standard_error
    _aligned(on, off, a, b)
    if a.dropna().lt(0).any() or b.dropna().lt(0).any():
        raise ValueError("Standard errors cannot be negative")
    pooled = np.sqrt((a.pow(2) + b.pow(2)) / 2)
    return dynamic_range(on, off) / pooled.where(pooled.gt(0))


def variant_robustness(reference: pd.Series, variants: pd.DataFrame) -> pd.Series:
    """STUB: variant mapping, measurement scale and robustness remain undefined."""
    raise NotImplementedError("Verify variant mapping, assay scale, and robustness definition")
