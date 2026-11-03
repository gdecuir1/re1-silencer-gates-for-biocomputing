"""STUBS: gate metrics require verified assay scales and experimental definitions."""

import pandas as pd


def dynamic_range(on_state: pd.Series, off_state: pd.Series) -> pd.Series:
    """STUB: quantify dynamic range on a verified measurement scale.

    Args:
        on_state: Aligned measurements in the on state.
        off_state: Aligned measurements in the off state.

    Returns:
        Element-indexed dynamic range once scientifically defined.

    Raises:
        NotImplementedError: Measurement scale and normalization are unspecified.
    """
    raise NotImplementedError("Verify on/off states, units, and normalization for dynamic range")


def off_state_leakiness(off_state: pd.Series, reference: pd.Series) -> pd.Series:
    """STUB: quantify residual expression relative to an agreed reference.

    Args:
        off_state: Off-state expression measurements.
        reference: Aligned reference measurements.

    Returns:
        Element-indexed leakiness once scientifically defined.

    Raises:
        NotImplementedError: Reference and normalization are unspecified.
    """
    raise NotImplementedError("Verify the off-state reference and normalization")


def separation(on_state: pd.DataFrame, off_state: pd.DataFrame) -> pd.Series:
    """STUB: quantify separation of replicated on/off measurements.

    Args:
        on_state: On-state replicate measurements.
        off_state: Off-state replicate measurements.

    Returns:
        Element-indexed separation once scientifically defined.

    Raises:
        NotImplementedError: Replicate structure and separation statistic are unspecified.
    """
    raise NotImplementedError("Verify replicate structure and the separation statistic")


def variant_robustness(reference: pd.Series, variants: pd.DataFrame) -> pd.Series:
    """STUB: quantify tolerance to sequence variants.

    Args:
        reference: Reference-element measurements.
        variants: Measured variants linked to their reference elements.

    Returns:
        Element-indexed robustness once scientifically defined.

    Raises:
        NotImplementedError: Variant set, scale, and aggregation are unspecified.
    """
    raise NotImplementedError("Verify variant mapping, assay scale, and robustness definition")
