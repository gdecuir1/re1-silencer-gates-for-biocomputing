"""Explicitly synthetic in-memory fixtures; no research data are used."""

import pandas as pd
import pytest

from src import config


@pytest.fixture
def tidy_table() -> pd.DataFrame:
    """Return a small synthetic table for validator tests only."""
    return pd.DataFrame(
        [
            {
                "element_id": f"synthetic_{index}",
                "sequence": base * config.ELEMENT_LENGTH,
                "cell_type": cell,
                "repression": float(index),
            }
            for index, base in enumerate(("A", "C"))
            for cell in config.CELL_TYPES
        ]
    )


@pytest.fixture
def intervals() -> pd.DataFrame:
    """Return two oriented synthetic 200 bp BED intervals."""
    return pd.DataFrame(
        {
            "element_id": ["synthetic_plus", "synthetic_minus"],
            "chrom": ["chrSynthetic", "chrSynthetic"],
            "start": [100, 400],
            "end": [300, 600],
            "strand": ["+", "-"],
        }
    )
