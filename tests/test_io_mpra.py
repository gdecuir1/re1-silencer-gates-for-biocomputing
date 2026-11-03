"""Validate source schemas and malformed inputs using synthetic measurements."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import config
from src.io_mpra import load_mpra, tidy_mpra, validate_mpra


def test_valid_tidy(tidy_table: pd.DataFrame) -> None:
    """Accept the complete synthetic element/cell product."""
    validate_mpra(tidy_table, expected_elements=2)


@pytest.mark.parametrize("column", ["element_id", "sequence", "cell_type", "repression"])
def test_missing_column(tidy_table: pd.DataFrame, column: str) -> None:
    """Name the missing canonical column."""
    with pytest.raises(ValueError, match=column):
        validate_mpra(tidy_table.drop(columns=column), expected_elements=2)


@pytest.mark.parametrize("length", [0, 150, 199, 201])
def test_bad_length(tidy_table: pd.DataFrame, length: int) -> None:
    """Reject non-200 bp sequences, including the incorrect 150 bp assumption."""
    tidy_table.loc[0, "sequence"] = "A" * length
    with pytest.raises(ValueError, match="200 bp|nonempty"):
        validate_mpra(tidy_table, expected_elements=2)


@pytest.mark.parametrize("value", [None, "not-a-number", np.nan, np.inf, -np.inf])
def test_invalid_activity(tidy_table: pd.DataFrame, value: object) -> None:
    """Reject missing, nonnumeric, and infinite measurements."""
    tidy_table["repression"] = tidy_table["repression"].astype(object)
    tidy_table.loc[0, "repression"] = value
    with pytest.raises(ValueError, match="repression"):
        validate_mpra(tidy_table, expected_elements=2)


def test_wrong_count(tidy_table: pd.DataFrame) -> None:
    """Use the real expected element count by default."""
    with pytest.raises(ValueError, match="8436.*2"):
        validate_mpra(tidy_table)


def test_duplicate_pairs(tidy_table: pd.DataFrame) -> None:
    """Reject repeated observations without choosing an aggregation method."""
    with pytest.raises(ValueError, match="Duplicate"):
        validate_mpra(pd.concat([tidy_table, tidy_table.iloc[:1]]), expected_elements=2)


def test_missing_cell_type(tidy_table: pd.DataFrame) -> None:
    """Require all four measurements for every element."""
    with pytest.raises(ValueError, match="all cell types"):
        validate_mpra(tidy_table.iloc[1:], expected_elements=2)


def test_unexpected_cell_type(tidy_table: pd.DataFrame) -> None:
    """Reject unrecognized cell-type labels."""
    tidy_table.loc[0, "cell_type"] = "synthetic_unknown"
    with pytest.raises(ValueError, match="cell_type"):
        validate_mpra(tidy_table, expected_elements=2)


def test_inconsistent_sequence(tidy_table: pd.DataFrame) -> None:
    """An element cannot change sequence between cell types."""
    tidy_table.loc[0, "sequence"] = "T" * 200
    with pytest.raises(ValueError, match="Inconsistent sequence"):
        validate_mpra(tidy_table, expected_elements=2)


def test_wide_load_preserves_identifiers_and_values(tmp_path: Path) -> None:
    """Load an explicitly mapped TSV without inferring signs or numeric IDs."""
    raw = pd.DataFrame(
        {
            "source_id": ["001", "002"],
            "source_dna": ["A" * 200, "T" * 200],
            **{f"measure_{cell}": [-2.5, 0.0] for cell in config.CELL_TYPES},
        }
    )
    path = tmp_path / "synthetic.tsv"
    raw.to_csv(path, sep="\t", index=False)
    result = load_mpra(
        path,
        columns={"element_id": "source_id", "sequence": "source_dna"},
        repression_columns={cell: f"measure_{cell}" for cell in config.CELL_TYPES},
        expected_elements=2,
    )
    assert len(result) == 8
    assert set(result["element_id"]) == {"001", "002"}
    assert set(result["repression"]) == {-2.5, 0.0}


def test_long_layout(tidy_table: pd.DataFrame) -> None:
    """The long-table path also validates the complete product."""
    result = tidy_mpra(
        tidy_table, columns={key: key for key in tidy_table}, layout="long", expected_elements=2
    )
    pd.testing.assert_frame_equal(result, tidy_table)


def test_unset_mapping() -> None:
    """Unfilled configuration fails clearly rather than guessing source headers."""
    with pytest.raises(ValueError, match="Fill column mapping.*element_id"):
        tidy_mpra(pd.DataFrame(), columns={"element_id": None, "sequence": None})


def test_missing_source_header() -> None:
    """Name the actual source column missing from a mapped file."""
    with pytest.raises(ValueError, match="absent_header"):
        tidy_mpra(
            pd.DataFrame(),
            columns={"element_id": "absent_header", "sequence": "dna"},
            repression_columns={cell: cell for cell in config.CELL_TYPES},
        )


def test_reused_source_mapping() -> None:
    """A source column cannot represent both ID and sequence."""
    with pytest.raises(ValueError, match="reuses"):
        tidy_mpra(pd.DataFrame(), columns={"element_id": "same", "sequence": "same"})
