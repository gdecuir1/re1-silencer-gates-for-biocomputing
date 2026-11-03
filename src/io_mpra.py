"""Load explicitly mapped MPRA tables without inferring scientific transformations."""

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from src import config
from src._utils import require_columns, validate_ids, validate_intervals, validate_sequences


def _mapping(mapping: Mapping[str, str | None], required: set[str]) -> dict[str, str]:
    """Check that a source-column mapping has been filled in.

    Args:
        mapping: Canonical names mapped to source headers.
        required: Canonical keys that must be present.

    Returns:
        Validated mapping with string values.
    """
    missing = sorted(required - set(mapping))
    unset = sorted(
        key for key, value in mapping.items() if not isinstance(value, str) or not value.strip()
    )
    if missing or unset:
        raise ValueError(f"Fill column mapping in config: missing keys={missing}; unset={unset}")
    result = dict(mapping)
    if len(set(result.values())) != len(result):
        raise ValueError("Column mapping reuses a source column")
    return result


def validate_mpra(
    frame: pd.DataFrame, *, expected_elements: int = config.EXPECTED_ELEMENT_COUNT
) -> None:
    """Validate a tidy table, including the complete element/cell-type product.

    Args:
        frame: Canonical table with element_id, sequence, cell_type, repression.
        expected_elements: Exact number of unique elements; default is 8,436.

    Raises:
        ValueError: Required columns, counts, lengths, values, or pairs are invalid.
    """
    require_columns(frame, ["element_id", "sequence", "cell_type", "repression"])
    validate_ids(frame["element_id"], unique=False)
    validate_sequences(frame["sequence"])
    count = frame["element_id"].nunique()
    if expected_elements < 1 or count != expected_elements:
        raise ValueError(f"Expected {expected_elements} elements; found {count}")
    if frame.duplicated(["element_id", "cell_type"]).any():
        raise ValueError("Duplicate (element_id, cell_type) pairs")
    if not frame["cell_type"].isin(config.CELL_TYPES).all():
        raise ValueError(f"cell_type must be one of {config.CELL_TYPES}; missing values forbidden")
    counts = frame.groupby("element_id")["cell_type"].nunique()
    if not counts.eq(len(config.CELL_TYPES)).all():
        raise ValueError(f"Each element must have all cell types: {config.CELL_TYPES}")
    activity = pd.to_numeric(frame["repression"], errors="coerce")
    if not np.isfinite(activity).all():
        raise ValueError("repression must contain finite numeric values; missing values forbidden")
    metadata = [column for column in frame if column not in ("cell_type", "repression")]
    for column in metadata:
        if frame.groupby("element_id")[column].nunique(dropna=False).gt(1).any():
            raise ValueError(f"Inconsistent {column} across cell types for an element")
    coordinates = {"chrom", "start", "end", "strand"}
    if coordinates.intersection(frame.columns):
        validate_intervals(frame.drop_duplicates("element_id"))


def tidy_mpra(
    raw: pd.DataFrame,
    *,
    columns: Mapping[str, str | None] | None = None,
    repression_columns: Mapping[str, str | None] | None = None,
    layout: str = config.MPRA_LAYOUT,
    expected_elements: int = config.EXPECTED_ELEMENT_COUNT,
) -> pd.DataFrame:
    """Convert wide or long source measurements using an explicit header mapping.

    Args:
        raw: Source table; no replicate averaging or sign changes are performed.
        columns: Canonical-to-source headers, defaulting to config.MPRA_COLUMNS.
        repression_columns: For wide input, cell type to source repression column.
        layout: Either 'wide' (one row/element) or 'long' (one row/element/cell).
        expected_elements: Exact expected unique element count.

    Returns:
        Validated tidy DataFrame, with numeric repression values as supplied.
    """
    if layout not in {"wide", "long"}:
        raise ValueError("layout must be 'wide' or 'long'")
    required = {"element_id", "sequence"}
    if layout == "long":
        required |= {"cell_type", "repression"}
    mapping = _mapping(config.MPRA_COLUMNS if columns is None else columns, required)
    if layout == "wide":
        if {"cell_type", "repression"}.intersection(mapping):
            raise ValueError("Wide columns mapping must not contain cell_type or repression")
        activity = _mapping(
            config.MPRA_REPRESSION_COLUMNS if repression_columns is None else repression_columns,
            set(config.CELL_TYPES),
        )
        if set(activity) != set(config.CELL_TYPES):
            raise ValueError(f"repression_columns must map exactly {config.CELL_TYPES}")
        if set(mapping.values()).intersection(activity.values()):
            raise ValueError("Metadata and repression mappings must use different source columns")
        require_columns(raw, [*mapping.values(), *activity.values()])
        metadata = raw[list(mapping.values())].rename(columns={v: k for k, v in mapping.items()})
        validate_ids(metadata["element_id"])
        tidy = pd.concat(
            [
                metadata.assign(cell_type=cell, repression=raw[activity[cell]].to_numpy())
                for cell in config.CELL_TYPES
            ],
            ignore_index=True,
        )
    else:
        require_columns(raw, mapping.values())
        tidy = raw[list(mapping.values())].rename(columns={v: k for k, v in mapping.items()}).copy()
    validate_mpra(tidy, expected_elements=expected_elements)
    tidy["repression"] = pd.to_numeric(tidy["repression"])
    return tidy.reset_index(drop=True)


def load_mpra(
    path: Path = config.MPRA_PATH,
    *,
    columns: Mapping[str, str | None] | None = None,
    repression_columns: Mapping[str, str | None] | None = None,
    layout: str = config.MPRA_LAYOUT,
    sep: str = config.MPRA_SEPARATOR,
    expected_elements: int = config.EXPECTED_ELEMENT_COUNT,
) -> pd.DataFrame:
    """Read and validate an explicitly mapped delimited MPRA file.

    Args:
        path: Local input file; nothing is downloaded.
        columns: Canonical-to-source metadata/long-table column mapping.
        repression_columns: Wide-table cell-type-to-source-column mapping.
        layout: 'wide' or 'long'.
        sep: File delimiter, defaulting to tab.
        expected_elements: Exact expected element count.

    Returns:
        One row per (element_id, cell_type).
    """
    mapping = config.MPRA_COLUMNS if columns is None else columns
    # Preserve numeric-looking identifiers and literal strings such as 'NA'.
    id_column = mapping.get("element_id")
    dtype = {id_column: str} if id_column else None
    raw = pd.read_csv(path, sep=sep, dtype=dtype, keep_default_na=False)
    return tidy_mpra(
        raw,
        columns=columns,
        repression_columns=repression_columns,
        layout=layout,
        expected_elements=expected_elements,
    )
