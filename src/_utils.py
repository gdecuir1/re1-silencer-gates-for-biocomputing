"""Shared validation and external-command mechanics."""

import subprocess
from collections.abc import Iterable, Mapping
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from src.config import ELEMENT_LENGTH

BED_COLUMNS = ["chrom", "start", "end", "element_id", "score", "strand"]


def require_columns(frame: pd.DataFrame, columns: Iterable[str]) -> None:
    """Raise a ValueError naming missing or duplicated column labels.

    Args:
        frame: Input table.
        columns: Required column labels.
    """
    if frame.columns.duplicated().any():
        raise ValueError("Duplicate column labels are not allowed")
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")


def validate_ids(values: pd.Series, *, unique: bool = True) -> None:
    """Validate nonempty, whitespace-free string identifiers.

    Args:
        values: Identifiers to check.
        unique: Whether duplicate identifiers are forbidden.
    """
    valid = values.map(
        lambda value: (
            isinstance(value, str) and bool(value) and not any(char.isspace() for char in value)
        )
    )
    if not valid.all():
        raise ValueError("element_id must contain nonempty strings without whitespace")
    if unique and values.duplicated().any():
        raise ValueError("Duplicate element_id values")


def validate_sequences(values: pd.Series, length: int | None = ELEMENT_LENGTH) -> None:
    """Check nonempty IUPAC DNA sequences and optional exact length.

    Args:
        values: Sequence strings.
        length: Required length, or None for variable-length target intervals.
    """
    valid = values.map(lambda value: isinstance(value, str) and bool(value))
    if not valid.all():
        raise ValueError("sequence must contain nonempty strings")
    if length is not None and not values.str.len().eq(length).all():
        raise ValueError(f"Every element sequence must be {length} bp")
    if not values.str.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+").all():
        raise ValueError("sequence contains invalid IUPAC DNA characters")


def validate_intervals(frame: pd.DataFrame, *, full_elements: bool = True) -> None:
    """Validate unique, oriented BED intervals without silently changing coordinates.

    Args:
        frame: Table with element_id, chrom, start, end, and strand.
        full_elements: Require 200 bp source intervals if True.
    """
    require_columns(frame, ["element_id", "chrom", "start", "end", "strand"])
    validate_ids(frame["element_id"])
    validate_ids(frame["chrom"], unique=False)
    for column in ("start", "end"):
        series = frame[column]
        if series.map(lambda value: isinstance(value, (bool, np.bool_))).any():
            raise ValueError(f"{column} must contain integer coordinates")
        numeric = pd.to_numeric(series, errors="coerce")
        if not np.isfinite(numeric).all() or not numeric.mod(1).eq(0).all():
            raise ValueError(f"{column} must contain finite integer coordinates")
        if not pd.api.types.is_numeric_dtype(series):
            raise ValueError(f"{column} must be numeric, not strings")
    if (frame["start"] < 0).any() or (frame["end"] <= frame["start"]).any():
        raise ValueError("Intervals require 0 <= start < end")
    if full_elements and not (frame["end"] - frame["start"]).eq(ELEMENT_LENGTH).all():
        raise ValueError(f"Source intervals must span the full {ELEMENT_LENGTH} bp element")
    if not frame["strand"].isin(["+", "-"]).all():
        raise ValueError("strand must be explicitly '+' or '-'")


def write_bed(frame: pd.DataFrame, path: Path) -> None:
    """Write validated intervals as BED6 with integer coordinates.

    Args:
        frame: Validated interval table.
        path: Output BED path.
    """
    bed = frame.assign(score=0).loc[:, BED_COLUMNS].copy()
    bed[["start", "end"]] = bed[["start", "end"]].astype("int64")
    bed.to_csv(path, sep="\t", header=False, index=False)


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA without discarding duplicate identifiers.

    Args:
        path: Input FASTA path.

    Returns:
        Ordered mapping of FASTA identifiers to DNA strings.
    """
    records: dict[str, str] = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in records:
            raise ValueError(f"Duplicate FASTA identifier: {record.id}")
        records[record.id] = str(record.seq)
    if not records:
        raise ValueError(f"No FASTA records in {path.name}")
    validate_ids(pd.Series(list(records)))
    validate_sequences(pd.Series(list(records.values())), length=None)
    return records


def write_fasta(records: Mapping[str, str], path: Path) -> None:
    """Write already validated FASTA records.

    Args:
        records: Identifier-to-sequence mapping.
        path: Output path.
    """
    with path.open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")


def run_command(args: list[str], *, stdout_path: Path | None = None) -> str:
    """Run an argv list without a shell and expose useful failure diagnostics.

    Args:
        args: Executable and its arguments.
        stdout_path: Stream potentially large standard output to this file.

    Returns:
        Standard output text.

    Raises:
        RuntimeError: Executable is missing or exits unsuccessfully.
    """
    try:
        with stdout_path.open("w") if stdout_path else nullcontext(subprocess.PIPE) as output:
            result = subprocess.run(
                args, check=True, stdout=output, stderr=subprocess.PIPE, text=True
            )
    except FileNotFoundError as exc:
        raise RuntimeError(f"Required executable is not on PATH: {args[0]}") from exc
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"{args[0]} failed ({exc.returncode}): {exc.stderr.strip()}") from exc
    return result.stdout or ""
