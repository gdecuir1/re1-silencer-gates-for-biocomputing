"""Lift full elements and extract oriented sequences, retaining mapping failures."""

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from src._utils import (
    BED_COLUMNS,
    run_command,
    validate_intervals,
    validate_sequences,
    write_bed,
)


@dataclass(frozen=True)
class LiftOverResult:
    """Mapped BED6 records and unmapped original BED6 records with reasons."""

    mapped: pd.DataFrame
    unmapped: pd.DataFrame


def _read_bed(path: Path, *, unmapped: bool = False) -> pd.DataFrame:
    """Read liftOver BED6 output, preserving failure reasons and original lines.

    Args:
        path: liftOver output file.
        unmapped: Whether to retain comments as failure reasons.

    Returns:
        BED6 table, optionally with reason and raw_record fields.
    """
    rows: list[list[str]] = []
    reason = ""
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            reason = line.lstrip("#").strip()
            continue
        fields = line.split("\t")
        if len(fields) != 6:
            raise ValueError(f"Expected BED6 from liftOver; received {len(fields)} fields")
        rows.append(fields + [reason, line] if unmapped else fields)
        reason = ""
    columns = BED_COLUMNS + (["reason", "raw_record"] if unmapped else [])
    frame = pd.DataFrame(rows, columns=columns)
    for column in ("start", "end"):
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    validate_intervals(frame, full_elements=unmapped)
    return frame


def lift_elements(
    elements: pd.DataFrame,
    chain: Path,
    *,
    min_match: float,
    executable: str = "liftOver",
) -> LiftOverResult:
    """Lift entire 200 bp source elements using an explicit mapping threshold.

    Args:
        elements: Unique BED intervals with element_id, chrom, start, end, strand.
        chain: Existing source-hg38-to-target chain file.
        min_match: Minimum mapped fraction; must be chosen by the caller.
        executable: UCSC liftOver executable.

    Returns:
        Both mapped and unmapped tables. Target intervals may differ in length.
        Unmapped records retain their source coordinates and tool failure reason.
        Multimapping is not enabled; no arbitrary target is selected.

    Raises:
        ValueError: Source intervals are invalid or output loses/duplicates records.
        RuntimeError: The executable fails, as distinct from an individual failed lift.
    """
    validate_intervals(elements)
    if not 0 < min_match <= 1:
        raise ValueError("min_match must be in (0, 1]")
    chain = Path(chain).resolve()
    if not chain.is_file():
        raise FileNotFoundError(chain)
    with TemporaryDirectory(prefix="re1-liftover-") as directory:
        temporary = Path(directory)
        source, mapped, unmapped = (
            temporary / name for name in ("source.bed", "mapped.bed", "unmapped.bed")
        )
        write_bed(elements, source)
        run_command(
            [
                executable,
                f"-minMatch={min_match}",
                str(source),
                str(chain),
                str(mapped),
                str(unmapped),
            ]
        )
        result = LiftOverResult(_read_bed(mapped), _read_bed(unmapped, unmapped=True))
    returned = pd.concat([result.mapped["element_id"], result.unmapped["element_id"]])
    if returned.duplicated().any() or set(returned) != set(elements["element_id"]):
        raise ValueError("liftOver output must account for every element exactly once")
    return result


def extract_sequences(
    intervals: pd.DataFrame, genome: Path, *, executable: str = "bedtools"
) -> pd.DataFrame:
    """Run bedtools getfasta with stranded extraction and check all records return.

    Args:
        intervals: Mapped BED intervals; target lengths need not be 200 bp.
        genome: Existing target FASTA. bedtools may create its adjacent FASTA index.
        executable: bedtools executable available on PATH.

    Returns:
        BED6 columns plus sequence, oriented according to strand.

    Raises:
        ValueError: bedtools silently skips a chromosome or truncates an interval.
    """
    validate_intervals(intervals, full_elements=False)
    genome = Path(genome).resolve()
    if not genome.is_file():
        raise FileNotFoundError(genome)
    if intervals.empty:
        return pd.DataFrame(columns=[*BED_COLUMNS, "sequence"])
    with TemporaryDirectory(prefix="re1-getfasta-") as directory:
        temporary = Path(directory)
        bed, output = temporary / "input.bed", temporary / "sequences.tsv"
        write_bed(intervals, bed)
        run_command(
            [executable, "getfasta", "-fi", str(genome), "-bed", str(bed), "-s", "-bedOut"],
            stdout_path=output,
        )
        if output.stat().st_size == 0:
            raise ValueError("bedtools returned no sequences for the supplied intervals")
        result = pd.read_csv(
            output, sep="\t", header=None, keep_default_na=False, dtype={0: str, 3: str}
        )
    if result.shape[1] != 7:
        raise ValueError("Expected BED6 plus sequence from bedtools getfasta")
    result.columns = [*BED_COLUMNS, "sequence"]
    validate_intervals(result, full_elements=False)
    if set(result["element_id"]) != set(intervals["element_id"]):
        raise ValueError("bedtools did not return every element; check chromosome names and bounds")
    expected = intervals.set_index("element_id").loc[result["element_id"]]
    for column in ("chrom", "start", "end", "strand"):
        if not (result[column].to_numpy() == expected[column].to_numpy()).all():
            raise ValueError(f"bedtools changed input {column}")
    validate_sequences(result["sequence"], length=None)
    if not result["sequence"].str.len().eq(result["end"] - result["start"]).all():
        raise ValueError("bedtools returned truncated sequences")
    return result
