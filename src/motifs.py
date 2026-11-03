"""FIMO scoring without p-value censoring or fabricated scores for missing hits."""

from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from src import config
from src._utils import read_fasta, require_columns, run_command, write_fasta

HIT_COLUMNS = ["motif_id", "sequence_name", "start", "stop", "strand", "score", "p-value"]


def _best_hits(table: pd.DataFrame) -> pd.DataFrame:
    """Choose maximum-score hits with deterministic ties.

    Args:
        table: FIMO hits.

    Returns:
        At most one hit per sequence, retaining one-based inclusive coordinates.
    """
    return table.sort_values(
        ["score", "p-value", "start", "stop", "strand"],
        ascending=[False, True, True, True, True],
        kind="stable",
    ).drop_duplicates("sequence_name")


def score_sequences(
    fasta: Path,
    motif_file: Path = config.MOTIF_PATH,
    *,
    executable: str = "fimo",
    motif_id: str = config.JASPAR_ID,
) -> pd.DataFrame:
    """Run FIMO and return the best hit for every input sequence.

    Args:
        fasta: Local DNA FASTA with unique identifiers; target lengths may vary.
        motif_file: Local MEME-format matrix containing the requested JASPAR ID.
        executable: FIMO executable available on PATH.
        motif_id: Exact MEME motif identifier to score.

    Returns:
        DataFrame keyed by element_id. FIMO start/stop remain one-based inclusive.
        Unscorable sequences have has_hit=False and missing scores/retains_motif.
        retains_motif uses score >= M_VALUE, not FIMO's p-value threshold.

    Raises:
        ValueError: Input matrix, sequences, or tool output are invalid.
        RuntimeError: FIMO cannot run.
    """
    records = read_fasta(Path(fasta))
    motif_file = Path(motif_file).resolve()
    motif_ids = {
        fields[1]
        for line in motif_file.read_text().splitlines()
        if len(fields := line.split()) >= 2 and fields[0] == "MOTIF"
    }
    if motif_id not in motif_ids:
        raise ValueError(f"MEME matrix missing motif {motif_id}")
    # Safe internal identifiers also prevent genomic-looking IDs being reinterpreted.
    aliases = {f"seq_{index}": name for index, name in enumerate(records)}
    best = pd.DataFrame(columns=HIT_COLUMNS)
    with TemporaryDirectory(prefix="re1-fimo-") as directory:
        temporary = Path(directory)
        input_path, output_path = temporary / "input.fa", temporary / "hits.tsv"
        write_fasta({alias: records[name] for alias, name in aliases.items()}, input_path)
        # Text mode avoids FIMO's max-stored-scores truncation. A threshold of 1
        # avoids censoring weak best scores. Stream to disk and reduce in chunks.
        run_command(
            [
                executable,
                "--text",
                "--thresh",
                "1",
                "--motif",
                motif_id,
                str(motif_file),
                str(input_path),
            ],
            stdout_path=output_path,
        )
        try:
            chunks = pd.read_csv(output_path, sep="\t", comment="#", chunksize=50_000)
            for chunk in chunks:
                require_columns(chunk, HIT_COLUMNS)
                if not chunk["sequence_name"].isin(aliases).all():
                    raise ValueError("FIMO returned unknown sequence identifiers")
                if not chunk["motif_id"].eq(motif_id).all():
                    raise ValueError("FIMO returned an unexpected motif identifier")
                for column in ("score", "p-value", "start", "stop"):
                    chunk[column] = pd.to_numeric(chunk[column], errors="raise")
                    if not np.isfinite(chunk[column]).all():
                        raise ValueError(f"FIMO returned nonfinite {column}")
                if not chunk.empty:
                    selected = _best_hits(chunk[HIT_COLUMNS])
                    best = _best_hits(selected if best.empty else pd.concat([best, selected]))
        except pd.errors.EmptyDataError as exc:
            raise ValueError("FIMO produced no TSV header") from exc
    best["sequence_name"] = best["sequence_name"].map(aliases)
    best = best.rename(columns={"sequence_name": "element_id", "p-value": "p_value"})
    result = pd.DataFrame({"element_id": list(records)}).merge(best, how="left", on="element_id")
    result["score"] = pd.to_numeric(result["score"])
    result["has_hit"] = result["score"].notna()
    result["retains_motif"] = result["score"].ge(config.M_VALUE).astype("boolean")
    result.loc[~result["has_hit"], "retains_motif"] = pd.NA
    return result
