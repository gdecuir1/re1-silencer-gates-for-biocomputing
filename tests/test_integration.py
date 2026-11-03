"""Installed-tool checks using temporary SYNTHETIC fixtures, never research data."""

import random
import shutil
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq

from src import folds, motifs, orthologs

pytestmark = pytest.mark.integration


def _tool(name: str) -> str:
    """Find an optional command-line executable or skip its integration check.

    Args:
        name: Executable name.

    Returns:
        Resolved executable path when installed.
    """
    executable = shutil.which(name)
    if executable is None:
        pytest.skip(f"{name} not installed/on PATH; unit tests remain data-free")
    return executable


def test_real_cdhit_and_folds(tmp_path: Path) -> None:
    """Cluster duplicate/reverse-complement synthetic sequences with real CD-HIT."""
    executable = _tool("cd-hit-est")
    rng = random.Random(11)
    sequence = "".join(rng.choices("ACGT", k=200))
    different = "".join(rng.choices("ACGT", k=200))
    source = pd.DataFrame(
        {
            "element_id": ["synthetic_a", "synthetic_duplicate", "synthetic_rc", "synthetic_other"],
            "chrom": ["chr1", "chr2", "chr2", "chr3"],
            "sequence": [sequence, sequence, str(Seq(sequence).reverse_complement()), different],
        }
    )
    result = folds.build_folds(
        source, output=tmp_path / "folds.csv", n_splits=2, executable=executable
    )
    related = result[result["element_id"] != "synthetic_other"]
    assert related["cluster_id"].nunique() == 1
    assert related["fold"].nunique() == 1
    assert result["fold"].nunique() == 2


def test_real_liftover(tmp_path: Path) -> None:
    """Lift a synthetic full interval through an identity chain, retaining failures."""
    executable = _tool("liftOver")
    chain = tmp_path / "synthetic.chain"
    chain.write_text("chain 1000 chrSynthetic 1000 + 0 1000 chrTarget 1000 + 0 1000 1\n1000\n\n")
    source = pd.DataFrame(
        {
            "element_id": ["synthetic_maps", "synthetic_fails"],
            "chrom": ["chrSynthetic", "chrAbsent"],
            "start": [100, 100],
            "end": [300, 300],
            "strand": ["+", "-"],
        }
    )
    result = orthologs.lift_elements(source, chain, min_match=0.95, executable=executable)
    assert result.mapped["chrom"].tolist() == ["chrTarget"]
    assert result.mapped["end"].tolist() == [300]
    assert result.unmapped["element_id"].tolist() == ["synthetic_fails"]
    assert result.unmapped.loc[0, "reason"]


def test_real_bedtools(tmp_path: Path, intervals: pd.DataFrame) -> None:
    """Extract both strands from a small synthetic genome using real bedtools."""
    executable = _tool("bedtools")
    sequence = "A" * 100 + "C" * 200 + "G" * 100 + "A" * 200 + "T" * 400
    genome = tmp_path / "synthetic.fa"
    genome.write_text(">chrSynthetic\n" + sequence + "\n")
    result = orthologs.extract_sequences(intervals, genome, executable=executable)
    assert result["sequence"].tolist() == ["C" * 200, "T" * 200]


def test_real_fimo(tmp_path: Path) -> None:
    """Check installed FIMO against an explicitly artificial matrix and sequences."""
    executable = _tool("fimo")
    matrix, fasta = tmp_path / "synthetic.meme", tmp_path / "synthetic.fa"
    matrix.write_text(
        "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
        "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF SYNTHETIC_TEST_ONLY\nletter-probability matrix: alength= 4 w= 21 nsites= 20 E= 0\n"
        + "0.97 0.01 0.01 0.01\n"
        * 21
    )
    fasta.write_text(">synthetic_hit\n" + "A" * 200 + "\n>synthetic_unknown\n" + "N" * 200 + "\n")
    result = motifs.score_sequences(
        fasta, matrix, motif_id="SYNTHETIC_TEST_ONLY", executable=executable
    )
    assert result.loc[0, "has_hit"]
    assert result.loc[0, "stop"] - result.loc[0, "start"] + 1 == 21
    assert not result.loc[1, "has_hit"]
