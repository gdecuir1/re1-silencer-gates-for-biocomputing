"""Compatibility exports for the unimplemented models in src.model_scaffold.

The fitted FIMO/phyloP rank-linear comparison is in src.heldout_validation.
This module preserves the original public imports without changing behavior.
"""

from src.model_scaffold import (
    baseline,
    evaluate_folds,
    gapped_kmer_svm,
    motif_cofactor,
    motif_conservation,
)

__all__ = [
    "baseline",
    "evaluate_folds",
    "gapped_kmer_svm",
    "motif_cofactor",
    "motif_conservation",
]
