"""Unimplemented scientific models and a shared held-out prediction evaluator."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from src._utils import require_columns, validate_ids


def baseline(motif_score: pd.Series, target: pd.Series, folds: pd.Series) -> pd.DataFrame:
    """STUB: fit a motif-score-only model within the supplied training folds.

    Args:
        motif_score: Element-indexed motif scores.
        target: Aligned measured repression values.
        folds: Aligned, precomputed held-out fold IDs.

    Returns:
        Held-out prediction table compatible with evaluate_folds.

    Raises:
        NotImplementedError: Model specification has not been scientifically verified.
    """
    raise NotImplementedError("Specify and validate the motif-only baseline before fitting")


def motif_conservation(features: pd.DataFrame, target: pd.Series, folds: pd.Series) -> pd.DataFrame:
    """STUB: fit motif plus conservation models using training-only preprocessing.

    Args:
        features: Aligned motif and conservation predictors.
        target: Measured repression.
        folds: Precomputed held-out fold IDs.

    Returns:
        Held-out prediction table compatible with evaluate_folds.

    Raises:
        NotImplementedError: Covariates and model specification remain to be verified.
    """
    raise NotImplementedError("Specify conservation features, missingness handling, and model")


def motif_cofactor(features: pd.DataFrame, target: pd.Series, folds: pd.Series) -> pd.DataFrame:
    """STUB: fit motif plus cofactor occupancy models.

    Args:
        features: Aligned motif scores and experimentally supported occupancy features.
        target: Measured repression.
        folds: Precomputed held-out fold IDs.

    Returns:
        Held-out prediction table compatible with evaluate_folds.

    Raises:
        NotImplementedError: Occupancy definition and model remain unspecified.
    """
    raise NotImplementedError("Verify cofactor sources and feature definitions before fitting")


def gapped_kmer_svm(
    fasta: Path,
    labels: pd.DataFrame,
    folds: pd.DataFrame,
    *,
    gkmtrain: str = "gkmtrain",
    gkmpredict: str = "gkmpredict",
) -> pd.DataFrame:
    """STUB: train/predict per fold via lsgkm subprocesses with nested model selection.

    Args:
        fasta: Local 200 bp element FASTA.
        labels: Element IDs and explicitly defined training targets.
        folds: Frozen chromosome/cluster fold assignments.
        gkmtrain: lsgkm training executable on PATH.
        gkmpredict: lsgkm prediction executable on PATH.

    Returns:
        Held-out prediction table compatible with evaluate_folds.

    Raises:
        NotImplementedError: Kernel parameters and label definitions need verification.
    """
    raise NotImplementedError("Verify lsgkm labels, parameters, and nested selection protocol")


def evaluate_folds(predictions: pd.DataFrame) -> pd.DataFrame:
    """Calculate Spearman and AUROC per cell type and held-out fold.

    Args:
        predictions: Columns element_id, cell_type, fold, observed, prediction,
            and positive (explicit binary ground truth, 0/1). Larger predictions
            must mean stronger repression and greater probability of positive.
            Each element/cell pair must occur once, in its held-out fold.

    Returns:
        fold, cell_type, n, spearman, auroc, and status fields. Undefined metrics
        are NaN with a reason; no positive-label cutoff is inferred.

    Raises:
        ValueError: Inputs are duplicated, missing, nonfinite, or not binary labels.
    """
    require_columns(
        predictions, ["element_id", "cell_type", "fold", "observed", "prediction", "positive"]
    )
    if predictions.empty:
        raise ValueError("No held-out predictions to evaluate")
    validate_ids(predictions["element_id"], unique=False)
    validate_ids(predictions["cell_type"], unique=False)
    if predictions["fold"].isna().any():
        raise ValueError("fold must not be missing")
    if predictions.duplicated(["element_id", "cell_type"]).any():
        raise ValueError("Duplicate held-out element/cell predictions")
    if predictions.groupby("element_id")["fold"].nunique().gt(1).any():
        raise ValueError("An element appears in multiple held-out folds")
    values = predictions.copy()
    for column in ("observed", "prediction", "positive"):
        values[column] = pd.to_numeric(values[column], errors="coerce")
        if not np.isfinite(values[column]).all():
            raise ValueError(f"{column} must contain finite numeric values")
    if not values["positive"].isin([0, 1]).all():
        raise ValueError("positive must contain explicit binary labels (0/1)")
    rows = []
    for (fold, cell), group in values.groupby(["fold", "cell_type"], sort=True):
        correlation_ok = (
            len(group) >= 2
            and group["observed"].nunique() > 1
            and group["prediction"].nunique() > 1
        )
        auc_ok = group["positive"].nunique() == 2
        rows.append(
            {
                "fold": fold,
                "cell_type": cell,
                "n": len(group),
                "spearman": float(spearmanr(group["observed"], group["prediction"]).statistic)
                if correlation_ok
                else np.nan,
                "auroc": float(roc_auc_score(group["positive"], group["prediction"]))
                if auc_ok
                else np.nan,
                "spearman_status": "ok" if correlation_ok else "insufficient_or_constant_values",
                "auroc_status": "ok" if auc_ok else "single_class",
            }
        )
    return pd.DataFrame(rows)
