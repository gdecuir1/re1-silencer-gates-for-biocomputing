"""Fixed two-model chromosome/sequence-grouped validation of gate outcomes.

Fit ordinary least squares to training-only empirical ranks. Test predictors
use the training empirical distribution; no test outcomes enter fitting or
preprocessing. Reuse cached, hash-verified complete cases and frozen targets.
"""

import json
import platform
import shutil
import subprocess

import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

from src import config, folds
from src.conservation_io import verified_elements
from src.cross_context import CONTEXTS
from src.dataset import digest, markdown

FOLDER = config.FOLDS_DIR
TABLES = config.TABLES_DIR
CASES = TABLES / "conservation_gate_preview_cases.csv"
FEATURES = ("fimo_score", "phyloP_motif_mean")
GROUP = ["outcome", "cell_type", "enhancer"]


def empirical_rank(training, values):
    """Training mid-distribution ranks: (count less + half count equal) / n.

    For training observations this is (average rank - 0.5) / n. Unseen test
    values are positioned using training observations alone, including tails.
    """
    training, values = np.asarray(training, float), np.asarray(values, float)
    if training.ndim != 1 or not len(training) or not np.isfinite(training).all():
        raise ValueError("Ranks require finite, nonempty training observations")
    if not np.isfinite(values).all():
        raise ValueError("Ranks do not impute missing values")
    ordered = np.sort(training)
    return (
        np.searchsorted(ordered, values, "left") + np.searchsorted(ordered, values, "right")
    ) / (2 * len(ordered))


def fit_predict(train, test, outcome, augmented):
    """OLS with intercept, one or two rank predictors, and training-ranked target."""
    columns = FEATURES if augmented else FEATURES[:1]
    design = np.column_stack(
        [np.ones(len(train))] + [empirical_rank(train[c], train[c]) for c in columns]
    )
    heldout = np.column_stack(
        [np.ones(len(test))] + [empirical_rank(train[c], test[c]) for c in columns]
    )
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("Rank-deficient training design; no substitute model is selected")
    coefficients = np.linalg.lstsq(
        design, empirical_rank(train[outcome], train[outcome]), rcond=None
    )[0]
    return heldout @ coefficients, dict(zip(("intercept", *columns), coefficients, strict=True))


def evaluate(cases, assignments):
    """Identical train/test complete cases for both models, shared folds across outcomes."""
    if assignments.element_id.duplicated().any():
        raise ValueError("Duplicate fold assignment")
    if set(assignments.fold) != set(range(config.N_FOLDS)):
        raise ValueError("Exactly five nonempty folds are required")
    for column in ("chrom", "cluster_id", "component_id"):
        if assignments.groupby(column).fold.nunique().gt(1).any():
            raise ValueError(f"Leakage: {column} crosses folds")
    if cases.duplicated([*GROUP, "element_id"]).any():
        raise ValueError("Duplicate element within endpoint")
    joined = cases.merge(assignments, on="element_id", how="left", validate="m:1")
    if joined.fold.isna().any():
        raise ValueError("Missing fold assignment")
    scores, predictions, coefficients = [], [], []
    for label, group in joined.groupby(GROUP, sort=True):
        outcome = label[0]
        tags = dict(zip(GROUP, label, strict=True))
        group = group.sort_values("element_id")
        if not np.isfinite(group[[*FEATURES, outcome]].to_numpy()).all():
            raise ValueError("Only complete cases may enter model fitting")
        for fold in range(config.N_FOLDS):
            train, test = group[group.fold != fold], group[group.fold == fold]
            if len(test) < 3 or test[outcome].nunique() < 2:
                raise ValueError("Held-out outcome is too small or constant for Spearman")
            row = dict(**tags, fold=fold, n_train=len(train), n_test=len(test))
            prediction = test[["element_id", "chrom", "cluster_id", "fold"]].copy()
            prediction["observed_outcome"] = test[outcome]
            for augmented, name in ((False, "baseline"), (True, "augmented")):
                pred, coef = fit_predict(train, test, outcome, augmented)
                if np.ptp(pred) == 0:
                    raise ValueError("Constant predictions; held-out Spearman undefined")
                row[f"{name}_spearman"] = float(spearmanr(test[outcome], pred).statistic)
                prediction[f"{name}_predicted_rank"] = pred
                coefficients.append(dict(**tags, fold=fold, model=name, **coef))
            row["delta"] = row["augmented_spearman"] - row["baseline_spearman"]
            row["improved"] = row["delta"] > 0
            scores.append(row)
            predictions.append(prediction.assign(**tags))
    return (
        pd.DataFrame(scores),
        pd.concat(predictions, ignore_index=True),
        pd.DataFrame(coefficients),
    )


def summarize(scores):
    """Show all-fold consistency and spread, without choosing a model or endpoint."""
    rows = []
    for label, group in scores.groupby(GROUP, sort=True):
        row = dict(zip(GROUP, label, strict=True))
        row.update(n_folds=len(group), folds_improved=int(group.improved.sum()))
        for metric in ("baseline_spearman", "augmented_spearman", "delta"):
            for stat, value in (
                ("mean", group[metric].mean()),
                ("median", group[metric].median()),
                ("min", group[metric].min()),
                ("max", group[metric].max()),
                ("sd", group[metric].std(ddof=1)),
            ):
                row[f"{metric}_{stat}"] = float(value)
        row["verdict"] = (
            "consistent improvement in every fold"
            if group.improved.all()
            else "no fold improved"
            if not group.improved.any()
            else "inconsistent improvement across folds"
        )
        rows.append(row)
    return pd.DataFrame(rows)


def read_cases():
    """Verify frozen preview provenance; reuse its exact complete-case eligibility."""
    manifest_path = config.PROCESSED_DIR / "conservation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if not manifest["sanity_gate_passed"] or manifest["m_value"] != config.M_VALUE:
        raise ValueError("Conservation gate/threshold disagrees with the preview")
    sources = [
        CASES,
        config.PROCESSED_DIR / "conservation_features.csv",
        config.PROCESSED_DIR / "gate_metrics.csv",
        TABLES / "gate_cross_context_candidates.csv",
    ]
    hashes = manifest["inputs"] | manifest["outputs"]
    for path in sources:
        if digest(path) != hashes[str(path.relative_to(config.ROOT))]:
            raise ValueError(f"Changed preview source: {path}")
    audit = pd.read_csv(CASES, low_memory=False)
    selected = audit[audit.cohort.isin(["above_m_value", "above_m_value_coverage_4_of_5"])].copy()
    if not selected.fimo_score.ge(config.M_VALUE).all():
        raise ValueError("Below-threshold rows in candidate cohort")
    cases = selected[selected.complete_case].copy()
    if not cases.phyloP_motif_coverage_status.eq("complete").all():
        raise ValueError("Complete-case flags disagree with motif coverage")
    expected = {("dynamic_range", cell, enhancer) for cell, enhancer in CONTEXTS}
    expected.add(("minimum_rank_percentile", "across_contexts", "five_passing_contexts"))
    if set(cases[GROUP].itertuples(index=False, name=None)) != expected:
        raise ValueError("Unexpected endpoint/context scope")
    dr = cases[cases.outcome.eq("dynamic_range")]
    gate = pd.read_csv(config.PROCESSED_DIR / "gate_metrics.csv", low_memory=False)
    actual = dr.merge(
        gate,
        on=["element_id", "cell_type", "enhancer"],
        how="left",
        validate="1:1",
        suffixes=("_preview", "_source"),
    )
    if not actual.reference_type.eq("scrambled").all() or not actual.population.eq("motif").all():
        raise ValueError("Preview contains a noncanonical or fallback reference")
    for name in ("dynamic_range", "fimo_score", "on_activity", "off_activity"):
        if not np.allclose(actual[name + "_preview"], actual[name + "_source"]):
            raise ValueError(f"Preview disagrees with gate source: {name}")
    if not np.allclose(dr.dynamic_range, dr.on_activity - dr.off_activity):
        raise ValueError("Scrambled ON minus intact OFF disagreement")
    robust = cases[cases.outcome.eq("minimum_rank_percentile")]
    if not robust.n_contexts_measured.ge(4).all():
        raise ValueError("Robustness requires at least four measured contexts")
    original = pd.read_csv(TABLES / "gate_cross_context_candidates.csv").set_index("element_id")
    for name in ("minimum_rank_percentile", "fimo_score", "n_contexts_measured"):
        if not np.array_equal(robust[name].to_numpy(), original.loc[robust.element_id, name]):
            raise ValueError(f"Frozen robustness target disagrees: {name}")
    coverage = (
        selected.groupby(GROUP, sort=True)
        .agg(n_source=("element_id", "size"), n_complete_cases=("complete_case", "sum"))
        .reset_index()
    )
    coverage["n_excluded"] = coverage.n_source - coverage.n_complete_cases
    return cases, coverage, [manifest_path, *sources]


def render(assignments, coverage, fold_counts, clusters, scores, summary):
    """Write the observed fold scores and an explicit endpoint-specific verdict."""
    verdicts = [
        f"- {r.outcome}, {r.cell_type}/{r.enhancer}: **{r.verdict}** "
        f"({r.folds_improved}/{r.n_folds} folds; mean Δρ={r.delta_mean:+.6f}, "
        f"range {r.delta_min:+.6f} to {r.delta_max:+.6f})."
        for r in summary.itertuples()
    ]
    consistent = int(summary.folds_improved.eq(config.N_FOLDS).sum())
    headline = (
        f"Conservation improves held-out prediction in every fold for {consistent} of "
        f"{len(summary)} outcomes. The other outcomes do not meet the every-fold "
        "consistency criterion."
        if consistent < len(summary)
        else "Conservation improves held-out prediction consistently across all requested outcomes."
    )
    spread = []
    for row in summary.to_dict("records"):
        for metric in ("baseline_spearman", "augmented_spearman", "delta"):
            spread.append(
                {
                    **{c: row[c] for c in GROUP},
                    "metric": metric,
                    **{s: row[f"{metric}_{s}"] for s in ("mean", "median", "min", "max", "sd")},
                }
            )
    lines = [
        "# Held-out chromosome validation of conservation",
        f"**Verdict: {headline}**",
        "\n".join(verdicts),
        "## Cohort and fixed folds",
        f"Canonical FIMO ≥{config.M_VALUE}, scrambled-paired only, in the five specified "
        f"passing contexts. There are **{len(assignments):,} distinct elements** in the union "
        "of complete cases across the six endpoints. All observations of an element use "
        "one shared fold. Each endpoint retains its own available complete cases; no "
        "intersection across all contexts is required and no missing outcome is imputed. "
        "Both models use exactly the same rows, including fully covered motif-window phyloP. "
        "No NegCtrl fallback is included.",
        markdown(coverage),
        f"CD-HIT-EST on full {config.ELEMENT_LENGTH}-bp sequences produced **{len(clusters):,} "
        f"clusters**, largest **{int(clusters.n_elements.max())} elements**; "
        f"**{int(clusters.n_elements.gt(1).sum())}** clusters contain multiple elements. "
        f"There are **{assignments.component_id.nunique()}** connected chromosome/cluster "
        "components. Components are greedily balanced by element count with deterministic "
        "ties. Chromosomes and clusters are never split across folds; fold IDs are zero-based.",
        f"Existing config, unchanged: global sequence identity {config.CDHIT_IDENTITY}, "
        f"coverage of both sequences {config.CDHIT_COVERAGE}, both strands, one thread, "
        "sorted element IDs, CD-HIT accurate greedy mode. Shared clusters link whole "
        "chromosomes before assignment. CD-HIT is representative-based; membership is not "
        "a guarantee that every distant relative or transposon-family member is grouped.",
        markdown(fold_counts),
        "[Fold assignments](../folds/heldout_assignments.csv) · "
        "[Cluster sizes](../folds/heldout_clusters.csv) · "
        "[Provenance](../folds/heldout_manifest.json)",
        "## Two fixed models and rank preprocessing",
        "Baseline: OLS with intercept and ranked FIMO alone. Augmented: OLS with intercept, "
        "ranked FIMO and ranked motif-window mean phyloP. The target is ranked within the "
        "training fold. No tuning, feature selection, extra models or pooled-enhancer fit.",
        "Ranks use the training empirical mid-distribution: "
        "(number of training values less than x + half the number equal to x) / training n. "
        "For training observations this equals (average-tie rank − 0.5) / n. Test predictors "
        "are mapped with that training distribution alone, with tails at 0 and 1. "
        "Spearman compares held-out predictions with observed held-out outcomes. Test outcome "
        "ranks enter evaluation only. This extends the in-sample rank-linear framing without "
        "letting test observations determine fitting or preprocessing.",
        "Dynamic range is scrambled ON minus intact OFF in log2 units, fit separately per "
        "context. Robustness is the existing minimum context rank percentile, one value per "
        "element with at least four measured contexts. Its original cohort-wide rank "
        "denominators are frozen as part of the target definition; the target is not "
        "recomputed or selected using validation scores.",
        "## All held-out folds",
        "Δρ = augmented − baseline held-out Spearman; positive means conservation helps. "
        "The expected balanced size is approximate: chromosome constraints and endpoint "
        "missingness determine the actual n_test below.",
        markdown(scores),
        "## Per-fold spread and consistency",
        "Means are unweighted across folds; SD is sample SD across folds. Shared training "
        "sets make folds dependent, so these spreads are descriptive, not confidence intervals.",
        markdown(summary[[*GROUP, "folds_improved", "n_folds", "verdict"]]),
        markdown(pd.DataFrame(spread)),
        "[All fold scores](heldout_scores.csv) · [Summary](heldout_summary.csv) · "
        "[Out-of-fold predictions](heldout_predictions.csv) · "
        "[Training coefficients](heldout_coefficients.csv)",
        "## Interpretation and reproducibility",
        "This is internal chromosome-held-out validation of an already observed association. "
        "The above-threshold cohort and passing contexts were chosen using the existing "
        "dataset before this validation; these results are not an independent external "
        "replication. No thresholds or folds were changed in response to model scores. "
        "A positive in-sample partial correlation alone does not establish predictive benefit; "
        "the fold deltas above are the relevant comparison.",
        "Run `make validate-models PYTHON=.venv/bin/python` with `cd-hit-est` on PATH. "
        "Only existing local sequences, conservation and gate outcomes are read; conservation "
        "is not re-extracted. The manifest records input and output SHA-256 hashes, clustering "
        "parameters and software versions. Fold assignments are stored under `results/folds/` "
        "for version control; nothing under `data/` is staged or committed. "
        "[Two-run receipt](../folds/heldout_reproducibility.json).",
    ]
    return "\n\n".join(lines) + "\n"


def main():
    cases, coverage, inputs = read_cases()
    elements = verified_elements()
    elements = elements[elements.element_id.isin(cases.element_id)].copy()
    if set(elements.element_id) != set(cases.element_id):
        raise ValueError("Missing verified sequence for validation cohort")
    executable = shutil.which("cd-hit-est")
    if executable is None:
        raise RuntimeError("Install the declared cd-hit dependency or put it on PATH")
    assignment_path = FOLDER / "heldout_assignments.csv"
    assignments = folds.build_folds(elements, output=assignment_path, executable=executable)
    scores, predictions, coefficients = evaluate(cases, assignments)
    summary = summarize(scores)
    clusters = (
        assignments.groupby("cluster_id", sort=True)
        .agg(
            n_elements=("element_id", "size"),
            fold=("fold", "first"),
            chromosomes=("chrom", lambda x: ",".join(sorted(set(x)))),
        )
        .reset_index()
    )
    fold_counts = (
        assignments.groupby("fold", sort=True)
        .agg(
            n_elements=("element_id", "size"),
            n_clusters=("cluster_id", "nunique"),
            n_components=("component_id", "nunique"),
            chromosomes=("chrom", lambda x: ",".join(sorted(set(x)))),
        )
        .reset_index()
    )
    outputs = {
        FOLDER / "heldout_clusters.csv": clusters,
        FOLDER / "heldout_counts.csv": fold_counts,
        TABLES / "heldout_scores.csv": scores,
        TABLES / "heldout_summary.csv": summary,
        TABLES / "heldout_predictions.csv": predictions,
        TABLES / "heldout_coefficients.csv": coefficients,
        TABLES / "heldout_coverage.csv": coverage,
    }
    for path, frame in outputs.items():
        frame.to_csv(path, index=False, float_format="%.17g")
    report = TABLES / "heldout_validation.md"
    report.write_text(render(assignments, coverage, fold_counts, clusters, scores, summary))
    help_result = subprocess.run([executable, "-h"], capture_output=True, text=True)
    version = next(
        x.strip()
        for x in (help_result.stdout + help_result.stderr).splitlines()
        if "CD-HIT version" in x
    )
    inputs.extend(
        [
            config.HG38_ELEMENTS,
            config.PROCESSED_DIR / "liftover_hg19_to_hg38/element_audit.csv",
            *[
                config.ROOT / "src" / name
                for name in (
                    "heldout_validation.py",
                    "folds.py",
                    "conservation_io.py",
                    "config.py",
                )
            ],
        ]
    )
    manifest = dict(
        m_value=config.M_VALUE,
        n_folds=config.N_FOLDS,
        n_elements=len(assignments),
        n_clusters=len(clusters),
        largest_cluster=int(clusters.n_elements.max()),
        cdhit_identity=config.CDHIT_IDENTITY,
        cdhit_coverage=config.CDHIT_COVERAGE,
        cdhit_version=version,
        cdhit_binary_sha256=digest(executable),
        models=["FIMO rank OLS", "FIMO + phyloP rank OLS"],
        ranking="training-only empirical mid-distribution; average ties; intercept",
        python=platform.python_version(),
        numpy=np.__version__,
        pandas=pd.__version__,
        scipy=scipy.__version__,
        remote_reads_performed=False,
        inputs={str(p.relative_to(config.ROOT)): digest(p) for p in inputs},
        outputs={
            str(p.relative_to(config.ROOT)): digest(p) for p in [assignment_path, *outputs, report]
        },
    )
    (FOLDER / "heldout_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(
        summary[
            [*GROUP, "folds_improved", "delta_mean", "delta_min", "delta_max", "verdict"]
        ].to_string(index=False)
    )
    print("Report: results/tables/heldout_validation.md")


if __name__ == "__main__":
    main()
