"""Build continuous MPRA observations without aggregating enhancer contexts.

Each row is one element × cell_type × enhancer. Repression is the median
activity of NegCtrl elements in that SAME cell/enhancer minus element activity.
Positive repression therefore means lower expression than the matched baseline,
even when raw log2FoldChange is positive. Missing activity remains missing.
Published padj is retained as metadata, never used to make binary silencer calls.
"""

import gzip
import hashlib
import json
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from src import config
from src.assembly_check import GEO_ARCHIVE

RAW = config.RAW_DIR / "scoring_validation"
VALIDATION = config.PROCESSED_DIR / "scoring_validation"
DATASET = config.PROCESSED_DIR / "analysis_dataset.csv"
CONTEXTS = config.PROCESSED_DIR / "analysis_dataset_contexts.csv"
MANIFEST = config.PROCESSED_DIR / "analysis_dataset_manifest.json"
REPORT = config.TABLES_DIR / "analysis_dataset.md"
POPULATIONS = ("motif", "noMotif", "NegCtrl")
CELLS = {"GM12878": "GM12878", "HepG2": "HepG2", "K562": "K562", "SKNSH": "SK-N-SH"}
KEY = ["element_id", "cell_type", "enhancer"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_source():
    """Read every selected source row; require exact GEO/author identity."""
    frames, inventory = [], []
    sources = json.loads((RAW / "sources.json").read_text())
    with tarfile.open(GEO_ARCHIVE) as archive:
        for original, cell in CELLS.items():
            names = [n for n in archive.getnames() if f"WGREST_{original}.txt.gz" in n]
            if len(names) != 1:
                raise ValueError(f"Expected exactly one source table for {cell}")
            with gzip.GzipFile(fileobj=archive.extractfile(names[0])) as handle:
                geo = pd.read_csv(handle, sep="\t", low_memory=False)
            path = RAW / f"RESTscreen_derived_{original}_results.run1.txt.gz"
            if digest(path) != sources[path.name]["sha256"]:
                raise ValueError(f"Author source hash mismatch: {path}")
            author = pd.read_csv(
                path, sep="\t", usecols=["ID", "project", "log2FoldChange", "padj"]
            )
            merged = geo.merge(
                author,
                left_on="id",
                right_on="ID",
                how="outer",
                validate="1:1",
                suffixes=("", "_author"),
                indicator=True,
            )
            same = merged.log2FoldChange.eq(merged.log2FoldChange_author) | (
                merged.log2FoldChange.isna() & merged.log2FoldChange_author.isna()
            )
            if not merged._merge.eq("both").all() or not same.all():
                raise ValueError(f"GEO and author activity/ID mismatch in {cell}")
            if not merged.project.eq(merged.project_author).all():
                raise ValueError(f"Population disagreement in {cell}")
            if not merged.id.eq(merged.enhancer + "^" + merged.silencer).all():
                raise ValueError("Construct ID disagrees with explicit source metadata")
            selected = merged.loc[
                merged.project.isin(POPULATIONS),
                ["id", "silencer", "project", "enhancer", "log2FoldChange", "padj"],
            ]
            selected = selected.rename(
                columns={
                    "id": "source_construct_id",
                    "silencer": "element_id",
                    "project": "population",
                    "log2FoldChange": "raw_activity",
                }
            )
            selected["cell_type"] = cell
            frames.append(selected)
            inventory.append(
                {
                    "cell_type": cell,
                    "source_rows_all_populations": len(geo),
                    "selected_source_rows": len(selected),
                    "excluded_other_populations": len(geo) - len(selected),
                }
            )
    source = pd.concat(frames, ignore_index=True)
    if source.duplicated(KEY).any():
        raise ValueError("Duplicate element/cell/enhancer in source")
    if source.groupby("element_id").population.nunique().ne(1).any():
        raise ValueError("Inconsistent element population")
    return source, pd.DataFrame(inventory)


def matched_repression(source):
    """Subtract observed activity from the context's observed NegCtrl median."""
    if source.duplicated(KEY).any() or source[KEY + ["population"]].isna().any().any():
        raise ValueError("Missing or duplicate source keys")
    if not np.isfinite(source.raw_activity.dropna()).all():
        raise ValueError("Nonfinite activity")
    if not source.padj.dropna().between(0, 1).all():
        raise ValueError("Invalid supplied adjusted p-value")
    controls = source[source.population.eq("NegCtrl")]
    baseline = (
        controls.groupby(["cell_type", "enhancer"])
        .raw_activity.agg(
            negctrl_n_source="size",
            negctrl_n_measured="count",
            negctrl_baseline="median",
            negctrl_sd="std",
            negctrl_q25=lambda x: x.quantile(0.25),
            negctrl_q75=lambda x: x.quantile(0.75),
            negctrl_p05=lambda x: x.quantile(0.05),
            negctrl_p95=lambda x: x.quantile(0.95),
        )
        .reset_index()
    )
    baseline["negctrl_iqr"] = baseline.negctrl_q75 - baseline.negctrl_q25
    result = source.merge(baseline, on=["cell_type", "enhancer"], how="left", validate="m:1")
    result["repression"] = result.negctrl_baseline - result.raw_activity
    for field in ["raw_activity", "padj", "negctrl_baseline", "repression"]:
        result[f"missing_{field}"] = result[field].isna()
    return result


def attach_metadata(frame, coordinates, audit, scores):
    """Retain unmapped/contrast/control rows; never invent coordinates or scores."""
    for name, table in [("coordinates", coordinates), ("lift audit", audit), ("scores", scores)]:
        if table.element_id.duplicated().any():
            raise ValueError(f"Duplicate element IDs in {name}")
    if not coordinates.working_build.eq("hg38").all():
        raise ValueError("Coordinates must be verified hg38")
    if not (coordinates.end - coordinates.start).eq(config.ELEMENT_LENGTH).all():
        raise ValueError("Verified intervals must be 200 bp")
    canonical_ids = set(frame.loc[frame.population.eq("motif"), "element_id"])
    if set(scores.element_id) != canonical_ids or not np.isfinite(scores.score).all():
        raise ValueError("Canonical scoring coverage changed or scores are missing/nonfinite")
    coords = coordinates[["element_id", "chrom", "start", "end", "strand"]].rename(
        columns={c: f"hg38_{c}" for c in ("chrom", "start", "end", "strand")}
    )
    result = frame.merge(coords, on="element_id", how="left", validate="m:1")
    result = result.merge(
        audit[["element_id", "bucket"]], on="element_id", how="left", validate="m:1"
    )
    result["coordinate_status"] = result.bucket.fillna("no_verified_interval_for_control")
    result["missing_hg38_coordinates"] = result.hg38_start.isna()
    if (result.population.ne("NegCtrl") & result.bucket.isna()).any():
        raise ValueError("Native element is missing from the exhaustive lift audit")
    if not result.bucket.eq("lifted+verified").eq(~result.missing_hg38_coordinates).all():
        raise ValueError("Coordinate availability disagrees with sequence verification")
    result["eligible_for_conservation"] = ~result.missing_hg38_coordinates
    for c in ("hg38_start", "hg38_end"):
        result[c] = result[c].astype("Int64")
    result = result.drop(columns="bucket").merge(
        scores[["element_id", "score"]].rename(columns={"score": "fimo_score"}),
        on="element_id",
        how="left",
        validate="m:1",
    )
    result["score_status"] = np.where(
        result.population.eq("motif"), "canonical_21bp_core", "not_applicable_no_canonical_site"
    )
    result["missing_fimo_score"] = result.fimo_score.isna()
    result["above_m_value"] = result.fimo_score.ge(config.M_VALUE).astype("boolean")
    result.loc[result.missing_fimo_score, "above_m_value"] = pd.NA
    return result


def assert_source_rows(source, result):
    """Assert equal row counts AND exact key sets, not a Cartesian reconstruction."""
    if len(result) != len(source):
        raise AssertionError(f"Source rows {len(source)} != dataset rows {len(result)}")
    if result.duplicated(KEY).any():
        raise AssertionError("Duplicate output observations")
    left = source[KEY].sort_values(KEY).reset_index(drop=True)
    right = result[KEY].sort_values(KEY).reset_index(drop=True)
    if not left.equals(right):
        raise AssertionError("Dataset keys differ from source despite row count")


def characterize(frame):
    """Rank within-context canonical-site correlations; no enhancer aggregation."""
    rows = []
    baseline_cols = [c for c in frame if c.startswith("negctrl_")]
    for (cell, enhancer), group in frame.groupby(["cell_type", "enhancer"], sort=True):
        motif = group[group.population.eq("motif")]
        contrast = group[group.population.eq("noMotif")]
        paired = motif.dropna(subset=["fimo_score", "repression"])
        usable = (
            len(paired) >= 3 and paired.fimo_score.nunique() > 1 and paired.repression.nunique() > 1
        )
        rho = float(spearmanr(paired.fimo_score, paired.repression).statistic) if usable else np.nan
        pearson = (
            float(pearsonr(paired.fimo_score, paired.repression).statistic) if usable else np.nan
        )
        row = {
            "cell_type": cell,
            "enhancer": enhancer,
            **group.iloc[0][baseline_cols].to_dict(),
            "motif_n_source": len(motif),
            "motif_n_measured": int(motif.raw_activity.notna().sum()),
            "noMotif_n_measured": int(contrast.raw_activity.notna().sum()),
            "motif_median_activity": motif.raw_activity.median(),
            "noMotif_median_activity": contrast.raw_activity.median(),
            "motif_median_repression": motif.repression.median(),
            "motif_vs_noMotif_gap": contrast.raw_activity.median() - motif.raw_activity.median(),
            "n_correlation_pairs": len(paired),
            "spearman_rho": rho,
            "pearson_r": pearson,
            "correlation_status": "measured" if usable else "insufficient_or_constant",
        }
        rows.append(row)
    result = pd.DataFrame(rows).sort_values(
        ["spearman_rho", "cell_type", "enhancer"], ascending=[False, True, True], na_position="last"
    )
    result.insert(0, "rank", range(1, len(result) + 1))
    result["recommendation"] = "do_not_prioritize"
    positive = result.spearman_rho.gt(0) & result.pearson_r.gt(0)
    result.loc[positive, "recommendation"] = "secondary_positive_context"
    # Predeclared descriptive selection: strongest concordant context per cell.
    primary = result[positive].groupby("cell_type", sort=False).head(1).index
    result.loc[primary, "recommendation"] = "primary_within_cell"
    return result.reset_index(drop=True)


def markdown(frame):
    def fmt(value):
        if pd.isna(value):
            return "NA"
        if isinstance(value, (float, np.floating)):
            return f"{value:.6f}"
        return str(value)

    return "\n".join(
        [
            "| " + " | ".join(frame.columns) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *(
                "| " + " | ".join(map(fmt, r)) + " |"
                for r in frame.itertuples(index=False, name=None)
            ),
        ]
    )


def render(frame, contexts, inventory, gate):
    population = (
        frame.groupby("population", sort=True)
        .agg(
            elements=("element_id", "nunique"),
            rows=("element_id", "size"),
            measured_activity=("raw_activity", "count"),
            missing_activity=("missing_raw_activity", "sum"),
            missing_padj=("missing_padj", "sum"),
            missing_hg38=("missing_hg38_coordinates", "sum"),
            missing_score=("missing_fimo_score", "sum"),
        )
        .reset_index()
    )
    example = frame[
        frame.population.eq("motif") & frame.raw_activity.gt(0) & frame.repression.gt(0)
    ].iloc[0]
    picks = contexts[contexts.recommendation.eq("primary_within_cell")]
    selection = "; ".join(
        f"{r.cell_type}/{r.enhancer} (ρ={r.spearman_rho:.3f})" for r in picks.itertuples()
    )
    ranking = [
        "rank",
        "cell_type",
        "enhancer",
        "negctrl_n_measured",
        "negctrl_baseline",
        "negctrl_sd",
        "negctrl_iqr",
        "motif_n_measured",
        "motif_vs_noMotif_gap",
        "n_correlation_pairs",
        "spearman_rho",
        "pearson_r",
        "recommendation",
    ]
    spread = [
        "rank",
        "cell_type",
        "enhancer",
        "negctrl_n_source",
        "negctrl_n_measured",
        "negctrl_q25",
        "negctrl_q75",
        "negctrl_p05",
        "negctrl_p95",
    ]
    coord_counts = (
        frame[["element_id", "coordinate_status"]]
        .drop_duplicates()
        .groupby("coordinate_status")
        .size()
        .rename("elements")
        .reset_index()
    )
    return (
        "\n\n".join(
            [
                "# Continuous MPRA analysis dataset",
                f"**Built {len(frame):,} rows at element × cell_type × enhancer resolution. "
                "No enhancer aggregation and no binary silencer calls.**",
                "## Authorized override and retained limitation",
                "gate_metrics_allowed = true. " + gate["override"]["reason"],
                "The original validation status remains recorded as FAIL; scoring and sign checks "
                "passed. "
                "The unresolved Figure 2B binary classification is not consumed here. Supplied padj"
                " is "
                "retained unchanged, not filtered or interpreted as a matched-enhancer repression "
                "FDR. "
                "This dataset authorizes continuous analysis; it does not demonstrate a working NOT"
                " gate. "
                "See [the scoring audit](scoring_validation.md).",
                "## Source accounting and assertions",
                "Only the explicitly requested motif, noMotif and NegCtrl populations are included."
                " "
                "Other WGREST populations (variants, scrambled sequences and other controls) are "
                "counted "
                "as outside this dataset's scope, not silently removed. Source missing measurements"
                " stay "
                "as rows. No missing context is synthesized; every existing source key appears "
                "exactly once.",
                markdown(inventory.assign(dataset_rows=inventory.selected_source_rows)),
                markdown(population),
                "Expected from selected source: "
                f"**{inventory.selected_source_rows.sum():,}** rows; "
                f"observed **{len(frame):,}**. Duplicate keys: **{frame.duplicated(KEY).sum()}**. "
                "`assert_source_rows` checks both the row count and exact sorted key sets; metadata"
                " joins "
                "are many-to-one and fail on duplicated metadata. All GEO activities, including "
                "missingness, "
                "are required to match the authors' tables before their padj column is joined.",
                "## Repression sign and example",
                "**repression = median(NegCtrl activity in the SAME cell and enhancer) − raw "
                "activity. "
                "Positive means greater repression relative to that context; negative means "
                "activity "
                "above its baseline.** Units are log2 activity difference. Missing activity or "
                "baseline "
                "gives missing repression, never zero. Baseline medians use measured NegCtrl "
                "elements "
                "only; SD uses ddof=1 and quantiles use linear interpolation. NegCtrl rows use the "
                "same "
                "full context baseline (not leave-one-out).",
                f"Observed example: `{example.element_id}`, "
                f"{example.cell_type}/{example.enhancer}: "
                f"baseline **{example.negctrl_baseline:.6f}** − raw activity "
                f"**{example.raw_activity:.6f}** = repression **{example.repression:.6f}**. "
                "The raw activity is positive, yet the element represses relative to its matched "
                "baseline.",
                "## Ranked contexts",
                "Ranking uses Spearman ρ between canonical-core FIMO score and continuous "
                "repression, "
                "separately in each cell/enhancer; Pearson r is a sensitivity check. Only paired "
                "finite "
                "measurements enter correlations; n is shown. Missingness, padj, and the 20.86 "
                "split do "
                "not filter these pairs. Expected contexts: 4 cells × 5 enhancers = 20; observed: "
                f"**{len(contexts)}**. The motif-versus-noMotif gap is "
                "median(noMotif raw activity) minus "
                "median(motif raw activity), equivalently motif minus noMotif median repression in "
                "that "
                "context. Positive gap means motif elements repress more. Correlations are "
                "unitless; "
                "baseline, spread and gaps are in log2 activity units.",
                markdown(contexts[ranking]),
                "### Additional baseline spread (same ranking)",
                markdown(contexts[spread]),
                "## Context recommendation",
                "Carry forward **" + selection + "** as the primary context in each cell. "
                "The explicit descriptive rule is the strongest positive Spearman association per "
                "cell, "
                "requiring positive Pearson association too. Retain other concordantly positive "
                "contexts "
                "as separate secondary analyses; do not average them into the primary measurements."
                " "
                "Nonpositive or undefined associations are not prioritized for REST-dependent "
                "metrics. "
                "The global ranking gives the order for follow-up across cells. This is "
                "prioritization "
                "on the observed assay, not held-out validation or a causal proof of "
                "REST activity.",
                "En02 has a near-zero NegCtrl baseline and was the source study's background "
                "control enhancer. Carry those high-correlation En02 contexts as sequence-response "
                "arms, not as evidence of a strong ON state. Also carry En19 separately in all "
                "four cells for active-enhancer follow-up, alongside GM12878/En09. "
                "En19 combines positive score associations with higher baseline expression; "
                "a gate design must evaluate both response and baseline rather than "
                "correlation alone.",
                "A population median gap and a within-motif score correlation answer different "
                "questions: "
                "a small or reversed overall gap need not imply no score response among canonical "
                "motifs. "
                "Baseline subtraction changes the interpretation and zero point; within one context"
                " it "
                "cannot manufacture correlation because it subtracts the same constant from every "
                "activity. "
                "Explicit on/off conditions and a measurement-scale definition are still required "
                "before "
                "interpreting dynamic range or leakiness as gate performance.",
                "## Coordinates, motif scoring and missingness",
                markdown(coord_counts),
                "Only sequence-verified hg38 intervals from elements_hg38.csv populate "
                "hg38_chrom/start/end/strand "
                "(zero-based, half-open). Native lift failures and sequence mismatches remain with "
                "missing "
                "hg38 coordinates and their audit status. NegCtrl source records lack complete "
                "intervals "
                "and are outside the native lift catalog; their coordinates are not reconstructed "
                "from "
                "IDs. All these rows still contribute measured activities and control baselines as "
                "appropriate.",
                "fimo_score is the validated MA0138.2 score of the canonical 21 bp site at "
                "zero-based "
                "inclusive offsets (90,110); above_m_value is score ≥20.86. All canonical scores "
                "must "
                "be present. For noMotif and NegCtrl there is no defined canonical site: score and "
                "boolean "
                "stay missing with score_status=not_applicable_no_canonical_site. No whole-element "
                "score "
                "is silently substituted. Missingness flags cover activity, padj, baseline, "
                "repression, "
                "coordinates and score. Missing scores do not mean below threshold.",
                "eligible_for_conservation marks verified hg38 intervals for a later phyloP join. "
                "No phyloP values are computed or imputed by this module; production track "
                "selection "
                "and extraction remain separate work.",
                "## Reproduction",
                "Run `make dataset PYTHON=.venv/bin/python`. This is offline and consumes the "
                "cached "
                "GEO archive, author tables, validated canonical scores, verified hg38 catalog and "
                "lift "
                "audit. It does not run upstream pipelines or change the threshold. Input and "
                "output "
                "SHA-256 hashes are in data/processed/analysis_dataset_manifest.json. "
                "The measured two-run comparison is [recorded "
                "separately](analysis_dataset_reproducibility.json). "
                "No data are staged or committed.",
            ]
        )
        + "\n"
    )


def main():
    gate_path = VALIDATION / "gate.json"
    gate = json.loads(gate_path.read_text())
    if not (
        gate.get("gate_metrics_allowed") is True
        and gate.get("override", {}).get("approved_for_continuous_analysis") is True
        and gate.get("scoring") == "PASS"
        and gate.get("sign") == "PASS"
    ):
        raise ValueError(
            "Continuous dataset requires the explicit user override and passed scoring/sign"
        )
    if gate.get("threshold") != config.M_VALUE or config.M_VALUE != 20.86:
        raise ValueError("Fixed threshold differs from validated threshold")
    score_path = VALIDATION / "canonical_scores.csv"
    score_manifest = json.loads((VALIDATION / "manifest.json").read_text())
    expected = score_manifest["outputs"][str(score_path.relative_to(config.ROOT))]
    if digest(score_path) != expected:
        raise ValueError("Validated canonical scores changed")
    if digest(GEO_ARCHIVE) != score_manifest["inputs"][str(GEO_ARCHIVE.relative_to(config.ROOT))]:
        raise ValueError("Validated GEO activity source changed")
    source, inventory = read_source()
    audit_path = config.PROCESSED_DIR / "liftover_hg19_to_hg38/element_audit.csv"
    result = attach_metadata(
        matched_repression(source),
        pd.read_csv(config.HG38_ELEMENTS),
        pd.read_csv(audit_path),
        pd.read_csv(score_path),
    )
    result = result.sort_values(KEY, kind="stable").reset_index(drop=True)
    result = result[[*KEY, "population", *[c for c in result if c not in [*KEY, "population"]]]]
    assert_source_rows(source, result)
    if not result.loc[result.missing_raw_activity, "repression"].isna().all():
        raise AssertionError("Missing activity was imputed")
    contexts = characterize(result)
    if len(contexts) != 20:
        raise AssertionError(f"Expected 20 contexts, observed {len(contexts)}")
    inputs = [
        GEO_ARCHIVE,
        gate_path,
        score_path,
        VALIDATION / "manifest.json",
        config.HG38_ELEMENTS,
        audit_path,
        RAW / "sources.json",
        Path(__file__),
        Path(config.__file__),
    ]
    inputs += [RAW / f"RESTscreen_derived_{cell}_results.run1.txt.gz" for cell in CELLS]
    hashes = {str(p.relative_to(config.ROOT)): digest(p) for p in inputs}
    report = render(result, contexts, inventory, gate)
    with TemporaryDirectory(prefix="dataset-", dir=config.PROCESSED_DIR) as directory:
        staging = Path(directory)
        result.to_csv(staging / DATASET.name, index=False, float_format="%.17g")
        contexts.to_csv(staging / CONTEXTS.name, index=False, float_format="%.17g")
        (staging / REPORT.name).write_text(report)
        outputs = {
            str(p.relative_to(config.ROOT)): digest(staging / p.name)
            for p in (DATASET, CONTEXTS, REPORT)
        }
        manifest = {
            "grain": KEY,
            "source_rows": len(source),
            "dataset_rows": len(result),
            "inputs": hashes,
            "outputs": outputs,
            "threshold": config.M_VALUE,
            "gate_metrics_allowed": True,
            "binary_silencer_calls_used": False,
        }
        (staging / MANIFEST.name).write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        for p in (DATASET, CONTEXTS, REPORT, MANIFEST):
            (staging / p.name).replace(p)
    print(f"Source rows = dataset rows = {len(result):,}; unique source keys preserved.")
    print(
        contexts[["rank", "cell_type", "enhancer", "spearman_rho", "recommendation"]].to_string(
            index=False
        )
    )
    print(f"Report: {REPORT.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
