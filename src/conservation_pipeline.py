"""Continuous conservation features and a mandatory pre-join sanity gate.

No activity aggregation across enhancers. Partial Spearman means Pearson
correlation of rank residuals after regressing each ranked variable on ranked
FIMO with an intercept, using identical complete cases for both correlations.
"""

import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from src import config, conservation_gate_preview
from src.conservation import motif_windows
from src.conservation_io import load_cached, verified_elements
from src.conservation_strata import CHECKS, corrected_checks, distribution_figure
from src.cross_context import CONTEXTS
from src.dataset import KEY, digest, markdown

REPORT = config.TABLES_DIR / "conservation.md"
FEATURES = config.PROCESSED_DIR / "conservation_features.csv"
JOINED = config.PROCESSED_DIR / "analysis_dataset_conservation.csv"
MANIFEST = config.PROCESSED_DIR / "conservation_manifest.json"
ORIGINAL = config.RESULTS_DIR / "archive/conservation_original_gate"


def window_summary(values):
    """Observed-base summaries with explicit denominators; unknown is never neutral."""
    if values.ndim != 2 or np.isinf(values).any():
        raise ValueError("Expected finite-or-missing per-base matrix")
    n = np.isfinite(values).sum(axis=1)
    total = values.shape[1]
    mean = np.divide(np.nansum(values, axis=1), n, out=np.full(len(values), np.nan), where=n > 0)
    maximum = np.max(np.where(np.isfinite(values), values, -np.inf), axis=1)
    maximum[n == 0] = np.nan
    positive = ((values > 0) & np.isfinite(values)).sum(axis=1)
    fraction = np.divide(positive, n, out=np.full(len(values), np.nan), where=n > 0)
    return dict(
        mean=mean,
        max=maximum,
        covered_bases=n,
        requested_bases=total,
        coverage_fraction=n / total,
        positive_bases=positive,
        positive_fraction_observed=fraction,
        coverage_status=np.where(n == 0, "no_data", np.where(n == total, "complete", "partial")),
    )


def summarize(elements, arrays, absent):
    """Derive motif, full-element and complementary flank statistics from exact bases."""
    windows = motif_windows(elements)
    result = elements[["element_id", "project", "chrom", "start", "end", "strand"]].copy()
    result = result.rename(columns={"project": "population"})
    result["motif_start"] = windows.start
    result["motif_end"] = windows.end
    left, right = config.MOTIF_SPAN
    for track in ("phyloP", "phastCons"):
        values = np.asarray(arrays[track], dtype=float).copy()
        if values.shape != (len(elements), config.ELEMENT_LENGTH):
            raise ValueError("Per-base array does not match the element catalog")
        if track == "phastCons" and ((values < 0) | (values > 1)).any():
            raise ValueError("Invalid phastCons probability")
        reverse = elements.strand.eq("-").to_numpy()
        values[reverse] = values[reverse, ::-1]
        regions = dict(
            element=values,
            motif=values[:, left : right + 1],
            flank=np.concatenate((values[:, :left], values[:, right + 1 :]), axis=1),
        )
        result[f"{track}_chromosome_absent"] = absent[track]
        for window, bases in regions.items():
            for stat, value in window_summary(bases).items():
                result[f"{track}_{window}_{stat}"] = value
    return result.sort_values("element_id").reset_index(drop=True)


def coverage_table(features):
    """Count complete, partial and absent windows per track and population."""
    rows = []
    for track in ("phyloP", "phastCons"):
        for window in ("motif", "element", "flank"):
            for population in ("all", "motif", "noMotif"):
                group = (
                    features
                    if population == "all"
                    else features[features.population.eq(population)]
                )
                status = group[f"{track}_{window}_coverage_status"]
                rows.append(
                    dict(
                        track=track,
                        window=window,
                        population=population,
                        n_elements=len(group),
                        any_values=int(status.ne("no_data").sum()),
                        complete=int(status.eq("complete").sum()),
                        partial=int(status.eq("partial").sum()),
                        no_data=int(status.eq("no_data").sum()),
                        chromosome_absent=int(group[f"{track}_chromosome_absent"].sum()),
                    )
                )
    return pd.DataFrame(rows)


def sanity_checks(features):
    """Predeclared direction checks on fully covered windows, evaluated for both tracks."""
    paired, distributions, verdicts = [], [], []
    for track in ("phyloP", "phastCons"):
        for population in ("motif", "noMotif"):
            group = features[features.population.eq(population)]
            full = group[f"{track}_element_coverage_status"].eq("complete")
            delta = group.loc[full, f"{track}_motif_mean"] - group.loc[full, f"{track}_flank_mean"]
            median = float(delta.median()) if len(delta) else np.nan
            positive = float(delta.gt(0).mean()) if len(delta) else np.nan
            paired.append(
                dict(
                    track=track,
                    population=population,
                    n=len(delta),
                    excluded_partial_or_missing=int((~full).sum()),
                    mean_difference=delta.mean(),
                    median_difference=median,
                    p05_difference=delta.quantile(0.05),
                    p95_difference=delta.quantile(0.95),
                    fraction_motif_greater=positive,
                    n_equal=int(delta.eq(0).sum()),
                )
            )
            if population == "motif":
                verdicts.append(
                    dict(
                        check="canonical motif versus own flanks",
                        track=track,
                        expected="median difference > 0 and fraction positive > 0.5",
                        observed_median_difference=median,
                        observed_fraction_positive=positive,
                        passed=bool(len(delta) and median > 0 and positive > 0.5),
                    )
                )
            for scope in ("complete", "any_observed"):
                status = group[f"{track}_motif_coverage_status"]
                selected = group.loc[
                    status.eq("complete") if scope == "complete" else status.ne("no_data"),
                    f"{track}_motif_mean",
                ]
                distributions.append(
                    dict(
                        track=track,
                        population=population,
                        scope=scope,
                        n=len(selected),
                        mean=selected.mean(),
                        sd=selected.std(),
                        p05=selected.quantile(0.05),
                        p25=selected.quantile(0.25),
                        median=selected.median(),
                        p75=selected.quantile(0.75),
                        p95=selected.quantile(0.95),
                    )
                )
        d = pd.DataFrame(distributions)
        current = d[d.track.eq(track) & d.scope.eq("complete")].set_index("population")
        difference = current.loc["motif", "median"] - current.loc["noMotif", "median"]
        verdicts.append(
            dict(
                check="canonical motif versus noMotif window",
                track=track,
                expected="canonical median minus noMotif median > 0",
                observed_median_difference=difference,
                observed_fraction_positive=np.nan,
                passed=bool(current.n.gt(0).all() and difference > 0),
            )
        )
    return pd.DataFrame(paired), pd.DataFrame(distributions), pd.DataFrame(verdicts)


def join_features(dataset, features, verdicts):
    """Left join only after every conservation sanity check passes; retain all keys."""
    if (
        len(verdicts) != 2
        or set(verdicts.check) != set(CHECKS)
        or not verdicts.track.eq("phyloP").all()
        or not verdicts.passed.all()
    ):
        raise ValueError("Conservation sanity gate failed; join prohibited")
    if dataset.duplicated(KEY).any() or features.element_id.duplicated().any():
        raise ValueError("Duplicate analysis key or conservation element")
    columns = ["element_id"] + [
        c for c in features if c.startswith(("phyloP_", "phastCons_element_"))
    ]
    result = dataset.merge(
        features[columns], on="element_id", how="left", validate="m:1", sort=False
    )
    if len(result) != len(dataset) or not result[KEY].equals(dataset[KEY]):
        raise AssertionError("Conservation join changed source row count or order")
    pd.testing.assert_frame_equal(result[dataset.columns], dataset)
    result["missing_conservation_coordinates"] = ~result.element_id.isin(features.element_id)
    for track in ("phyloP", "phastCons"):
        for window in ("motif", "element", "flank") if track == "phyloP" else ("element",):
            result[f"missing_{track}_{window}_mean"] = result[f"{track}_{window}_mean"].isna()
    return result


def partial_spearman(x, y, z):
    """Rank residual correlation, controlling a single covariate with an intercept."""
    values = np.column_stack([x, y, z])
    if len(values) < 4 or not np.isfinite(values).all():
        return np.nan
    ranks = np.column_stack([rankdata(values[:, i], method="average") for i in range(3)])
    if any(np.ptp(ranks[:, i]) == 0 for i in range(3)):
        return np.nan
    design = np.column_stack([np.ones(len(values)), ranks[:, 2]])
    residuals = ranks[:, :2] - design @ np.linalg.lstsq(design, ranks[:, :2], rcond=None)[0]
    if np.any(np.linalg.norm(residuals, axis=0) < 1e-10):
        return np.nan
    return float(np.corrcoef(residuals.T)[0, 1])


def preview(joined):
    """Descriptive canonical-element correlations per passing context, same-case raw/partial."""
    rows = []
    for cell, enhancer in CONTEXTS:
        g = joined[
            joined.population.eq("motif") & joined.cell_type.eq(cell) & joined.enhancer.eq(enhancer)
        ]
        full = g.phyloP_motif_coverage_status.eq("complete")
        valid = g.loc[full].dropna(subset=["phyloP_motif_mean", "repression", "fimo_score"])
        raw = (
            float(spearmanr(valid.phyloP_motif_mean, valid.repression).statistic)
            if len(valid) > 2
            and valid.phyloP_motif_mean.nunique() > 1
            and valid.repression.nunique() > 1
            else np.nan
        )
        rows.append(
            dict(
                cell_type=cell,
                enhancer=enhancer,
                n_canonical=len(g),
                n_complete_cases=len(valid),
                n_excluded=len(g) - len(valid),
                missing_coordinates=int(g.missing_conservation_coordinates.sum()),
                incomplete_or_missing_motif=int((~full).sum()),
                missing_repression=int(g.repression.isna().sum()),
                missing_fimo=int(g.fimo_score.isna().sum()),
                spearman_phyloP_repression=raw,
                partial_spearman_given_fimo=partial_spearman(
                    valid.phyloP_motif_mean, valid.repression, valid.fimo_score
                ),
            )
        )
    return pd.DataFrame(rows)


def render(
    features,
    coverage,
    paired,
    distributions,
    verdicts,
    identities,
    dataset,
    joined,
    correlations,
    original_verdicts,
    observed,
    strata,
):
    """Report the authorized check-design correction without revising extracted values."""
    passed = bool(verdicts.passed.all())
    primary = strata[strata.scope.eq("complete")]
    canonical = strata[strata.scope.eq("complete") & strata.stratum.ne("noMotif")]
    weak_n = int(
        canonical.loc[canonical.stratum.eq("canonical_below_m_value"), "source_elements"].iloc[0]
    )
    canonical_n = int(canonical.source_elements.sum())
    lines = [
        "# Continuous mammalian conservation",
        f"**Corrected sanity gate: {'PASS' if passed else 'FAIL'}.** "
        f"Reused the existing windows for {len(features):,} verified hg38 elements. "
        "No conservation was re-extracted; feature values, coverage and the m-value are unchanged.",
        "## Why the gate was revised",
        "This is an explicitly user-authorized correction to check design after observi"
        "ng the original "
        "failure, not a correction to the conservation data or an independently preregi"
        "stered test. "
        "The original failing gate is retained below and in the "
        "[archived report](../archive/conservation_original_gate/conservation.md).",
        "**Primary motif feature: phyloP.** phyloP evaluates individual alignment columns. "
        "**Secondary feature: full-element phastCons.** phastCons uses a hidden Markov model and "
        "neighboring columns to infer conserved elements, so a motif and its flanks can share the "
        "same conserved block. Its 21-base-versus-flank contrast is therefore removed f"
        "rom the gate; "
        "phastCons remains a secondary 200-base feature. This is a scale distinction, not a fixed "
        "claim that every phastCons block is hundreds of bases long. "
        "[UCSC method description](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g="
        "cons470way).",
        "**Population check:** compare above-threshold canonical sites (FIMO ≥20.86) with noMotif "
        "using mean phyloP in fully covered motif windows. The directional median crite"
        "rion remains "
        ">0; only the canonical comparison population changes. Below-threshold sites re"
        "main visible "
        "as a third stratum. The existing within-element phyloP criterion remains media"
        "n difference "
        ">0 and fraction positive >0.5 across all fully covered canonical elements.",
        f"The measured source inventory has **{weak_n:,}/{canonical_n:,} "
        f"({100 * weak_n / canonical_n:.2f}%)** canonical elements below 20.86. "
        "This is the full-library "
        "fraction, not the roughly 74% chr8-pilot fraction. Transposon enrichment is a proposed "
        "biological explanation, not something measured by this analysis; low score alone does "
        "not establish evolutionary age or absence of function.",
        "## Corrected checks",
        markdown(verdicts),
        "### Score-stratified motif-window phyloP distributions",
        "Each source element is counted once, irrespective of its measured cell/enhancer contexts. "
        "Primary comparisons require all 21 bases. Exclusions for absent verified coordinates, "
        "partial windows and no-data windows are explicit. All source canonical FIMO scores must "
        "be finite; missing scores trigger an error rather than assignment below the threshold.",
        markdown(primary[["stratum", "n", "mean", "median", "p05", "p25", "p75", "p95"]]),
        markdown(
            primary[
                ["stratum", "source_elements", "missing_coordinates", "partial", "no_data", "n"]
            ]
        ),
        "[Full per-element stratum values, including exclusions](conservation_stratum_v"
        "alues.csv) · "
        "[Complete and any-observed distribution summaries](conservation_strata.csv). "
        "The per-element CSV contains every observed value, so these are not just selec"
        "ted quantiles.",
        "![PhyloP distributions by FIMO stratum](../figures/conservation_strata.png)",
        "[PDF](../figures/conservation_strata.pdf) · [SVG](../figures/conservation_strata.svg). "
        "The histogram uses shared 0.25-score-unit bins and includes pooled canonical elements; "
        "the empirical CDF includes every complete observation. Neither fits a mixture or smooths "
        "the data to manufacture modes. At this displayed resolution, the pooled canonical "
        "histogram has a dominant peak near zero and a long right tail, without two clearly "
        "separated peaks. Strong sites have a broader conserved tail. These are descriptive "
        "observations, not a formal multimodality test; bimodality is not established here.",
        "### Any-observed coverage sensitivity (not used to change the gate)",
        markdown(strata[strata.scope.eq("any_observed")]),
        "## Original failing gate: preserved",
        markdown(original_verdicts),
        "The original checks required both tracks to pass the paired comparison and the whole-"
        "canonical population median comparison. Their failure is unchanged. The phastCons "
        "motif-scale requirement and the pooled canonical median requirement are superseded by "
        "the authorized criteria above; they are not erased or relabelled as passes.",
        "### Original pooled motif-window distributions",
        markdown(distributions),
        "### Within-element comparisons",
        markdown(paired),
        "The canonical phyloP paired comparison is unchanged and remains part of the corrected "
        "gate. The phastCons and noMotif paired rows are retained as descriptive diagnostics.",
        "## Secondary full-element phastCons",
    ]
    if passed:
        partial = correlations.partial_spearman_given_fimo
        if len(partial) == len(CONTEXTS) and partial.notna().all() and partial.abs().lt(0.1).all():
            lines.insert(
                2,
                "**Headline descriptive result: after controlling for FIMO, "
                "conservation–repression "
                f"correlations are near zero in all {len(partial)} contexts "
                f"(partial Spearman {partial.min():.6f} to {partial.max():.6f}).** "
                "These data show little residual monotonic association in this preview; "
                "incremental predictive value remains a question for held-out modeling.",
            )
    full = features[features.phastCons_element_coverage_status.eq("complete")]
    secondary = (
        full.groupby("population")
        .phastCons_element_mean.agg(
            n="count",
            mean="mean",
            median="median",
            p05=lambda x: x.quantile(0.05),
            p95=lambda x: x.quantile(0.95),
        )
        .reset_index()
    )
    lines += [
        markdown(secondary),
        "These distributions describe fully covered 200-base elements and do not gate the "
        "motif analysis. All extracted phastCons windows remain in the feature CSV for audit; "
        "the joined modeling dataset includes phastCons **full-element** columns only.",
        "## Coverage and windows",
        markdown(coverage),
        "Motif offsets (90, 110) are zero-based inclusive: 21 bases, with 179 complementary "
        "flank bases in each 200-base element. Reverse-strand windows are oriented before "
        "slicing. Means, maxima and positive fractions use observed bases, with explicit "
        "covered/requested counts and complete/partial/no-data status. No-data values stay "
        "missing; negative phyloP values are retained. The gate and preview use complete "
        "windows. [Unchanged feature CSV](../../data/processed/conservation_features.csv).",
    ]
    if passed:
        lines += [
            "## Join assertion",
            f"**PASS: {len(dataset):,} source rows → {len(joined):,} joined rows.** "
            "The many-to-one element-ID join preserves element × cell × enhancer keys, order, "
            "all original columns and missingness. No enhancers are aggregated or activitie"
            "s imputed.",
            "Rows without verified coordinates: "
            f"**{joined.missing_conservation_coordinates.sum():,}**; "
            f"rows without a phyloP motif mean: **{joined.missing_phyloP_motif_mean.sum():,}**. "
            "Controls and unverified native coordinates remain present with missing features.",
            "[Joined dataset](../../data/processed/analysis_dataset_conservation.csv).",
            "## Descriptive preview — not modeling",
            "All canonical elements, above and below 20.86, are eligible in each of the five fixed "
            "passing contexts. Repression is the same-cell/same-enhancer NegCtrl median minus "
            "activity; positive means greater repression. This is not scrambled dynamic range "
            "or cross-context robustness. Raw and adjusted correlations use the same complete "
            "cases: fully covered phyloP motif window, measured repression and FIMO. "
            "Missingness counts overlap.",
            "Partial Spearman is Pearson correlation between residual ranks of phyloP and "
            "repression after separately regressing each on ranked FIMO with an intercept; "
            "ties use average ranks. This is residual rank association, not causal evidence "
            "or held-out predictive improvement.",
            markdown(correlations),
        ]
        for row in correlations.itertuples(index=False):
            value = row.partial_spearman_given_fimo
            interpretation = (
                "unavailable"
                if not np.isfinite(value)
                else "near zero (the existing descriptive |ρ| < 0.10 rule)"
                if abs(value) < 0.1
                else "positive residual association"
                if value > 0
                else "negative residual association"
            )
            lines.append(
                f"**{row.cell_type}/{row.enhancer}:** raw ρ = "
                f"{row.spearman_phyloP_repression:.6f}; after FIMO, **ρ = {value:.6f}** — "
                f"{interpretation}."
            )
    else:
        lines += [
            "## Join and descriptive preview: NOT RUN",
            "The corrected gate failed. No join or partial correlations were computed.",
        ]
    tracks = pd.DataFrame([dict(track=k, **v) for k, v in identities.items()])
    lines += [
        "## Track provenance and reproducibility",
        f"[{config.CONSERVATION_TRACK}]({config.CONSERVATION_TRACK_SOURCE}): "
        f"**{config.CONSERVATION_ALIGNED_GENOMES} assemblies from "
        f"{config.CONSERVATION_SPECIES_COUNT} distinct species**, as documented by UCSC. "
        "The earlier metadata incorrectly called 470 a species count. Correcting this label "
        "does not change the files or their values. Both tracks and their exact URLs remain fixed.",
        markdown(tracks),
        "The original extraction used verified HTTP byte ranges after HTTPS timeouts, "
        "never a genome-wide bigWig download. This revision is completely offline: "
        "array/receipt hashes, coordinate hashes and the original extraction manifest "
        "are checked locally. Remote identities shown above are the recorded extraction "
        "identities, not freshly queried headers. Missing or changed caches fail; there is "
        "no automatic remote fallback.",
        "Run `make conservation PYTHON=.venv/bin/python`. "
        "[Original extraction manifest](../archive/conservation_original_gate/conservat"
        "ion_manifest.json) · "
        "[Two-run receipt](conservation_reproducibility.json). "
        "The manifest records the user-authorized revision, original failed verdict, fixed "
        "m-value and unchanged feature hash. No data are staged or committed.",
    ]
    return "\n\n".join(lines) + "\n"


def main():
    elements = verified_elements()
    regions_path = config.TABLES_DIR / "gate_operating_regions.csv"
    regions = pd.read_csv(regions_path)
    passing = set(
        regions.loc[regions.usable, ["cell_type", "enhancer"]].itertuples(index=False, name=None)
    )
    if passing != set(CONTEXTS):
        raise ValueError("Passing contexts changed; review the descriptive preview scope")
    dataset_path = config.PROCESSED_DIR / "analysis_dataset.csv"
    dataset_manifest = config.PROCESSED_DIR / "analysis_dataset_manifest.json"
    expected = json.loads(dataset_manifest.read_text())["outputs"][
        str(dataset_path.relative_to(config.ROOT))
    ]
    if digest(dataset_path) != expected:
        raise ValueError("Original analysis dataset failed its recorded hash")
    dataset = pd.read_csv(dataset_path, low_memory=False)
    arrays, absent, identities, snapshots = load_cached(
        elements, ORIGINAL / "conservation_manifest.json"
    )
    features = summarize(elements, arrays, absent)
    original_manifest = json.loads((ORIGINAL / "conservation_manifest.json").read_text())
    feature_hash = hashlib.sha256(
        features.to_csv(index=False, float_format="%.17g").encode()
    ).hexdigest()
    if feature_hash != original_manifest["outputs"][str(FEATURES.relative_to(config.ROOT))]:
        raise ValueError("Recomputed cached features differ from the original extraction")
    coverage = coverage_table(features)
    paired, distributions, original_verdicts = sanity_checks(features)
    verdicts, observed, strata = corrected_checks(features, dataset, paired)
    passed = bool(verdicts.passed.all())
    joined = join_features(dataset, features, verdicts) if passed else None
    correlations = preview(joined) if passed else pd.DataFrame()
    gate_inputs = []
    gate_table = gate_cases = gate_summary = None
    if passed:
        gate_manifest_path = config.PROCESSED_DIR / "gate_metrics_manifest.json"
        gate_manifest = json.loads(gate_manifest_path.read_text())
        gate_path = config.PROCESSED_DIR / "gate_metrics.csv"
        candidate_path = config.TABLES_DIR / "gate_cross_context_candidates.csv"
        for path in (gate_path, candidate_path):
            if digest(path) != gate_manifest["outputs"][str(path.relative_to(config.ROOT))]:
                raise ValueError(f"Gate outcome source hash mismatch: {path}")
        gate_table, gate_cases, gate_summary = conservation_gate_preview.analyze(
            features,
            pd.read_csv(gate_path, low_memory=False),
            pd.read_csv(candidate_path),
            partial_spearman,
        )
        gate_inputs = [gate_manifest_path, gate_path, candidate_path]

    text = render(
        features,
        coverage,
        paired,
        distributions,
        verdicts,
        identities,
        dataset,
        joined,
        correlations,
        original_verdicts,
        observed,
        strata,
    )
    if passed:
        gate_text = conservation_gate_preview.render(gate_table, gate_summary)
        primary = gate_table[
            gate_table.outcome.eq("dynamic_range") & gate_table.cohort.eq("above_m_value")
        ].partial_spearman_given_fimo
        robustness = gate_table[gate_table.outcome.eq("minimum_rank_percentile")].iloc[0]
        gate_headline = (
            "**Gate-outcome preview (existing above-threshold candidates):** "
            f"FIMO-adjusted phyloP correlations are {primary.min():.6f}–{primary.max():.6f} "
            f"for scrambled-pair dynamic range and {robustness.partial_spearman_given_fimo:.6f} "
            "for cross-context minimum percentile. These endpoints differ from "
            "NegCtrl-relative repression; the all-score scrambled comparison "
            "is reported separately."
        )
        text = text.replace(
            "## Why the gate was revised", gate_headline + "\n\n## Why the gate was revised"
        )
        text = text.replace(
            "## Track provenance and reproducibility",
            gate_text + "\n\n## Track provenance and reproducibility",
        )
    outputs = {
        FEATURES: features,
        config.TABLES_DIR / "conservation_strata.csv": strata,
        config.TABLES_DIR / "conservation_stratum_values.csv": observed,
        config.TABLES_DIR / "conservation_original_checks.csv": original_verdicts,
        config.TABLES_DIR / "conservation_coverage.csv": coverage,
        config.TABLES_DIR / "conservation_paired.csv": paired,
        config.TABLES_DIR / "conservation_populations.csv": distributions,
        config.TABLES_DIR / "conservation_sanity_checks.csv": verdicts,
    }
    if passed:
        outputs[JOINED] = joined
        outputs[config.TABLES_DIR / "conservation_gate_preview.csv"] = gate_table
        outputs[config.TABLES_DIR / "conservation_gate_preview_cases.csv"] = gate_cases
        outputs[config.TABLES_DIR / "conservation_preview.csv"] = correlations
    inputs = [
        config.HG38_ELEMENTS,
        ORIGINAL / "conservation_manifest.json",
        ORIGINAL / "conservation_sanity_checks.csv",
        dataset_path,
        dataset_manifest,
        regions_path,
        config.PROCESSED_DIR / "liftover_hg19_to_hg38/element_audit.csv",
        *[
            config.ROOT / "src" / name
            for name in (
                "conservation.py",
                "conservation_io.py",
                "conservation_pipeline.py",
                "conservation_strata.py",
                "conservation_gate_preview.py",
                "config.py",
                "dataset.py",
                "cross_context.py",
            )
        ],
        *snapshots,
        *gate_inputs,
    ]
    with TemporaryDirectory(prefix="conservation-publish-", dir=config.PROCESSED_DIR) as temp:
        staging = Path(temp)
        staged = {}
        for path, table in outputs.items():
            table.to_csv(staging / path.name, index=False, float_format="%.17g")
            staged[path] = staging / path.name
        (staging / REPORT.name).write_text(text)
        staged[REPORT] = staging / REPORT.name
        for path in distribution_figure(observed, staging):
            staged[config.FIGURES_DIR / path.name] = path
        manifest = dict(
            gate_version="phyloP_score_stratified_v2",
            gate_outcome_preview=gate_summary,
            gate_revision_reason=(
                "User-authorized check-design correction: phyloP motif-scale primary; "
                "above-m-value canonical population comparison; phastCons secondary at full-ele"
                "ment scale"
            ),
            original_gate_passed=bool(original_verdicts.passed.all()),
            remote_reads_performed=False,
            feature_values_unchanged=True,
            feature_sha256=feature_hash,
            m_value=config.M_VALUE,
            tracks=identities,
            track_name=config.CONSERVATION_TRACK,
            aligned_genomes=config.CONSERVATION_ALIGNED_GENOMES,
            distinct_species=config.CONSERVATION_SPECIES_COUNT,
            motif_span_zero_based_inclusive=list(config.MOTIF_SPAN),
            element_length=config.ELEMENT_LENGTH,
            sanity_gate_passed=passed,
            join_performed=passed,
            source_rows=len(dataset),
            output_rows=len(joined) if passed else None,
            inputs={str(p.relative_to(config.ROOT)): digest(p) for p in inputs},
            outputs={str(p.relative_to(config.ROOT)): digest(q) for p, q in staged.items()},
        )
        (staging / MANIFEST.name).write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        staged[MANIFEST] = staging / MANIFEST.name
        for path, stored in staged.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            stored.replace(path)
    print(verdicts.to_string(index=False))
    print(
        correlations.to_string(index=False)
        if passed
        else "STOP: sanity gate failed; no join or preview"
    )
    print(f"Report: {REPORT.relative_to(config.ROOT)}")
    if not passed:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
