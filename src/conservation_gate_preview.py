"""Offline descriptive conservation associations with scrambled-pair gate outcomes."""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src import config
from src.cross_context import CONTEXTS
from src.dataset import KEY, markdown


def correlate(cases, outcome, partial):
    """Identical complete cases for raw/adjusted coefficients; explicit exclusion flags."""
    cases = cases.copy()
    cases["missing_conservation"] = cases.phyloP_motif_mean.isna()
    cases["incomplete_motif"] = ~cases.phyloP_motif_coverage_status.eq("complete")
    cases["missing_outcome"] = cases[outcome].isna()
    cases["missing_fimo"] = cases.fimo_score.isna()
    cases["missing_pair"] = (
        cases[["on_activity", "off_activity"]].isna().any(axis=1)
        if outcome == "dynamic_range"
        else False
    )
    flags = [
        "missing_conservation",
        "incomplete_motif",
        "missing_outcome",
        "missing_fimo",
        "missing_pair",
    ]
    cases["complete_case"] = ~cases[flags].any(axis=1)
    for column in ("phyloP_motif_mean", outcome, "fimo_score"):
        if not np.isfinite(cases[column].dropna()).all():
            raise ValueError(f"Nonfinite {column} in gate preview")
    valid = cases[cases.complete_case]
    x, y, z = (valid[c] for c in ("phyloP_motif_mean", outcome, "fimo_score"))

    def rho(a, b):
        return (
            float(spearmanr(a, b).statistic)
            if len(a) > 2 and a.nunique() > 1 and b.nunique() > 1
            else np.nan
        )

    raw, adjusted, fimo = rho(x, y), partial(x, y, z), rho(z, y)
    result = dict(
        n_source=len(cases),
        n_complete_cases=len(valid),
        n_excluded=len(cases) - len(valid),
        **{name: int(cases[name].sum()) for name in flags},
        spearman_phyloP_outcome=raw,
        partial_spearman_given_fimo=adjusted,
        spearman_fimo_outcome_same_cases=fimo,
        fraction_residual_rank_variance=adjusted**2,
        incremental_rank_r_squared=(1 - fimo**2) * adjusted**2,
    )
    return result, cases


def analyze(features, frame, candidates, partial):
    """Never rerank after conservation exclusions; one robustness observation per element."""
    if features.element_id.duplicated().any() or frame.duplicated(KEY).any():
        raise ValueError("Duplicate conservation element or gate observation")
    if candidates.element_id.duplicated().any():
        raise ValueError("Robustness requires one value per element")
    if (
        not candidates.n_contexts_measured.ge(4).all()
        or not candidates.fimo_score.ge(config.M_VALUE).all()
    ):
        raise ValueError("Unexpected cross-context candidate eligibility")
    scr = frame[frame.reference_type.eq("scrambled") & frame.population.eq("motif")].copy()
    if not set(candidates.element_id) <= set(scr.element_id):
        raise ValueError("Robustness candidates must have scrambled references")
    scores = scr[["element_id", "fimo_score"]].drop_duplicates()
    if scores.element_id.duplicated().any():
        raise ValueError("Inconsistent motif scores across gate contexts")
    joined_scores = candidates[["element_id", "fimo_score"]].merge(
        scores, on="element_id", validate="1:1"
    )
    if not np.array_equal(joined_scores.fimo_score_x, joined_scores.fimo_score_y, equal_nan=True):
        raise ValueError("FIMO disagreement between candidate and gate source")
    cols = ["element_id", "phyloP_motif_mean", "phyloP_motif_coverage_status"]
    scr = scr.merge(features[cols], on="element_id", how="left", validate="m:1")
    rows, audits = [], []
    for cohort in ("above_m_value", "all_canonical_scores"):
        for cell, enhancer in CONTEXTS:
            group = scr[scr.cell_type.eq(cell) & scr.enhancer.eq(enhancer)]
            if cohort == "above_m_value":
                group = group[group.fimo_score.ge(config.M_VALUE)]
            used = group.dropna(subset=["dynamic_range", "on_activity", "off_activity"])
            if not np.allclose(used.dynamic_range, used.on_activity - used.off_activity):
                raise ValueError("Dynamic range disagrees with scrambled ON minus intact OFF")
            keep = [
                *KEY,
                "fimo_score",
                "on_activity",
                "off_activity",
                "dynamic_range",
                "phyloP_motif_mean",
                "phyloP_motif_coverage_status",
            ]
            result, cases = correlate(group[keep], "dynamic_range", partial)
            labels = dict(outcome="dynamic_range", cohort=cohort, cell_type=cell, enhancer=enhancer)
            rows.append(dict(**labels, **result))
            audits.append(cases.assign(**labels))
    robust = candidates.merge(features[cols], on="element_id", how="left", validate="1:1")
    keep = [
        "element_id",
        "fimo_score",
        "n_contexts_measured",
        "minimum_rank_percentile",
        "phyloP_motif_mean",
        "phyloP_motif_coverage_status",
    ]
    result, cases = correlate(robust[keep], "minimum_rank_percentile", partial)
    labels = dict(
        outcome="minimum_rank_percentile",
        cohort="above_m_value_coverage_4_of_5",
        cell_type="across_contexts",
        enhancer="five_passing_contexts",
    )
    rows.append(dict(**labels, **result))
    audits.append(cases.assign(**labels))
    baseline = float(spearmanr(candidates.minimum_rank_percentile, candidates.fimo_score).statistic)
    summary = dict(original_candidate_n=len(candidates), original_fimo_robustness_spearman=baseline)
    return pd.DataFrame(rows), pd.concat(audits, ignore_index=True), summary


def render(table, summary):
    """Separate candidate endpoint results from the NegCtrl repression preview."""
    dynamic = table[table.outcome.eq("dynamic_range")]
    primary = dynamic[dynamic.cohort.eq("above_m_value")]
    robust = table[table.outcome.eq("minimum_rank_percentile")].iloc[0]
    columns = [
        "cell_type",
        "enhancer",
        "n_source",
        "n_complete_cases",
        "n_excluded",
        "spearman_phyloP_outcome",
        "partial_spearman_given_fimo",
    ]
    adjusted = robust.partial_spearman_given_fimo
    description = (
        "unavailable"
        if not np.isfinite(adjusted)
        else "near zero under the existing descriptive |ρ| <0.10 rule"
        if abs(adjusted) < 0.1
        else "a positive residual rank association"
        if adjusted > 0
        else "a negative residual rank association"
    )
    lines = [
        "## Conservation versus gate outcomes — descriptive preview",
        f"**Cross-context robustness: raw Spearman {robust.spearman_phyloP_outcome:.6f}; "
        f"FIMO-adjusted Spearman {adjusted:.6f}, n={int(robust.n_complete_cases):,} "
        f"— {description}.**",
        "Primary dynamic-range analysis follows the candidate rankings: canonical FIMO ≥20.86, "
        "scrambled references only, in the same five passing contexts. A separate all-score "
        "scrambled cohort matches the wider score scope of the existing repression preview. "
        "NegCtrl fallback is excluded from both. Dynamic range is scrambled ON minus intact OFF "
        "in log2 units; both states must be measured.",
        "The same complete-case rule and partial-correlation implementation as the repression "
        "preview are used: fully covered 21-base phyloP window, measured outcome and finite FIMO. "
        "Raw and partial coefficients use identical rows. Partial Spearman is Pearson correlation "
        "of residual average-tie ranks after regressing each ranked variable on ranked FIMO "
        "plus an intercept. Missing measurements remain missing.",
        "### Primary per-context dynamic range (above m-value)",
        markdown(primary[columns]),
        "### Cross-context minimum rank percentile (one row per element)",
        "Use the existing above-threshold, scrambled-only candidates measured in at least four "
        "of five contexts. The saved minimum percentile is retained: ranks and their original "
        "context denominators are **not recomputed** after conservation exclusions. Missing "
        "contexts are not assigned a worst rank.",
        markdown(table[table.outcome.eq("minimum_rank_percentile")][columns]),
        f"The original {summary['original_candidate_n']:,} candidates give FIMO-versus-robustness "
        f"Spearman **{summary['original_fimo_robustness_spearman']:.6f}**. In the "
        f"{int(robust.n_complete_cases):,} conservation-complete cases that coefficient is "
        f"**{robust.spearman_fimo_outcome_same_cases:.6f}**; "
        "this is the matched cohort for adjustment.",
        "### Does conservation account for remaining variation?",
        f"For a descriptive linear regression on **ranks**, squared partial correlation is "
        f"**{100 * robust.fraction_residual_rank_variance:.4f}%** of the rank variance remaining "
        "after FIMO. Adding ranked phyloP increases in-sample rank R² by "
        f"**{100 * robust.incremental_rank_r_squared:.4f} percentage points**. "
        "These quantities concern ranks, not raw activity variance or held-out predictive gain. "
        "The original 0.58 Spearman coefficient does not establish that two thirds of biological "
        "activity variance is unexplained.",
        (
            "Conservation contributes very little residual monotonic information about robustness "
            "in this preview; this is not evidence of exactly zero association."
            if np.isfinite(adjusted) and abs(adjusted) < 0.1
            else "A residual association is observed; its predictive usefulness "
            "must be tested on held-out data."
            if np.isfinite(adjusted)
            else "Residual association cannot be estimated from the available cases."
        ),
        "### All-score scrambled dynamic-range comparison",
        markdown(dynamic[dynamic.cohort.eq("all_canonical_scores")][columns]),
        "### Coverage and exclusions",
        "Source counts include designed scrambled references whose activity may be missing. "
        "Missingness reasons overlap; n_excluded is the union, not their sum. "
        "Incomplete motif coverage includes partial, no-data and unavailable coordinate windows.",
        markdown(
            table[
                [
                    "outcome",
                    "cohort",
                    "cell_type",
                    "enhancer",
                    "n_source",
                    "n_complete_cases",
                    "missing_conservation",
                    "incomplete_motif",
                    "missing_outcome",
                    "missing_fimo",
                    "missing_pair",
                ]
            ]
        ),
        "[Full coefficients and rank-variance calculations](conservation_gate_preview.csv) · "
        "[Per-element inclusion/exclusion audit](conservation_gate_preview_cases.csv).",
    ]
    return "\n\n".join(lines)
