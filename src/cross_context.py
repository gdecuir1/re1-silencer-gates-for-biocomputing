"""Worst-context percentile prioritization of measured scrambled pairs.

Raw effects remain context-specific. Cross-context summaries describe consistency;
missing contexts are never imputed, and no metric composite is constructed.
"""

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src import config
from src.dataset import KEY, markdown

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Explicit study scope requested by the user; do not add scrambled-only screen passes.
CONTEXTS = (
    ("HepG2", "En19"),
    ("K562", "En09"),
    ("K562", "En11"),
    ("K562", "En19"),
    ("SK-N-SH", "En19"),
)
MIN_CONTEXTS = 4


def reference_rankings(ranked, reference):
    """Recompute ordinal metric ranks within context AND reference type."""
    result = ranked[ranked.reference_type.eq(reference)].copy()
    for metric, ascending in (("dynamic_range", False), ("separation", False), ("leakiness", True)):
        ordered = result.dropna(subset=[metric]).sort_values(
            ["cell_type", "enhancer", metric, "element_id"],
            ascending=[True, True, ascending, True],
            kind="stable",
        )
        ranks = ordered.groupby(["cell_type", "enhancer"], sort=False).cumcount() + 1
        result[f"rank_{metric}"] = ranks.reindex(result.index).astype("Int64")
    return result.sort_values(["cell_type", "enhancer", "rank_dynamic_range"], kind="stable")


def analyze(frame, contexts=CONTEXTS, min_contexts=MIN_CONTEXTS):
    """Average-tie ascending percentile = 100 * rank / measured context n.

    Population: canonical FIMO >= M_VALUE with designed scrambled references.
    Rank on minimum measured percentile, then median, coverage, and stable ID.
    Full audit retains elements with no measured pair in the selected contexts.
    """
    if frame.duplicated(KEY).any():
        raise ValueError("Duplicate element/cell/enhancer observations")
    if len(set(contexts)) != len(contexts) or not 1 <= min_contexts <= len(contexts):
        raise ValueError("Invalid context list or minimum coverage")
    eligible = frame[frame.reference_type.eq("scrambled") & frame.fimo_score.ge(config.M_VALUE)]
    if eligible.groupby("element_id").fimo_score.nunique().gt(1).any():
        raise ValueError("An element has inconsistent motif scores")
    ids = sorted(eligible.element_id.unique())
    audit = (
        eligible[["element_id", "fimo_score"]]
        .drop_duplicates()
        .set_index("element_id")
        .reindex(ids)
    )
    raw_cols, pct_cols, off_cols, denominators = [], [], [], []
    for cell, enhancer in contexts:
        tag = f"{cell}__{enhancer}"
        g = eligible[eligible.cell_type.eq(cell) & eligible.enhancer.eq(enhancer)].copy()
        valid = g[["on_activity", "off_activity", "dynamic_range"]].notna().all(axis=1)
        measured = g.loc[valid].set_index("element_id")
        if not np.isfinite(measured[["on_activity", "off_activity", "dynamic_range"]]).all().all():
            raise ValueError("Infinite paired measurement")
        if not np.allclose(measured.dynamic_range, measured.on_activity - measured.off_activity):
            raise ValueError("Dynamic range disagrees with the measured pair")
        raw, pct, off = f"dynamic_range__{tag}", f"percentile__{tag}", f"off_activity__{tag}"
        audit[raw] = measured.dynamic_range.reindex(ids)
        audit[pct] = (100 * measured.dynamic_range.rank(method="average", pct=True)).reindex(ids)
        audit[off] = measured.off_activity.reindex(ids)
        raw_cols.append(raw)
        pct_cols.append(pct)
        off_cols.append(off)
        denominators.append(
            dict(cell_type=cell, enhancer=enhancer, n_measured_scrambled=len(measured))
        )
    audit["n_contexts_measured"] = audit[raw_cols].notna().sum(axis=1)
    audit["median_rank_percentile"] = audit[pct_cols].median(axis=1)
    audit["minimum_rank_percentile"] = audit[pct_cols].min(axis=1)
    audit["dynamic_range_min"] = audit[raw_cols].min(axis=1)
    audit["dynamic_range_max"] = audit[raw_cols].max(axis=1)
    audit["dynamic_range_spread"] = audit.dynamic_range_max - audit.dynamic_range_min
    audit["dynamic_range_sd"] = audit[raw_cols].std(axis=1, ddof=1)
    audit["median_off_activity"] = audit[off_cols].median(axis=1)
    audit["max_off_activity"] = audit[off_cols].max(axis=1)
    audit["coverage_qualified"] = audit.n_contexts_measured.ge(min_contexts)
    audit["missing_contexts"] = audit.apply(
        lambda row: ", ".join(
            f"{cell}/{enh}"
            for (cell, enh), col in zip(contexts, raw_cols, strict=True)
            if pd.isna(row[col])
        ),
        axis=1,
    )
    audit = audit.reset_index()
    candidates = (
        audit[audit.coverage_qualified]
        .sort_values(
            [
                "minimum_rank_percentile",
                "median_rank_percentile",
                "n_contexts_measured",
                "element_id",
            ],
            ascending=[False, False, False, True],
            kind="stable",
        )
        .copy()
    )
    candidates.insert(0, "cross_context_rank", range(1, len(candidates) + 1))
    ranked_ids = candidates.set_index("element_id").cross_context_rank
    audit["cross_context_rank"] = audit.element_id.map(ranked_ids).astype("Int64")
    coverage = (
        audit.groupby("n_contexts_measured").size().reindex(range(len(contexts) + 1), fill_value=0)
    )
    coverage = coverage.rename("elements").rename_axis("n_contexts_measured").reset_index()
    valid = (
        len(candidates) > 2
        and candidates.minimum_rank_percentile.nunique() > 1
        and candidates.fimo_score.nunique() > 1
    )
    rho = (
        float(spearmanr(candidates.minimum_rank_percentile, candidates.fimo_score).statistic)
        if valid
        else None
    )
    summary = dict(
        n_above_threshold=frame[frame.fimo_score.ge(config.M_VALUE)].element_id.nunique(),
        n_scrambled_eligible=len(audit),
        n_qualified=len(candidates),
        n_excluded_coverage=int((~audit.coverage_qualified).sum()),
        min_contexts=min_contexts,
        contexts=[f"{a}/{b}" for a, b in contexts],
        correlation_n=len(candidates),
        spearman_min_percentile_fimo=rho,
        percentile_definition=(
            "100 * ascending average-tie rank / n measured scrambled pairs in context"
        ),
        dynamic_range_spread_definition="max minus min over measured contexts; log2 units",
        dynamic_range_sd_definition="sample SD, ddof=1, over measured contexts; log2 units",
    )
    return candidates, audit, coverage, pd.DataFrame(denominators), summary


def report(candidates, coverage, denominators, summary):
    """Report coverage, worst-context ranking, raw effects and OFF levels separately."""
    top = candidates.head(25)
    lines = [
        "## Cross-context candidate selection: weakest measured percentile",
        "Primary analysis uses **scrambled-paired, FIMO ≥20.86** observations only, in "
        + ", ".join(summary["contexts"])
        + ". The operating screen and these five contexts are "
        "held fixed; no additional contexts are selected to improve a candidate's score.",
        "Within each context, rank percentile = **100 × ascending average-tie dynamic-range "
        "rank / measured scrambled-pair count**. Higher is better; the maximum possible is 100, "
        "the lowest untied value is 100/n, and ties share their average percentile. Missing pairs "
        "are excluded from that context's denominator and remain NA for the element. "
        "Percentiles are computed using all measured eligible scrambled pairs in that context "
        "before applying the cross-context coverage filter.",
        markdown(denominators),
        f"Of **{summary['n_above_threshold']:,}** above-threshold canonical elements, "
        f"**{summary['n_scrambled_eligible']:,}** have a designed scrambled reference; "
        f"**{summary['n_above_threshold'] - summary['n_scrambled_eligible']:,}** fallback elements "
        "are outside the primary cross-context population. "
        f"**{summary['n_qualified']:,} qualify** with measurements in at least four of five "
        f"contexts; **{summary['n_excluded_coverage']:,} are excluded for coverage**, including "
        "elements with no measured pair in these contexts.",
        markdown(coverage),
        "Candidates are ranked primarily by **minimum measured percentile**, then median "
        "percentile, number of measured contexts, and element ID. This is a worst-context "
        "criterion applied to one metric, not a composite of dynamic range, leakiness and "
        "separation. A four-context candidate has an untested fifth context; no worst rank "
        "or activity is imputed for it. Counts and missing context names remain visible.",
        "Raw magnitudes are reported separately: **spread = max(dynamic range) − min(dynamic "
        "range)** and **SD = sample standard deviation (ddof=1)** over measured contexts, "
        "both in log2 units. These describe cross-context variability, not replicate "
        "uncertainty. A stable rank does not guarantee a stable effect magnitude.",
        "[Full ranked candidates](gate_cross_context_candidates.csv) · "
        "[All eligible elements and coverage exclusions](gate_cross_context_coverage.csv).",
        "### Motif score and robustness",
    ]
    rho = summary["spearman_min_percentile_fimo"]
    if rho is None:
        lines.append(
            "Correlation is unavailable because there are too few or constant observations."
        )
    else:
        strength = "strong" if abs(rho) >= 0.6 else "moderate" if abs(rho) >= 0.3 else "weak"
        lines.append(
            f"Spearman correlation of minimum percentile with FIMO score: "
            f"**ρ = {rho:.6f}, n = {summary['correlation_n']:,}** qualified elements. "
            f"This is a {strength} rank association (descriptive |ρ| bands: weak <0.30, "
            "moderate 0.30–<0.60, strong ≥0.60). "
            + (
                "Motif strength tracks robustness, but correlation alone does not establish "
                "that it captures all useful candidate variation. "
                if strength == "strong"
                else "Motif strength only partially tracks robustness. Additional predictive "
                "value requires a held-out comparison; this association does not establish "
                "independence. "
            )
            + "Conservation has not been evaluated here; this correlation cannot establish "
            "whether conservation adds or duplicates the information."
        )
    columns = [
        "cross_context_rank",
        "element_id",
        "n_contexts_measured",
        "minimum_rank_percentile",
        "median_rank_percentile",
        "dynamic_range_spread",
        "dynamic_range_sd",
        "fimo_score",
        "max_off_activity",
        "missing_contexts",
    ]
    lines += [
        "### Top 25 cross-context candidates",
        markdown(
            top[columns].rename(columns={"max_off_activity": "maximum leakiness (= OFF, log2)"})
        ),
        "### Top candidates: raw dynamic range in every context (log2 units)",
    ]
    for prefix in ("dynamic_range", "off_activity"):
        if prefix == "off_activity":
            lines.append(
                "### Top candidates: leakiness = OFF activity in every context (log2 units)"
            )
        cols = [f"{prefix}__{cell}__{enh}" for cell, enh in CONTEXTS]
        lines.append(
            markdown(
                top[["cross_context_rank", "element_id", *cols]].rename(
                    columns={
                        col: f"{cell}/{enh}"
                        for col, (cell, enh) in zip(cols, CONTEXTS, strict=True)
                    }
                )
            )
        )
    if len(top):
        lines.append(
            f"Among the displayed candidates, maximum measured OFF activity ranges "
            f"from **{top.max_off_activity.min():.2f} to {top.max_off_activity.max():.2f}** "
            "log2 units. Large, consistently ranked dynamic range can coexist with "
            "substantial residual OFF activity; leakiness is shown but does not enter "
            "the robustness ranking."
        )
    lines += [
        "### Dynamic range across contexts: top 10 candidates",
        "![Dynamic range across contexts](../figures/gate_cross_context_top10.png)",
        "[Vector PDF](../figures/gate_cross_context_top10.pdf) · "
        "[Vector SVG](../figures/gate_cross_context_top10.svg). Each line is one candidate "
        "in worst-percentile rank order; missing contexts break the line and have no point. "
        "The dashed line marks the existing log2 range screen, not a fitted threshold.",
    ]
    return "\n\n".join(lines)


def figure(candidates, directory):
    """Fixed-order raw effects; preserve NaN gaps; deterministic print exports."""
    top = candidates.head(10)
    plt.rcParams.update({"font.size": 11, "svg.hashsalt": "re1-cross-context"})
    fig, ax = plt.subplots(figsize=(12, 8.5))
    colors = plt.get_cmap("tab10").colors
    markers = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h")
    columns = [f"dynamic_range__{cell}__{enh}" for cell, enh in CONTEXTS]
    for i, row in enumerate(top.itertuples(index=False)):
        values = top.iloc[i][columns].to_numpy(dtype=float)
        ax.plot(
            range(len(CONTEXTS)),
            values,
            marker=markers[i],
            color=colors[i],
            linewidth=1.6,
            markersize=6,
            label=f"{row.cross_context_rank}. {row.element_id}",
        )
    ax.axhline(1, color="#777777", linestyle="--", linewidth=1)
    ax.set_xticks(range(len(CONTEXTS)), [f"{a}\n{b}" for a, b in CONTEXTS])
    ax.set_ylabel("Dynamic range: scrambled ON − intact OFF (log2 units)")
    ax.set_xlabel("Cell type / enhancer context")
    ax.set_title(
        "Cross-context consistency of the top 10 candidates", loc="left", fontsize=16, pad=28
    )
    ax.text(
        0,
        1.025,
        "Ranked by minimum within-context percentile; ≥4/5 measured scrambled pairs",
        transform=ax.transAxes,
        fontsize=10,
    )
    ax.grid(axis="y", alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    if len(top):
        ax.legend(
            loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2, fontsize=8.5, frameon=False
        )
    else:
        ax.text(0.5, 0.5, "No coverage-qualified candidates", transform=ax.transAxes, ha="center")
    fig.subplots_adjust(left=0.09, right=0.98, top=0.87, bottom=0.32)
    paths = []
    for ext in ("png", "pdf", "svg"):
        path = directory / f"gate_cross_context_top10.{ext}"
        metadata = {"Creator": "RE1 gate metrics"}
        if ext == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        if ext == "svg":
            metadata.update(Date=None)
        fig.savefig(path, dpi=300, metadata=metadata)
        paths.append(path)
    plt.close(fig)
    return paths
