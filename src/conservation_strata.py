"""User-corrected conservation gate and complete score-stratified distributions."""

import matplotlib
import numpy as np
import pandas as pd

from src import config

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STRATA = ("canonical_above_m_value", "canonical_below_m_value", "noMotif")
LABELS = ("Canonical FIMO ≥20.86", "Canonical FIMO <20.86", "No canonical motif")
CHECKS = ("canonical motif versus own flanks", "above-threshold canonical versus noMotif")


def corrected_checks(features, dataset, paired):
    """Keep the original phyloP paired check; stratify the population median check at M_VALUE."""
    native = dataset[dataset.population.isin(["motif", "noMotif"])]
    metadata = native[["element_id", "population", "fimo_score"]].drop_duplicates()
    if metadata.element_id.duplicated().any():
        raise ValueError("Inconsistent source population or FIMO score across contexts")
    motif = metadata.population.eq("motif")
    if not np.isfinite(metadata.loc[motif, "fimo_score"]).all():
        raise ValueError("Canonical FIMO score missing or nonfinite; cannot assign a stratum")
    metadata["stratum"] = np.where(
        ~motif, "noMotif", np.where(metadata.fimo_score.ge(config.M_VALUE), STRATA[0], STRATA[1])
    )
    cols = [
        "element_id",
        "population",
        "phyloP_motif_mean",
        "phyloP_motif_coverage_status",
        "phyloP_motif_covered_bases",
        "phastCons_element_mean",
        "phastCons_element_coverage_status",
    ]
    if features.element_id.duplicated().any():
        raise ValueError("Duplicate feature element")
    observed = metadata.merge(
        features[cols], on=["element_id", "population"], how="left", validate="1:1", indicator=True
    )
    if (
        len(observed) != len(metadata)
        or not set(features.element_id) <= set(metadata.element_id)
        or observed["_merge"].eq("both").sum() != len(features)
    ):
        raise ValueError("Conservation stratum membership changed")
    observed["missing_verified_coordinates"] = observed["_merge"].eq("left_only")
    observed = observed.drop(columns="_merge").sort_values(["stratum", "element_id"])
    observed["phyloP_motif_coverage_status"] = observed.phyloP_motif_coverage_status.fillna(
        "no_verified_coordinates"
    )
    summaries = []
    for stratum in STRATA:
        group = observed[observed.stratum.eq(stratum)]
        for scope in ("complete", "any_observed"):
            status = group.phyloP_motif_coverage_status
            use = (
                status.eq("complete")
                if scope == "complete"
                else status.isin(["complete", "partial"])
            )
            values = group.loc[use, "phyloP_motif_mean"]
            summaries.append(
                dict(
                    stratum=stratum,
                    scope=scope,
                    source_elements=len(group),
                    missing_coordinates=int(group.missing_verified_coordinates.sum()),
                    partial=int(status.eq("partial").sum()),
                    no_data=int(status.eq("no_data").sum()),
                    n=len(values),
                    mean=values.mean(),
                    sd=values.std(),
                    min=values.min(),
                    p01=values.quantile(0.01),
                    p05=values.quantile(0.05),
                    p25=values.quantile(0.25),
                    median=values.median(),
                    p75=values.quantile(0.75),
                    p95=values.quantile(0.95),
                    p99=values.quantile(0.99),
                    max=values.max(),
                )
            )
    summary = pd.DataFrame(summaries)
    primary = summary[summary.scope.eq("complete")].set_index("stratum")
    old = paired[paired.track.eq("phyloP") & paired.population.eq("motif")].iloc[0]
    difference = primary.loc[STRATA[0], "median"] - primary.loc["noMotif", "median"]
    verdicts = pd.DataFrame(
        [
            dict(
                check=CHECKS[0],
                track="phyloP",
                expected="median difference > 0 and fraction positive > 0.5",
                observed_median_difference=old.median_difference,
                observed_fraction_positive=old.fraction_motif_greater,
                passed=bool(
                    old.n > 0 and old.median_difference > 0 and old.fraction_motif_greater > 0.5
                ),
            ),
            dict(
                check=CHECKS[1],
                track="phyloP",
                expected="above-threshold median minus noMotif median > 0",
                observed_median_difference=difference,
                observed_fraction_positive=np.nan,
                passed=bool(
                    primary.loc[[STRATA[0], "noMotif"], "n"].gt(0).all() and difference > 0
                ),
            ),
        ]
    )
    return verdicts, observed.reset_index(drop=True), summary


def distribution_figure(observed, directory):
    """Shared fixed-width histogram and exact empirical CDF; no fitted mixture or KDE."""
    primary = observed[observed.phyloP_motif_coverage_status.eq("complete")]
    values = primary.phyloP_motif_mean.to_numpy()
    if not len(values):
        raise ValueError("No complete motif windows for distribution plot")
    edges = np.arange(np.floor(values.min() * 4) / 4, np.ceil(values.max() * 4) / 4 + 0.251, 0.25)
    plt.rcParams.update({"font.size": 11, "svg.hashsalt": "re1-conservation-strata"})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = ("#087E8B", "#D2691E", "#6B5B95")
    for stratum, label, color in zip(STRATA, LABELS, colors, strict=True):
        x = np.sort(primary.loc[primary.stratum.eq(stratum), "phyloP_motif_mean"].to_numpy())
        if not len(x):
            continue
        axes[0].hist(
            x,
            bins=edges,
            density=True,
            histtype="step",
            linewidth=1.8,
            color=color,
            label=f"{label} (n={len(x):,})",
        )
        axes[1].step(
            np.r_[x[0], x],
            np.r_[0, np.arange(1, len(x) + 1) / len(x)],
            where="post",
            linewidth=1.8,
            color=color,
        )
    pooled = primary.loc[primary.population.eq("motif"), "phyloP_motif_mean"]
    axes[0].hist(
        pooled,
        bins=edges,
        density=True,
        histtype="step",
        color="#333333",
        linestyle="--",
        linewidth=1.2,
        label="All canonical (pooled)",
    )
    axes[0].set_ylabel("Probability density (per phyloP score unit)")
    axes[1].set_ylabel("Cumulative fraction of elements")
    axes[0].set_title("Observed score-stratified distributions")
    axes[1].set_title("Full empirical distributions")
    for ax in axes:
        ax.set_xlabel("Mean phyloP in 21-base motif window")
        ax.axvline(0, color="#888888", linewidth=0.8, linestyle=":")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=0.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False, fontsize=10)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.27, wspace=0.28)
    paths = []
    for extension in ("png", "pdf", "svg"):
        path = directory / f"conservation_strata.{extension}"
        metadata = {"Creator": "RE1 conservation analysis"}
        if extension == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        if extension == "svg":
            metadata.update(Date=None)
        fig.savefig(path, dpi=300, metadata=metadata)
        paths.append(path)
    plt.close(fig)
    return paths
