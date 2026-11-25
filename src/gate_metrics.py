"""Assess operating regions for canonical silencer/enhancer pairs, offline.

Never average enhancers. ON is the element's own designed scrambled counterpart;
only elements without a designed partner use the same-context NegCtrl median.
Missing scrambled measurements do not trigger fallback. OFF is intact activity.
These are reporter sequence states, not measured REST occupancy or a circuit test.
"""

import hashlib
import json
import os
import shutil
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.transforms import Bbox
from scipy.stats import spearmanr

from src import config, cross_context, metrics
from src._utils import write_fasta
from src.assembly_check import GEO_ARCHIVE, TABLE_S2, TABLE_S2_MEMBER
from src.dataset import (
    CELLS,
    KEY,
    RAW,
    VALIDATION,
    assert_source_rows,
    digest,
    markdown,
    read_source,
)
from src.motifs import score_sequences

matplotlib.use("Agg")
import matplotlib.pyplot as plt

POPULATIONS = ("motif", "NegCtrl", "scrambled", "scrambled,SNP")
# Descriptive screening cutoffs specified before inspecting these paired results.
MIN_ON = 1.0
MIN_DYNAMIC_RANGE = 1.0
BOOTSTRAPS = 2000
BOOTSTRAP_SEED = 1729
TOP_N = 25
REPORT = config.TABLES_DIR / "gate_metrics.md"
MANIFEST = config.PROCESSED_DIR / "gate_metrics_manifest.json"
PAIRS = config.PROCESSED_DIR / "gate_metrics.csv"


def pair_catalog(source, oligos):
    """Exact Ref/Scr identifier pairing plus independent published sequence checks.

    No chromosome coordinates are inferred from names. The observed designed edit
    extends through zero-based position 111; the scored core stays (90,110).
    """
    native = sorted(source.loc[source.population.eq("motif"), "element_id"].unique())
    scr = source[source.population.isin(("scrambled", "scrambled,SNP"))]
    if scr.groupby("element_id").population.nunique().gt(1).any():
        raise ValueError("A scrambled ID has conflicting populations")
    scr_ids = set(scr.element_id)
    if any(not x.endswith("_Ref") for x in native) or any(not x.endswith("_Scr") for x in scr_ids):
        raise ValueError("Unexpected published Ref/Scr identifier convention")
    if {x.removesuffix("_Scr") + "_Ref" for x in scr_ids} - set(native):
        raise ValueError("Scrambled construct has no exact native counterpart")
    aliases = {}
    for r in oligos[oligos.library.eq("WG_REST")].itertuples():
        tokens = r.oligo_name[1:-1].split(";") if r.oligo_name.startswith("(") else [r.oligo_name]
        for token in tokens:
            aliases.setdefault(token, []).append(r.sequence.upper())
    rows, cores = [], {}
    left, right = config.MOTIF_SPAN
    if (left, right) != (90, 110):
        raise ValueError("Gate pairing requires the validated canonical core")
    for element in native:
        candidate = element.removesuffix("_Ref") + "_Scr"
        exists = candidate in scr_ids
        row = dict(
            element_id=element,
            scrambled_id=candidate if exists else None,
            has_designed_partner=exists,
            pair_status="no_designed_partner",
            changed_bases=np.nan,
            edited_adjacent_base=pd.NA,
        )
        if exists:
            a, b = aliases.get(element, []), aliases.get(candidate, [])
            if len(a) != 1 or len(b) != 1:
                raise ValueError(f"Missing or ambiguous published oligo: {element}")
            a, b = a[0], b[0]
            if len(a) != config.ELEMENT_LENGTH or len(b) != config.ELEMENT_LENGTH:
                raise ValueError("Paired oligos must have the full element length")
            if set(a + b) - set("ACGT"):
                raise ValueError("Ambiguous paired oligo sequence")
            # The extra adjacent base is a measured library property, not a new motif offset.
            if a[:left] != b[:left] or a[right + 2 :] != b[right + 2 :]:
                raise ValueError(f"Changes outside the observed motif-edit region: {element}")
            if a[left : right + 1] == b[left : right + 1]:
                raise ValueError("A scrambled partner has an unchanged canonical core")
            row.update(
                pair_status="verified_sequence_pair",
                changed_bases=sum(x != y for x, y in zip(a, b, strict=True)),
                edited_adjacent_base=a[right + 1] != b[right + 1],
            )
            cores[candidate] = b[left : right + 1]
        rows.append(row)
    return pd.DataFrame(rows), cores


def control_baselines(source):
    """Same-context median plus bootstrap SE over observed NegCtrl constructs.

    This SE measures uncertainty of the control-population median, including
    between-control heterogeneity; it is not a biological replicate SD. Sorting
    inputs and seeding each context makes the bootstrap independent of row order.
    """
    rows = []
    for (cell, enhancer), g in source[source.population.eq("NegCtrl")].groupby(
        ["cell_type", "enhancer"], sort=True
    ):
        values = np.sort(g.raw_activity.dropna().to_numpy())
        median = np.median(values) if len(values) else np.nan
        se = np.nan
        if len(values) > 1:
            token = f"{BOOTSTRAP_SEED}:{cell}:{enhancer}".encode()
            seed = int.from_bytes(hashlib.sha256(token).digest()[:8], "big")
            rng = np.random.default_rng(seed)
            draws = rng.choice(values, size=(BOOTSTRAPS, len(values)), replace=True)
            se = float(np.median(draws, axis=1).std(ddof=1))
        rows.append(
            dict(
                cell_type=cell,
                enhancer=enhancer,
                negctrl_n=len(values),
                negctrl_baseline=median,
                negctrl_median_se=se,
            )
        )
    return pd.DataFrame(rows)


def assemble(source, catalog, baseline, scores):
    """Retain exactly all canonical source keys and explicit reference/missingness."""
    if source.duplicated(KEY).any():
        raise ValueError("Duplicate source observation")
    for column in ("raw_activity", "lfcSE"):
        if not np.isfinite(source[column].dropna()).all():
            raise ValueError(f"Nonfinite {column}")
    if source.lfcSE.dropna().lt(0).any():
        raise ValueError("Negative source standard error")
    canonical = source[source.population.eq("motif")].copy()
    result = canonical.rename(
        columns={"raw_activity": "off_activity", "lfcSE": "off_se", "padj": "off_padj"}
    )
    result = result.merge(catalog, on="element_id", how="left", validate="m:1")
    if result.has_designed_partner.isna().any():
        raise ValueError("Pair catalog does not cover canonical source elements")
    scrambled = source[source.population.isin(("scrambled", "scrambled,SNP"))].rename(
        columns={
            "element_id": "scrambled_id",
            "raw_activity": "scrambled_activity",
            "lfcSE": "scrambled_se",
            "padj": "scrambled_padj",
            "population": "scrambled_population",
            "source_construct_id": "scrambled_construct_id",
        }
    )
    result = result.merge(
        scrambled,
        on=["scrambled_id", "cell_type", "enhancer"],
        how="left",
        validate="m:1",
        indicator="scrambled_join",
    )
    result = result.merge(baseline, on=["cell_type", "enhancer"], how="left", validate="m:1")
    result["reference_type"] = np.where(result.has_designed_partner, "scrambled", "NegCtrl")
    result["reference_reason"] = np.select(
        [
            ~result.has_designed_partner,
            result.scrambled_join.ne("both"),
            result.scrambled_activity.isna(),
        ],
        ["no_designed_partner", "partner_context_absent", "partner_activity_missing"],
        default="measured_scrambled_partner",
    )
    result["on_activity"] = result.scrambled_activity.where(
        result.has_designed_partner, result.negctrl_baseline
    )
    result["on_se"] = result.scrambled_se.where(
        result.has_designed_partner, result.negctrl_median_se
    )
    result["on_uncertainty_basis"] = np.where(
        result.has_designed_partner, "published_lfcSE", "bootstrap_NegCtrl_median_SE"
    )
    result["dynamic_range"] = metrics.dynamic_range(result.on_activity, result.off_activity)
    result["leakiness"] = metrics.off_state_leakiness(result.off_activity, result.on_activity)
    result["pooled_se"] = np.sqrt((result.on_se.pow(2) + result.off_se.pow(2)) / 2)
    on = pd.DataFrame({"activity": result.on_activity, "standard_error": result.on_se})
    off = pd.DataFrame({"activity": result.off_activity, "standard_error": result.off_se})
    result["separation"] = metrics.separation(on, off)
    result["separation_basis"] = "RMS_standard_error_not_replicate_SD"
    result["separation_status"] = np.select(
        [
            result.on_activity.isna() | result.off_activity.isna(),
            result.on_se.isna() | result.off_se.isna(),
            result.pooled_se.eq(0),
        ],
        ["missing_activity", "missing_standard_error", "zero_pooled_standard_error"],
        default="measured",
    )
    if set(scores.element_id) != set(canonical.element_id):
        raise ValueError("Canonical score IDs do not exactly match the source")
    if not np.isfinite(scores.score).all():
        raise ValueError("Missing canonical motif scores")
    result = result.merge(
        scores[["element_id", "score"]].rename(columns={"score": "fimo_score"}),
        on="element_id",
        how="left",
        validate="m:1",
    )
    result["above_m_value"] = result.fimo_score.ge(config.M_VALUE)
    for field in (
        "on_activity",
        "off_activity",
        "on_se",
        "off_se",
        "dynamic_range",
        "leakiness",
        "separation",
    ):
        result[f"missing_{field}"] = result[field].isna()
    assert_source_rows(canonical, result)
    return result.sort_values(KEY, kind="stable").reset_index(drop=True)


def coverage_table(frame):
    """Separate designed, measured, paired, fallback, and above-threshold coverage."""
    rows = []
    for (cell, enhancer), g in frame.groupby(["cell_type", "enhancer"], sort=True):
        scr = g.reference_type.eq("scrambled")
        rows.append(
            dict(
                cell_type=cell,
                enhancer=enhancer,
                canonical_rows=len(g),
                designed_scrambled=int(scr.sum()),
                partner_context_present=int(g.scrambled_join.eq("both").sum()),
                scrambled_on_measured=int((scr & g.on_activity.notna()).sum()),
                scrambled_both_measured=int((scr & g.dynamic_range.notna()).sum()),
                fallback_assigned=int((~scr).sum()),
                fallback_both_measured=int((~scr & g.dynamic_range.notna()).sum()),
                missing_on=int(g.on_activity.isna().sum()),
                missing_off=int(g.off_activity.isna().sum()),
                missing_either=int(g.dynamic_range.isna().sum()),
                above_threshold=int(g.above_m_value.sum()),
                above_scrambled_both=int((g.above_m_value & scr & g.dynamic_range.notna()).sum()),
                above_fallback_both=int((g.above_m_value & ~scr & g.dynamic_range.notna()).sum()),
            )
        )
    return pd.DataFrame(rows)


def operating_regions(frame):
    """Context medians use the same above-threshold, complete ON/OFF pairs.

    Each reference stratum is also reported; a mixed-reference median is never
    presented as if all observations used a matched scrambled sequence.
    """
    rows, strata = [], []
    for (cell, enhancer), g in frame.groupby(["cell_type", "enhancer"], sort=True):
        high = g[g.above_m_value]
        paired = high.dropna(subset=["on_activity", "off_activity"])
        row = dict(
            cell_type=cell,
            enhancer=enhancer,
            n_above=len(high),
            n_pairs=len(paired),
            n_missing_pair=len(high) - len(paired),
            median_on=paired.on_activity.median(),
            median_off=paired.off_activity.median(),
            median_dynamic_range=paired.dynamic_range.median(),
            n_separation=int(paired.separation.notna().sum()),
        )
        for reference in ("scrambled", "NegCtrl"):
            subset = paired[paired.reference_type.eq(reference)]
            row[f"n_{reference}"] = len(subset)
            row[f"median_on_{reference}"] = subset.on_activity.median()
            row[f"median_dr_{reference}"] = subset.dynamic_range.median()
            strata.append(
                dict(
                    cell_type=cell,
                    enhancer=enhancer,
                    reference_type=reference,
                    n_pairs=len(subset),
                    median_on=subset.on_activity.median(),
                    median_dynamic_range=subset.dynamic_range.median(),
                )
            )
        rows.append(row)
    regions = pd.DataFrame(rows)
    regions["usable"] = regions.median_on.ge(MIN_ON) & regions.median_dynamic_range.ge(
        MIN_DYNAMIC_RANGE
    )
    regions["verdict"] = np.select(
        [
            regions.n_pairs.eq(0),
            regions.usable,
            regions.median_on.lt(MIN_ON) & regions.median_dynamic_range.lt(MIN_DYNAMIC_RANGE),
            regions.median_on.lt(MIN_ON),
        ],
        ["missing_pairs", "usable_screen", "low_ON_and_low_range", "low_ON"],
        default="low_dynamic_range",
    )
    regions["scrambled_only_usable"] = regions.median_on_scrambled.ge(
        MIN_ON
    ) & regions.median_dr_scrambled.ge(MIN_DYNAMIC_RANGE)
    return regions, pd.DataFrame(strata)


def rank_pairs(frame, regions):
    """Rank within each usable context; report overlap on both all and common rows."""
    selection = regions.loc[regions.usable, ["cell_type", "enhancer"]]
    ranked = frame[frame.above_m_value & frame.dynamic_range.notna()].merge(
        selection, on=["cell_type", "enhancer"], how="inner", validate="m:1"
    )
    # No context means an empty, explicitly reported ranking rather than invented candidates.
    pieces, agreements = [], []
    for (cell, enhancer), group in ranked.groupby(["cell_type", "enhancer"], sort=True):
        group = group.sort_values(["dynamic_range", "element_id"], ascending=[False, True]).copy()
        group["rank_dynamic_range"] = range(1, len(group) + 1)
        for metric, ascending in (("separation", False), ("leakiness", True)):
            ordered = group.dropna(subset=[metric]).sort_values(
                [metric, "element_id"], ascending=[ascending, True]
            )
            ranks = pd.Series(range(1, len(ordered) + 1), index=ordered.index)
            group[f"rank_{metric}"] = ranks.reindex(group.index).astype("Int64")
        pieces.append(group)
        for reference in ("all_explicit_references", "scrambled", "NegCtrl"):
            subset = (
                group
                if reference == "all_explicit_references"
                else group[group.reference_type.eq(reference)]
            )
            common = subset.dropna(subset=["separation"])
            dr = set(
                subset.sort_values(["dynamic_range", "element_id"], ascending=[False, True])
                .head(TOP_N)
                .element_id
            )
            sep = set(
                common.sort_values(["separation", "element_id"], ascending=[False, True])
                .head(TOP_N)
                .element_id
            )
            leak = set(subset.sort_values(["leakiness", "element_id"]).head(TOP_N).element_id)
            common_dr = set(
                common.sort_values(["dynamic_range", "element_id"], ascending=[False, True])
                .head(TOP_N)
                .element_id
            )
            denominator = min(len(dr), len(sep))
            rho = (
                float(spearmanr(common.dynamic_range, common.separation).statistic)
                if (
                    len(common) > 2
                    and common.dynamic_range.nunique() > 1
                    and common.separation.nunique() > 1
                )
                else np.nan
            )
            n_overlap = len(dr & sep)
            agreements.append(
                dict(
                    cell_type=cell,
                    enhancer=enhancer,
                    reference_type=reference,
                    n_dynamic_range=len(subset),
                    n_separation=len(common),
                    top_dr_n=len(dr),
                    top_separation_n=len(sep),
                    overlap_dr_separation=n_overlap,
                    overlap_fraction=n_overlap / denominator if denominator else np.nan,
                    common_eligible_overlap=len(common_dr & sep),
                    overlap_dr_low_leakiness=len(dr & leak),
                    spearman_dr_separation=rho,
                    agreement="unavailable"
                    if not denominator
                    else "low_top_rank_agreement"
                    if n_overlap / denominator < 0.5
                    else "at_least_half_overlap",
                )
            )
    if pieces:
        return pd.concat(pieces, ignore_index=True), pd.DataFrame(agreements)
    for name in ("rank_dynamic_range", "rank_leakiness", "rank_separation"):
        ranked[name] = pd.Series(dtype="Int64")
    return ranked, pd.DataFrame(columns=["cell_type", "enhancer", "reference_type", "agreement"])


def sensitivity(regions):
    rows = []
    for on in (0.0, 0.5, 1.0, 1.5, 2.0):
        for dr in (0.0, 0.5, 1.0, 1.5, 2.0):
            # At the origin require positivity; elsewhere use inclusive explicit cutoffs.
            keep = (regions.median_on.gt(on) if on == 0 else regions.median_on.ge(on)) & (
                regions.median_dynamic_range.gt(dr)
                if dr == 0
                else regions.median_dynamic_range.ge(dr)
            )
            labels = [f"{r.cell_type}/{r.enhancer}" for r in regions[keep].itertuples()]
            rows.append(
                dict(
                    min_on=on,
                    min_dynamic_range=dr,
                    n_usable=len(labels),
                    contexts=", ".join(labels),
                )
            )
    return pd.DataFrame(rows)


def operating_figure(regions, directory):
    """Print-sized labelled context plane, with deterministic vector/raster output."""
    plt.rcParams.update({"font.size": 11, "svg.hashsalt": "re1-operating-region"})
    fig, ax = plt.subplots(figsize=(11, 8))
    colors = dict(zip(CELLS.values(), ("#0072B2", "#D55E00", "#009E73", "#7A5195"), strict=True))
    markers = dict(
        zip(("En02", "En09", "En11", "En19", "En21"), ("o", "s", "^", "D", "P"), strict=True)
    )
    xmax = max(2.5, regions.median_on.max() + 1.1)
    xmin = min(-0.6, regions.median_on.min() - 0.8)
    ymin = min(-0.4, regions.median_dynamic_range.min() - 0.6)
    ymax = max(2.0, regions.median_dynamic_range.max() + 0.7)
    ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
    ax.fill_between([MIN_ON, xmax], MIN_DYNAMIC_RANGE, ymax, color="#e3f2e5", zorder=0)
    ax.axvline(MIN_ON, color="#555555", linestyle="--", linewidth=1)
    ax.axhline(MIN_DYNAMIC_RANGE, color="#555555", linestyle="--", linewidth=1)
    ax.axhline(0, color="#aaaaaa", linewidth=0.7)
    annotations = []
    for r in regions.itertuples():
        if not np.isfinite(r.median_on) or not np.isfinite(r.median_dynamic_range):
            continue
        ax.scatter(
            r.median_on,
            r.median_dynamic_range,
            c=colors[r.cell_type],
            marker=markers[r.enhancer],
            s=72,
            edgecolor="black" if r.usable else "white",
            linewidth=1,
            zorder=4,
        )
        annotations.append(
            ax.annotate(
                f"{r.cell_type}/{r.enhancer}",
                (r.median_on, r.median_dynamic_range),
                xytext=(7, 7),
                textcoords="offset points",
                color=colors[r.cell_type],
                fontsize=9,
                zorder=5,
                arrowprops={"arrowstyle": "-", "color": colors[r.cell_type], "lw": 0.5},
            )
        )
    ax.set_xlabel("Median ON activity (log2 normalized reporter activity)", fontsize=12)
    ax.set_ylabel("Median dynamic range: ON − OFF (log2 units)", fontsize=12)
    ax.set_title("Operating regions of silencer–enhancer pairs", fontsize=16, loc="left", pad=30)
    ax.text(
        0,
        1.02,
        f"Canonical FIMO ≥ {config.M_VALUE}; complete pairs; references explicitly recorded",
        transform=ax.transAxes,
        fontsize=10,
    )
    ax.text(
        0.98,
        0.97,
        "Usable screen: ON ≥ 1 and range ≥ 1",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="#245731",
        fontsize=10,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.15)
    usable = (
        ", ".join(f"{r.cell_type}/{r.enhancer}" for r in regions[regions.usable].itertuples())
        or "none"
    )
    import textwrap

    fig.text(0.09, 0.025, "\n".join(textwrap.wrap("Usable contexts: " + usable, 120)), fontsize=9)
    fig.subplots_adjust(left=0.10, right=0.97, bottom=0.17, top=0.87)
    # Deterministic vertical label repulsion in display coordinates.
    for _ in range(60):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        moved = False
        for i, label in enumerate(annotations):
            # Text-only bounds avoid treating connector lines as label collisions.
            box = label._get_layout(renderer)[0].translated(
                *label.get_transform().transform(label.get_position())
            )
            obstacles = [
                other._get_layout(renderer)[0].translated(
                    *other.get_transform().transform(other.get_position())
                )
                for other in annotations[:i]
            ]
            radius = 6 * fig.dpi / 72
            for other in annotations:
                px, py = ax.transData.transform(other.xy)
                obstacles.append(Bbox.from_bounds(px - radius, py - radius, 2 * radius, 2 * radius))
            if any(box.expanded(1.03, 1.1).overlaps(otherbox) for otherbox in obstacles):
                x, y = label.get_position()
                label.set_position((x, y + 6))
                moved = True
        if not moved:
            break
    paths = []
    for ext in ("png", "pdf", "svg"):
        path = directory / f"gate_operating_region.{ext}"
        metadata = {"Creator": "RE1 gate metrics"}
        if ext == "pdf":
            metadata.update(CreationDate=None, ModDate=None)
        if ext == "svg":
            metadata.update(Date=None)
        fig.savefig(path, dpi=300, metadata=metadata)
        paths.append(path)
    plt.close(fig)
    return paths


def render(
    frame,
    catalog,
    coverage,
    regions,
    strata,
    ranked,
    agreement,
    cutoffs,
    primary,
    fallback,
    cross_report,
):
    """Report actual counts, reference choices, criteria and uncertainty semantics."""
    n = len(catalog)
    n_scr = int(catalog.has_designed_partner.sum())
    observed_any = frame[
        frame.reference_type.eq("scrambled") & frame.dynamic_range.notna()
    ].element_id.nunique()
    usable = regions[regions.usable]
    selected = ", ".join(f"{r.cell_type}/{r.enhancer}" for r in usable.itertuples()) or "none"
    scope = regions[
        [
            "cell_type",
            "enhancer",
            "n_above",
            "n_pairs",
            "n_missing_pair",
            "n_scrambled",
            "n_NegCtrl",
            "median_on",
            "median_off",
            "median_dynamic_range",
            "verdict",
            "scrambled_only_usable",
        ]
    ]
    state_counts = (
        frame.groupby(["reference_type", "reference_reason"], observed=True)
        .size()
        .rename("rows")
        .reset_index()
    )
    errors = frame.groupby("separation_status").size().rename("rows").reset_index()
    dr_total = int(frame.dynamic_range.notna().sum())
    gate = json.loads((VALIDATION / "gate.json").read_text())
    lines = [
        "# Silencer–enhancer operating regions",
        f"**Usable contexts under the declared screen: {selected}.** This is a reporter "
        "operating-region screen, with individual sequence pairs ranked within each context.",
        "## State definitions and source integrity",
        f"The canonical source contributes **{len(frame):,}** element × cell × enhancer rows "
        f"from **{n:,}** elements. Exact source key sets and row counts are asserted. "
        "Enhancers are never aggregated. GEO activity and published lfcSE must match the author "
        "tables exactly, including missingness; author padj is carried "
        "without significance filtering. "
        "The canonical rows are also checked against analysis_dataset.csv. Cached source hashes "
        "and the validated canonical score file are checked before calculation.",
        f"gate_metrics_allowed = {str(gate['gate_metrics_allowed']).lower()}; the explicit "
        "continuous-analysis override remains in force. Published Figure 2B binary silencer labels "
        "are unverified and unused. The original gate finding is preserved.",
        "ON is the observed activity of the exact same element's scrambled-motif construct; OFF "
        "is intact-element activity, in the SAME cell and enhancer. These names describe sequence "
        "states designed to alter REST binding; "
        "no REST occupancy was measured in this calculation. "
        "An element without a designed scrambled partner uses its same-context NegCtrl median "
        "as an explicitly marked baseline proxy. A designed partner with missing activity or a "
        "missing context stays missing; it is never replaced by NegCtrl. No activity is imputed.",
        "## Scrambled coverage and pairing",
        f"**{n_scr:,}/{n:,} ({100 * n_scr / n:.2f}%)** canonical elements "
        "have a designed scrambled "
        f"partner; **{n - n_scr:,}** do not. **{observed_any:,}** "
        "have at least one measured ON/OFF "
        "pair. Coverage is substantial but incomplete; the coverage table distinguishes library "
        "design from measured pairs and the NegCtrl fallback. The source populations `scrambled` "
        "and `scrambled,SNP` are both included; neither is silently omitted.",
        "Pairing uses the exact published `_Ref` → `_Scr` identifier relation, with independent "
        "Table S2 oligo/alias verification, not coordinate-neighbor inference. "
        "All paired sequences "
        "have equal full-element length and identical flanks "
        "outside zero-based inclusive positions "
        "90–111. The scored canonical core remains **(90,110)**. "
        f"**{int(catalog.edited_adjacent_base.fillna(False).sum()):,}** pairs also change adjacent "
        "position 111; these are motif-region perturbations, not strictly core-only substitutions. "
        f"All **{n_scr:,}** scrambled canonical cores were rescored using the shared MA0138.2 "
        f"FIMO procedure: **{int(catalog.scrambled_core_score.ge(config.M_VALUE).sum())}** reach "
        f"20.86; maximum score **{catalog.scrambled_core_score.max():.5f}**. "
        "This supports disruption "
        "of the intended core, not absence of every possible motif elsewhere in the element.",
        markdown(coverage),
        markdown(state_counts),
        "Full sources, the exact pairing catalog and all canonical rows are retained in "
        "`data/processed/gate_metrics_sources.csv`, `data/processed/gate_pair_catalog.csv` and "
        "`data/processed/gate_metrics.csv`. The per-row `reference_type`, `reference_reason`, "
        "`scrambled_id`, source construct IDs and uncertainty basis expose every reference choice.",
        "## Metrics and units",
        "- **dynamic_range = ON − OFF**, in log2 activity units. Positive values mean lower "
        "intact-element activity than its chosen ON reference.\n"
        "- **leakiness = OFF activity**, in log2 normalized reporter-activity units. Lower is "
        "a tighter OFF state; there is no clipping, floor, exponentiation "
        "or baseline normalization.\n"
        "- **separation = dynamic_range / sqrt((SE_ON² + SE_OFF²)/2)**, dimensionless. This is "
        "explicitly **uncertainty-standardized separation**, not a replicate-SD effect size or "
        "Cohen's d. The RMS gives equal weight to the two state variances and stays in log2 units.",
        "For scrambled pairs, SE_ON and SE_OFF are the published DESeq log2FoldChange standard "
        "errors (`lfcSE`), verified against GEO. Raw RNA/DNA count replicates exist, but replicate "
        "activity on the fitted normalized scale and the pair covariance are not supplied. "
        "Treating raw counts, or SE × sqrt(replicate count), as a fitted-activity SD would invent "
        "a spread estimate. Thus separation ranks precision-standardized effects; it is neither "
        "a paired significance test nor a claim about replicate distribution overlap. A "
        "biological replicate-SD separation remains unavailable from these fitted summaries.",
        f"For NegCtrl fallback only, SE_ON is the sample SD (ddof=1) of **{BOOTSTRAPS:,}** "
        f"bootstrap medians of measured controls in that context, with seed **{BOOTSTRAP_SEED}** "
        "combined deterministically with the context key. This estimates uncertainty of the "
        "control-population median, including control heterogeneity, not the uncertainty of an "
        "individual scrambled construct. The bootstrap does not fill missing activities. "
        "Fallback and scrambled precision rankings are therefore also compared separately.",
        "Missing ON or OFF gives missing dynamic range, leakiness and separation. Missing "
        "standard error or zero pooled SE makes separation missing while preserving calculable "
        "dynamic range and leakiness. Negative dynamic range and negative OFF values are retained.",
        f"Observed computable dynamic ranges/leakiness values: **{dr_total:,}**; missing "
        f"pairs: **{len(frame) - dr_total:,}**. Computable separations: "
        f"**{int(frame.separation.notna().sum()):,}**.",
        markdown(errors),
        "## Operating-region screen",
        "The explicit screening rule, specified before inspecting paired outcomes, is "
        f"**median ON ≥{MIN_ON:.2f} log2 activity "
        f"and median dynamic range ≥{MIN_DYNAMIC_RANGE:.2f} "
        "log2 units**. These correspond to a twofold normalized-activity level relative to the "
        "source assay's shared zero and a twofold ON/OFF ratio. They are engineering screening "
        "choices, not experimentally calibrated requirements. Absolute expression or transcript "
        "copy number cannot be inferred from this normalized assay scale.",
        f"All summaries use canonical FIMO **≥{config.M_VALUE}** only. ON and dynamic-range "
        "medians use the identical complete-pair cohort within each context; median(ON−OFF) "
        "is computed per pair, not as a difference of independent medians. n_pairs and missing "
        "pairs are reported. The headline combines the explicitly chosen references; separate "
        "reference strata below show whether this changes the conclusion. A context verdict "
        "describes its median response and is not an assertion that every element passes.",
        markdown(scope),
        "![Labelled operating-region plot](../figures/gate_operating_region.png)",
        "[Vector PDF](../figures/gate_operating_region.pdf) · "
        "[Vector SVG](../figures/gate_operating_region.svg). The shaded upper-right region is "
        "the declared usable screen; black marker outlines identify passing contexts.",
        "### Reference-specific medians (no enhancer aggregation)",
        markdown(strata),
        "### Verdict by context",
    ]
    for r in regions.itertuples():
        reason = {
            "usable_screen": "adequate ON and dynamic range under the declared screen",
            "low_ON": "insufficient ON level",
            "low_dynamic_range": "insufficient repression range",
            "low_ON_and_low_range": "both ON level and repression range are insufficient",
            "missing_pairs": "no measured ON/OFF pairs",
        }[r.verdict]
        lines.append(
            f"- **{r.cell_type}/{r.enhancer}: {r.verdict}** — {reason}; "
            f"median ON {r.median_on:.2f}, range {r.median_dynamic_range:.2f}. "
            f"Scrambled-only screen: {'passes' if r.scrambled_only_usable else 'does not pass'}."
        )
    marginal = regions[
        regions.usable
        & (
            (regions.median_on - MIN_ON).lt(0.05)
            | (regions.median_dynamic_range - MIN_DYNAMIC_RANGE).lt(0.05)
        )
    ]
    for r in marginal.itertuples():
        lines.append(
            f"**Borderline: {r.cell_type}/{r.enhancer}** has median ON {r.median_on:.6f} "
            f"and range {r.median_dynamic_range:.6f}; it lies within 0.05 log2 units of "
            "a screening boundary. Treat its passing label as cutoff-sensitive."
        )
    for r in regions[regions.usable.ne(regions.scrambled_only_usable)].itertuples():
        lines.append(
            f"**Reference-dependent: {r.cell_type}/{r.enhancer}** has all-reference "
            f"range {r.median_dynamic_range:.6f} versus scrambled-only range "
            f"{r.median_dr_scrambled:.6f}. The screen verdict changes with reference "
            "composition. Retain the scrambled-only result as a separate follow-up cohort."
        )
    lines += [
        "Low_ON and low_dynamic_range refer to failure of the declared cutoffs, not proof "
        "of literally zero expression or zero repression. En21 must be judged from the paired, "
        "above-threshold measurements here: the motif-versus-noMotif population gap used in the "
        "dataset report is a different contrast and cannot decide this operating-region question.",
        "### Sensitivity to screening choices",
        markdown(cutoffs),
        "The zero-cutoff row uses strictly positive ON and range. Other rows use inclusive "
        "cutoffs. Contexts are not relabelled by searching this grid for a preferred result.",
        "## Per-element rankings within usable contexts",
        "**Headline rankings use measured scrambled pairs only.** The "
        "[scrambled-only ranking](gate_metrics_scrambled_rankings.csv) and "
        "[separate fallback ranking](gate_metrics_fallback_rankings.csv) recompute all ranks "
        "within context and reference type. The [full CSV](gate_metrics_rankings.csv) retains "
        "both groups and its existing columns; rank columns now refer to their own reference "
        "stratum, not a competition between the two references. All activity and metric values "
        "are unchanged. Ties use element ID for deterministic ordinal ranks.",
        "NegCtrl fallback rows share a context-wide ON value, so their dynamic range is "
        "baseline minus OFF rather than an element-specific scrambled perturbation; the "
        "two reference groups are not directly comparable.",
        "Report tables retain only **leakiness (= OFF activity, log2)** and omit the duplicate "
        "off_activity column. Both names remain in the CSV schema for compatibility. "
        "No composite score is constructed; individual values can be negative or below "
        "screen cutoffs even within a passing context.",
        "### Rank agreement",
        "Top sets are taken independently for each metric within each context. The table reports "
        "dynamic-range versus separation overlap, its fraction of the smaller available top set, "
        "and agreement when restricting both rankings to rows with calculable separation. "
        "Dynamic-range versus lowest-leakiness overlap is also shown. Fewer than 25 eligible "
        "rows uses all available rows; unavailable comparisons are labelled, never scored zero. "
        "An overlap below half is explicitly called low agreement; this descriptive rule does "
        "not authorize a composite score.",
        markdown(agreement[agreement.reference_type.eq("scrambled")])
        if len(agreement)
        else "No contexts pass: rankings and rank agreement are unavailable, not omitted.",
    ]
    if len(agreement):
        all_refs = agreement[agreement.reference_type.eq("scrambled")]
        low = all_refs[all_refs.agreement.eq("low_top_rank_agreement")]
        low_labels = ", ".join(f"{r.cell_type}/{r.enhancer}" for r in low.itertuples()) or "none"
        lines.append(
            f"**Low scrambled-only top-rank agreement contexts: {low_labels}.** Keep the metrics "
            "separate: large effects, low OFF activity and precise effects need not "
            "select the same sequences."
        )
    columns = [
        "rank_dynamic_range",
        "element_id",
        "on_activity",
        "dynamic_range",
        "leakiness",
        "separation",
        "rank_separation",
        "fimo_score",
    ]
    for title, table in (
        ("Primary scrambled-paired top 25", primary),
        ("Separate NegCtrl fallback top 25", fallback),
    ):
        lines.append("## " + title)
        for (cell, enhancer), g in table.groupby(["cell_type", "enhancer"], sort=True):
            lines += [
                f"### {cell}/{enhancer}: top {min(TOP_N, len(g))} by dynamic range",
                markdown(
                    g.head(TOP_N)[columns].rename(
                        columns={"leakiness": "leakiness (= OFF activity, log2)"}
                    )
                ),
            ]
    lines += [
        cross_report,
        "## REST response versus a usable gate",
        "The strongest score–repression correlations occur in En02/En11 contexts, yet the "
        "passing screen is concentrated in En19/En09 (with K562/En11 also passing), because "
        "En02 has almost no ON state to switch: a context can be highly REST-responsive "
        "and still be unusable as a gate under the declared screen.",
    ]
    lines += [
        "## Reproduction and boundaries",
        "Run `make gate-metrics PYTHON=.venv/bin/python` with FIMO on PATH. The command is "
        "offline. It consumes the existing continuous dataset, validated scores, GEO archive, "
        "author results and Table S2. It does not rewrite the analysis dataset, classification "
        "results or gate override. It writes the gate-specific source/pair tables, full rankings, "
        "context summaries, report and deterministic PNG/PDF/SVG files. Input, output and FIMO "
        "hashes are recorded in `data/processed/gate_metrics_manifest.json`. "
        "The [two-run receipt](gate_metrics_reproducibility.json) records actual byte comparisons.",
        "The assay is episomal and these are cell lines. Screen-positive sequence pairs are "
        "experimental candidates, not demonstrated circuits or verified occupancy switches. "
        "NegCtrl fallback measures repression relative to a population baseline, not a "
        "sequence-matched perturbation. Keep it distinguishable in any downstream selection. "
        "No data artifacts are staged or committed.",
    ]
    return "\n\n".join(lines) + "\n"


def main():
    gate_path = VALIDATION / "gate.json"
    gate = json.loads(gate_path.read_text())
    if not (
        gate.get("gate_metrics_allowed")
        and gate.get("scoring") == "PASS"
        and gate.get("sign") == "PASS"
    ):
        raise ValueError("Gate analysis requires passed scoring/sign and authorization")
    if config.M_VALUE != 20.86 or gate.get("threshold") != config.M_VALUE:
        raise ValueError("Threshold changed from the validated value")
    dataset_path = config.PROCESSED_DIR / "analysis_dataset.csv"
    dataset_manifest = config.PROCESSED_DIR / "analysis_dataset_manifest.json"
    dm = json.loads(dataset_manifest.read_text())
    if digest(dataset_path) != dm["outputs"][str(dataset_path.relative_to(config.ROOT))]:
        raise ValueError("Analysis dataset differs from its recorded output hash")
    score_path = VALIDATION / "canonical_scores.csv"
    score_manifest = VALIDATION / "manifest.json"
    sm = json.loads(score_manifest.read_text())
    for path, group in (
        (score_path, "outputs"),
        (GEO_ARCHIVE, "inputs"),
        (TABLE_S2, "inputs"),
        (config.MOTIF_PATH, "inputs"),
    ):
        if digest(path) != sm[group][str(path.relative_to(config.ROOT))]:
            raise ValueError(f"Validated source has changed: {path}")
    source, _ = read_source(populations=POPULATIONS, include_uncertainty=True)
    source = source.sort_values(KEY, kind="stable").reset_index(drop=True)
    with zipfile.ZipFile(TABLE_S2) as archive:
        oligos = pd.read_csv(archive.open(TABLE_S2_MEMBER), sep="\t")
    catalog, cores = pair_catalog(source, oligos)
    with TemporaryDirectory(prefix="gate-scrambled-") as tmp:
        path = Path(tmp) / "cores.fa"
        write_fasta(cores, path)
        scr_scores = score_sequences(path)
    if not scr_scores.has_hit.all() or scr_scores.retains_motif.any():
        raise ValueError(
            "Scrambled canonical core is unscorable or still above threshold; review ON definition"
        )
    catalog = catalog.merge(
        scr_scores[["element_id", "score"]].rename(
            columns={"element_id": "scrambled_id", "score": "scrambled_core_score"}
        ),
        on="scrambled_id",
        how="left",
        validate="m:1",
    )
    frame = assemble(source, catalog, control_baselines(source), pd.read_csv(score_path))
    # Match the production canonical activities exactly, retaining missingness.
    existing = pd.read_csv(dataset_path, float_precision="round_trip", low_memory=False)
    existing = existing[existing.population.eq("motif")].sort_values(KEY).reset_index(drop=True)
    assert_source_rows(existing, frame)
    np.testing.assert_array_equal(existing.raw_activity.to_numpy(), frame.off_activity.to_numpy())
    np.testing.assert_array_equal(
        existing.negctrl_baseline.to_numpy(), frame.negctrl_baseline.to_numpy()
    )
    np.testing.assert_array_equal(existing.fimo_score.to_numpy(), frame.fimo_score.to_numpy())
    coverage = coverage_table(frame)
    regions, strata = operating_regions(frame)
    if len(regions) != 20:
        raise AssertionError("Expected the full 20-context operating table")
    ranked, agreement = rank_pairs(frame, regions)
    primary = cross_context.reference_rankings(ranked, "scrambled")
    fallback = cross_context.reference_rankings(ranked, "NegCtrl")
    ranked = pd.concat([primary, fallback], ignore_index=True).sort_values(
        ["cell_type", "enhancer", "reference_type", "rank_dynamic_range"], kind="stable"
    )
    if set(
        regions.loc[regions.usable, ["cell_type", "enhancer"]].itertuples(index=False, name=None)
    ) != set(cross_context.CONTEXTS):
        raise ValueError("Passing contexts changed; review the specified cross-context scope")
    candidates, candidate_audit, candidate_coverage, denominators, cross_summary = (
        cross_context.analyze(frame)
    )
    cross_report = cross_context.report(candidates, candidate_coverage, denominators, cross_summary)
    cutoffs = sensitivity(regions)
    text = render(
        frame,
        catalog,
        coverage,
        regions,
        strata,
        ranked,
        agreement,
        cutoffs,
        primary,
        fallback,
        cross_report,
    )
    inputs = [
        dataset_path,
        dataset_manifest,
        gate_path,
        score_path,
        score_manifest,
        GEO_ARCHIVE,
        TABLE_S2,
        config.MOTIF_PATH,
        RAW / "sources.json",
        Path(__file__),
        Path(metrics.__file__),
        Path(cross_context.__file__),
        config.ROOT / "src/dataset.py",
        config.ROOT / "src/config.py",
        config.ROOT / "src/motifs.py",
    ]
    inputs += [RAW / f"RESTscreen_derived_{cell}_results.run1.txt.gz" for cell in CELLS]
    output_tables = {
        PAIRS: frame,
        config.PROCESSED_DIR / "gate_metrics_sources.csv": source,
        config.PROCESSED_DIR / "gate_pair_catalog.csv": catalog,
        config.TABLES_DIR / "gate_metrics_coverage.csv": coverage,
        config.TABLES_DIR / "gate_operating_regions.csv": regions,
        config.TABLES_DIR / "gate_metrics_reference_strata.csv": strata,
        config.TABLES_DIR / "gate_metrics_rankings.csv": ranked,
        config.TABLES_DIR / "gate_metrics_scrambled_rankings.csv": primary,
        config.TABLES_DIR / "gate_metrics_fallback_rankings.csv": fallback,
        config.TABLES_DIR / "gate_cross_context_candidates.csv": candidates,
        config.TABLES_DIR / "gate_cross_context_coverage.csv": candidate_audit,
        config.TABLES_DIR / "gate_cross_context_coverage_counts.csv": candidate_coverage,
        config.TABLES_DIR / "gate_cross_context_denominators.csv": denominators,
        config.TABLES_DIR / "gate_metrics_rank_agreement.csv": agreement,
        config.TABLES_DIR / "gate_metrics_sensitivity.csv": cutoffs,
    }
    with TemporaryDirectory(prefix="gate-publish-", dir=config.PROCESSED_DIR) as tmp:
        staging = Path(tmp)
        destinations = {}
        for path, table in output_tables.items():
            table.to_csv(staging / path.name, index=False, float_format="%.17g")
            destinations[path] = staging / path.name
        (staging / REPORT.name).write_text(text)
        destinations[REPORT] = staging / REPORT.name
        for path in operating_figure(regions, staging) + cross_context.figure(candidates, staging):
            destinations[config.FIGURES_DIR / path.name] = path
        executable = Path(shutil.which("fimo"))
        manifest = dict(
            inputs={str(p.relative_to(config.ROOT)): digest(p) for p in inputs},
            outputs={str(p.relative_to(config.ROOT)): digest(q) for p, q in destinations.items()},
            tools={str(executable): digest(executable)},
            cross_context=cross_summary,
            source_rows=len(source),
            canonical_rows=len(frame),
            operating_rule={"min_on": MIN_ON, "min_dynamic_range": MIN_DYNAMIC_RANGE},
            spread="RMS published lfcSE; fallback ON bootstrap median SE, not replicate SD",
            bootstrap={"resamples": BOOTSTRAPS, "seed": BOOTSTRAP_SEED},
        )
        (staging / MANIFEST.name).write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        destinations[MANIFEST] = staging / MANIFEST.name
        for path, temp in destinations.items():
            path.parent.mkdir(exist_ok=True, parents=True)
            os.replace(temp, path)
    print(
        f"Canonical source rows preserved: {len(frame):,}; "
        f"scrambled partners: {catalog.has_designed_partner.sum():,}/{len(catalog):,}"
    )
    print(
        regions[
            ["cell_type", "enhancer", "n_pairs", "median_on", "median_dynamic_range", "verdict"]
        ].to_string(index=False)
    )
    print(
        f"Cross-context qualified: {len(candidates)}; "
        f"excluded for coverage: {cross_summary['n_excluded_coverage']}; "
        f"Spearman with FIMO: {cross_summary['spearman_min_percentile_fimo']}"
    )
    print(f"Report: {REPORT.relative_to(config.ROOT)}")


if __name__ == "__main__":
    main()
