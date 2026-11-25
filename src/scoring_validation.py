"""Offline scoring/sign audit; honor an explicit continuous-analysis override.

Run from the repository root with ``python -m src.scoring_validation``.
The available padj values test activity against the shared background, not the
matched-enhancer silencer null. Diagnostic calls must never become gate labels.
"""

import gzip
import hashlib
import json
import shutil
import tarfile

import matplotlib
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from src import config
from src._utils import write_fasta
from src.assembly_check import GEO_ARCHIVE, TABLE_S2
from src.lift import published_elements
from src.motifs import score_sequences

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAW = config.RAW_DIR / "scoring_validation"
OUT = config.PROCESSED_DIR / "scoring_validation"
REPORT = config.TABLES_DIR / "scoring_validation.md"
CELLS = {"GM12878": "GM12878", "HepG2": "HepG2", "K562": "K562", "SKNSH": "SK-N-SH"}
FDR = 0.01
EXPECTED_REPRESSORS = {"motif": 2657, "noMotif": 486}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def aggregate(frame):
    """One element per cell; median activity and any observed diagnostic call.

    Missing tests stay missing. A false value means no detected call in the
    available contexts, not proof of no repression in unmeasured contexts.
    """
    frame = frame.copy()
    if not np.isfinite(frame.log2FoldChange.dropna()).all():
        raise ValueError("Nonfinite activity")
    if not frame.padj.dropna().between(0, 1).all():
        raise ValueError("Adjusted p-values must lie in [0,1]")
    valid = frame.log2FoldChange.notna() & frame.padj.notna()
    frame["diagnostic_call"] = (frame.log2FoldChange.lt(0) & frame.padj.le(FDR)).astype("boolean")
    frame.loc[~valid, "diagnostic_call"] = pd.NA
    result = frame.groupby(["element_id", "project"], as_index=False).agg(
        activity=("log2FoldChange", "median"),
        n_measured_contexts=("log2FoldChange", "count"),
        n_testable_contexts=("diagnostic_call", "count"),
        diagnostic_call=("diagnostic_call", "any"),
    )
    result.loc[result.n_testable_contexts.eq(0), "diagnostic_call"] = pd.NA
    return result


def sign_check(summary):
    """Refuse an inverted or unavailable canonical-versus-contrast ordering."""
    for cell in CELLS.values():
        sub = summary[summary.cell_type.eq(cell)].set_index("population")
        if not sub.loc["motif", "median"] < sub.loc["noMotif", "median"]:
            raise ValueError(f"SIGN GATE FAILED: canonical median is not below contrast in {cell}")


def stats(frame):
    measured = frame[frame.activity.notna()]
    calls = measured.diagnostic_call.dropna()
    n, k = len(calls), int(calls.sum())
    return {
        "n_elements": len(frame),
        "n_measured": len(measured),
        "n_testable": n,
        "n_diagnostic_calls": k,
        "median_activity": measured.activity.median(),
        "diagnostic_fraction": k / n if n else np.nan,
        "published_repressing_fraction": np.nan,
        "n_incomplete_contexts": int(measured.n_testable_contexts.lt(5).sum()),
    }


def bins(frame, edges):
    """Include empty bins and every measured element; left-closed, right-open."""
    rows = []
    for left, right in zip(edges[:-1], edges[1:], strict=True):
        sub = frame[frame.score.ge(left) & frame.score.lt(right)]
        rows.append({"bin_left": left, "bin_right": right, **stats(sub)})
    if sum(row["n_elements"] for row in rows) != len(frame):
        raise ValueError("Score bins lost or duplicated elements")
    return pd.DataFrame(rows)


def md(frame):
    """Markdown without an optional tabulate dependency; missing is explicit."""

    def fmt(x):
        if pd.isna(x):
            return "NA"
        if isinstance(x, (float, np.floating)):
            return f"{x:.6f}"
        return str(x)

    return "\n".join(
        [
            "| " + " | ".join(frame.columns) + " |",
            "| " + " | ".join(["---"] * len(frame.columns)) + " |",
            *(
                "| " + " | ".join(map(fmt, row)) + " |"
                for row in frame.itertuples(index=False, name=None)
            ),
        ]
    )


def plot(table, splits):
    """Show available calls honestly; never title them verified silencer calls."""
    plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    for ax, cell in zip(axes.flat, CELLS.values(), strict=True):
        sub = table[table.cell_type.eq(cell)].copy()
        x = (sub.bin_left + sub.bin_right) / 2
        n, k = sub.n_testable.to_numpy(), sub.n_diagnostic_calls.to_numpy()
        valid = n > 0
        p = np.divide(k, n, out=np.full(len(n), np.nan), where=valid)
        # Descriptive Wilson intervals; elements may share sequence/context.
        z = 1.959963984540054
        nv = n[valid]
        center = (p[valid] + z * z / (2 * nv)) / (1 + z * z / nv)
        half = (
            z * np.sqrt(p[valid] * (1 - p[valid]) / nv + z * z / (4 * nv * nv)) / (1 + z * z / nv)
        )
        ax.fill_between(x[valid], center - half, center + half, color="#207DA2", alpha=0.16)
        ax.plot(x, p, color="#207DA2", lw=1.5)
        ax.scatter(x[valid], p[valid], s=15 + np.sqrt(nv), color="#207DA2", zorder=3)
        ax.axvline(config.M_VALUE, ls="--", color="#C74835", lw=1.5)
        ax.set_title(cell, loc="left", fontweight="bold")
        ax.set_ylim(-0.025, 1.025)
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        split = splits[splits.cell_type.eq(cell)].set_index("group")
        ax.text(
            0.03,
            0.96,
            f"Below: {split.loc['below', 'diagnostic_fraction']:.1%}\n"
            f"At/above: {split.loc['at_or_above', 'diagnostic_fraction']:.1%}",
            transform=ax.transAxes,
            va="top",
            fontsize=10,
        )
    fig.suptitle("REST score and available activity/FDR calls", fontsize=16, fontweight="bold")
    fig.supxlabel("Canonical-core FIMO score (log-odds score units; bin width 1.0)", y=0.12)
    fig.supylabel("Fraction with negative activity and published padj ≤ 0.01", x=0.025)
    fig.text(
        0.10,
        0.066,
        "Diagnostic only: published silencer classification is NOT reproduced.",
        fontsize=11,
        color="#A32C23",
        fontweight="bold",
    )
    fig.text(
        0.10,
        0.037,
        "Dashed line: fixed 20.86. Any measured enhancer context per element. "
        "Shading: 95% Wilson intervals.",
        fontsize=9,
    )
    fig.text(
        0.10,
        0.016,
        "Point size reflects tested elements per bin. Missing bins are not zero. "
        "Binary silencer labels remain unverified.",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.21, top=0.90, wspace=0.16, hspace=0.22)
    paths = []
    for ext in ("png", "pdf", "svg"):
        path = config.FIGURES_DIR / f"scoring_validation.{ext}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path)
    plt.close(fig)
    return paths


def main():
    OUT.mkdir(exist_ok=True)
    gate_path = OUT / "gate.json"
    previous_gate = json.loads(gate_path.read_text()) if gate_path.exists() else {}
    override = previous_gate.get("override", {})
    continuous_allowed = override.get("approved_for_continuous_analysis") is True
    gate_path.write_text(
        json.dumps({"status": "INCOMPLETE", "gate_metrics_allowed": False, "override": override})
        + "\n"
    )
    inputs = [GEO_ARCHIVE, TABLE_S2, config.MOTIF_PATH, config.ROOT / __file__]
    sources = json.loads((RAW / "sources.json").read_text())
    for name, meta in sources.items():
        path = RAW / name
        if sha(path) != meta["sha256"]:
            raise ValueError(f"Source hash mismatch: {path}")
        inputs.append(path)
    inputs.extend(RAW.glob("*.R"))
    inputs.extend(RAW.glob("*.ipynb"))
    inputs.extend([RAW / "author_README.md", RAW / "sources.json"])
    elements = published_elements()
    canonical = elements[elements.project.eq("motif")].copy()
    if len(canonical) != 8436 or config.MOTIF_SPAN != (90, 110) or config.M_VALUE != 20.86:
        raise ValueError("Published population, coordinate or fixed threshold contract changed")
    fasta = OUT / "canonical_cores.fa"
    start, stop = config.MOTIF_SPAN
    write_fasta(
        dict(
            zip(
                canonical.element_id,
                canonical.sequence.str.upper().str.slice(start, stop + 1),
                strict=True,
            )
        ),
        fasta,
    )
    scores = score_sequences(fasta)
    scores.to_csv(OUT / "canonical_scores.csv", index=False)
    reference = pd.read_csv(RAW / "ref_all_fimo.txt", sep="\t")
    reference = reference[
        reference["#pattern name"].eq(config.JASPAR_ID)
        & reference.start.eq(91)
        & reference.stop.eq(111)
    ]
    reference = reference.sort_values("score", ascending=False).drop_duplicates("sequence name")
    comparison = scores.merge(
        reference[["sequence name", "score"]],
        left_on="element_id",
        right_on="sequence name",
        how="left",
        suffixes=("_ours", "_published"),
        validate="1:1",
    )
    comparison["difference"] = comparison.score_ours - comparison.score_published
    comparison.to_csv(OUT / "published_score_comparison.csv", index=False)
    compared = comparison[comparison.score_published.notna()]
    if not compared.difference.eq(0).all() or scores.score.isna().any():
        raise ValueError("SCORING GATE FAILED: score mismatch or unscorable canonical core")

    sign, context_sign, audit, all_elements, contexts, agreement = [], [], [], [], [], []
    unions = {project: set() for project in EXPECTED_REPRESSORS}
    with tarfile.open(GEO_ARCHIVE) as archive:
        for source_cell, cell in CELLS.items():
            names = [n for n in archive.getnames() if f"WGREST_{source_cell}.txt.gz" in n]
            if len(names) != 1:
                raise ValueError(f"Expected one GEO member for {cell}")
            with gzip.GzipFile(fileobj=archive.extractfile(names[0])) as handle:
                geo = pd.read_csv(handle, sep="\t", low_memory=False)
            author = pd.read_csv(
                RAW / f"RESTscreen_derived_{source_cell}_results.run1.txt.gz",
                sep="\t",
                low_memory=False,
            )
            joined = geo.merge(
                author[["ID", "log2FoldChange", "padj", "pvalue", "stat"]],
                left_on="id",
                right_on="ID",
                how="outer",
                validate="1:1",
                suffixes=("", "_author"),
                indicator=True,
            )
            same = joined.log2FoldChange.eq(joined.log2FoldChange_author) | (
                joined.log2FoldChange.isna() & joined.log2FoldChange_author.isna()
            )
            if not joined._merge.eq("both").all() or not same.all():
                raise ValueError(f"GEO/author activity identity failed: {cell}")
            if not joined.id.eq(joined.enhancer + "^" + joined.silencer).all():
                raise ValueError("Construct identifiers do not match explicit metadata")
            if joined.duplicated(["silencer", "enhancer"]).any():
                raise ValueError("Duplicate element/enhancer observations")
            joined = joined.rename(columns={"silencer": "element_id"})
            agreement.append(
                {
                    "cell_type": cell,
                    "expected_rows": 139325,
                    "observed_rows": len(joined),
                    "identical_activities_including_NA": int(same.sum()),
                }
            )
            agg = aggregate(joined)
            agg["cell_type"] = cell
            for population in ("all", "motif", "noMotif", "NegCtrl"):
                group = agg if population == "all" else agg[agg.project.eq(population)]
                values = group.activity.dropna()
                sign.append(
                    {
                        "cell_type": cell,
                        "population": population,
                        "n_elements": len(group),
                        "n_measured": len(values),
                        "n_missing": int(group.activity.isna().sum()),
                        "median": values.median(),
                        "p05": values.quantile(0.05),
                        "p95": values.quantile(0.95),
                    }
                )
            for enhancer, group in joined.groupby("enhancer"):
                row = {"cell_type": cell, "enhancer": enhancer}
                for population in ("motif", "noMotif", "NegCtrl"):
                    row[population + "_median"] = group.loc[
                        group.project.eq(population), "log2FoldChange"
                    ].median()
                context_sign.append(row)
            for project in unions:
                selected = joined[
                    joined.project.eq(project) & joined.log2FoldChange.lt(0) & joined.padj.le(FDR)
                ]
                unions[project].update(selected.element_id)
            native = agg[agg.project.eq("motif")].merge(
                scores[["element_id", "score"]],
                on="element_id",
                how="outer",
                validate="1:1",
                indicator=True,
            )
            native["join_reason"] = np.select(
                [
                    native._merge.eq("left_only"),
                    native._merge.eq("right_only"),
                    native.score.isna(),
                    native.activity.isna(),
                    native.n_testable_contexts.eq(0),
                ],
                [
                    "no_score_record",
                    "no_activity_record",
                    "unscorable",
                    "all_activity_missing",
                    "activity_present_no_FDR",
                ],
                default="joined",
            )
            audit.append(
                {
                    "cell_type": cell,
                    "expected_canonical": 8436,
                    "observed_canonical": len(native),
                    "scores_present": int(native.score.notna().sum()),
                    "measured_with_score": int(
                        (native.score.notna() & native.activity.notna()).sum()
                    ),
                    "all_activity_missing": int(native.activity.isna().sum()),
                    "unmatched_ids": int(native._merge.ne("both").sum()),
                    "measured_without_FDR": int(
                        (native.activity.notna() & native.n_testable_contexts.eq(0)).sum()
                    ),
                    "missing_activity_contexts": int(
                        joined.loc[joined.project.eq("motif"), "log2FoldChange"].isna().sum()
                    ),
                    "measured_contexts_without_FDR": int(
                        (
                            joined.project.eq("motif")
                            & joined.log2FoldChange.notna()
                            & joined.padj.isna()
                        ).sum()
                    ),
                }
            )
            all_elements.append(native.drop(columns="_merge"))
            contexts.append(
                joined[joined.project.eq("motif")][
                    ["element_id", "enhancer", "log2FoldChange", "padj", "pvalue", "stat"]
                ].assign(cell_type=cell)
            )
    sign = pd.DataFrame(sign)
    sign.to_csv(OUT / "sign_summary.csv", index=False)
    sign_check(sign)
    all_elements = pd.concat(all_elements, ignore_index=True)
    all_elements.to_csv(OUT / "element_audit.csv", index=False)
    pd.concat(contexts).to_csv(OUT / "construct_audit.csv", index=False)
    pd.DataFrame(context_sign).to_csv(OUT / "sign_by_enhancer.csv", index=False)
    measured = all_elements[all_elements.activity.notna() & all_elements.score.notna()]
    edges = scores.score.min() + np.arange(
        int(np.floor(scores.score.max() - scores.score.min())) + 2
    )
    binned, splits, nearby = [], [], []
    for cell in CELLS.values():
        group = measured[measured.cell_type.eq(cell)]
        binned.append(bins(group, edges).assign(cell_type=cell))
        for label, mask in [
            ("below", group.score.lt(config.M_VALUE)),
            ("at_or_above", group.score.ge(config.M_VALUE)),
        ]:
            splits.append({"cell_type": cell, "group": label, **stats(group[mask])})
        low = group[group.score.ge(config.M_VALUE - 2) & group.score.lt(config.M_VALUE)]
        high = group[group.score.ge(config.M_VALUE) & group.score.lt(config.M_VALUE + 2)]
        a, b = stats(low), stats(high)
        table = [
            [b["n_diagnostic_calls"], b["n_testable"] - b["n_diagnostic_calls"]],
            [a["n_diagnostic_calls"], a["n_testable"] - a["n_diagnostic_calls"]],
        ]
        nearby.append(
            {
                "cell_type": cell,
                "n_below_2_units": a["n_testable"],
                "fraction_below": a["diagnostic_fraction"],
                "n_above_2_units": b["n_testable"],
                "fraction_above": b["diagnostic_fraction"],
                "difference_pp": 100 * (b["diagnostic_fraction"] - a["diagnostic_fraction"]),
                "fisher_two_sided_p": fisher_exact(table).pvalue,
            }
        )
    binned, splits, nearby = pd.concat(binned), pd.DataFrame(splits), pd.DataFrame(nearby)
    for name, table in [
        ("bins", binned),
        ("threshold_split", splits),
        ("threshold_neighborhood", nearby),
        ("join_summary", pd.DataFrame(audit)),
        ("geo_author_agreement", pd.DataFrame(agreement)),
    ]:
        table.to_csv(OUT / f"{name}.csv", index=False)
    figures = plot(binned, splits)
    observed = pd.DataFrame(
        [
            {
                "population": p,
                "expected_published_repressors": expected,
                "observed_negative_padj_calls": len(unions[p]),
            }
            for p, expected in EXPECTED_REPRESSORS.items()
        ]
    )
    n_chr8 = int(compared.element_id.str.startswith("REST_chr8:").sum())
    pilot_path = config.PROCESSED_DIR / "pilot_mouse_evidence/calls.csv"
    inputs.append(pilot_path)
    pilot_ids = pd.read_csv(pilot_path, usecols=["element_id"]).element_id
    n_pilot = int(compared.element_id.isin(pilot_ids).sum())
    lines = [
        "# Motif scoring and activity-sign validation",
        "**Correctness gate: FAIL / STOP. No gate metrics were computed or authorized by this"
        " audit.**",
        "**Verdict: the threshold does not establish the requested published-FDR repressor "
        "separation "
        "in our data. The available activity-test labels do not reproduce the published "
        "silencer set.**",
        "Scoring and the direction of activity pass. The blocking issue is the "
        "repression-label/null "
        "definition, not evidence for a score-scale mismatch. Diagnostic enrichment below is not a "
        "validated silencer classification.",
        "## Expected versus observed",
        f"- Canonical cores: expected 8,436; observed {len(scores):,}; "
        f"missing scores {scores.score.isna().sum()}.",
        "- Published score agreement: expected exact agreement; "
        f"observed {len(compared):,}/{len(compared):,} "
        f"exact matches (maximum absolute difference {compared.difference.abs().max():g}), "
        "including "
        f"{n_chr8} native chr8 sites ({n_pilot} in the current hg38 chr8 pilot). "
        f"{len(scores) - len(compared):,} recomputed sites "
        "have no published canonical-position "
        "hit for a direct score comparison; their scores are measured by uncensored FIMO, not"
        " filled with zero.",
        f"- Fixed threshold: expected 20.86; observed configuration {config.M_VALUE}; "
        "not fitted or adjusted.",
        "- Sign: expected canonical median < no-motif median in each cell; observed PASS in "
        "all four.",
        "- Repression labels: expected published Figure 2B totals; observed mismatch:",
        md(observed),
        "## Source and sign evidence",
        "All four GEO WGREST tables were loaded, including controls, scrambled constructs and"
        " variants. "
        "The authors' processed result tables supply padj, which GEO/Table S6 omit. Their "
        "activities "
        "match GEO exactly, including missingness:",
        md(pd.DataFrame(agreement)),
        "The author pipeline explicitly contrasts cellular RNA with DNA, applies a "
        "shared-background "
        "normalization based on En02 plus negative-control sequences, and reruns the DESeq "
        "Wald test. "
        "Lower log2FoldChange means lower reporter activity. It does NOT mean that every "
        "repressing "
        "element must have log2FoldChange < 0: active enhancer baselines remain above zero. "
        "Repression relative to a matched enhancer requires its own reference comparison. "
        "No signs, zeros or normalization offsets were changed in this audit.",
        "Element-level activity below is the median of available enhancer-context "
        "log2FoldChange values "
        "within each cell, followed by the requested population statistics. Enhancers are distinct "
        "contexts, not biological replicates. All means/medians and percentiles use observed "
        "values only; "
        "all includes every WGREST population (not the separate benchmark/noncanonical libraries). "
        "Units: log2 normalized RNA/DNA activity. Percentiles use linear interpolation.",
        md(sign),
        "### Enhancer-context sign check (no pooling across cell types)",
        md(pd.DataFrame(context_sign)),
        "## Why the published repression FDR is not verified",
        "The article reports 1% FDR repression in at least one cell/enhancer combination, but the "
        "retrieved results' padj belongs to a shared-background activity test. The inspected "
        "public "
        "analysis notebooks do not supply the Figure 2B matched-enhancer silencer-call "
        "procedure or "
        "its per-element labels. Neither Table S6 nor GEO includes those labels. Using "
        "negative activity "
        "and supplied padj ≤0.01 is an auditable diagnostic, but does not recover the "
        "reported totals "
        "and is not substituted for the missing repression labels. No new FDR calculation or "
        "threshold tuning was used to force agreement.",
        "The published 20.86 is an average breakpoint in native-versus-scrambled motif "
        "contribution "
        "across 18 responsive cell/enhancer combinations (reported range 19.0–21.9). It is not a "
        "guaranteed discontinuity in the binary repression fraction after enhancer aggregation. "
        "No step in that different endpoint would not, by itself, contradict exact per-site "
        "score agreement.",
        "## Joins and missingness",
        "Scores use all published native oligos, including elements excluded by hg38 conversion; "
        "this avoids conditioning an assay reproduction on genome lift success. Joins use "
        "exact native "
        "element IDs and the existing explicit Table S2 alias resolver, never coordinates "
        "guessed from IDs.",
        md(pd.DataFrame(audit)),
        "The 493 elements without activity lack fitted values in every source enhancer context. "
        "The source tables do not supply a per-row rejection reason; the paper describes "
        "barcode/DNA-count filters, but an individual cause cannot be assigned "
        "from missingness alone.",
        "Every canonical element/cell pair and exclusion reason is retained in element_audit.csv. "
        "Missing activity in all five contexts is excluded only from measured-activity summaries; "
        "it is not a non-repressor. Missing adjusted p-values are excluded from call denominators. "
        "An element with at least one available FDR test contributes once; any observed negative "
        "significant context gives a diagnostic call. Incomplete contexts are counted, and a "
        "negative "
        "diagnostic result is only 'not detected in available contexts'. This any-context summary "
        "does not establish an element-level 1% FDR after combining tests.",
        "## Fixed above/below split: diagnostic calls only",
        "The requested **published repressing fraction is NA (not verified)** throughout. "
        "diagnostic_fraction is the explicitly different available quantity: fraction of tested "
        "elements with log2FoldChange <0 and supplied padj ≤0.01 in at least one enhancer context. "
        "n_measured and n_testable state both denominators; no missing value is scored as zero.",
        md(splits),
        "### Shape near 20.86",
        "The fixed local comparison is [18.86,20.86) versus [20.86,22.86); these windows were not "
        "searched for a favorable result. The Fisher tests are descriptive, two-sided and "
        "unadjusted; "
        "they do not identify a breakpoint or repair the missing silencer labels.",
        md(nearby),
        "**Observed shape:** GM12878 stays near a low background fraction across the score range; "
        "HepG2 and SK-N-SH rise broadly with increasing score, and K562 rises more modestly "
        "at high "
        "scores. There is no clear abrupt step exactly at 20.86 in this diagnostic endpoint. "
        "Across the fixed neighboring windows, all four unadjusted Fisher p-values exceed 0.05; "
        "this does not prove absence of a biological transition. The broad above/below enrichment "
        "in three cells is not a substitute for validating the published repression calls.",
        "## Bins (width 1.0, anchored to the observed minimum)",
        f"Score range: {scores.score.min():.5f} to {scores.score.max():.5f}. "
        "Intervals are [left,right); the final interval extends beyond the maximum. "
        "Empty bins remain present with n=0 and fraction/median NA. n_elements counts measured "
        "canonical elements with a score, without counting enhancer contexts as separate elements.",
    ]
    for cell in CELLS.values():
        lines += [f"### {cell}", md(binned[binned.cell_type.eq(cell)].drop(columns="cell_type"))]
    lines += [
        "## Figure and reproduction",
        "![Activity/FDR diagnostic by score](../figures/scoring_validation.png)",
        "PNG (300 dpi), PDF and SVG are saved in results/figures/. The figure deliberately labels "
        "the available endpoint and the failed gate rather than claiming published repressor "
        "fractions.",
        "Run `PATH=/opt/miniconda3/envs/re1-gates/bin:$PATH MPLCONFIGDIR=/tmp/re1-matplotlib "
        "python -m src.scoring_validation` from the repository root with the "
        "cached inputs. "
        "Exit code 2 is intentional: the diagnostic artifacts are published, but "
        "gate_metrics_allowed "
        "remains false in data/processed/scoring_validation/gate.json. No gate metrics, models or "
        "candidate ranking are run. Input/output SHA-256 hashes are in manifest.json.",
        "**Required resolution:** obtain or reconstruct and independently validate the exact "
        "Figure 2B repression null, multiplicity scope and per-element calls, including enhancer "
        "handling. Then rerun this threshold analysis with those labels. Do not reverse signs or "
        "change 20.86 to resolve a label mismatch.",
        "## Sources",
        "- [Mouri et al., Figure 2 and STAR "
        "Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC9903721/).",
        "- [Published Table "
        "S6](https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc7.zip).",
        "- [Author processed "
        "data](https://github.com/tewhey-lab/duoREST/tree/main/processed_data).",
        "- [Author motif-score "
        "selection](https://github.com/tewhey-lab/duoREST/blob/main/README.md).",
        "- [DUOmodel.R, normalization and "
        "contrasts](https://github.com/tewhey-lab/MPRAduo/blob/master/DUOmodel.R).",
        "- [Author figure-generation "
        "notebook](https://github.com/tewhey-lab/duoREST/blob/main/figure_generation/duoREST_fig_generation_R.ipynb).",
    ]
    if continuous_allowed:
        lines.insert(
            1,
            "**User override for continuous analysis ("
            + override["date"]
            + "):** "
            + override["reason"]
            + "\n\nThe binary-label audit below is retained. "
            "Its blocking verdict and exit-code description apply to the original "
            "binary-label requirement; the recorded override permits continuous analysis "
            "and a successful audit rerun exits 0. Scoring/sign failures still stop execution.",
        )
    REPORT.write_text("\n\n".join(lines) + "\n")
    gate = {
        "status": "FAIL",
        "gate_metrics_allowed": False,
        "scoring": "PASS",
        "sign": "PASS",
        "published_repression_labels": "UNVERIFIED",
        "observed_union": {p: len(v) for p, v in unions.items()},
        "expected_union": EXPECTED_REPRESSORS,
        "threshold": config.M_VALUE,
    }
    if continuous_allowed:
        gate.update(
            original_status="FAIL",
            status="OVERRIDDEN_FOR_CONTINUOUS_ANALYSIS",
            gate_metrics_allowed=True,
            override=override,
        )
    gate_path.write_text(json.dumps(gate, indent=2) + "\n")
    outputs = [p for p in OUT.iterdir() if p.name != "manifest.json"] + [REPORT, *figures]
    executable = shutil.which("fimo")
    manifest = {
        "inputs": {str(p.relative_to(config.ROOT)): sha(p) for p in inputs},
        "outputs": {str(p.relative_to(config.ROOT)): sha(p) for p in outputs},
        "fimo": {"path": executable, "sha256": sha(config.ROOT / executable)},
        "gate": gate,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(sign.to_string(index=False))
    print(splits.to_string(index=False))
    print(nearby.to_string(index=False))
    if continuous_allowed:
        print("Continuous analysis authorized by user override; binary labels remain unverified.")
        return 0
    print("FAIL / STOP: published repression labels are unverified; gate metrics remain blocked.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
