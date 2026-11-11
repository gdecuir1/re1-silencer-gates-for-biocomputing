"""Flank-anchored synteny rescue, independent of the target assembly.

Source flanks are adjacent to, and exclude, the full element. BED strands are
deliberately '+' so target strand describes the mapping of source genomic order,
not library-oligo orientation. No existing mapping or FIMO output is rewritten.
"""

import argparse
import json
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from src import config, motifs
from src._utils import read_fasta, require_columns, run_command, validate_intervals, write_bed
from src.classify import OrthologState
from src.downloads import local_md5
from src.lift import lift_intervals, run_tool
from src.pilot import parse_sizes

FLANK_SWEEP = (500, 1000, 2000)
MAX_GAP_BP = 10_000
WINDOW_BP = 101
# An operational span guard, not a biological deletion cutoff.
STATE_MAP = {
    "RESCUED": OrthologState.RETAINED,
    "ABSENT_CONFIRMED": OrthologState.NO_MOTIF,
    "NO_SYNTENY": OrthologState.FAILED_LIFT,
}
REASONS = (
    "source_flank_out_of_bounds",
    "both_flanks_unmapped",
    "left_flank_unmapped",
    "right_flank_unmapped",
    "different_target_chromosomes",
    "different_target_orientations",
    "order_not_preserved",
    "negative_gap",
    "gap_too_large",
    "target_out_of_bounds",
)


def absent_elements(sidecar: pd.DataFrame, source: pd.DataFrame, assembly: str) -> pd.DataFrame:
    """Select exactly the old ABSENT IDs from verified source-coordinate columns."""
    require_columns(sidecar, ["element_id", "call", "source_build", "target_build"])
    require_columns(source, ["element_id", "working_build"])
    if not sidecar.element_id.is_unique or not source.element_id.is_unique:
        raise ValueError("Duplicate element IDs in source or sidecar")
    if not sidecar.target_build.eq(assembly).all():
        raise ValueError("Sidecar target_build differs from the requested assembly")
    if not sidecar.source_build.eq(config.GENOME_BUILD).all():
        raise ValueError("Sidecar source_build differs from the working assembly")
    if not source.working_build.eq(config.GENOME_BUILD).all():
        raise ValueError("Use verified working-build source elements")
    ids = set(sidecar.loc[sidecar.call.eq("ABSENT"), "element_id"])
    result = source.loc[source.element_id.isin(ids)].copy()
    if not ids or set(result.element_id) != ids:
        raise ValueError("Every ABSENT ID must exist exactly once in the source elements")
    validate_intervals(result)
    # Detect stale sidecars, without ever reconstructing coordinates from an ID.
    old = sidecar.set_index("element_id")
    for row in result.itertuples(index=False):
        for key in ("chrom", "start", "end", "strand"):
            column = f"{config.GENOME_BUILD}_element_{key}"
            if column not in old or old.loc[row.element_id, column] != getattr(row, key):
                raise ValueError(f"Source/sidecar coordinate mismatch: {row.element_id}: {key}")
    return result.sort_values("element_id", kind="stable").reset_index(drop=True)


def build_flanks(
    elements: pd.DataFrame, sizes: dict[str, int], flank: int
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Build exact left/right flanks; reject boundary cases rather than shorten anchors."""
    validate_intervals(elements)
    if isinstance(flank, bool) or not isinstance(flank, int) or flank <= 0:
        raise ValueError("Flank must be a positive integer")
    left, right, excluded = [], [], set()
    for row in elements.sort_values("element_id").itertuples(index=False):
        if row.chrom not in sizes or row.start < flank or row.end + flank > sizes[row.chrom]:
            excluded.add(row.element_id)
            continue
        base = dict(element_id=row.element_id, chrom=row.chrom, strand="+")
        left.append(base | dict(start=row.start - flank, end=row.start))
        right.append(base | dict(start=row.end, end=row.end + flank))
    columns = ["element_id", "chrom", "start", "end", "strand"]
    return pd.DataFrame(left, columns=columns), pd.DataFrame(right, columns=columns), excluded


def anchor_gap(
    left: dict | None, right: dict | None, sizes: dict[str, int], max_gap: int = MAX_GAP_BP
) -> dict:
    """Check both anchors, chromosome, orientation, order, gap, then reference bounds.

    One mutually exclusive first-failure reason is reported per pair. Reverse
    mappings must reverse genomic order: the human-left anchor is target-right.
    A zero gap is a candidate complete deletion; it is not an extraction error.
    """
    result = dict(
        valid_pair=False,
        rejection_reason="",
        gap_chrom="",
        gap_start="",
        gap_end="",
        gap_strand="",
        gap_width="",
    )
    if max_gap < 0:
        raise ValueError("Maximum gap must be nonnegative")
    if left is None or right is None:
        reason = (
            "both_flanks_unmapped"
            if left is None and right is None
            else "left_flank_unmapped"
            if left is None
            else "right_flank_unmapped"
        )
    elif left["chrom"] != right["chrom"]:
        reason = "different_target_chromosomes"
    elif left["strand"] != right["strand"] or left["strand"] not in {"+", "-"}:
        reason = "different_target_orientations"
    else:
        low, high = (left, right) if left["strand"] == "+" else (right, left)
        if low["start"] >= high["start"] or low["end"] > high["end"]:
            reason = "order_not_preserved"
        else:
            start, end = low["end"], high["start"]
            result.update(
                gap_chrom=left["chrom"],
                gap_start=start,
                gap_end=end,
                gap_strand=left["strand"],
                gap_width=end - start,
            )
            if end < start:
                reason = "negative_gap"
            elif end - start > max_gap:
                reason = "gap_too_large"
            elif any(
                r["chrom"] not in sizes or not 0 <= r["start"] < r["end"] <= sizes[r["chrom"]]
                for r in (left, right)
            ):
                reason = "target_out_of_bounds"
            else:
                reason = ""
    result.update(valid_pair=not reason, rejection_reason=reason)
    return result


def run_anchors(
    elements: pd.DataFrame,
    source_sizes: dict[str, int],
    target_sizes: dict[str, int],
    assembly: str,
    chain: Path,
    directory: Path,
    flank: int,
    max_gap: int = MAX_GAP_BP,
) -> pd.DataFrame:
    """Lift each flank separately at the shared threshold and retain all failure records."""
    left, right, excluded = build_flanks(elements, source_sizes, flank)
    tracks = {}
    for label, frame in (("left", left), ("right", right)):
        track_dir = directory / label
        track_dir.mkdir(parents=True, exist_ok=True)
        bed = track_dir / "input.bed"
        write_bed(frame, bed)
        tracks[label] = (
            lift_intervals(
                bed,
                config.GENOME_BUILD,
                assembly,
                chain,
                track_dir,
                min_match=config.MINMATCH_CROSS_SPECIES,
            )
            if len(frame)
            else dict(mapped={}, unmapped={})
        )
    rows = []
    for element in elements.itertuples(index=False):
        name = element.element_id
        lrow, rrow = (tracks[k]["mapped"].get(name) for k in ("left", "right"))
        pair = anchor_gap(lrow, rrow, target_sizes, max_gap)
        if name in excluded:
            pair.update(valid_pair=False, rejection_reason="source_flank_out_of_bounds")
        row = dict(
            element_id=name,
            assembly=assembly,
            flank_bp=flank,
            min_match=config.MINMATCH_CROSS_SPECIES,
            max_gap_bp=max_gap,
            source_chrom=element.chrom,
            source_start=element.start,
            source_end=element.end,
        )
        for label, mapped in (("left", lrow), ("right", rrow)):
            row[f"{label}_lifted"] = mapped is not None
            row[f"{label}_failure_reason"] = (
                tracks[label]["unmapped"].get(name, {}).get("failure_reason", "")
            )
            for key in ("chrom", "start", "end", "strand"):
                row[f"{label}_{key}"] = mapped[key] if mapped else ""
        rows.append(row | pair)
    result = pd.DataFrame(rows).sort_values("element_id").reset_index(drop=True)
    assert result.element_id.is_unique and set(result.element_id) == set(elements.element_id)
    assert result.valid_pair.sum() + result.rejection_reason.ne("").sum() == len(elements)
    return result


def sweep_summary(frames: dict[int, pd.DataFrame]) -> pd.DataFrame:
    """Summarize geometry, never use motif scores to choose the flank length."""
    rows = []
    for flank, frame in sorted(frames.items()):
        widths = pd.to_numeric(frame.loc[frame.valid_pair, "gap_width"])
        rows.append(
            dict(
                flank_bp=flank,
                input_count=len(frame),
                valid_pairs=len(widths),
                median_gap_bp=widths.median() if len(widths) else None,
                min_gap_bp=widths.min() if len(widths) else None,
                max_gap_bp=widths.max() if len(widths) else None,
                **{r: int(frame.rejection_reason.eq(r).sum()) for r in REASONS},
            )
        )
    return pd.DataFrame(rows)


def scan_intervals(frame: pd.DataFrame, sizes: dict[str, int], mode: str) -> pd.DataFrame:
    """Use the full gap or an exact 101 bp midpoint window, on genomic plus strand.

    Midpoints round down for half-base ties. Short gaps require extending into
    the anchors; record that explicitly. Refuse chromosome-end clipping so the
    comparator always has the promised length. FIMO scans both strands.
    """
    if mode not in {"whole", "center"}:
        raise ValueError("Unknown scan mode")
    rows = []
    for row in frame.loc[frame.valid_pair].itertuples(index=False):
        start, end = int(row.gap_start), int(row.gap_end)
        if mode == "center":
            start = (start + end - WINDOW_BP) // 2
            end = start + WINDOW_BP
        if not 0 <= start <= end <= sizes[row.gap_chrom]:
            raise ValueError(f"Cannot extract an unclipped {mode} window: {row.element_id}")
        rows.append(
            dict(
                element_id=row.element_id,
                chrom=row.gap_chrom,
                start=start,
                end=end,
                strand="+",
                scan_bp=end - start,
                extends_outside_gap=start < row.gap_start or end > row.gap_end,
            )
        )
    return pd.DataFrame(
        rows,
        columns=["element_id", "chrom", "start", "end", "strand", "scan_bp", "extends_outside_gap"],
    )


def score_intervals(
    intervals: pd.DataFrame, genome: Path, matrix: Path, directory: Path
) -> pd.DataFrame:
    """Extract verified DNA and call the existing scorer, retaining missing scores.

    A sub-motif-length gap cannot contain a full motif (score stays missing).
    A negative call containing ambiguous DNA, or an unexplained missing hit,
    fails publication rather than being recoded as a loss or lack of synteny.
    """
    directory.mkdir(parents=True, exist_ok=True)
    nonempty = intervals.loc[intervals.scan_bp.gt(0)]
    records, hits = {}, pd.DataFrame()
    if len(nonempty):
        bed, fasta = directory / "windows.bed", directory / "windows.fa"
        write_bed(nonempty, bed)
        tool = shutil.which("twoBitToFa")
        if tool is None:
            raise FileNotFoundError("twoBitToFa is not on PATH")
        run_tool(
            [tool, f"-bed={bed}", "-noMask", str(genome), str(fasta)], directory / "twoBitToFa.log"
        )
        records = read_fasta(fasta)
        if set(records) != set(nonempty.element_id):
            raise ValueError("Extracted FASTA IDs differ from the accepted scan intervals")
        for row in nonempty.itertuples(index=False):
            if len(records[row.element_id]) != row.scan_bp:
                raise ValueError(f"Wrong extracted length: {row.element_id}")
        hits = motifs.score_sequences(fasta, matrix)
        if not hits.element_id.is_unique or set(hits.element_id) != set(nonempty.element_id):
            raise ValueError("FIMO must return exactly one row per extracted interval")
        hits = hits.set_index("element_id")
    motif_width = config.MOTIF_SPAN[1] - config.MOTIF_SPAN[0] + 1
    rows = []
    for interval in intervals.to_dict("records"):
        name = interval["element_id"]
        seq = records.get(name, "").upper()
        hit = hits.loc[name] if name in hits.index else None
        has_hit = hit is not None and bool(hit.has_hit)
        ambiguous = sum(base not in "ACGT" for base in seq)
        short = interval["scan_bp"] < motif_width
        if has_hit:
            present = bool(hit.retains_motif)
            if present != (float(hit.score) >= config.M_VALUE):
                raise ValueError("Scorer threshold differs from the configured m-value")
            if not present and ambiguous:
                raise ValueError(f"Unscorable negative with ambiguous DNA: {name}")
        elif short:
            present = False
        else:
            raise ValueError(f"Unscorable anchored interval: {name}; never infer loss")
        interval.update(
            status="present" if present else "absent",
            has_hit=has_hit,
            ambiguous_bp=ambiguous,
            score_basis="FIMO" if has_hit else "shorter_than_motif",
            score=hit.score if has_hit else None,
            p_value=hit.p_value if has_hit else None,
            fimo_start=int(hit.start) if has_hit else None,
            fimo_stop=int(hit.stop) if has_hit else None,
            fimo_strand=hit.strand if has_hit else "",
            hit_start=interval["start"] + int(hit.start) - 1 if has_hit else None,
            hit_end=interval["start"] + int(hit.stop) if has_hit else None,
        )
        rows.append(interval)
    columns = list(intervals.columns) + [
        "status",
        "has_hit",
        "ambiguous_bp",
        "score_basis",
        "score",
        "p_value",
        "fimo_start",
        "fimo_stop",
        "fimo_strand",
        "hit_start",
        "hit_end",
    ]
    return pd.DataFrame(rows, columns=columns)


def resolve(frame: pd.DataFrame, scans: dict[str, pd.DataFrame], expected: int) -> pd.DataFrame:
    """Resolve the old ABSENT bucket and explicitly adapt to classify.py states."""
    result = frame.copy()
    for mode in ("whole", "center"):
        scores = scans[mode]
        if not scores.element_id.is_unique or set(scores.element_id) != set(
            result.loc[result.valid_pair, "element_id"]
        ):
            raise ValueError("Both scans must cover exactly the valid anchor pairs")
        scores = scores.rename(columns={c: f"{mode}_{c}" for c in scores if c != "element_id"})
        result = result.merge(scores, on="element_id", how="left", validate="one_to_one")
        status = result[f"{mode}_status"]
        if not status.loc[result.valid_pair].isin(["present", "absent"]).all():
            raise ValueError("An unscorable anchor pair cannot be assigned a loss")
        result[f"{mode}_resolution"] = [
            "NO_SYNTENY" if not valid else "RESCUED" if call == "present" else "ABSENT_CONFIRMED"
            for valid, call in zip(result.valid_pair, status, strict=True)
        ]
        result[f"{mode}_ortholog_state"] = result[f"{mode}_resolution"].map(STATE_MAP)
        assert len(result) == expected and result.element_id.is_unique
        assert (
            result[f"{mode}_resolution"].value_counts().reindex(STATE_MAP, fill_value=0).sum()
            == expected
        )
        assert (
            result.loc[~result.valid_pair, f"{mode}_ortholog_state"]
            .eq(OrthologState.FAILED_LIFT)
            .all()
        )
    result["search_space_sensitive"] = result.whole_resolution.ne(result.center_resolution)
    return result.sort_values("element_id", kind="stable").reset_index(drop=True)


def markdown(frame: pd.DataFrame) -> str:
    """Render small deterministic tables without an optional tabulate dependency."""

    def cell(value: object) -> str:
        if pd.isna(value):
            return "—"
        return str(value).replace("|", "\\|").replace("\n", " ")

    rows = [
        "| " + " | ".join(map(str, frame.columns)) + " |",
        "| " + " | ".join(["---"] * len(frame.columns)) + " |",
    ]
    rows.extend(
        "| " + " | ".join(cell(v) for v in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    )
    return "\n".join(rows)


def render(
    result: pd.DataFrame,
    summary: pd.DataFrame,
    baseline: pd.DataFrame,
    assembly: str,
    flank: int,
    max_gap: int,
    fingerprints: dict,
    prefix: Path,
) -> str:
    """Report measured counts, short-gap caveats, width strata, and explicit state mapping."""
    valid = result.loc[result.valid_pair].copy()
    changed = int(valid.search_space_sensitive.sum())
    any_positive = int((valid.whole_status.eq("present") | valid.center_status.eq("present")).sum())
    # Operational sensitivity rule fixed before inspecting scores. With few positives,
    # even one disagreement can dominate the inferred rescue rate.
    small = changed == 0 or (
        changed / max(len(valid), 1) <= 0.05 and changed / max(any_positive, 1) <= 0.10
    )
    counts = (
        pd.DataFrame(
            {
                mode: result[f"{mode}_resolution"].value_counts().reindex(STATE_MAP, fill_value=0)
                for mode in ("whole", "center")
            }
        )
        .rename_axis("state")
        .reset_index()
    )
    base = baseline.loc[baseline.call.eq("MAPPED")]
    widths = base.window_end - base.window_start
    strata = []
    for low, high in ((0, 100), (101, 200), (201, 500), (501, 1000), (1001, max_gap)):
        subset = valid.loc[pd.to_numeric(valid.gap_width).between(low, high)]
        strata.append(
            dict(
                gap_bp=f"{low}–{high}",
                n=len(subset),
                whole_present=int(subset.whole_status.eq("present").sum()),
                whole_absent=int(subset.whole_status.eq("absent").sum()),
                center_present=int(subset.center_status.eq("present").sum()),
                center_absent=int(subset.center_status.eq("absent").sum()),
                changed=int(subset.search_space_sensitive.sum()),
            )
        )
    rejection = summary[["flank_bp", *REASONS]].set_index("flank_bp").T.reset_index()
    rejection = rejection.rename(columns={"index": "first_failure"})
    recommended = int(
        summary.sort_values(["valid_pairs", "flank_bp"], ascending=[False, True]).iloc[0].flank_bp
    )
    medians = ", ".join(
        f"{int(r.flank_bp)} bp: {r.median_gap_bp:g} bp"
        for r in summary.itertuples()
        if pd.notna(r.median_gap_bp)
    )
    lines = [
        f"# Flank-anchored rescue: {config.GENOME_BUILD} → {assembly}",
        "",
        f"**{len(result)} previously ABSENT elements; {len(valid)} valid pairs at F={flank} bp.** "
        f"Whole-gap calls: {int(result.whole_resolution.eq('RESCUED').sum())} RESCUED, "
        f"{int(result.whole_resolution.eq('ABSENT_CONFIRMED').sum())} ABSENT_CONFIRMED, "
        f"{int(result.whole_resolution.eq('NO_SYNTENY').sum())} NO_SYNTENY. "
        "NO_SYNTENY is unresolved evidence, never motif loss.",
        "",
        "## MAPPED FIMO baseline",
        "",
        f"{len(base)} MAPPED elements: **{int(base.motif_status.eq('present').sum())} present, "
        f"{int(base.motif_status.eq('absent').sum())} absent**, "
        f"{int(base.motif_status.eq('no_hit').sum())} no-hit. "
        f"The other {len(baseline) - len(base)} sidecar rows are unscored.",
        "",
        f"The baseline uses lifted motif ±40 bp: nominally 101 bp, actually "
        f"{int(widths.min())}–{int(widths.max())} bp because lifted motif widths vary. "
        "Its existing mapping, extraction and scoring logic is unchanged. "
        "Only safe publication of its completed score table/report was changed.",
        "",
        "## Flank sweep (geometry selected before scoring)",
        "",
        markdown(
            summary[
                [
                    "flank_bp",
                    "input_count",
                    "valid_pairs",
                    "median_gap_bp",
                    "min_gap_bp",
                    "max_gap_bp",
                ]
            ]
        ),
        "",
        f"Recommend **F={recommended} bp** from this sweep; selected F={flank} bp "
        f"(`config.RESCUE_FLANK_BP` defaults to {config.RESCUE_FLANK_BP}): "
        "maximize valid-pair coverage, choosing the smallest flank on a tie. "
        f"Observed median widths: {medians}. No monotone improvement or gap widening "
        "is supported by this sweep. Human inner flank edges stay at the element boundaries; "
        "longer flanks change alignment coverage and endpoints, not the source gap. "
        "Each sweep has a different accepted subset; its medians are not a paired comparison.",
        "",
        f"All flanks use configured minMatch={config.MINMATCH_CROSS_SPECIES}; "
        "full 200 bp source elements are excluded from both anchors. "
        "Both anchors must return one mapping (liftOver without -multiple; duplicate IDs "
        "are rejected), share chromosome and "
        "orientation, preserve order including inversions, and leave a nonnegative gap "
        f"≤{max_gap:,} bp. This maximum is a prespecified operational guard, not a validated "
        "biological cutoff. Coordinates are zero-based, half-open. Source boundary flanks "
        "are rejected rather than shortened. Reasons are mutually exclusive first failures.",
        "",
        "## Rejection reasons at every flank size",
        "",
        markdown(rejection),
        "",
        "Raw left/right liftOver reasons and coordinates are retained for every element "
        "in each flank CSV. Missing anchors can reflect divergence or ambiguous mapping; "
        "they do not demonstrate a deletion.",
        "",
        "## Resolution and search-space sensitivity",
        "",
        markdown(counts),
        "",
        f"Both columns sum to {len(result)}; asserted in code. "
        f"**{changed}/{len(valid)} anchored calls change** "
        f"({100 * changed / max(len(valid), 1):.1f}%; "
        f"{changed}/{any_positive} of elements positive "
        "in either scan).",
        "",
        (
            "The observed change is small; recommend whole-gap calls while retaining the "
            "midpoint sensitivity check. This small pilot does not calibrate a false-positive rate."
            if small
            else "**The rescue rate is search-space dependent.** Do not use the "
            "whole-gap rescue rate as directly comparable to MAPPED. Carry both analyses; "
            "use midpoint calls for the matched-length sensitivity analysis, remembering that "
            "a midpoint can miss a genuine off-center motif. Resolve discordant sites before "
            "making a biological loss claim."
        ),
        "",
        "The operational definition of small is ≤5% of anchored calls AND ≤10% of "
        "any-positive calls changing (or zero changes), fixed before examining scores.",
        "",
        "Both scans use `motifs.score_sequences`, the same MA0138.2 matrix, both strands, "
        f"FIMO `--text --thresh 1`, and best score ≥{config.M_VALUE}. Ties use p-value, "
        "start, stop and strand in the shared helper. Larger windows have more opportunities "
        "to exceed the same threshold; no multiple-testing adjustment is added to either path.",
        "",
        f"Midpoint windows are exactly {WINDOW_BP} bp, with half-base ties rounded down. "
        f"Accepted gaps shorter than {WINDOW_BP} bp: "
        f"{int(valid.center_extends_outside_gap.sum())}. For these, the "
        "midpoint windows extend into the anchors, explicitly flagged "
        "per element. They are not a strict subset of the whole-gap scan. "
        "No chromosome-end clipping is allowed. A sub-21 bp whole gap cannot fit this motif; "
        "its score remains missing, with basis `shorter_than_motif`. An unexplained missing "
        "FIMO hit or ambiguous negative sequence stops publication, never becomes a loss.",
        "",
        "### Width strata at the selected flank",
        "",
        markdown(pd.DataFrame(strata)),
        "",
        "### Every accepted interval and both best scores",
        "",
        markdown(
            valid[
                [
                    "element_id",
                    "gap_chrom",
                    "gap_start",
                    "gap_end",
                    "gap_strand",
                    "gap_width",
                    "whole_score",
                    "whole_status",
                    "center_score",
                    "center_status",
                    "center_extends_outside_gap",
                ]
            ]
        ),
        "",
        "## Downstream meaning",
        "",
        "Existing `classify.OrthologState` can absorb all three states; no enum extension "
        "is needed. The output supplies **separate** `whole_ortholog_state` and "
        "`center_ortholog_state` columns, plus `search_space_sensitive`; consumers must "
        "choose an analysis explicitly. Existing MAPPED outputs are not rewritten.",
        "",
        markdown(
            pd.DataFrame(
                [{"rescue_state": k, "classify_state": str(v)} for k, v in STATE_MAP.items()]
            )
        ),
        "",
        "`failed_lift` preserves unresolved evidence in classification and never enters "
        "`no_motif_species`. `lifted_unscorable` remains appropriate for future unscorable "
        "anchored DNA; this run refuses publication if encountered rather than forcing it "
        "into the requested three-state partition. `ABSENT_CONFIRMED` means no full motif "
        "above the cutoff in the anchored region, not experimental proof of functional loss. "
        "Synteny is inferred from these anchors; independent chains could still mispair. "
        "A retained PWM match does not demonstrate REST binding.",
        "",
        "## Outputs and reproduction",
        "",
        f"`{prefix}.resolved.csv` contains all {len(result)} original ABSENT IDs, both scans, "
        "coordinates, widths, best hits, states and failure reasons. "
        f"`{prefix}.flank{{500,1000,2000}}.csv` and `{prefix}.sweep.csv` retain the sweep; "
        f"`{prefix}.work/` retains BEDs, lifted/unmapped records, FASTAs and command logs.",
        "",
        "Run `make pilot-rescue` with liftOver, twoBitInfo, twoBitToFa and FIMO on PATH. "
        "The generic CLI accepts assembly, chain, genome, prefix, baseline and report paths. "
        "It requires cached source/target chromosome sizes and target publisher MD5 metadata "
        "(overridable CLI paths); it downloads nothing. It checks the genome MD5 and "
        "twoBitInfo sizes, baseline IDs/statuses, and baseline report input hashes. "
        "Published CSVs/report are staged after successful computation. No timestamps or "
        "random sampling enter the output. Tool hashes below pin executable versions.",
        "",
        markdown(pd.DataFrame([{"input_or_tool": k, "MD5": v} for k, v in fingerprints.items()])),
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    """Run the parameterized sweep, dual scoring, and explicit three-state resolution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--source-elements", type=Path, default=config.HG38_ELEMENTS)
    parser.add_argument("--max-gap-bp", type=int, default=MAX_GAP_BP)
    parser.add_argument("--expected-absent", type=int, required=True)
    parser.add_argument("--baseline-scores", type=Path, required=True)
    parser.add_argument("--baseline-report", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--flank-bp", type=int, choices=FLANK_SWEEP, default=config.RESCUE_FLANK_BP)
    parser.add_argument("--source-sizes", type=Path)
    parser.add_argument("--target-sizes", type=Path)
    parser.add_argument("--genome-md5", type=Path)
    args = parser.parse_args()
    sidecar = pd.read_csv(args.sidecar)
    elements = absent_elements(sidecar, pd.read_csv(args.source_elements), args.assembly)
    assert len(elements) == args.expected_absent, "Unexpected ABSENT input count"
    if not args.genome.is_file():
        raise FileNotFoundError(args.genome)
    metadata = config.RAW_DIR / "genomes" / "verification"
    source_path = args.source_sizes or metadata / f"{config.GENOME_BUILD}.chrom.sizes"
    target_path = args.target_sizes or metadata / f"{args.assembly}.chrom.sizes"
    md5_path = args.genome_md5 or metadata / f"{args.assembly}.md5sum.txt"
    source_sizes, target_sizes = (
        parse_sizes(source_path.read_text()),
        parse_sizes(target_path.read_text()),
    )
    fingerprints = {
        str(path): local_md5(path)
        for path in (
            args.sidecar,
            args.source_elements,
            args.chain,
            args.genome,
            config.MOTIF_PATH,
            args.baseline_scores,
            args.baseline_report,
            source_path,
            target_path,
            md5_path,
        )
    }
    published = {
        row.split()[1].lstrip("*"): row.split()[0]
        for row in md5_path.read_text().splitlines()
        if row.strip()
    }
    if published.get(f"{args.assembly}.2bit") != fingerprints[str(args.genome)]:
        raise ValueError("Genome does not match the target assembly publisher MD5")
    if parse_sizes(run_command(["twoBitInfo", str(args.genome), "stdout"])) != target_sizes:
        raise ValueError("Genome and target chromosome sizes disagree")
    baseline = pd.read_csv(args.baseline_scores)
    require_columns(
        baseline, ["element_id", "call", "motif_status", "score", "window_start", "window_end"]
    )
    if not baseline.element_id.is_unique or set(baseline.element_id) != set(sidecar.element_id):
        raise ValueError("Baseline IDs differ from the sidecar")
    if (
        not baseline.set_index("element_id")
        .call.sort_index()
        .equals(sidecar.set_index("element_id").call.sort_index())
    ):
        raise ValueError("Baseline mapping calls differ from the sidecar")
    mapped = baseline.loc[baseline.call.eq("MAPPED")]
    if mapped.empty or not mapped.motif_status.isin(["present", "absent"]).all():
        raise ValueError("Baseline MAPPED rows must all have actual FIMO scores")
    if (
        mapped.score.isna().any()
        or not mapped.score.ge(config.M_VALUE).eq(mapped.motif_status.eq("present")).all()
    ):
        raise ValueError("Baseline scores do not use the current m-value")
    evidence = args.baseline_report.read_text()
    for path in (args.sidecar, args.genome, config.MOTIF_PATH):
        if fingerprints[str(path)] not in evidence:
            raise ValueError("Baseline report input fingerprints do not match this rescue run")
    for name in ("liftOver", "twoBitInfo", "twoBitToFa", "fimo"):
        executable = shutil.which(name)
        if executable is None:
            raise FileNotFoundError(f"Missing tool: {name}")
        fingerprints[f"tool:{executable}"] = local_md5(Path(executable))
    prefix = args.output_prefix
    frames = {}
    for flank in FLANK_SWEEP:
        frames[flank] = run_anchors(
            elements,
            source_sizes,
            target_sizes,
            args.assembly,
            args.chain,
            prefix.parent / f"{prefix.name}.work" / f"flank_{flank}",
            flank,
            args.max_gap_bp,
        )
    summary = sweep_summary(frames)
    print(summary.to_string(index=False), flush=True)
    selected = frames[args.flank_bp]
    scans = {
        mode: score_intervals(
            scan_intervals(selected, target_sizes, mode),
            args.genome,
            config.MOTIF_PATH,
            prefix.parent / f"{prefix.name}.work" / mode,
        )
        for mode in ("whole", "center")
    }
    result = resolve(selected, scans, args.expected_absent)
    report = render(
        result,
        summary,
        baseline,
        args.assembly,
        args.flank_bp,
        args.max_gap_bp,
        fingerprints,
        prefix,
    )
    args.report.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".rescue-publish-", dir=prefix.parent) as temporary:
        staged = []
        for suffix, frame in [(f"flank{f}", x) for f, x in frames.items()] + [
            ("sweep", summary),
            ("resolved", result),
        ]:
            path = Path(temporary) / f"{suffix}.csv"
            frame.to_csv(path, index=False)
            staged.append((path, Path(f"{prefix}.{suffix}.csv")))
        report_path = Path(temporary) / "report.md"
        report_path.write_text(report)
        for source, destination in staged + [(report_path, args.report)]:
            os.replace(source, destination)
    print(json.dumps(result.whole_resolution.value_counts().to_dict(), sort_keys=True), flush=True)
    print(f"Report: {args.report}", flush=True)


if __name__ == "__main__":
    main()
