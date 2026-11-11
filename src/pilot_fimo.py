"""Score lifted mouse motif windows with FIMO and record presence at M_VALUE.

The repeatable entry point is ``python -m src.pilot_fimo`` or ``make pilot-fimo``.
It reads ``orthologs.mm39.csv`` and does not rebuild the lift. Only MAPPED
elements are scanned. The scan is the lifted motif interval expanded by
FLANK_BP on each genomic side, not the full extracted element. UNDETERMINED
and ABSENT stay in the table with no score.

FIMO coordinates are one-based and inclusive on the strand-oriented window.
Genomic coordinates in the output stay zero-based and half-open. A score at
or above ``config.M_VALUE`` is present. A hit below that threshold is absent.
No hit is ``no_hit``, not a fabricated zero. Re-running replaces the table and
the draft report, including blank browser-observation cells.
"""

import os
import random
import shutil
from pathlib import Path
from statistics import median

import pandas as pd

from src import config, pilot, pilot_orthologs
from src.downloads import download_small_file, local_md5
from src.lift import run_tool
from src.motifs import score_sequences

FLANK_BP = 40
SPOTCHECK_COUNT = 20
OUTPUT = pilot_orthologs.ROOT / "fimo"
SCORES = pilot_orthologs.ROOT / "fimo_scores.csv"
REPORT = config.TABLES_DIR / "pilot_mouse_fimo.md"
GENOME = pilot.GENOMES / "mm39.2bit"
STATUSES = ("present", "absent", "no_hit", "unscored")
CHROM_SIZES_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chrom.sizes"


def _integer(value: object, label: str) -> int:
    """Require one finite integer coordinate and reject booleans and partial numbers."""
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer coordinate")
    number = pd.to_numeric(value, errors="coerce")
    if pd.isna(number) or float(number) != int(number):
        raise ValueError(f"{label} must be an integer coordinate")
    return int(number)


def scoring_windows(
    audit: pd.DataFrame, sizes: dict[str, int], flank: int = FLANK_BP
) -> pd.DataFrame:
    """Expand each MAPPED motif interval by ``flank`` bp, clipping at chromosome ends.

    Non-MAPPED rows keep empty window fields. Padding is genomic, not strand-relative,
    matching the ortholog extraction. A clipped side is recorded rather than hidden.
    """
    if flank < 0:
        raise ValueError("Flank must be nonnegative")
    if not audit["element_id"].is_unique:
        raise ValueError("Duplicate element_id")
    if not set(audit["call"]).issubset(set(pilot_orthologs.CALLS)):
        raise ValueError("Unexpected mapping call")
    rows = []
    for row in audit.itertuples(index=False):
        blank = dict(
            element_id=row.element_id,
            window_chrom="",
            window_start=pd.NA,
            window_end=pd.NA,
            window_strand="",
            flank_left=pd.NA,
            flank_right=pd.NA,
            flank_clipped=pd.NA,
        )
        if row.call != "MAPPED":
            rows.append(blank)
            continue
        chrom = row.mm39_motif_chrom
        start = _integer(row.mm39_motif_start, f"{row.element_id} motif start")
        end = _integer(row.mm39_motif_end, f"{row.element_id} motif end")
        strand = row.mm39_motif_strand
        if chrom not in sizes or strand not in {"+", "-"} or not 0 <= start < end <= sizes[chrom]:
            raise ValueError(f"Lifted motif lies outside mm39 or has no strand: {row.element_id}")
        left, right = max(0, start - flank), min(sizes[chrom], end + flank)
        rows.append(
            dict(
                element_id=row.element_id,
                window_chrom=chrom,
                window_start=left,
                window_end=right,
                window_strand=strand,
                flank_left=start - left,
                flank_right=right - end,
                flank_clipped=start - left != flank or right - end != flank,
            )
        )
    windows = pd.DataFrame(rows)
    windows["window_start"] = windows["window_start"].astype("Int64")
    windows["window_end"] = windows["window_end"].astype("Int64")
    windows["flank_left"] = windows["flank_left"].astype("Int64")
    windows["flank_right"] = windows["flank_right"].astype("Int64")
    return windows


def hit_interval(
    window_start: int, window_end: int, strand: str, fimo_start: int, fimo_stop: int
) -> tuple[int, int]:
    """Convert one FIMO hit on the oriented window to a genomic half-open interval.

    Plus-strand position 1 is ``window_start``. Minus-strand position 1 is the
    last genomic base, because twoBitToFa reverse-complements that window.
    """
    length = window_end - window_start
    if strand not in {"+", "-"} or not 1 <= fimo_start <= fimo_stop <= length:
        raise ValueError("FIMO hit is outside the scored window or has no strand")
    if strand == "+":
        return window_start + fimo_start - 1, window_start + fimo_stop
    return window_end - fimo_stop, window_end - fimo_start + 1


def _status(call: str, has_hit: object, retains: object) -> str:
    """Map one row to present, absent, no_hit, or unscored without inventing a score."""
    if call != "MAPPED":
        return "unscored"
    if pd.isna(has_hit):
        raise ValueError("MAPPED element was not scored")
    if not bool(has_hit):
        return "no_hit"
    if pd.isna(retains):
        raise ValueError("A FIMO hit is missing its threshold call")
    return "present" if bool(retains) else "absent"


def assemble(audit: pd.DataFrame, windows: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """Join mapping calls, scan windows, and best FIMO hits. One row per pilot element."""
    mapped = set(audit.loc[audit["call"].eq("MAPPED"), "element_id"])
    if set(hits["element_id"]) != mapped:
        raise ValueError("FIMO results must cover every MAPPED element and no others")
    window_ids = set(windows["element_id"])
    same_ids = windows["element_id"].is_unique and window_ids == set(audit["element_id"])
    if not same_ids:
        raise ValueError("Scan windows must match the sidecar IDs")
    columns = [
        "element_id",
        "score",
        "p_value",
        "has_hit",
        "retains_motif",
        "start",
        "stop",
        "strand",
    ]
    slim = hits.loc[:, columns].rename(
        columns={"start": "fimo_start", "stop": "fimo_stop", "strand": "fimo_strand"}
    )
    merged = audit.merge(windows, on="element_id", how="left")
    merged = merged.merge(slim, on="element_id", how="left")
    merged["motif_status"] = [
        _status(row.call, row.has_hit, row.retains_motif) for row in merged.itertuples(index=False)
    ]
    hit_start, hit_end = [], []
    for row in merged.itertuples(index=False):
        if row.motif_status not in {"present", "absent"}:
            hit_start.append(pd.NA)
            hit_end.append(pd.NA)
            continue
        genomic_start, genomic_end = hit_interval(
            int(row.window_start),
            int(row.window_end),
            row.window_strand,
            _integer(row.fimo_start, f"{row.element_id} FIMO start"),
            _integer(row.fimo_stop, f"{row.element_id} FIMO stop"),
        )
        hit_start.append(genomic_start)
        hit_end.append(genomic_end)
    merged["hit_start"] = pd.Series(hit_start, dtype="Int64")
    merged["hit_end"] = pd.Series(hit_end, dtype="Int64")
    counts = merged["motif_status"].value_counts()
    if int(counts.reindex(STATUSES).fillna(0).sum()) != len(merged):
        raise AssertionError("Motif statuses do not partition the sidecar")
    if int(merged["call"].eq("MAPPED").sum()) != int(merged["motif_status"].ne("unscored").sum()):
        raise AssertionError("Every MAPPED element must be scored, and only those")
    return merged.sort_values("element_id", kind="stable").reset_index(drop=True)


def choose_spotcheck(
    element_ids: list[str], count: int = SPOTCHECK_COUNT, seed: int = pilot.SAMPLE_SEED
) -> list[str]:
    """Choose a fixed sample of scored elements. The order in the report is sorted."""
    ordered = sorted(element_ids)
    if count < 1 or len(ordered) < count:
        raise ValueError(
            f"Need at least {count} scored elements to spot-check; found {len(ordered)}"
        )
    return sorted(random.Random(seed).sample(ordered, count))


def _interval(chrom: object, start: object, end: object, strand: object) -> str:
    """Format one zero-based half-open interval, or a blank when it was not scored."""
    if pd.isna(start) or pd.isna(end) or chrom == "":
        return ""
    return f"{chrom}:{int(start)}-{int(end)} ({strand})"


def _browser_position(chrom: object, start: object, end: object) -> str:
    """Format a zero-based half-open interval as a one-based UCSC position."""
    if pd.isna(start) or pd.isna(end) or pd.isna(chrom) or chrom == "":
        return ""
    return f"{chrom}:{int(start) + 1}-{int(end)}"


def _read_strand(window_strand: object, fimo_strand: object) -> str:
    """Genomic strand of the hit. A minus-strand window flips FIMO's strand."""
    if pd.isna(window_strand) or pd.isna(fimo_strand) or window_strand == "":
        return ""
    return "+" if window_strand == fimo_strand else "-"


def _score_text(value: object) -> str:
    """Render a FIMO score with fixed precision. Missing scores stay blank."""
    return "" if pd.isna(value) else f"{float(value):.4f}"


def ucsc_url(chrom: str, start: int, end: int) -> str:
    """Link the scan window in the UCSC browser. The position is one-based inclusive."""
    return (
        "https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm39"
        f"&position={chrom}:{start + 1}-{end}"
    )


def render(frame: pd.DataFrame, spotcheck: list[str], fingerprints: dict[str, str]) -> str:
    """Write the draft report. Browser observations are intentionally blank."""
    counts = frame["motif_status"].value_counts()
    scored = frame.loc[frame["motif_status"].ne("unscored")]
    finite = scored.loc[scored["score"].notna(), "score"]
    if finite.empty:
        score_line = "No MAPPED element produced a finite FIMO score."
    else:
        score_line = (
            f"Finite scores (present and absent only): minimum **{float(finite.min()):.4f}**, "
            f"median **{float(median(finite)):.4f}**, maximum **{float(finite.max()):.4f}**."
        )
    lengths = (scored["window_end"] - scored["window_start"]).astype("int64")
    clipped = int(scored["flank_clipped"].eq(True).sum())
    lines = [
        "# Mouse pilot FIMO scores",
        "",
        "## Scoring contract",
        "",
        f"Input: **{len(frame)}** pilot elements from `orthologs.mm39.csv`. "
        f"Only **MAPPED** elements are scanned. The scan is the lifted mm39 motif "
        f"interval expanded by **{FLANK_BP} bp** on each genomic side, clipped at a "
        "chromosome end when the flank does not fit. The full extracted element "
        "sequence is not scanned.",
        "",
        f"Matrix: `{config.JASPAR_ID}` at `{config.MOTIF_PATH.name}`. "
        f"A best hit with score **>= {config.M_VALUE}** is **present**. "
        "A hit below that threshold is **absent**. No hit is **no_hit**, not a score of zero. "
        "UNDETERMINED and ABSENT stay **unscored**. FIMO is run with `--text --thresh 1` "
        "so a weak best score is kept, and the retained hit is the maximum score.",
        "",
        "Window coordinates are zero-based and half-open. FIMO start/stop are one-based "
        "and inclusive on the strand-oriented window. `hit_start`/`hit_end` convert that "
        "hit back to a genomic half-open interval. A minus-strand window is reverse-complemented "
        "before scoring, so position 1 of that sequence is the last genomic base.",
        "",
        "## Calls",
        "",
        "| Status | Elements | Meaning |",
        "| --- | ---: | --- |",
        f"| present | {int(counts.get('present', 0))} | Best score >= {config.M_VALUE} |",
        f"| absent | {int(counts.get('absent', 0))} | Best hit is below {config.M_VALUE} |",
        f"| no_hit | {int(counts.get('no_hit', 0))} | MAPPED window with no FIMO hit |",
        f"| unscored | {int(counts.get('unscored', 0))} | "
        "UNDETERMINED or ABSENT; not a motif-loss call |",
        "",
        f"Sum: **{int(counts.reindex(STATUSES).fillna(0).sum())} / {len(frame)}**. "
        f"Scored elements: **{len(scored)}**.",
        "",
        score_line,
        "",
        f"Scan-window length: minimum **{int(lengths.min())} bp**, "
        f"maximum **{int(lengths.max())} bp**. "
        f"Windows clipped at a chromosome end: **{clipped}**.",
        "",
        "## Spot check",
        "",
        f"**{len(spotcheck)}** MAPPED elements, chosen with seed `{pilot.SAMPLE_SEED}` "
        "and then sorted. Intervals below are zero-based and half-open. The link opens "
        "the scan window in the UCSC mm39 browser; that position is one-based. "
        "Best hit is the one-based UCSC position of the top-scoring FIMO hit and links "
        "to it. Read is the genomic strand to compare against the RE1 consensus. "
        "The observation column is left blank for a manual check.",
        "",
        "| Element | Motif interval | Scan window | Best hit | Read | Status | Score "
        "| Browser observation |",
        "| --- | --- | --- | --- | :---: | --- | ---: | --- |",
    ]
    chosen = frame.set_index("element_id").loc[spotcheck]
    for element_id, row in chosen.iterrows():
        motif = _interval(
            row.mm39_motif_chrom, row.mm39_motif_start, row.mm39_motif_end, row.mm39_motif_strand
        )
        window = _interval(row.window_chrom, row.window_start, row.window_end, row.window_strand)
        link = ucsc_url(row.window_chrom, int(row.window_start), int(row.window_end))
        hit = _browser_position(row.window_chrom, row.hit_start, row.hit_end)
        hit_cell = (
            f"[{hit}]({ucsc_url(row.window_chrom, int(row.hit_start), int(row.hit_end))})"
            if hit
            else ""
        )
        read = _read_strand(row.window_strand, row.fimo_strand)
        lines.append(
            f"| {element_id} | {motif} | [{window}]({link}) | {hit_cell} | {read} | "
            f"{row.motif_status} | {_score_text(row.score)} | |"
        )
    lines += [
        "",
        "## Boundaries",
        "",
        "A MAPPED window means the motif interval could be scored. **present** means the "
        f"best FIMO score on that window reached {config.M_VALUE}. It does not prove REST "
        "binding in mouse. **unscored** is missing alignment, not a negative score. "
        "**ABSENT** remains an alignment gap, not a demonstrated deletion. This table "
        "does not estimate repression or compare models.",
        "",
        "## Reproduction",
        "",
        "Run `make pilot-fimo` after `orthologs.mm39.csv` exists. The command is offline. "
        "It does not rebuild the lift. Re-running replaces this report and "
        "`fimo_scores.csv`, including any browser notes added to the draft.",
        "",
        "| Input | MD5 |",
        "| --- | --- |",
    ]
    for name, digest in fingerprints.items():
        lines.append(f"| {name} | `{digest}` |")
    lines.append("")
    return "\n".join(lines)


def _extract(windows: pd.DataFrame, genome: Path, directory: Path) -> Path:
    """Extract strand-oriented scan windows with twoBitToFa. The FASTA is not the element."""
    mapped = windows.loc[windows["window_chrom"].ne("")]
    records = {
        row.element_id: dict(
            chrom=row.window_chrom,
            start=int(row.window_start),
            end=int(row.window_end),
            strand=row.window_strand,
        )
        for row in mapped.itertuples(index=False)
    }
    directory.mkdir(parents=True, exist_ok=True)
    bed = directory / "windows.mm39.bed"
    fasta = directory / "windows.mm39.fa"
    pilot_orthologs.write_records(records, bed)
    tool = shutil.which("twoBitToFa")
    if not tool:
        raise FileNotFoundError("Missing extraction tool: twoBitToFa")
    run_tool(
        [tool, f"-bed={bed}", "-noMask", str(genome), str(fasta)],
        directory / "twoBitToFa.log",
    )
    return fasta


def score_sidecar(
    audit: pd.DataFrame, sizes: dict[str, int], genome: Path, motif_file: Path, directory: Path
) -> tuple[pd.DataFrame, str]:
    """Scan MAPPED windows and return the full table plus the draft report."""
    windows = scoring_windows(audit, sizes)
    fasta = _extract(windows, genome, directory)
    hits = score_sequences(fasta, motif_file)
    frame = assemble(audit, windows, hits)
    spotcheck = choose_spotcheck(frame.loc[frame["call"].eq("MAPPED"), "element_id"].tolist())
    fingerprints = {
        "ortholog_sidecar": local_md5(pilot_orthologs.SIDECAR),
        "JASPAR_matrix": local_md5(motif_file),
        "mm39_2bit": local_md5(genome),
    }
    return frame, render(frame, spotcheck, fingerprints)


def _output_columns(frame: pd.DataFrame) -> pd.DataFrame:
    """Keep the published columns in a stable order. Missing scores stay empty."""
    columns = [
        "element_id",
        "call",
        "call_reason",
        "motif_status",
        "score",
        "p_value",
        "mm39_motif_chrom",
        "mm39_motif_start",
        "mm39_motif_end",
        "mm39_motif_strand",
        "window_chrom",
        "window_start",
        "window_end",
        "window_strand",
        "flank_left",
        "flank_right",
        "flank_clipped",
        "fimo_start",
        "fimo_stop",
        "fimo_strand",
        "hit_start",
        "hit_end",
    ]
    return frame.loc[:, columns]


def main() -> None:
    """Score the existing mouse sidecar. Fail before publishing a partial table."""
    if not pilot_orthologs.SIDECAR.is_file():
        raise FileNotFoundError("Missing orthologs.mm39.csv; run make pilot-orthologs first")
    SCORES.unlink(missing_ok=True)
    REPORT.unlink(missing_ok=True)
    try:
        audit = pd.read_csv(pilot_orthologs.SIDECAR)
        if len(audit) != pilot_orthologs.EXPECTED_PILOT_COUNT:
            raise ValueError(
                f"Expected {pilot_orthologs.EXPECTED_PILOT_COUNT} sidecar rows; found {len(audit)}"
            )
        sizes_path = download_small_file(CHROM_SIZES_URL, pilot.METADATA / "mm39.chrom.sizes")
        sizes = pilot.parse_sizes(sizes_path.read_text())
        if not GENOME.is_file():
            raise FileNotFoundError(f"Missing mouse genome: {GENOME}")
        frame, text = score_sidecar(audit, sizes, GENOME, config.MOTIF_PATH, OUTPUT)
        SCORES.parent.mkdir(parents=True, exist_ok=True)
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        temporary = SCORES.with_suffix(".csv.tmp")
        _output_columns(frame).to_csv(temporary, index=False)
        os.replace(temporary, SCORES)
        REPORT.write_text(text)
        print(text, flush=True)
    except (ValueError, OSError, RuntimeError, AssertionError) as exc:
        failure = (
            "# Mouse pilot FIMO scores: FAILED\n\n"
            f"```text\n{type(exc).__name__}: {exc}\n```\n\n"
            "No score table was published.\n"
        )
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(failure)
        print(failure, flush=True)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
