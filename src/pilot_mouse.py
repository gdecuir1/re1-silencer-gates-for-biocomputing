"""Measure the chr8 hg38->mm39 minMatch tradeoff, then retain a complete mouse audit.

Only hg38 and mm39 are verified. This module never downloads a genome, touches
other species, or interprets a mapped interval as a retained REST motif. Width
stability is a useful diagnostic, not a measurement of orthology precision.
All four runs differ only in minMatch; unique mapping behavior stays unchanged.
"""

import argparse
from collections import Counter
from pathlib import Path
from statistics import median

import pandas as pd

from src import config, pilot
from src.downloads import local_md5
from src.lift import lift_intervals, parse_bed

SWEEP = (0.95, 0.5, 0.2, 0.1)
OUTPUT = config.PROCESSED_DIR / "pilot_mouse"
REPORT = config.TABLES_DIR / "pilot_mouse_lift.md"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz"


def summarize(result: dict, input_count: int) -> dict:
    """Count all mapped records, including abnormal widths, with no silent filtering."""
    widths = [record["end"] - record["start"] for record in result["mapped"].values()]
    count = len(widths)
    failed = len(result["unmapped"])
    assert count + failed == input_count, "Mapped + unmapped must equal the pilot count"
    return dict(
        min_match=result["min_match"],
        lifted=count,
        input_count=input_count,
        lift_rate=count / input_count,
        median_width=median(widths) if widths else None,
        width_deviation_gt20pct=sum(abs(width - 200) > 40 for width in widths),
        width_within20pct=sum(160 <= width <= 240 for width in widths),
        unmapped=failed,
    )


def validate_target(result: dict, sizes: dict[str, int]) -> None:
    """Reject out-of-reference coordinates without imposing a 200 bp target width."""
    for name, row in result["mapped"].items():
        if row["chrom"] not in sizes or not 0 <= row["start"] < row["end"] <= sizes[row["chrom"]]:
            raise ValueError(f"Mapped interval outside verified mm39 reference: {name}: {row}")


def write_audit(source: dict, result: dict, path: Path) -> None:
    """Keep one record per input ID, including unchanged failed coordinates/reasons."""
    rows = []
    for name, original in sorted(source.items()):
        mapped = result["mapped"].get(name)
        width = mapped["end"] - mapped["start"] if mapped else None
        rows.append(
            dict(
                element_id=name,
                source_build=result["source_build"],
                target_build=result["target_build"],
                min_match=result["min_match"],
                hg38_chrom=original["chrom"],
                hg38_start=original["start"],
                hg38_end=original["end"],
                hg38_strand=original["strand"],
                status="lifted" if mapped else "failed-to-lift",
                target_chrom=mapped["chrom"] if mapped else "",
                target_start=mapped["start"] if mapped else "",
                target_end=mapped["end"] if mapped else "",
                target_strand=mapped["strand"] if mapped else "",
                target_width=width if mapped else "",
                width_deviation_gt20pct=abs(width - 200) > 40 if mapped else "",
                failure_reason="" if mapped else result["unmapped"][name]["failure_reason"],
            )
        )
    assert len(rows) == len(source) and len({row["element_id"] for row in rows}) == len(source)
    pd.DataFrame(rows).to_csv(path, index=False)


def render(
    genomes: list[dict], observations: list[dict], source: dict, final=None, *, error: str = ""
) -> str:
    """Produce a deterministic observed-versus-expected report; failures remain explicit."""
    lines = [
        "# Pilot mouse lift: hg38 to mm39",
        "",
        "## Required genome verification",
        "",
        "This step requires only hg38 and mm39. Each uses the unchanged twoBitInfo, "
        "per-contig sizes, total bases, stable-file, and publisher MD5 checks. "
        "Other species are not inspected.",
        "",
        "| Assembly | Sequences observed / UCSC | Bases observed / UCSC | "
        "MD5 observed / UCSC | Result |",
        "| --- | ---: | ---: | --- | --- |",
    ]
    for row in genomes:
        lines.append(
            f"| {row['assembly']} | {row['sequence_count']} / {row['expected_sequence_count']} | "
            f"{row['total_bases']} / {row['expected_total_bases']} | "
            f"{row['observed_md5']} / {row['expected_md5']} | {row['status']} |"
        )
        for issue in row["errors"]:
            lines.append(f"\n{row['assembly']} FAIL: {issue}\n")
    lines += [
        "",
        "## Pilot and minMatch sweep",
        "",
        f"Input: **{len(source)} canonical {config.PILOT_CHROM} elements**; "
        "expected 200 bp hg38 intervals, no contrast elements. "
        "The pilot BED is checked against the verified hg38 CSV before lifting.",
        "",
        "| minMatch | Lifted / input | Lift rate | Median target width (bp) | "
        "Width outside 160–240 bp (>20% change) | Within 160–240 bp |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in observations:
        med = "not estimable" if row["median_width"] is None else f"{row['median_width']:g}"
        lines.append(
            f"| {row['min_match']:g} | {row['lifted']} / {row['input_count']} | "
            f"{row['lift_rate']:.2%} | {med} | {row['width_deviation_gt20pct']} | "
            f"{row['width_within20pct']} |"
        )
    lines += [
        "",
        "Only -minMatch changes; -multiple is never enabled. All lifted intervals "
        "are counted, including those outside the width band. Exactly 160 or 240 bp "
        "is a 20% change and is not counted as greater than 20%. minMatch measures "
        "the fraction of source bases remapped, not nucleotide identity. "
        "[UCSC explanation](https://groups.google.com/a/soe.ucsc.edu/g/genome/c/YbyOe7AEXj4).",
        "",
        "## Recommendation and final run",
        "",
    ]
    if final:
        chosen = summarize(final, len(source))
        lines += [
            f"Chosen **MINMATCH_CROSS_SPECIES = {final['min_match']:g}**.",
            "",
            config.MINMATCH_CROSS_SPECIES_RATIONALE,
            "",
            "This is a mouse-pilot starting policy, not a validated optimum for every "
            "species. Width stability does not establish orthology or motif retention; "
            "future classification must retain mapping uncertainty.",
            "",
            f"Final lift: **{chosen['lifted']}/{len(source)} ({chosen['lift_rate']:.2%})**; "
            f"failed: **{chosen['unmapped']}**. Sum equals input: **PASS**.",
            "",
            "| Unmapped reason (verbatim category) | Elements |",
            "| --- | ---: |",
        ]
        reasons = Counter(row["failure_reason"] for row in final["unmapped"].values())
        lines += [f"| {reason} | {count} |" for reason, count in sorted(reasons.items())] or [
            "| none | 0 |"
        ]
        lines += [
            "",
            f"Mapped BED: `{final['mapped_path'].relative_to(config.ROOT)}`.",
            f"Unmapped BED: `{final['unmapped_path'].relative_to(config.ROOT)}`.",
            "All inputs, mapped coordinates, width flags, and failure reasons: "
            "`data/processed/pilot_mouse/final/element_audit.csv`. "
            "No length filtering or sequence/motif classification is performed here.",
        ]
    else:
        lines += [
            "Final run: **NOT RUN**. Review the measured sweep before setting the config value."
        ]
    if error:
        lines += ["", "## FAILURE", "", "```text", error, "```"]
    if observations:
        lines += [
            "",
            "## Provenance and reproduction",
            "",
            f"Chain: [{config.CHAIN_PATHS['mm39'].name}]({CHAIN_URL}); "
            f"local MD5 `{local_md5(config.CHAIN_PATHS['mm39'])}`.",
            f"Pilot BED MD5: `{local_md5(pilot.BED)}`.",
            "Download provenance is appended to data/README.md. Each sweep directory "
            "retains mapped/unmapped BEDs, a log containing the exact command, and "
            "a full per-element audit. `sweep.csv` preserves all measured summaries.",
            "",
            "Run `make pilot-mouse` offline. It rebuilds the hg38 pilot, verifies hg38/mm39, "
            "runs all four thresholds, and repeats the configured threshold for the final outputs. "
            "The final outputs are asserted identical to that sweep entry. No downloads occur.",
        ]
    text = "\n".join(lines) + "\n"
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text)
    print(text, flush=True)
    return text


def main() -> None:
    """Run the small sweep first; --sweep-only deliberately publishes no final choice."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-only", action="store_true")
    args = parser.parse_args()
    # Avoid stale final outputs on a failed dependency check or configuration change.
    final_dir = OUTPUT / "final"
    for name in ("lifted.mm39.bed", "unmapped.hg38.bed", "element_audit.csv", "liftOver.log"):
        (final_dir / name).unlink(missing_ok=True)
    genomes = pilot.verify_step("mouse")
    observations, source, final = [], {}, None
    try:
        if any(row["status"] != "PASS" for row in genomes):
            raise ValueError("Mouse dependencies failed; no lift was run")
        OUTPUT.mkdir(parents=True, exist_ok=True)
        expected = pilot.choose_pilot(pd.read_csv(config.HG38_ELEMENTS, keep_default_na=False))
        expected_path = OUTPUT / "expected_pilot.hg38.bed"
        pilot.write_bed(expected, expected_path)
        if not pilot.BED.is_file() or pilot.BED.read_bytes() != expected_path.read_bytes():
            raise ValueError(
                "Pilot BED is absent/stale or differs from canonical hg38 source; run make pilot"
            )
        source = parse_bed(pilot.BED)
        sizes = pilot.parse_sizes((pilot.METADATA / "mm39.chrom.sizes").read_text())
        runs = {}
        for value in SWEEP:
            directory = OUTPUT / f"minmatch_{value:g}"
            result = lift_intervals(
                pilot.BED, "hg38", "mm39", config.CHAIN_PATHS["mm39"], directory, min_match=value
            )
            validate_target(result, sizes)
            write_audit(source, result, directory / "element_audit.csv")
            observations.append(summarize(result, len(source)))
            runs[value] = result
        pd.DataFrame(observations).to_csv(OUTPUT / "sweep.csv", index=False)
        if not args.sweep_only:
            choice = config.MINMATCH_CROSS_SPECIES
            if choice not in SWEEP or not config.MINMATCH_CROSS_SPECIES_RATIONALE:
                raise ValueError(
                    "Record a measured sweep recommendation and rationale in config first"
                )
            final = lift_intervals(
                pilot.BED, "hg38", "mm39", config.CHAIN_PATHS["mm39"], final_dir, min_match=choice
            )
            validate_target(final, sizes)
            assert final["mapped"] == runs[choice]["mapped"]
            assert final["unmapped"] == runs[choice]["unmapped"]
            write_audit(source, final, final_dir / "element_audit.csv")
    except (ValueError, OSError, RuntimeError, AssertionError) as exc:
        render(genomes, observations, source, error=str(exc))
        raise SystemExit(1) from exc
    render(genomes, observations, source, final)


if __name__ == "__main__":
    main()
