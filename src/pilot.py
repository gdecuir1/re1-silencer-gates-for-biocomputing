"""Verify the genomes needed by each step before publishing a canonical chr8 pilot.

``make pilot`` is offline: it never completes or repairs a genome download.
``python -m src.pilot --fetch-metadata`` may fetch only small UCSC chromosome-size
and checksum files, plus the official twoBitInfo executable if absent. Every
such download uses the existing append-only provenance helper.

Genome verification is a hard gate. A missing, changing, truncated, wrong-build,
or checksum-mismatched REQUIRED reference prevents that step. Unused genomes
are neither read nor repaired.
Published chromosome sizes include every sequence in the .2bit (alternate and
unplaced sequences as well as chromosomes); a primary-assembly-only total would
be the wrong comparator. MD5 is checked in addition to sizes because a readable
index or plausible total alone cannot establish that all DNA bytes arrived.
"""

import argparse
import json
import platform
import random
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

from src import config
from src.downloads import download_small_file, local_md5

ASSEMBLIES = ("hg38", "rheMac10", "mm39", "canFam6", "bosTau9", "monDom5")
BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"
GENOMES = config.RAW_DIR / "genomes"
METADATA = GENOMES / "verification"
LOCAL_INFO = config.RAW_DIR / "tools" / "twoBitInfo"
REPORT = config.TABLES_DIR / "pilot_setup.md"
BED = config.PROCESSED_DIR / f"pilot_{config.PILOT_CHROM}.bed"
EVIDENCE = config.PROCESSED_DIR / "pilot_setup"
SAMPLE_SEED = 20261108
SAMPLE_COUNT = 10


def info_executable() -> str | None:
    """Prefer the activated environment; otherwise use the locally recorded utility."""
    return shutil.which("twoBitInfo") or (str(LOCAL_INFO) if LOCAL_INFO.is_file() else None)


def fetch_metadata(assemblies: tuple[str, ...] = ("hg38",)) -> list[str]:
    """Acquire bounded metadata only; surface every failed fetch without inventing data."""
    errors = []
    for assembly in assemblies:
        for suffix, remote in (
            ("chrom.sizes", f"{assembly}.chrom.sizes"),
            ("md5sum.txt", "md5sum.txt"),
        ):
            try:
                download_small_file(
                    f"{BASE}/{assembly}/bigZips/{remote}",
                    METADATA / f"{assembly}.{suffix}",
                    fetch=True,
                )
            except (ValueError, OSError, RuntimeError) as exc:
                errors.append(f"{assembly}/{remote}: {exc}")
    if not info_executable():
        systems = {
            ("Darwin", "arm64"): "macOSX.arm64",
            ("Darwin", "x86_64"): "macOSX.x86_64",
            ("Linux", "x86_64"): "linux.x86_64",
        }
        system = systems.get((platform.system(), platform.machine()))
        if not system:
            errors.append("No supported UCSC twoBitInfo download for this platform; install it")
        else:
            try:
                download_small_file(
                    f"https://hgdownload.soe.ucsc.edu/admin/exe/{system}/twoBitInfo",
                    LOCAL_INFO,
                    fetch=True,
                    max_bytes=20_000_000,
                )
                LOCAL_INFO.chmod(0o755)
            except (ValueError, OSError, RuntimeError) as exc:
                errors.append(f"twoBitInfo: {exc}")
    return errors


def parse_sizes(text: str) -> dict[str, int]:
    """Require unique names and positive integer lengths; reject HTML/error responses."""
    sizes = {}
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 2 or fields[0] in sizes:
            raise ValueError(f"Invalid or duplicate chromosome-size row: {line!r}")
        length = int(fields[1])
        if length <= 0:
            raise ValueError(f"Nonpositive sequence length: {line!r}")
        sizes[fields[0]] = length
    if not sizes:
        raise ValueError("Empty chromosome-size table")
    return sizes


def verify_genome(assembly: str, executable: str | None) -> dict:
    """Measure one file and compare every contig and its whole-file publisher MD5.

    Continue across assemblies so the report names every failure. Size/mtime are
    checked before and after reading: a file actively being downloaded cannot
    pass using a mixture of observations from different moments.
    """
    path = GENOMES / f"{assembly}.2bit"
    result = dict(
        assembly=assembly,
        file_bytes=None,
        sequence_count=None,
        total_bases=None,
        expected_sequence_count=None,
        expected_total_bases=None,
        observed_md5="",
        expected_md5="",
        status="FAIL",
        errors=[],
        sizes_url=f"{BASE}/{assembly}/bigZips/{assembly}.chrom.sizes",
    )
    expected, observed = None, None
    for suffix in ("chrom.sizes", "md5sum.txt"):
        try:
            metadata_path = METADATA / f"{assembly}.{suffix}"
            remote = f"{assembly}.chrom.sizes" if suffix == "chrom.sizes" else "md5sum.txt"
            # Offline cache reads still check their download provenance/digest.
            download_small_file(f"{BASE}/{assembly}/bigZips/{remote}", metadata_path)
            text = metadata_path.read_text()
            if suffix == "chrom.sizes":
                expected = parse_sizes(text)
                result.update(
                    expected_sequence_count=len(expected),
                    expected_total_bases=sum(expected.values()),
                )
            else:
                matches = [
                    line.split()[0]
                    for line in text.splitlines()
                    if len(line.split()) == 2
                    and line.split()[1].lstrip("*").removeprefix("./") == f"{assembly}.2bit"
                ]
                if len(matches) != 1 or len(matches[0]) != 32:
                    raise ValueError(f"No unique publisher MD5 for {assembly}.2bit")
                result["expected_md5"] = matches[0].lower()
        except (ValueError, OSError, RuntimeError) as exc:
            result["errors"].append(str(exc))
    if not path.is_file():
        result["errors"].append(f"Missing file: {path}")
        return result
    before = path.stat()
    result["file_bytes"] = before.st_size
    if not before.st_size:
        result["errors"].append("File is empty")
        return result
    if not executable:
        result["errors"].append("twoBitInfo is missing; install it or use --fetch-metadata")
    else:
        try:
            command = [executable, str(path), "stdout"]
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            EVIDENCE.mkdir(parents=True, exist_ok=True)
            (EVIDENCE / f"{assembly}.twoBitInfo.txt").write_text(completed.stdout)
            (EVIDENCE / f"{assembly}.twoBitInfo.log").write_text(
                json.dumps(command) + f"\nexit={completed.returncode}\n" + completed.stderr
            )
            if completed.returncode:
                raise RuntimeError(
                    f"twoBitInfo exit {completed.returncode}: {completed.stderr.strip()}"
                )
            observed = parse_sizes(completed.stdout)
            result.update(sequence_count=len(observed), total_bases=sum(observed.values()))
            if expected is not None and observed != expected:
                result["errors"].append("Sequence names/lengths differ from UCSC chrom.sizes")
        except (ValueError, OSError, RuntimeError) as exc:
            result["errors"].append(str(exc))
    result["observed_md5"] = local_md5(path)
    if result["expected_md5"] and result["observed_md5"] != result["expected_md5"]:
        result["errors"].append(
            "Whole-file MD5 differs from UCSC; incomplete or incorrect download"
        )
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        result["errors"].append("File changed during verification; download may still be active")
    if not result["errors"] and observed is not None and expected is not None:
        result["status"] = "PASS"
    return result


# Dependencies name the actual consumers. Week 3 also uses the hg38 source;
# all five targets are required there, but this step is not run by the pilot.
STEP_ASSEMBLIES = {
    "pilot": ("hg38",),
    "mouse": ("hg38", "mm39"),
    "multi_species": ASSEMBLIES,
}


def verify_step(step: str) -> list[dict]:
    """Apply the unchanged integrity checks only to a named step's dependencies.

    Return every required assembly's observation, including failures. Callers must
    stop if any status is not PASS. Unrequired paths are never inspected.
    """
    if step not in STEP_ASSEMBLIES:
        raise ValueError(f"Unknown genome-verification step: {step}")
    return [verify_genome(assembly, info_executable()) for assembly in STEP_ASSEMBLIES[step]]


def choose_pilot(frame: pd.DataFrame) -> pd.DataFrame:
    """Select only exact project='motif' in the working hg38 coordinate columns."""
    required = {"element_id", "project", "working_build", "chrom", "start", "end", "sequence"}
    if missing := required - set(frame):
        raise ValueError(f"Missing columns: {sorted(missing)}")
    if not frame.working_build.eq("hg38").all():
        raise ValueError("Input working_build must be hg38 for every row")
    pilot = frame.loc[frame.project.eq("motif") & frame.chrom.eq(config.PILOT_CHROM)].copy()
    if not 150 <= len(pilot) <= 600:
        alternatives = frame.loc[frame.project.eq("motif")].groupby("chrom").size()
        alternatives = alternatives.loc[alternatives.between(150, 600)]
        ranked = sorted(alternatives.items(), key=lambda item: (abs(item[1] - 350), item[0]))
        suggestion = (
            f"Consider {ranked[0][0]} ({ranked[0][1]} canonical elements)"
            if ranked
            else "No chromosome meets 150–600"
        )
        raise ValueError(
            f"Observed {len(pilot)} canonical {config.PILOT_CHROM} elements; expected roughly "
            f"300–400, allowed 150–600. {suggestion}; no alternative selected automatically"
        )
    if pilot.element_id.isna().any() or pilot.element_id.duplicated().any():
        raise ValueError("Pilot IDs must be present and unique")
    if not pilot.element_id.str.fullmatch(r"[^\s]+").all():
        raise ValueError("Pilot IDs must be valid whitespace-free BED/FASTA names")
    for column in ("start", "end"):
        values = pd.to_numeric(pilot[column], errors="raise")
        if not (np.isfinite(values).all() and values.eq(np.floor(values)).all()):
            raise ValueError(f"Invalid integer coordinates in {column}; refusing to round")
        pilot[column] = values.astype("int64")
    widths = pilot.end - pilot.start
    invalid = widths.ne(200) | pilot.start.lt(0)
    if invalid.any():
        raise ValueError(
            "Expected 200 bp nonnegative intervals; invalid records:\n"
            + pilot.loc[invalid, ["element_id", "start", "end"]]
            .assign(width=widths[invalid])
            .to_csv(index=False)
        )
    assert widths.eq(config.ELEMENT_LENGTH).all(), "Pilot interval widths must equal 200"
    if "strand" in pilot and not pilot.strand.isin(["+", "-"]).all():
        raise ValueError("Pilot strand must be '+' or '-'")
    if not pilot.sequence.str.fullmatch("[ACGT]{200}").all():
        raise ValueError("Expected uppercase, sequence-verified 200 bp library cores")
    return pilot.sort_values(["start", "end", "element_id"], kind="stable").reset_index(drop=True)


def write_bed(frame: pd.DataFrame, path: Path) -> None:
    """Write BED6 if strand exists, otherwise BED4; never re-center an interval."""
    assert (frame.end - frame.start).eq(200).all(), "BED requires exactly 200 bp"
    with path.open("w") as handle:
        for row in frame.itertuples(index=False):
            suffix = f"\t0\t{row.strand}" if "strand" in frame else ""
            handle.write(f"{row.chrom}\t{row.start}\t{row.end}\t{row.element_id}{suffix}\n")


def regression(pilot: pd.DataFrame, candidate: Path) -> list[dict]:
    """Extract ten seeded random BED rows and require oriented, exact sequence identity.

    Sample from the WRITTEN BED, not a separately reconstructed set of coordinates.
    twoBitToFa respects BED6 strand; do not accept a reverse-complement rescue here,
    since that could conceal a strand-writing error. Preserve sample IDs and DNA.
    """
    executable = shutil.which("twoBitToFa")
    if not executable:
        raise FileNotFoundError("twoBitToFa is not on PATH")
    lines = candidate.read_text().splitlines()
    assert len(lines) == len(pilot), "Written BED count differs from selected pilot"
    indices = sorted(random.Random(SAMPLE_SEED).sample(range(len(lines)), SAMPLE_COUNT))
    sampled = EVIDENCE / "sample.bed"
    sampled.write_text("\n".join(lines[i] for i in indices) + "\n")
    fasta = EVIDENCE / "sample.fa"
    command = [executable, f"-bed={sampled}", "-noMask", str(config.HG38_TWOBIT), str(fasta)]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    (EVIDENCE / "twoBitToFa.log").write_text(
        json.dumps(command) + "\n" + result.stdout + result.stderr
    )
    if result.returncode:
        raise RuntimeError(f"twoBitToFa exit {result.returncode}: {result.stderr.strip()}")
    records = list(SeqIO.parse(fasta, "fasta"))
    sequences = {record.id: str(record.seq) for record in records}
    sample = pilot.iloc[indices]
    if len(records) != SAMPLE_COUNT or set(sequences) != set(sample.element_id):
        raise ValueError("Extracted FASTA IDs/count do not match the ten selected BED records")
    checks = [
        dict(
            element_id=row.element_id,
            chrom=row.chrom,
            start=row.start,
            end=row.end,
            strand=getattr(row, "strand", "+"),
            length=len(sequences[row.element_id]),
            exact=sequences[row.element_id] == row.sequence,
        )
        for row in sample.itertuples(index=False)
    ]
    pd.DataFrame(checks).to_csv(EVIDENCE / "regression.csv", index=False)
    return checks


def report(genomes: list[dict], *, pilot=None, checks=None, errors=()) -> str:
    """Write a deterministic report with missing observations explicitly marked unmeasured."""

    def show(value):
        return "not measured" if value is None else f"{value:,}"

    lines = [
        "# Week 2 pilot setup",
        "",
        "## Genome gate",
        "",
        "Required assembly for pilot selection, BED writing, and regression: hg38 only. "
        "The mouse lift separately requires hg38 and mm39; the full multi-species step "
        "requires hg38 and all five targets. Unrequired genomes are not inspected. "
        "Counts and totals "
        "include all sequences in UCSC's published .chrom.sizes, including alternate/unplaced "
        "sequences. Exact per-contig equality and whole-file publisher MD5 are also required.",
        "",
        "| Assembly | File bytes | Observed sequences | UCSC sequences | "
        "Observed bases | UCSC bases | Result |",
        "| --- | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in genomes:
        fields = [
            row["assembly"],
            *[
                show(row[key])
                for key in (
                    "file_bytes",
                    "sequence_count",
                    "expected_sequence_count",
                    "total_bases",
                    "expected_total_bases",
                )
            ],
            row["status"],
        ]
        lines.append("| " + " | ".join(fields) + " |")
    lines += ["", "## Integrity details", ""]
    for row in genomes:
        lines += [
            f"### {row['assembly']}",
            "",
            f"Published sizes: {row['sizes_url']}",
            f"Observed MD5: `{row['observed_md5'] or 'not measured'}`; "
            f"expected: `{row['expected_md5'] or 'unavailable'}`.",
            "",
        ]
        lines += [f"- FAIL: {error}" for error in row["errors"]] or [
            "- PASS: nonzero, stable file; contigs and publisher MD5 match."
        ]
        lines.append("")
    lines += [
        "## Pilot and regression",
        "",
        f"Configured chromosome: **{config.PILOT_CHROM}**. "
        "Chosen by the project request for a single-chromosome pilot before scaling; "
        "only canonical project='motif' elements qualify. Contrast elements are excluded.",
        "",
    ]
    if pilot is None:
        lines += [
            "Pilot selection: **NOT RUN**. Pilot BED: **not published**. "
            "Ten-element regression: **NOT RUN**. See failures above/below."
        ]
    else:
        lines += [
            f"Pilot count: **{len(pilot)}**; expected roughly 300–400; allowed 150–600.",
            "Width distribution: "
            f"{(pilot.end - pilot.start).value_counts().sort_index().to_dict()}; expected 200 bp.",
            f"Random sample: {SAMPLE_COUNT} elements, fixed seed {SAMPLE_SEED}; "
            "exact oriented sequence equality required, with -noMask.",
            "",
            "| Element | hg38 interval | Strand | Extracted bp | Exact match |",
            "| --- | --- | --- | ---: | --- |",
        ]
        for row in checks or []:
            lines.append(
                f"| {row['element_id']} | {row['chrom']}:{row['start']}-{row['end']} | "
                f"{row['strand']} | {row['length']} | {'PASS' if row['exact'] else 'FAIL'} |"
            )
        passed = sum(row["exact"] for row in checks or [])
        lines += [
            "",
            f"Regression: **{passed}/{SAMPLE_COUNT}** exact matches. "
            f"Pilot BED: **{'published' if not errors else 'not published'}**.",
        ]
    if errors:
        lines += ["", "## Blocking errors", "", *[f"- {error}" for error in errors]]
    lines += [
        "",
        "## Reproduction",
        "",
        "Run `make pilot` offline with a verified hg38 genome. "
        "Use `python -m src.pilot --fetch-metadata` only to acquire missing small UCSC metadata "
        "and twoBitInfo. This pipeline never downloads/repairs a genome or changes hg19 data. "
        "Fixed sample seed, stable sorting, and absence of report timestamps make output "
        "identical for identical inputs/tools. "
        "An actively changing download is not a stable input.",
        "",
    ]
    text = "\n".join(lines)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(text)
    print(text, flush=True)
    return text


def main() -> None:
    """Verify first, fail closed, and publish the pilot only after ten exact matches."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch-metadata", action="store_true")
    args = parser.parse_args()
    BED.unlink(missing_ok=True)  # Never leave an old consumable pilot after a failed rerun.
    errors = fetch_metadata(STEP_ASSEMBLIES["pilot"]) if args.fetch_metadata else []
    genomes = verify_step("pilot")
    if errors or any(row["status"] != "PASS" for row in genomes):
        report(genomes, errors=errors)
        raise SystemExit(1)
    pilot, checks = None, None
    candidate = EVIDENCE / "candidate.bed"
    try:
        frame = pd.read_csv(config.HG38_ELEMENTS, keep_default_na=False)
        pilot = choose_pilot(frame)
        EVIDENCE.mkdir(parents=True, exist_ok=True)
        write_bed(pilot, candidate)
        checks = regression(pilot, candidate)
        if not all(row["exact"] for row in checks):
            raise ValueError(
                "Pilot BED regression failed: inspect coordinates/strand; do not change the lift"
            )
        BED.parent.mkdir(parents=True, exist_ok=True)
        candidate.replace(BED)
    except (OSError, ValueError, RuntimeError, AssertionError) as exc:
        errors.append(str(exc))
    finally:
        candidate.unlink(missing_ok=True)
    report(genomes, pilot=pilot, checks=checks, errors=errors)
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
