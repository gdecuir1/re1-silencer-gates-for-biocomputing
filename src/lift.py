"""Lift native MPRA elements once, then admit only exact hg38 sequence matches.

The default input is the native motif/noMotif population from the four GEO cell
files, joined to the published WG_REST oligos in Table S2. Parenthesized Table S2
names are explicit semicolon-separated aliases, not coordinates to be guessed.
No activity is aggregated and no reference sequence is substituted for an oligo.

All coordinates are zero-based, half-open. liftOver runs at its installed default
-minMatch (0.95), without -multiple. Both output files, the input BED, extracted
FASTA, complete per-element audit, and tool diagnostics remain on disk. Only the
verified output is intended for downstream hg38 analysis. This module is offline.
"""

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from src import config
from src.assembly_check import GEO_ARCHIVE, TABLE_S2, TABLE_S2_MEMBER, load_evidence
from src.downloads import local_md5

AUDIT_DIR = config.PROCESSED_DIR / "liftover_hg19_to_hg38"
REPORT = config.TABLES_DIR / "liftover_hg19_to_hg38.md"
BUCKETS = ("lifted+verified", "lifted+mismatch", "failed-to-lift")
MISMATCH_TYPES = (
    "length_differs",
    "one_base_difference",
    "multiple_base_differences",
    "high_divergence_or_unrelated",
    "ambiguous_bases",
    "adapter_unresolved",
)
# Operational reporting/stop rules, NOT biological thresholds or liftOver flags.
# Fix them before examining outcomes. Never tune them to recover more elements.
MIN_SUCCESS_RATE = 0.95
LARGE_BIAS_PP = 5.0


def published_elements() -> pd.DataFrame:
    """Resolve every native element to exactly one explicitly named library oligo.

    load_evidence already checks agreement of native metadata across all four
    cell types. Its exact-name join intentionally omits compound aliases; here
    we expand only the publisher's '(ID;ID;...)' syntax. Ambiguity is fatal even
    if two conflicting rows happen to have identical DNA.
    """
    elements, _ = load_evidence()
    with zipfile.ZipFile(TABLE_S2) as archive:
        oligos = pd.read_csv(archive.open(TABLE_S2_MEMBER), sep="\t")
    aliases = {}
    for row in oligos.loc[oligos.library.eq("WG_REST")].itertuples(index=False):
        name = row.oligo_name
        tokens = name[1:-1].split(";") if name.startswith("(") and name.endswith(")") else [name]
        for token in tokens:
            aliases.setdefault(token, []).append((row.sequence, name))
    records = []
    for row in elements.itertuples(index=False):
        matches = aliases.get(row.silencer, [])
        if len(matches) != 1:
            raise ValueError(
                f"{row.silencer}: expected one explicit oligo alias, got {len(matches)}"
            )
        sequence, name = matches[0]
        records.append(
            dict(
                element_id=row.silencer,
                chrom=f"chr{row.chromosome}",
                start=row.start,
                end=row.end,
                project=row.project,
                sequence=sequence,
                source_oligo_name=name,
                sequence_join="exact_name" if name == row.silencer else "explicit_alias",
                source_build=config.SOURCE_BUILD,
            )
        )
    return pd.DataFrame(records)


def validate_input(frame: pd.DataFrame) -> pd.DataFrame:
    """Reject invalid IDs, ambiguous builds, missing DNA, or non-200 bp intervals.

    Invalid widths are listed with their IDs in the error, never centered, padded,
    rounded, or trimmed. The CLI copies this diagnostic to its failure report.
    A missing strand becomes '+' solely to track mapping orientation; it is not
    evidence that the library oligo was synthesized on the genomic plus strand.
    """
    frame = frame.copy()
    required = {"element_id", "chrom", "start", "end", "sequence", "source_build"}
    if missing := required - set(frame.columns):
        raise ValueError(f"Missing required input columns: {sorted(missing)}")
    if frame.empty or frame[list(required)].isna().any().any():
        raise ValueError("Input must be nonempty with all required metadata and library sequences")
    if set(frame.source_build) != {config.SOURCE_BUILD}:
        raise ValueError(f"Expected source_build=hg19; observed {sorted(set(frame.source_build))}")
    if frame.element_id.duplicated().any() or not frame.element_id.str.fullmatch(r"[^\s]+\Z").all():
        raise ValueError("Element IDs must be unique, nonempty, whitespace-free BED/FASTA names")
    if not frame.chrom.str.fullmatch(r"chr[^\s]+").all():
        raise ValueError("Expected UCSC chromosome names (chr1, chrX, ...)")
    for col in ("start", "end"):
        numeric = pd.to_numeric(frame[col], errors="raise")
        if not (np.isfinite(numeric).all() and numeric.eq(np.floor(numeric)).all()):
            raise ValueError(f"{col} contains non-integer coordinates; refusing to round")
        frame[col] = numeric.astype("int64")
    widths = frame.end - frame.start
    invalid = widths.ne(config.ELEMENT_LENGTH) | frame.start.lt(0)
    if invalid.any():
        observed = frame.loc[invalid, ["element_id", "start", "end"]].assign(width=widths[invalid])
        raise ValueError(
            "Expected nonnegative 200 bp intervals; invalid records:\n"
            + observed.to_csv(index=False)
        )
    frame["strand_supplied"] = "strand" in frame.columns
    if "strand" not in frame:
        frame["strand"] = "+"
    if not frame.strand.isin(["+", "-"]).all():
        raise ValueError("Supplied strand must be '+' or '-' for every element")
    for col, default in (
        ("project", "unspecified"),
        ("source_oligo_name", ""),
        ("sequence_join", "provided_input"),
    ):
        if col not in frame:
            frame[col] = default
    frame["sequence"] = frame.sequence.str.upper()
    if not frame.sequence.str.fullmatch("[ACGTN]+").all():
        raise ValueError(
            "Library sequences must be DNA (ACGTN); whitespace is not silently removed"
        )
    return frame.sort_values("element_id", kind="stable").reset_index(drop=True)


def prepare_cores(frame: pd.DataFrame, adapter5: str = "", adapter3: str = "") -> pd.DataFrame:
    """Strip only explicit or unambiguously shared flanks, independently of hg38.

    For automatic inference, group long oligos by length and require at least two
    distinct sequences. Enumerate all ways to leave 200 bp; every nonempty flank
    must be at least eight bases and identical across the group. Exactly one
    split must qualify. Ambiguous/singleton long oligos remain untrimmed and fail
    verification. Explicit --adapter-5p/--adapter-3p can resolve a known library
    design. Never search the reference for a convenient matching subwindow.
    """
    frame = frame.copy()
    adapter5, adapter3 = adapter5.upper(), adapter3.upper()
    if any(base not in "ACGT" for base in adapter5 + adapter3):
        raise ValueError("Adapters must contain only ACGT")
    inferred = {}
    if not (adapter5 or adapter3):
        for length, group in frame.groupby(frame.sequence.str.len()):
            seqs = list(group.sequence.unique())
            if length <= 200 or len(seqs) < 2:
                continue
            candidates = []
            extra = length - 200
            for left in range(extra + 1):
                right = extra - left
                if (0 < left < 8) or (0 < right < 8):
                    continue
                prefixes = {seq[:left] for seq in seqs}
                suffixes = {seq[-right:] if right else "" for seq in seqs}
                if len(prefixes) == len(suffixes) == 1:
                    candidates.append((next(iter(prefixes)), next(iter(suffixes))))
            if len(candidates) == 1:
                inferred[length] = candidates[0]
    cores, prefixes, suffixes, statuses = [], [], [], []
    for sequence in frame.sequence:
        prefix, suffix, status = "", "", "none_present"
        if len(sequence) > 200:
            prefix, suffix = (
                (adapter5, adapter3)
                if adapter5 or adapter3
                else inferred.get(len(sequence), ("", ""))
            )
            valid = (
                bool(prefix or suffix)
                and sequence.startswith(prefix)
                and sequence.endswith(suffix)
                and len(sequence) - len(prefix) - len(suffix) == 200
            )
            if valid:
                status = "explicit_adapters" if adapter5 or adapter3 else "shared_flanks_inferred"
            else:
                prefix, suffix, status = "", "", "unresolved_long_oligo"
        stop = len(sequence) - len(suffix)
        cores.append(sequence[len(prefix) : stop])
        prefixes.append(prefix)
        suffixes.append(suffix)
        statuses.append(status)
    frame["library_sequence"] = frame.sequence
    frame["sequence"] = cores
    frame["adapter_5p"] = prefixes
    frame["adapter_3p"] = suffixes
    frame["adapter_status"] = statuses
    return frame


def run_tool(argv: list[str], log: Path) -> None:
    """Keep stdout/stderr and fail explicitly on a nonzero tool exit status."""
    result = subprocess.run(argv, capture_output=True, text=True, check=False)
    log.write_text(
        "ARGV: " + json.dumps(argv) + "\nSTDOUT:\n" + result.stdout + "\nSTDERR:\n" + result.stderr
    )
    if result.returncode:
        raise RuntimeError(f"{argv[0]} exited {result.returncode}: {result.stderr.strip()}")


def parse_bed(path: Path, *, unmapped: bool = False) -> dict[str, dict]:
    """Parse BED6 and attach each liftOver '#reason' to its following record.

    Duplicate mappings, missing reasons, and malformed records are errors rather
    than silent deduplication. Reasons are retained verbatim except the '#' prefix.
    """
    records, reason = {}, None
    for line in path.read_text().splitlines():
        if not line:
            continue
        if line.startswith("#"):
            reason = line[1:].strip()
            continue
        fields = line.split("\t")
        if len(fields) != 6 or fields[5] not in {"+", "-"}:
            raise ValueError(f"Malformed BED6 in {path}: {line}")
        chrom, start, end, name, score, strand = fields
        if name in records:
            raise ValueError(f"Multiple output records for {name} in {path}")
        if unmapped and not reason:
            raise ValueError(f"Unmapped record lacks a liftOver failure reason: {name}")
        records[name] = dict(
            chrom=chrom,
            start=int(start),
            end=int(end),
            strand=strand,
            failure_reason=reason if unmapped else "",
        )
        reason = None
    return records


def lift_intervals(
    source_bed: Path,
    source_build: str,
    target_build: str,
    chain: Path,
    output_dir: Path,
    *,
    min_match: float | None = None,
    executable: str = "liftOver",
) -> dict:
    """Map BED6 intervals between declared builds, preserving both output files.

    This mechanical wrapper does not impose same-species success-rate or exact-DNA
    requirements on other species. The hg19->hg38 caller retains those stricter
    checks. min_match=None omits the flag, preserving liftOver's default; an
    explicit float controls only minMatch, never -multiple or any other option.
    The caller supplies the directional chain and verifies its reference genomes.

    Every input ID must appear once in exactly one output; duplicate/missing IDs,
    changed unmapped records, and invalid target intervals fail loudly. An unmapped
    record is returned with the tool's reason, not discarded. Output paths encode
    source/target builds; adding a target only requires a new call/chain entry.
    """
    if any(not re.fullmatch(r"[A-Za-z0-9_]+", build) for build in (source_build, target_build)):
        raise ValueError("Build names must be nonempty alphanumeric identifiers")
    if source_build == target_build:
        raise ValueError("Source and target builds must differ")
    if min_match is not None and not 0 < min_match <= 1:
        raise ValueError("min_match must be in (0, 1]")
    tool = shutil.which(executable)
    if not tool:
        raise FileNotFoundError(f"Required executable is not on PATH: {executable}")
    if not chain.is_file():
        raise FileNotFoundError(f"Missing directional chain: {chain}")
    original = parse_bed(source_bed)
    if not original:
        raise ValueError("Cannot lift an empty BED")
    for name, record in original.items():
        if record["start"] < 0 or record["end"] <= record["start"]:
            raise ValueError(f"Invalid source interval: {name}: {record}")
    output_dir.mkdir(parents=True, exist_ok=True)
    mapped_path = output_dir / f"lifted.{target_build}.bed"
    unmapped_path = output_dir / f"unmapped.{source_build}.bed"
    options = [] if min_match is None else [f"-minMatch={min_match:g}"]
    command = [tool, *options, str(source_bed), str(chain), str(mapped_path), str(unmapped_path)]
    run_tool(command, output_dir / "liftOver.log")
    mapped, unmapped = parse_bed(mapped_path), parse_bed(unmapped_path, unmapped=True)
    if set(mapped) & set(unmapped) or set(mapped) | set(unmapped) != set(original):
        raise AssertionError("Mapped and unmapped IDs must partition the input exactly")
    for name, record in unmapped.items():
        if any(record[key] != original[name][key] for key in ("chrom", "start", "end", "strand")):
            raise AssertionError(f"Unmapped source interval changed: {name}")
    for name, record in mapped.items():
        if record["start"] < 0 or record["end"] <= record["start"]:
            raise ValueError(f"Invalid target interval: {name}: {record}")
    return dict(
        source_build=source_build,
        target_build=target_build,
        min_match=min_match,
        mapped=mapped,
        unmapped=unmapped,
        mapped_path=mapped_path,
        unmapped_path=unmapped_path,
        command=command,
    )


def verify_sequence(core: str, reference: str) -> dict:
    """Compare every base in both orientations; classify failures diagnostically.

    The minimum Hamming distance is computed only for equal lengths. A difference
    of >=25% is labeled high_divergence_or_unrelated; this is an explicit reporting
    convention, NOT a homology test. Distances of 2..49/200 get their own category.
    N bases never establish a verified match, even when both strings contain N.
    """
    core, reference = core.upper(), reference.upper()
    result = dict(
        verified=False,
        mismatch_type="",
        best_orientation="",
        differences="",
        differing_positions_1based="",
        extracted_length=len(reference),
    )
    if len(core) != 200 or len(reference) != 200:
        return result | {"mismatch_type": "length_differs"}
    reverse = str(Seq(reference).reverse_complement())
    distances = [
        sum(a != b for a, b in zip(core, seq, strict=True)) for seq in (reference, reverse)
    ]
    index = int(distances[1] < distances[0])
    distance = distances[index]
    result.update(
        best_orientation="+" if index == 0 else "-",
        differences=distance,
        differing_positions_1based=",".join(
            str(i)
            for i, (a, b) in enumerate(zip(core, (reference, reverse)[index], strict=True), start=1)
            if a != b
        ),
    )
    if any(base not in "ACGT" for base in core + reference):
        result["mismatch_type"] = "ambiguous_bases"
    elif distance == 0:
        result["verified"] = True
    elif distance == 1:
        result["mismatch_type"] = "one_base_difference"
    elif distance < 50:
        result["mismatch_type"] = "multiple_base_differences"
    else:
        result["mismatch_type"] = "high_divergence_or_unrelated"
    return result


def repeat_overlap(elements: pd.DataFrame, path: Path) -> tuple[np.ndarray, int]:
    """Intersect source intervals with hg19 rmsk using >=1 bp, ignoring strand.

    UCSC rmsk.sql defines columns 5/6/7 (zero-based) as genoName/genoStart/genoEnd.
    Scan the entire gzip in chunks. Within each chromosome/chunk, sorting repeat
    starts and taking cumulative maximum ends handles nested/overlapping repeats.
    Binary search asks whether ANY repeat begins before an element ends AND ends
    after it starts. Endpoint touching is not overlap. Each input contributes once
    to the inside/outside denominator regardless of how many repeats overlap it.
    """
    with gzip.open(path, "rt") as handle:
        if len(handle.readline().rstrip("\n").split("\t")) != 17:
            raise ValueError("Expected the 17-column UCSC hg19 rmsk table")
    overlap = np.zeros(len(elements), dtype=bool)
    groups = {chrom: np.asarray(index) for chrom, index in elements.groupby("chrom").groups.items()}
    total = 0
    for chunk in pd.read_csv(
        path, sep="\t", header=None, usecols=[5, 6, 7], names=list(range(17)), chunksize=200_000
    ):
        total += len(chunk)
        for chrom, repeats in chunk.groupby(5):
            if chrom not in groups:
                continue
            intervals = repeats.sort_values(6)
            starts = intervals[6].to_numpy(dtype=np.int64)
            ends = intervals[7].to_numpy(dtype=np.int64)
            if np.any(starts < 0) or np.any(ends <= starts):
                raise ValueError(f"Invalid RepeatMasker interval on {chrom}")
            ends = np.maximum.accumulate(ends)
            indices = groups[chrom]
            query = elements.loc[indices]
            positions = np.searchsorted(starts, query.end.to_numpy(), side="left") - 1
            hit = (positions >= 0) & (ends[np.maximum(positions, 0)] > query.start.to_numpy())
            overlap[indices] |= hit
    if total == 0:
        raise ValueError("RepeatMasker track is empty; outside-repeat labels would be invalid")
    return overlap, total


def bias_rows(audit: pd.DataFrame) -> tuple[list[str], str]:
    """Report failures and all exclusion separately; do not equate small with absent bias."""
    lines = [
        "| Population | Repeat overlap | Input | Failed lift | Failure rate | Total excluded |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ]
    rates = {}
    for population, frame in [("all", audit), *list(audit.groupby("project", sort=True))]:
        for inside in (True, False):
            subset = frame.loc[frame.in_repeat.eq(inside)]
            failures = int(subset.bucket.eq("failed-to-lift").sum())
            excluded = int(subset.bucket.ne("lifted+verified").sum())
            rate = failures / len(subset) if len(subset) else None
            label = "inside" if inside else "outside"
            display = f"{100 * rate:.4f}%" if rate is not None else "not estimable"
            lines.append(
                f"| {population} | {label} | {len(subset)} | {failures} | {display} | {excluded} |"
            )
            if population == "all":
                rates[inside] = rate
    if any(rate is None for rate in rates.values()):
        conclusion = "Repeat-bias difference is not estimable: one denominator is zero."
    else:
        inside, outside = rates[True], rates[False]
        difference = 100 * (inside - outside)
        ratio = f"{inside / outside:.3f}" if outside else "undefined (outside failure rate is zero)"
        magnitude = "LARGE" if abs(difference) >= LARGE_BIAS_PP else "not large"
        conclusion = (
            f"Inside minus outside failure rate: {difference:+.4f} percentage points; "
            f"risk ratio: {ratio}. The absolute difference is **{magnitude}** under the "
            f"predeclared descriptive rule of >= {LARGE_BIAS_PP:g} percentage points. "
            "This is not evidence that evolutionary-class bias is absent: repeats are "
            "a proxy, not a primate-specific label; inspect both rates and exclusions."
        )
    return lines, conclusion


def render_report(
    audit: pd.DataFrame, repeat_rows: int, fingerprints: dict, commands: list[list[str]]
) -> str:
    """Render deterministic Markdown: no wall-clock timestamp or runtime duration."""
    n = len(audit)
    counts = audit.bucket.value_counts()
    lifted_lengths = (
        audit.loc[audit.bucket.ne("failed-to-lift"), "extracted_length"]
        .value_counts()
        .sort_index()
        .to_dict()
    )
    lines = [
        "# hg19 to hg38 MPRA conversion and sequence verification",
        "",
        f"**Verified hg38 elements: {counts.get('lifted+verified', 0):,} / {n:,}.**",
        "",
        "## Input, coordinate conventions, and sequence provenance",
        "",
        "Source hg19; working build hg38. Coordinates are zero-based, half-open. "
        "The native motif and noMotif sets are kept separately; enhancer contexts and "
        "cell types share element coordinates and are not duplicate element records.",
        "",
        "Input interval lengths, expected 200: "
        f"{(audit.hg19_end - audit.hg19_start).value_counts().sort_index().to_dict()}.",
        "Library sequence lengths: "
        f"{audit.library_sequence.str.len().value_counts().sort_index().to_dict()}.",
        "Core lengths after trimming: "
        f"{audit.sequence.str.len().value_counts().sort_index().to_dict()}.",
        f"Sequence joins: {audit.sequence_join.value_counts().sort_index().to_dict()}.",
        "Explicit aliases are tokens inside Table S2's parenthesized, semicolon-separated "
        "oligo names; no coordinate-based guessing or reference-based sequence imputation.",
        "",
        "## Exhaustive accounting",
        "",
        "| Population | Input | lifted+verified | lifted+mismatch | failed-to-lift | Sum check |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for project, subset in [("all", audit), *list(audit.groupby("project", sort=True))]:
        observed = [int(subset.bucket.eq(bucket).sum()) for bucket in BUCKETS]
        assert sum(observed) == len(subset), "Bucket counts do not sum to input"
        lines.append(
            f"| {project} | {len(subset)} | " + " | ".join(map(str, observed)) + " | PASS |"
        )
    lines += [
        "",
        "Every input ID occurs exactly once in the audit. Only lifted+verified IDs "
        "are exported to data/processed/elements_hg38.csv.",
        "",
        "## Mapping and sequence checks",
        "",
        "Lifted interval length distribution (including failed sequence checks): "
        f"{lifted_lengths}. "
        "Expected for verification: 200 bp.",
        "Verified comparison orientations relative to extracted DNA: "
        f"{audit.loc[audit.verified, 'best_orientation'].value_counts().sort_index().to_dict()}.",
        "liftOver used its default minMatch=0.95; no -minMatch override or -multiple. "
        f"Stop guard: mapping success below {MIN_SUCCESS_RATE:.0%} stops publication. "
        f"Observed mapping success: {audit.bucket.ne('failed-to-lift').mean():.4%}.",
        "twoBitToFa used -noMask and BED6 strand. Library/core DNA was uppercased. "
        "Both extracted orientations were compared by exact full-length equality; "
        "N bases cannot pass. Exported sequence is the library core; exported strand "
        "describes that sequence relative to hg38, after accounting for extraction orientation.",
        "",
        "| Mismatch type | Count |",
        "| --- | ---: |",
    ]
    for kind in MISMATCH_TYPES:
        lines.append(f"| {kind} | {int(audit.mismatch_type.eq(kind).sum())} |")
    lines += [
        "",
        "Equal-length failures: one base, 2–49 bases, or >=50/200 bases "
        "(high_divergence_or_unrelated). The last label is descriptive, not proof of "
        "unrelated ancestry. Length changes are never repaired. Detailed distances, "
        "positions, extracted DNA, and original coordinates are in the audit CSV.",
        "",
        "| liftOver failure reason | Count |",
        "| --- | ---: |",
    ]
    reasons = audit.loc[audit.bucket.eq("failed-to-lift"), "failure_reason"].value_counts()
    lines += [f"| {reason} | {count} |" for reason, count in reasons.items()] or ["| none | 0 |"]
    lines += [
        "",
        "## Adapter handling",
        "",
        "Adapter inference uses shared library flanks independently of hg38. Ambiguous "
        "flanks are not trimmed. Exactly 200 bp oligos require no trimming.",
        "",
        "| Status | 5-prime stripped | 3-prime stripped | Count |",
        "| --- | --- | --- | ---: |",
    ]
    for (status, left, right), count in (
        audit.groupby(["adapter_status", "adapter_5p", "adapter_3p"], dropna=False).size().items()
    ):
        lines.append(f"| {status} | {left or '(none)'} | {right or '(none)'} | {count} |")
    repeat_table, conclusion = bias_rows(audit)
    lines += [
        "",
        "## Repeat-associated lift failure",
        "",
        f"Scanned {repeat_rows:,} hg19 RepeatMasker rows. Inside means >=1 bp overlap "
        "between the original hg19 element and any RepeatMasker annotation, any strand. "
        "All classes (including simple/low-complexity repeats) count. Denominator: all "
        "input elements, including failures. Each element is counted once.",
        "",
        *repeat_table,
        "",
        conclusion,
        "",
        "## Reproduction and input identity",
        "",
        "Run `make lift` with liftOver and twoBitToFa on PATH. This is offline. "
        "The report, BEDs, FASTA, audit CSV, and verified CSV are deterministic for "
        "identical inputs and executables. Download timestamps live only in data/README.md.",
        "",
        "| Input/tool | MD5 |",
        "| --- | --- |",
    ]
    lines += [f"| {name} | `{digest}` |" for name, digest in sorted(fingerprints.items())]
    lines += [
        "",
        "Commands (argv; no shell interpretation):",
        "",
        "```text",
        *[json.dumps(command) for command in commands],
        "```",
        "",
        "[UCSC chain source](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) · "
        "[RepeatMasker schema]"
        "(https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.sql) · "
        "[hg38 publisher checksums](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt)",
        "",
    ]
    return "\n".join(lines)


def convert(
    elements: pd.DataFrame,
    chain: Path,
    genome: Path,
    repeats: Path,
    output_dir: Path,
    verified_path: Path,
    report_path: Path,
    *,
    liftover: str = "liftOver",
    twobit: str = "twoBitToFa",
    adapter5: str = "",
    adapter3: str = "",
    extra_inputs: dict | None = None,
) -> pd.DataFrame:
    """Run the complete conversion and return the exhaustive per-element audit.

    Output files use stable names and sorted input IDs. A previous verified CSV
    is removed before a new attempt so a failing rerun cannot leave a stale file
    looking like the result of that run. The new CSV is published only after all
    accounting, sequence, repeat, and report-generation checks succeed.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    verified_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    verified_path.unlink(missing_ok=True)
    elements = prepare_cores(validate_input(elements), adapter5, adapter3)
    executables = []
    for executable in (liftover, twobit):
        resolved = shutil.which(executable)
        if not resolved:
            raise FileNotFoundError(f"Required executable is not on PATH: {executable}")
        executables.append(resolved)
    fingerprints = dict(extra_inputs or {})
    for label, path in (
        ("hg19_to_hg38_chain", chain),
        ("hg38_2bit", genome),
        ("hg19_RepeatMasker", repeats),
        ("liftOver_executable", Path(executables[0])),
        ("twoBitToFa_executable", Path(executables[1])),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")
        fingerprints[label] = local_md5(path)
    if genome.resolve() == config.HG38_TWOBIT.resolve():
        manifest = genome.parent / "hg38.md5sum.txt"
        checksums = [
            line.split()[0]
            for line in manifest.read_text().splitlines()
            if len(line.split()) == 2 and line.split()[1].lstrip("*") == "hg38.2bit"
        ]
        if len(checksums) != 1 or checksums[0] != fingerprints["hg38_2bit"]:
            raise ValueError("hg38.2bit does not match its required publisher MD5 manifest")
    source_bed, mapped_bed, unmapped_bed = [
        output_dir / name for name in ("input.hg19.bed", "lifted.hg38.bed", "unmapped.hg19.bed")
    ]
    elements.to_csv(output_dir / "input_elements.csv", index=False)
    with source_bed.open("w") as handle:
        for row in elements.itertuples(index=False):
            handle.write(
                f"{row.chrom}\t{row.start}\t{row.end}\t{row.element_id}\t0\t{row.strand}\n"
            )
    conversion = lift_intervals(
        source_bed,
        config.SOURCE_BUILD,
        config.WORKING_BUILD,
        chain,
        output_dir,
        executable=executables[0],
    )
    commands = [conversion["command"]]
    mapped, failed = conversion["mapped"], conversion["unmapped"]
    input_ids = set(elements.element_id)
    if set(mapped) & set(failed) or set(mapped) | set(failed) != input_ids:
        raise AssertionError("Mapped and unmapped IDs must partition the input exactly")
    # An unmapped record is biological evidence. Check that it still describes the
    # original interval, then retain its complete BED record and reason in the audit.
    for row in elements.itertuples(index=False):
        if row.element_id in failed:
            record = failed[row.element_id]
            if (record["chrom"], record["start"], record["end"], record["strand"]) != (
                row.chrom,
                row.start,
                row.end,
                row.strand,
            ):
                raise AssertionError(f"Unmapped interval changed: {row.element_id}")
    rate = len(mapped) / len(elements)
    print(
        f"liftOver: input={len(elements)}, mapped={len(mapped)}, unmapped={len(failed)}, "
        f"success={rate:.4%}; stop guard={MIN_SUCCESS_RATE:.0%}",
        flush=True,
    )
    if rate < MIN_SUCCESS_RATE:
        reasons = pd.Series([row["failure_reason"] for row in failed.values()]).value_counts()
        raise RuntimeError(
            f"LOW LIFT SUCCESS: {len(mapped)}/{len(elements)} ({rate:.4%}) "
            f"< {MIN_SUCCESS_RATE:.0%}; stopped without changing minMatch. "
            f"Unmapped={len(failed)}; reasons={reasons.to_dict()}. "
            "Both BED outputs are preserved; sequence verification was not run."
        )
    fasta = output_dir / "lifted.hg38.fa"
    commands.append([executables[1], f"-bed={mapped_bed}", "-noMask", str(genome), str(fasta)])
    run_tool(commands[1], output_dir / "twoBitToFa.log")
    extracted = {}
    for record in SeqIO.parse(fasta, "fasta"):
        if record.id in extracted:
            raise ValueError(f"Duplicate FASTA name: {record.id}")
        extracted[record.id] = str(record.seq).upper()
    if set(extracted) != set(mapped):
        raise AssertionError("Extracted FASTA IDs do not exactly equal mapped BED IDs")
    rows = []
    for row in elements.to_dict("records"):
        element_id = row["element_id"]
        original = {f"hg19_{col}": row.pop(col) for col in ("chrom", "start", "end", "strand")}
        result = (
            row
            | original
            | dict(
                working_build=config.WORKING_BUILD,
                chrom="",
                start="",
                end="",
                strand="",
                mapping_strand="",
                failure_reason="",
                mismatch_type="",
                best_orientation="",
                differences="",
                differing_positions_1based="",
                extracted_length="",
                extracted_sequence="",
                verified=False,
            )
        )
        if element_id in failed:
            result.update(
                bucket="failed-to-lift", failure_reason=failed[element_id]["failure_reason"]
            )
        else:
            location = mapped[element_id]
            reference = extracted[element_id]
            if location["end"] <= location["start"] or location["start"] < 0:
                raise ValueError(f"Invalid lifted interval: {element_id}: {location}")
            if len(reference) != location["end"] - location["start"]:
                raise AssertionError(f"FASTA/BED length disagreement: {element_id}")
            check = verify_sequence(row["sequence"], reference)
            if row["adapter_status"] == "unresolved_long_oligo":
                check.update(verified=False, mismatch_type="adapter_unresolved")
            orientation = location["strand"]
            if check["best_orientation"] == "-":
                orientation = "+" if orientation == "-" else "-"
            result.update(location)
            result.update(check)
            result.update(
                mapping_strand=location["strand"],
                strand=orientation,
                extracted_sequence=reference,
                bucket="lifted+verified" if check["verified"] else "lifted+mismatch",
            )
        rows.append(result)
    audit = pd.DataFrame(rows)
    # These assertions are the acceptance gate, not an optimistic report label.
    assert audit.element_id.is_unique and set(audit.element_id) == input_ids
    assert audit.bucket.isin(BUCKETS).all()
    assert sum(int(audit.bucket.eq(bucket).sum()) for bucket in BUCKETS) == len(elements)
    print("Intersecting all source elements with hg19 RepeatMasker...", flush=True)
    audit["in_repeat"], repeat_count = repeat_overlap(elements, repeats)
    verified = audit.loc[audit.bucket.eq("lifted+verified")].copy()
    assert verified.verified.all() and verified.sequence.str.len().eq(200).all()
    assert (verified.end - verified.start).eq(200).all()
    assert all(
        verify_sequence(row.sequence, row.extracted_sequence)["verified"]
        for row in verified.itertuples(index=False)
    )
    audit.to_csv(output_dir / "element_audit.csv", index=False)
    columns = [
        "element_id",
        "project",
        "chrom",
        "start",
        "end",
        "strand",
        "sequence",
        "source_build",
        "working_build",
        "hg19_chrom",
        "hg19_start",
        "hg19_end",
        "hg19_strand",
        "strand_supplied",
        "source_oligo_name",
        "sequence_join",
        "in_repeat",
    ]
    text = render_report(audit, repeat_count, fingerprints, commands)
    temporary = verified_path.with_suffix(".csv.tmp")
    verified[columns].to_csv(temporary, index=False)
    os.replace(temporary, verified_path)
    report_path.write_text(text)
    print(text, flush=True)
    return audit


def main() -> None:
    """Expose an offline default workflow and a strictly declared custom hg19 CSV."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, help="CSV with explicit source_build=hg19 and library DNA"
    )
    parser.add_argument("--chain", type=Path, default=config.SOURCE_TO_WORKING_CHAIN)
    parser.add_argument("--genome", type=Path, default=config.HG38_TWOBIT)
    parser.add_argument("--repeats", type=Path, default=config.HG19_REPEATS)
    parser.add_argument("--liftover", default="liftOver")
    parser.add_argument("--twobit", default="twoBitToFa")
    parser.add_argument("--adapter-5p", default="")
    parser.add_argument("--adapter-3p", default="")
    args = parser.parse_args()
    # Invalidate the previous consumable CSV even if loading the new source fails.
    config.HG38_ELEMENTS.unlink(missing_ok=True)
    try:
        if args.input:
            elements = pd.read_csv(args.input, keep_default_na=False)
            fingerprints = {"custom_input": local_md5(args.input)}
        else:
            elements = published_elements()
            fingerprints = {"GEO_archive": local_md5(GEO_ARCHIVE), "Table_S2": local_md5(TABLE_S2)}
        convert(
            elements,
            args.chain,
            args.genome,
            args.repeats,
            AUDIT_DIR,
            config.HG38_ELEMENTS,
            REPORT,
            liftover=args.liftover,
            twobit=args.twobit,
            adapter5=args.adapter_5p,
            adapter3=args.adapter_3p,
            extra_inputs=fingerprints,
        )
    except (ValueError, OSError, RuntimeError, AssertionError) as exc:
        failure = (
            "# hg19 to hg38 conversion: FAILED\n\n"
            f"```text\n{type(exc).__name__}: {exc}\n```\n\n"
            "No verified CSV was published. Retained intermediates are diagnostic only.\n"
        )
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(failure)
        print(failure, flush=True)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
