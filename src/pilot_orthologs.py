"""Separate motif-window mapping from element mapping before extracting mouse DNA.

A full element can pass minMatch while its motif window fails. Such an element
must remain UNDETERMINED, not become a motif loss. The converse also occurs in
this pilot: the user explicitly approved treating motif-only lifts as
UNDETERMINED too, with a distinct reason. All 390 IDs remain in the sidecar.

The configured MOTIF_SPAN is zero-based with BOTH endpoints inclusive. Windows
are oriented with the library sequence, then written as standard half-open BED6.
Both tracks use the unchanged MINMATCH_CROSS_SPECIES=0.5. Only jointly mapped,
concordant intervals are extracted; no motif scoring or evolutionary inference
is performed by this module. Genome verification requires only hg38 and mm39.
"""

import json
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from src import config, pilot, pilot_mouse
from src.downloads import local_md5
from src.lift import lift_intervals, parse_bed, run_tool

ROOT = pilot_mouse.OUTPUT
WINDOW_DIR = ROOT / "motif_window"
EXTRACTION_DIR = ROOT / "ortholog_extraction"
FASTA = ROOT / "orthologs.mm39.fa"
SIDECAR = ROOT / "orthologs.mm39.csv"
REPORT = config.TABLES_DIR / "pilot_mouse_orthologs.md"
PADDING = 15
EXPECTED_PILOT_COUNT = 390
CALLS = ("MAPPED", "UNDETERMINED", "ABSENT")


def motif_windows(elements: dict, span: tuple[int, int] = config.MOTIF_SPAN) -> dict:
    """Translate inclusive, zero-based oligo offsets into strand-aware genomic BED.

    For plus-strand [s,e), inclusive [a,b] becomes [s+a,s+b+1).
    For minus-strand [s,e), it becomes [e-b-1,e-a). In both cases the
    21 bp window is carried with the same ID and strand as the full element.
    No coordinate is inferred from the historically named element ID.
    """
    a, b = span
    if not (0 <= a <= b < config.ELEMENT_LENGTH and b - a + 1 == 21):
        raise ValueError(f"Expected an inclusive 21 bp motif span within 200 bp; observed {span}")
    windows = {}
    for name, row in elements.items():
        if row["end"] - row["start"] != 200 or row["start"] < 0:
            raise ValueError(f"Expected a 200 bp nonnegative source interval: {name}: {row}")
        if row["strand"] not in ("+", "-"):
            raise ValueError(f"Invalid source strand for {name}: {row['strand']}")
        start, end = (
            (row["start"] + a, row["start"] + b + 1)
            if row["strand"] == "+"
            else (row["end"] - b - 1, row["end"] - a)
        )
        assert end - start == 21
        assert row["start"] <= start < end <= row["end"]
        windows[name] = dict(chrom=row["chrom"], start=start, end=end, strand=row["strand"])
    return windows


def write_records(records: dict, path: Path) -> None:
    """Write stable BED6 rows with original element names and no width correction."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, row in sorted(records.items()):
            handle.write(
                f"{row['chrom']}\t{row['start']}\t{row['end']}\t{name}\t0\t{row['strand']}\n"
            )


def require_partition(source: dict, mapped: dict, unmapped: dict, label: str) -> None:
    """Assert each track accounts for every ID once, including failures with reasons."""
    if set(mapped) & set(unmapped) or set(mapped) | set(unmapped) != set(source):
        raise AssertionError(f"{label} outputs must partition all input IDs")
    for name, row in unmapped.items():
        if not row["failure_reason"]:
            raise ValueError(f"{label} failure lacks a reason: {name}")
        if any(row[key] != source[name][key] for key in ("chrom", "start", "end", "strand")):
            raise AssertionError(f"{label} failed interval changed: {name}")


def coordinate_fields(prefix: str, row: dict | None) -> dict:
    """Preserve missing coordinates as empty fields, not zeros or false overlap claims."""
    result = {
        f"{prefix}_{key}": row[key] if row else "" for key in ("chrom", "start", "end", "strand")
    }
    result[f"{prefix}_width"] = row["end"] - row["start"] if row else ""
    return result


def classify_tracks(
    source: dict, windows: dict, element_track: dict, motif_track: dict
) -> pd.DataFrame:
    """Join both tracks and apply the user-approved exhaustive three-way policy.

    Both lift -> MAPPED; exactly one lifts -> UNDETERMINED; neither -> ABSENT.
    The two UNDETERMINED cases remain distinct in call_reason. ABSENT is an
    alignment-based candidate loss, not proof that a biological motif was lost.
    A jointly lifted window outside its parent interval (or on a different
    chromosome/strand) fails publication instead of extracting an unrelated site.
    """
    require_partition(source, element_track["mapped"], element_track["unmapped"], "element")
    require_partition(windows, motif_track["mapped"], motif_track["unmapped"], "motif")
    if set(source) != set(windows):
        raise AssertionError("Motif-window IDs differ from element IDs")
    rows = []
    for name, original in sorted(source.items()):
        element = element_track["mapped"].get(name)
        motif = motif_track["mapped"].get(name)
        if element and motif:
            if (
                element["chrom"] != motif["chrom"]
                or element["strand"] != motif["strand"]
                or not element["start"] <= motif["start"] < motif["end"] <= element["end"]
            ):
                raise ValueError(f"Discordant jointly lifted tracks for {name}: {element}; {motif}")
            call, reason = "MAPPED", "both_lifted"
        elif element:
            call, reason = "UNDETERMINED", "element_only_motif_unmapped"
        elif motif:
            call, reason = "UNDETERMINED", "motif_only_element_unmapped"
        else:
            call, reason = "ABSENT", "neither_lifted"
        width = element["end"] - element["start"] if element else None
        rows.append(
            dict(
                element_id=name,
                call=call,
                call_reason=reason,
                element_lifted=bool(element),
                motif_window_lifted=bool(motif),
                include_in_fasta=call == "MAPPED",
                exclude_mouse_downstream=call == "UNDETERMINED",
                source_build="hg38",
                target_build="mm39",
                min_match=config.MINMATCH_CROSS_SPECIES,
                **coordinate_fields("hg38_element", original),
                **coordinate_fields("hg38_motif", windows[name]),
                **coordinate_fields("mm39_element", element),
                **coordinate_fields("mm39_motif", motif),
                width_outside_160_240=not 160 <= width <= 240 if element else "",
                element_failure_reason=""
                if element
                else element_track["unmapped"][name]["failure_reason"],
                motif_failure_reason=""
                if motif
                else motif_track["unmapped"][name]["failure_reason"],
            )
        )
    audit = pd.DataFrame(rows)
    assert audit.element_id.is_unique and set(audit.element_id) == set(source)
    assert audit.call.isin(CALLS).all()
    assert sum(int(audit.call.eq(call).sum()) for call in CALLS) == len(source)
    return audit


def padded_records(
    audit: pd.DataFrame, sizes: dict[str, int], padding: int = PADDING
) -> tuple[dict, pd.DataFrame]:
    """Pad the mapped ELEMENT interval, recording any clipping at chromosome ends.

    Padding is genomic left/right (15 bp each regardless of strand). At a chromosome
    boundary only available bases can be extracted; the sidecar records actual
    padding on each side so shortened context is never concealed.
    """
    if padding < 0:
        raise ValueError("Padding must be nonnegative")
    audit = audit.copy()
    for key in (
        "padded_chrom",
        "padded_start",
        "padded_end",
        "padded_strand",
        "padding_left",
        "padding_right",
        "padding_clipped",
        "extracted_length",
        "sequence_under_50bp",
        "sequence_n_bases",
    ):
        audit[key] = pd.Series([""] * len(audit), index=audit.index, dtype=object)
    records = {}
    for index, row in audit.loc[audit.call.eq("MAPPED")].iterrows():
        chrom, start, end, strand = [
            row[f"mm39_element_{key}"] for key in ("chrom", "start", "end", "strand")
        ]
        if chrom not in sizes or not 0 <= start < end <= sizes[chrom]:
            raise ValueError(f"Element lies outside the verified mouse genome: {row.element_id}")
        left, right = max(0, start - padding), min(sizes[chrom], end + padding)
        records[row.element_id] = dict(chrom=chrom, start=left, end=right, strand=strand)
        values = dict(
            padded_chrom=chrom,
            padded_start=left,
            padded_end=right,
            padded_strand=strand,
            padding_left=start - left,
            padding_right=right - end,
            padding_clipped=(start - left != padding or right - end != padding),
        )
        for key, value in values.items():
            audit.at[index, key] = value
    return records, audit


def extract(
    audit: pd.DataFrame,
    sizes: dict[str, int],
    genome: Path,
    directory: Path,
    candidate: Path,
    *,
    executable: str = "twoBitToFa",
) -> pd.DataFrame:
    """Extract exactly the MAPPED IDs and validate headers, alphabet, and every length.

    A short sequence is flagged, never dropped or reclassified. Non-MAPPED rows
    retain empty sequence measurements and never appear as fabricated empty FASTA
    records. twoBitToFa handles orientation from BED6 and -noMask removes casing.
    """
    records, audit = padded_records(audit, sizes)
    bed = directory / "padded_elements.mm39.bed"
    write_records(records, bed)
    if records:
        tool = shutil.which(executable)
        if not tool:
            raise FileNotFoundError(f"Missing extraction tool: {executable}")
        run_tool(
            [tool, f"-bed={bed}", "-noMask", str(genome), str(candidate)],
            directory / "twoBitToFa.log",
        )
    else:
        candidate.write_text("")
        (directory / "twoBitToFa.log").write_text(
            "NOT RUN: no MAPPED elements; empty FASTA expected.\n"
        )
    sequences = {}
    for row in SeqIO.parse(candidate, "fasta"):
        if row.id in sequences:
            raise ValueError(f"Duplicate extracted FASTA ID: {row.id}")
        dna = str(row.seq)
        if row.description != row.id or set(dna) - set("ACGTN"):
            raise ValueError(f"Expected bare element ID and uppercase DNA: {row.id}")
        sequences[row.id] = dna
    assert set(sequences) == set(records), "FASTA IDs must equal exactly the MAPPED IDs"
    for index, row in audit.loc[audit.call.eq("MAPPED")].iterrows():
        dna = sequences[row.element_id]
        expected = row.padded_end - row.padded_start
        if len(dna) != expected:
            raise AssertionError(f"{row.element_id}: extracted {len(dna)} bp; expected {expected}")
        audit.at[index, "extracted_length"] = len(dna)
        audit.at[index, "sequence_under_50bp"] = len(dna) < 50
        audit.at[index, "sequence_n_bases"] = dna.count("N")
    assert len(sequences) == int(audit.call.eq("MAPPED").sum())
    return audit


def render(audit: pd.DataFrame, genomes: list[dict], fingerprints: dict) -> str:
    """Render actual counts and the complete observed sequence-length distribution."""
    mapped = audit.loc[audit.call.eq("MAPPED")]
    counts = audit.call.value_counts()
    distorted = int(mapped.width_outside_160_240.eq(True).sum())
    short = int(mapped.sequence_under_50bp.eq(True).sum())
    clipped = int(mapped.padding_clipped.eq(True).sum())
    a, b = config.MOTIF_SPAN
    lines = [
        "# Pilot mouse ortholog sequences and motif-window coverage",
        "",
        "## Coordinate and mapping contract",
        "",
        f"Input: **{len(audit)} canonical {config.PILOT_CHROM} elements**, expected 390; "
        + f"element width 200 bp. Both tracks use **minMatch={config.MINMATCH_CROSS_SPECIES:g}**, "
        + "with no -multiple flag. No thresholds were lowered.",
        "",
        f"Configured MOTIF_SPAN=({a}, {b}) is **zero-based, both endpoints inclusive**. "
        + f"For plus-strand elements [s,e), motif BED is [s+{a},s+{b + 1}); "
        + f"for minus-strand elements it is [e-{b + 1},e-{a}). All input motif windows "
        + "are asserted to be 21 bp. The authors' one-based FIMO positions 91–111 "
        + "become zero-based inclusive 90–110. Published oligo/score regression fixtures "
        + "are checked separately by make pilot-validate.",
        "",
        "## Three-way calls",
        "",
        "| Call | Observed | Meaning |",
        "| --- | ---: | --- |",
    ]
    meanings = {
        "MAPPED": "Both tracks lifted; motif lies within its element on the same chromosome/strand",
        "UNDETERMINED": "Exactly one track lifted; exclude mouse downstream, never call a loss",
        "ABSENT": "Neither track lifted; candidate alignment loss, not proven biological loss",
    }
    for call in CALLS:
        lines.append(f"| {call} | {counts.get(call, 0)} | {meanings[call]} |")
    lines += [
        "",
        f"Sum: **{sum(counts)} / {len(audit)}**, expected **390 / 390**: **PASS**. "
        + "Every input ID occurs once in the sidecar.",
        "",
        "| Track combination | Elements |",
        "| --- | ---: |",
    ]
    for reason in (
        "both_lifted",
        "element_only_motif_unmapped",
        "motif_only_element_unmapped",
        "neither_lifted",
    ):
        lines.append(f"| {reason} | {int(audit.call_reason.eq(reason).sum())} |")
    lines += [
        "",
        "The user explicitly approved assigning motif-only lifts to UNDETERMINED. "
        + "They retain a separate reason and their successfully mapped motif coordinates. "
        + "No UNDETERMINED element is in the FASTA or silently recoded as ABSENT.",
        "",
        f"Element-track successes: {int(audit.element_lifted.sum())}; motif-track successes: "
        + f"{int(audit.motif_window_lifted.sum())}. "
        + "Joint-track chromosome/strand/containment checks: "
        + f"{len(mapped)}/{len(mapped)} PASS.",
        "",
        "A window passing minMatch 0.5 need not align every one of its 21 source bases. "
        + "This report distinguishes mapping availability; "
        + "it does not score motifs or establish motif loss.",
        "",
        "## Width-band breakdown within MAPPED",
        "",
        "| Original lifted element width | Count |",
        "| --- | ---: |",
        f"| Within 160–240 bp inclusive | {len(mapped) - distorted} |",
        f"| Outside 160–240 bp (>20% deviation from 200) | {distorted} |",
        "",
        "Width flags refer to the unpadded element, not the extracted sequence. "
        + "All MAPPED intervals are retained, including those outside the band.",
        "",
        "## Sequence extraction and length sanity checks",
        "",
        f"FASTA records: **{len(mapped)}**; expected MAPPED count **{len(mapped)}**: **PASS**. "
        + f"Sidecar rows: **{len(audit)}**, expected **390**: **PASS**.",
        "",
        f"Padding requested: {PADDING} bp on each side of the lifted ELEMENT interval. "
        + f"Chromosome-boundary clipping: **{clipped}** elements; "
        + "actual left/right padding is recorded. "
        + "Extraction uses BED6 strand and `twoBitToFa -bed=... -noMask mm39.2bit`. "
        + "Each observed length must equal padded_end - padded_start.",
        "",
        f"Sequences under 50 bp: **{short}** "
        + "(flagged as too short for motif plus context; never dropped). "
        + f"Sequences containing N: **{int(mapped.sequence_n_bases.gt(0).sum())}**.",
        "",
        "| Extracted length (bp) | Observed count |",
        "| ---: | ---: |",
    ]
    for length, count in sorted(Counter(mapped.extracted_length).items()):
        lines.append(f"| {length} | {count} |")
    if short:
        lines += [
            "",
            "Short IDs: "
            + ", ".join(mapped.loc[mapped.sequence_under_50bp.eq(True), "element_id"]),
        ]
    lines += [
        "",
        "## Required genomes and provenance",
        "",
        "| Assembly | Status | Observed / expected sequences | Observed / expected bases |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in genomes:
        lines.append(
            f"| {row['assembly']} | {row['status']} | {row['sequence_count']} / "
            + f"{row['expected_sequence_count']} | "
            + f"{row['total_bases']} / {row['expected_total_bases']} |"
        )
    lines += [
        "",
        "Both required references also passed their unchanged publisher-MD5 "
        + "and stable-file checks. "
        + "No other species genome is inspected.",
        "",
        "| Input | Local MD5 |",
        "| --- | --- |",
    ]
    lines += [f"| {name} | `{value}` |" for name, value in sorted(fingerprints.items())]
    lines += [
        "",
        "## Files and reproduction",
        "",
        "- `data/processed/pilot_mouse/orthologs.mm39.fa`: "
        "exactly MAPPED IDs; no fabricated sequences.",
        "- `data/processed/pilot_mouse/orthologs.mm39.csv`: one row per pilot ID, "
        "calls/reasons, both tracks, width flags, padded coordinates, and sequence measurements.",
        "- `data/processed/pilot_mouse/motif_window/`: "
        "original 21 bp BED, both liftOver outputs, and command log.",
        "- `data/processed/pilot_mouse/ortholog_extraction/`: "
        "replayed element track, padded BED, input identities, and extraction log.",
        "",
        "Run `make pilot-orthologs` offline. It rebuilds the pilot and element lift, verifies only "
        + "hg38/mm39, lifts the motif track at the same threshold, and publishes the FASTA/sidecar "
        + "only after exhaustive accounting and extraction checks. "
        + "A separate replay of the element "
        + "track is asserted identical to the existing final outputs, "
        + "preventing a stale input join. "
        + "No timestamps or randomized ordering enter these outputs.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    """Publish only after both-track accounting, containment, and extraction checks pass."""
    FASTA.unlink(missing_ok=True)
    SIDECAR.unlink(missing_ok=True)
    candidate = EXTRACTION_DIR / "orthologs.mm39.fa.tmp"
    candidate_csv = EXTRACTION_DIR / "orthologs.mm39.csv.tmp"
    try:
        genomes = pilot.verify_step("mouse")
        if any(row["status"] != "PASS" for row in genomes):
            raise ValueError("Required hg38/mm39 genome gate failed: " + json.dumps(genomes))
        source = parse_bed(pilot.BED)
        assert len(source) == EXPECTED_PILOT_COUNT, (
            f"Expected 390 pilot IDs; observed {len(source)}"
        )
        expected = pilot.choose_pilot(pd.read_csv(config.HG38_ELEMENTS, keep_default_na=False))
        EXTRACTION_DIR.mkdir(parents=True, exist_ok=True)
        expected_path = EXTRACTION_DIR / "expected_pilot.hg38.bed"
        pilot.write_bed(expected, expected_path)
        if expected_path.read_bytes() != pilot.BED.read_bytes():
            raise ValueError("Pilot differs from the canonical hg38 inputs; rebuild make pilot")
        element = lift_intervals(
            pilot.BED,
            "hg38",
            "mm39",
            config.CHAIN_PATHS["mm39"],
            EXTRACTION_DIR / "element_track",
            min_match=config.MINMATCH_CROSS_SPECIES,
        )
        for key, name in [("mapped", "lifted.mm39.bed"), ("unmapped", "unmapped.hg38.bed")]:
            previous = parse_bed(ROOT / "final" / name, unmapped=key == "unmapped")
            assert element[key] == previous, (
                f"Element {key} output is stale; rerun make pilot-mouse"
            )
        windows = motif_windows(source)
        window_bed = WINDOW_DIR / "input.hg38.bed"
        write_records(windows, window_bed)
        motif = lift_intervals(
            window_bed,
            "hg38",
            "mm39",
            config.CHAIN_PATHS["mm39"],
            WINDOW_DIR,
            min_match=config.MINMATCH_CROSS_SPECIES,
        )
        sizes = pilot.parse_sizes((pilot.METADATA / "mm39.chrom.sizes").read_text())
        pilot_mouse.validate_target(element, sizes)
        pilot_mouse.validate_target(motif, sizes)
        audit = classify_tracks(source, windows, element, motif)
        assert len(audit) == 390 and sum(int(audit.call.eq(c).sum()) for c in CALLS) == 390
        audit = extract(audit, sizes, pilot.GENOMES / "mm39.2bit", EXTRACTION_DIR, candidate)
        audit.to_csv(candidate_csv, index=False)
        fingerprints = {
            "pilot BED": local_md5(pilot.BED),
            "hg38 element CSV": local_md5(config.HG38_ELEMENTS),
            "hg38->mm39 chain": local_md5(config.CHAIN_PATHS["mm39"]),
            "element mapped BED": local_md5(element["mapped_path"]),
            "element unmapped BED": local_md5(element["unmapped_path"]),
            "motif input BED": local_md5(window_bed),
            "motif mapped BED": local_md5(motif["mapped_path"]),
            "motif unmapped BED": local_md5(motif["unmapped_path"]),
        }
        (EXTRACTION_DIR / "inputs.json").write_text(
            json.dumps(fingerprints, sort_keys=True, indent=2) + "\n"
        )
        text = render(audit, genomes, fingerprints)
        candidate.replace(FASTA)
        candidate_csv.replace(SIDECAR)
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(text)
        print(text, flush=True)
    except (OSError, ValueError, RuntimeError, AssertionError) as exc:
        FASTA.unlink(missing_ok=True)
        SIDECAR.unlink(missing_ok=True)
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        text = (
            "# Pilot mouse ortholog extraction: FAILED\n\n```text\n"
            f"{type(exc).__name__}: {exc}\n```\n\n"
            "No FASTA/sidecar published. Diagnostic track outputs are retained.\n"
        )
        REPORT.write_text(text)
        print(text, flush=True)
        raise SystemExit(1) from exc
    finally:
        candidate.unlink(missing_ok=True)
        candidate_csv.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
