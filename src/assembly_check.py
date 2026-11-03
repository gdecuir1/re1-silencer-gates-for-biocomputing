"""Resolve GEO's coordinate assembly by matching published full-length oligos.

The evidence is independent of motif strength: Table S2 provides the synthesized
200 bp sequences, GEO supplies the interval numbers, and UCSC supplies reference
DNA at those exact numbers in hg19 and hg38. No liftover, aligner, statistical model,
or motif threshold is used. Both DNA orientations are tested by exact equality.

By default this script is offline. ``python -m src.assembly_check --fetch`` permits
only the small Table S2 archive and bounded 202 bp UCSC sequence requests. Complete
genomes are never requested. Every response is cached with download provenance.
"""

import argparse
import gzip
import json
import tarfile
import time
import zipfile
from datetime import UTC, datetime

import pandas as pd
from Bio.Seq import Seq

from src import config
from src.downloads import download_small_file, local_md5

EVIDENCE_DIR = config.RAW_DIR / "assembly_check"
TABLE_S2 = EVIDENCE_DIR / "Table_S2_publisher.zip"
TABLE_S2_URL = "https://ars.els-cdn.com/content/image/1-s2.0-S2666979X22001926-mmc3.zip"
GEO_ARCHIVE = config.RAW_DIR / "geo" / "GSE196171" / "GSE196171_RAW.tar"
TABLE_S2_MEMBER = "xgen_100234_mmc03.tsv"
REPORT_PATH = config.TABLES_DIR / "coordinate_assembly.md"
COMPARISONS_PATH = config.TABLES_DIR / "coordinate_assembly_comparisons.csv"


def load_evidence() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join native GEO intervals to the appropriate published oligo library.

    Returns:
        Joined element table and per-cell coordinate-consistency summary. Rows
        without an exact oligo-name match remain visible, never guessed or dropped.

    Raises:
        ValueError: Archives, schemas, interval definitions, or cell tables disagree.
    """
    if not zipfile.is_zipfile(TABLE_S2):
        raise ValueError("Table S2 response is not a ZIP; check for an HTML challenge")
    with zipfile.ZipFile(TABLE_S2) as archive:
        with archive.open(TABLE_S2_MEMBER) as handle:
            oligos = pd.read_csv(handle, sep="\t", low_memory=False)
    # The same oligo name can occur in different experimental libraries. Filtering
    # by the explicit library label avoids a misleading many-to-many merge.
    oligos = oligos.loc[oligos["library"].eq("WG_REST")].copy()
    if oligos["oligo_name"].duplicated().any():
        raise ValueError("Duplicate oligo names within the WG_REST library")
    columns = ["silencer", "chromosome", "start", "end", "project"]
    reference = None
    summaries = []
    with tarfile.open(GEO_ARCHIVE) as archive:
        members = sorted(member for member in archive.getnames() if "counts_WGREST_" in member)
        if len(members) != 4:
            raise ValueError(f"Expected four WGREST tables; observed {len(members)}")
        for member in members:
            with archive.extractfile(member) as compressed:
                with gzip.open(compressed, "rt") as handle:
                    frame = pd.read_csv(handle, sep="\t", low_memory=False)
            subset = frame.loc[frame["project"].isin(["motif", "noMotif"]), columns]
            # Collapsing identical interval metadata is NOT aggregating activity.
            intervals = subset.drop_duplicates().sort_values("silencer").reset_index(drop=True)
            if intervals["silencer"].duplicated().any() or intervals.isna().any().any():
                raise ValueError(f"Ambiguous or missing interval metadata in {member}")
            if not (intervals["end"] - intervals["start"]).eq(200).all():
                raise ValueError(f"Non-200 bp intervals in {member}")
            if reference is None:
                reference = intervals
            else:
                pd.testing.assert_frame_equal(reference, intervals)
            summaries.append(
                {"file": member, "elements": len(intervals), "coordinates_identical": True}
            )
    joined = reference.merge(
        oligos, left_on="silencer", right_on="oligo_name", how="left", validate="one_to_one"
    )
    matched = joined["sequence"].notna()
    if not joined.loc[matched, "sequence"].str.fullmatch("[ACGTacgt]{200}").all():
        raise ValueError("Matched native oligos must contain exactly 200 unambiguous DNA bases")
    return joined, pd.DataFrame(summaries)


def select_samples(elements: pd.DataFrame) -> pd.DataFrame:
    """Select one directly matched native oligo per chromosome and project.

    Args:
        elements: Joined GEO/Table S2 table, including unresolved sequence names.

    Returns:
        Deterministic sample chosen by chromosome and ascending interval position,
        before inspecting either reference assembly or any expression measurement.
    """
    available = elements.loc[elements["sequence"].notna()].copy()
    return (
        available.sort_values(["project", "chromosome", "start", "silencer"])
        .drop_duplicates(["project", "chromosome"])
        .reset_index(drop=True)
    )


def exact_orientation(expected: str, reference: str) -> str:
    """Compare a published oligo against reference DNA in both orientations.

    Args:
        expected: Published synthesized oligo sequence.
        reference: Reference sequence extracted at the candidate coordinates.

    Returns:
        '+', '-', 'both' for a palindrome, or 'none' for no exact full-length match.
    """
    expected, reference = expected.upper(), reference.upper()
    forward = expected == reference
    reverse = expected == str(Seq(reference).reverse_complement())
    return "both" if forward and reverse else "+" if forward else "-" if reverse else "none"


def compare_sample(row: pd.Series, assembly: str, *, fetch: bool) -> dict[str, object]:
    """Compare one complete oligo at unchanged and off-by-one coordinate starts.

    Args:
        row: Published sequence and its original GEO interval metadata.
        assembly: Candidate reference, restricted to hg19 or hg38.
        fetch: Permit fetching absent small cache entries.

    Returns:
        Exact-match results for the original interval and both adjacent 200 bp
        windows. These adjacent tests diagnose origin errors; they never silently
        change the original coordinates or decide a match at another genomic locus.
    """
    if assembly not in {"hg19", "hg38"}:
        raise ValueError("This audit compares only hg19 and hg38")
    chrom = "chr" + str(row["chromosome"])
    start, end = int(row["start"]), int(row["end"])
    if start < 1:
        raise ValueError("Cannot request a left diagnostic flank before coordinate zero")
    url = (
        f"https://api.genome.ucsc.edu/getData/sequence?genome={assembly};chrom={chrom};"
        f"start={start - 1};end={end + 1}"
    )
    path = EVIDENCE_DIR / "ucsc" / f"{assembly}_{chrom}_{start - 1}_{end + 1}.json"
    if fetch and not path.exists():
        # Respect UCSC's recommendation of no more than one request per second.
        time.sleep(1)
    download_small_file(url, path, fetch=fetch, max_bytes=20_000)
    payload = json.loads(path.read_text())
    expected_metadata = {"genome": assembly, "chrom": chrom, "start": start - 1, "end": end + 1}
    for key, expected in expected_metadata.items():
        if payload.get(key) != expected:
            raise ValueError(
                f"UCSC response mismatch: {key} expected {expected}, got {payload.get(key)}"
            )
    dna = payload.get("dna", "").upper()
    if len(dna) != 202 or any(base not in "ACGTN" for base in dna):
        raise ValueError("UCSC must return exactly the requested 202 DNA bases")
    return {
        "element_id": row["silencer"],
        "project": row["project"],
        "chrom": chrom,
        "start": start,
        "end": end,
        "assembly": assembly,
        "exact_orientation": exact_orientation(row["sequence"], dna[1:-1]),
        "minus_one_orientation": exact_orientation(row["sequence"], dna[:-2]),
        "plus_one_orientation": exact_orientation(row["sequence"], dna[2:]),
        "url": url,
        "response_md5": local_md5(path),
    }


def write_report(
    elements: pd.DataFrame,
    samples: pd.DataFrame,
    comparisons: pd.DataFrame,
    consistency: pd.DataFrame,
) -> str:
    """Write all observed results and distinguish sample evidence from census checks.

    Args:
        elements: All canonical and no-motif GEO interval definitions.
        samples: Deterministically selected sequence-bearing oligos.
        comparisons: Every original and shifted comparison, including failures.
        consistency: Per-cell interval agreement results.

    Returns:
        The full Markdown report text for printing to standard output.
    """
    rows = [
        "# Coordinate assembly evidence",
        "",
        f"Generated UTC: {datetime.now(UTC).isoformat(timespec='seconds')}",
        "",
        "## Assembly discrimination",
        "",
        "Exact full-200-bp sequence equality, allowing either orientation; no liftover, "
        "alignment tolerance, FIMO cutoff, or expression filtering.",
        "",
        "| Candidate | Tested | Exact original-coordinate matches | Start −1 | Start +1 |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for assembly, group in comparisons.groupby("assembly", sort=True):
        counts = [
            int(group[column].ne("none").sum())
            for column in ["exact_orientation", "minus_one_orientation", "plus_one_orientation"]
        ]
        rows.append(f"| {assembly} | {len(group)} | {counts[0]} | {counts[1]} | {counts[2]} |")
    winners = [
        assembly
        for assembly, group in comparisons.groupby("assembly")
        if len(group) == len(samples) and group["exact_orientation"].ne("none").all()
    ]
    conclusion = (
        f"The tested native intervals support **{winners[0]}**."
        if len(winners) == 1
        else "**Assembly unresolved: no unique all-matching candidate.**"
    )
    rows += [
        "",
        conclusion,
        "",
        "This is empirical assembly discrimination for the sampled native intervals, "
        "not a claim that every record or variant in the study was validated.",
        "",
        "## Source and population checks",
        "",
        f"Table S2: [{TABLE_S2.name}]({TABLE_S2_URL}); local MD5 `{local_md5(TABLE_S2)}`.",
        f"GEO archive local MD5: `{local_md5(GEO_ARCHIVE)}`.",
        "",
        "| Project | GEO elements | Expected | Direct sequence-name matches | Unmatched |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for project, group in elements.groupby("project"):
        count = int(group["sequence"].notna().sum())
        expected = 8436 if project == "motif" else 4430
        rows.append(f"| {project} | {len(group)} | {expected} | {count} | {len(group) - count} |")
    lengths = elements.loc[elements["sequence"].notna(), "sequence"].str.len().value_counts()
    rows += [
        "",
        f"Observed sequence-length distribution for direct matches: {lengths.to_dict()}.",
        "Unmatched names are not repaired, imputed, or counted as reference mismatches.",
        "",
        "| GEO table | Native intervals | Coordinates identical across cell types |",
        "| --- | ---: | --- |",
    ]
    for row in consistency.itertuples(index=False):
        rows.append(f"| {row.file} | {row.elements} | {row.coordinates_identical} |")
    rows += [
        "",
        "## Deterministic sample and individual observations",
        "",
        "One directly name-matched sequence per available project/chromosome group, "
        "chosen by ascending genomic start and then ID, before querying references. "
        "Table S2 is restricted to library WG_REST to avoid cross-library duplicate names. "
        "Identical enhancer-context metadata are deduplicated without using activity values.",
        "",
        "| Element | Project | Candidate | Original | Start −1 | Start +1 |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for row in comparisons.itertuples(index=False):
        rows.append(
            f"| {row.element_id} | {row.project} | {row.assembly} | "
            f"{row.exact_orientation} | {row.minus_one_orientation} | "
            f"{row.plus_one_orientation} |"
        )
    rows += [
        "",
        "`+`/`-` denote exact forward/reverse-complement equality; `none` means "
        "no exact 200 bp match. Coordinates in the CSV remain the unchanged GEO values.",
        "",
        "## Reproduction and limits",
        "",
        "Run `python -m src.assembly_check` for an offline replay. Use `--fetch` only "
        "to obtain missing small evidence files. All UCSC responses are cached with "
        "URL, UTC timestamp, byte count, and local MD5 in data/README.md. "
        "No whole genomes are downloaded. The comparison CSV records query URLs and hashes.",
        "",
        "The queried UCSC intervals include one diagnostic flanking base at each end "
        "(202 bp total). The original extraction uses zero-based, half-open [start, end). "
        "See [UCSC API coordinate documentation](https://genome.ucsc.edu/goldenPath/help/api.html).",
        "",
        "The PMC Table S2 endpoint returned an HTML browser challenge (1,817 bytes), "
        "not a ZIP. That failed response is retained as data/raw/assembly_check/mmc3.zip; "
        "the actual ZIP was obtained from the publisher URL above. The publisher article "
        "page returned HTTP 403 to the browsing tool; its supplement URL succeeded.",
        "",
    ]
    text = "\n".join(rows)
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    comparisons.to_csv(COMPARISONS_PATH, index=False)
    REPORT_PATH.write_text(text)
    return text


def main() -> None:
    """Run the explicitly scoped assembly audit, offline unless --fetch is supplied."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Fetch missing small evidence files")
    args = parser.parse_args()
    download_small_file(TABLE_S2_URL, TABLE_S2, fetch=args.fetch)
    elements, consistency = load_evidence()
    samples = select_samples(elements)
    print(
        f"Testing {len(samples)} native oligos against two assemblies (202 bp per request).",
        flush=True,
    )
    observations = []
    for index, (_, row) in enumerate(samples.iterrows(), start=1):
        for assembly in ("hg19", "hg38"):
            observations.append(compare_sample(row, assembly, fetch=args.fetch))
        print(f"Compared {index}/{len(samples)}: {row['silencer']}", flush=True)
    comparisons = pd.DataFrame(observations)
    print(write_report(elements, samples, comparisons, consistency))
    complete_matches = comparisons.groupby("assembly")["exact_orientation"].agg(
        lambda states: states.ne("none").all()
    )
    if complete_matches.sum() != 1:
        raise SystemExit("Assembly unresolved: review the written evidence report")


if __name__ == "__main__":
    main()
