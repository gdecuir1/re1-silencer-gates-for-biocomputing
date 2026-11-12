"""Score the human starting site and publish uncertainty-aware pilot comparisons.

The raw lift and FIMO tables remain measurements. Only a strong human canonical
site with uniquely reproduced, same-chain, intact motif correspondence enters
the direct comparison. Gapped alignments are review cases, not presumed errors
or automatically certified losses. Flank-rescued regional hits stay unresolved
at the motif level. No evolutionary direction or experimental binding is inferred.
"""

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from Bio.Seq import Seq

from src import config, motifs, pilot, pilot_fimo, pilot_orthologs, rescue
from src._utils import (
    read_fasta,
    run_command,
    validate_intervals,
    write_bed,
    write_fasta,
)
from src.alignment_evidence import audit_mappings
from src.classify import OrthologState
from src.downloads import local_md5

OUTPUT = config.PROCESSED_DIR / "pilot_mouse_evidence"
REPORT = config.TABLES_DIR / "pilot_mouse_evidence.md"
REFERENCE = config.ROOT / "references" / "re1_coordinate_validation.json"
SNAPSHOT = config.RESULTS_DIR / "archive" / "pilot_before_motif_coordinate_fix_20261111_45402c8"
RESCUE_PREFIX = config.PROCESSED_DIR / "pilot_mouse_rescue" / "rescue"
COMPARISONS = (
    "HUMAN_WEAK",
    "SOURCE_UNSCORABLE",
    "UNRESOLVED",
    "RETAINED_CANDIDATE",
    "NONRETENTION_CANDIDATE",
)


def extract_intervals(frame: pd.DataFrame, genome: Path, directory: Path) -> dict[str, str]:
    """Extract oriented intervals and require exact IDs and lengths."""
    directory.mkdir(parents=True, exist_ok=True)
    validate_intervals(frame, full_elements=False)
    bed, fasta = directory / "input.bed", directory / "input.fa"
    write_bed(frame, bed)
    if frame.empty:
        fasta.write_text("")
        return {}
    run_command(["twoBitToFa", f"-bed={bed}", "-noMask", str(genome), str(fasta)])
    records = read_fasta(fasta)
    if set(records) != set(frame.element_id):
        raise ValueError("Extracted IDs differ from requested intervals")
    for row in frame.itertuples(index=False):
        if len(records[row.element_id]) != row.end - row.start:
            raise ValueError(f"Wrong sequence length for {row.element_id}")
    return records


def score_records(records: dict[str, str], directory: Path) -> pd.DataFrame:
    """Use the common scorer; retain weak scores and make unscorable DNA explicit."""
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "sequences.fa"
    write_fasta(records, fasta)
    if not records:
        return pd.DataFrame(columns=["element_id", "score", "status"])
    hits = motifs.score_sequences(fasta, config.MOTIF_PATH)
    if not hits.element_id.is_unique or set(hits.element_id) != set(records):
        raise ValueError("Scorer must cover exactly the supplied sequences")
    hits["sequence_bp"] = hits.element_id.map(lambda name: len(records[name]))
    hits["ambiguous_bp"] = hits.element_id.map(
        lambda name: sum(base.upper() not in "ACGT" for base in records[name])
    )
    hits["status"] = [
        "unscorable"
        if not row.has_hit or (row.ambiguous_bp and not row.retains_motif)
        else "present"
        if row.retains_motif
        else "below_threshold"
        for row in hits.itertuples(index=False)
    ]
    return hits.sort_values("element_id").reset_index(drop=True)


def check_reference(directory: Path) -> pd.DataFrame:
    """Independently pin the span to published oligos, matched bases and FIMO scores."""
    reference = json.loads(REFERENCE.read_text())
    a, b = config.MOTIF_SPAN
    records = {}
    for row in reference["records"]:
        if (a, b) != (row["fimo_start"] - 1, row["fimo_stop"] - 1):
            raise ValueError("Configured motif span differs from published FIMO coordinates")
        core = row["sequence"][a : b + 1]
        expected = row["matched_sequence"]
        if row["fimo_strand"] == "-":
            expected = str(Seq(expected).reverse_complement())
        if core != expected:
            raise ValueError("Published motif bases disagree with the configured oligo slice")
        records[row["element_id"]] = core
    observed = score_records(records, directory)
    expected = pd.DataFrame(reference["records"])[["element_id", "score"]]
    compared = observed.merge(
        expected, on="element_id", suffixes=("", "_published"), validate="one_to_one"
    )
    if not compared.score.sub(compared.score_published).abs().le(0.00005).all():
        raise ValueError("Current FIMO scores disagree with the published score scale")
    return compared


def human_scores(
    source: pd.DataFrame, sidecar: pd.DataFrame, genome: Path, directory: Path
) -> pd.DataFrame:
    """Verify human reference DNA, then separately score core, 101 bp window and element."""
    if not source.element_id.is_unique or not sidecar.element_id.is_unique:
        raise ValueError("Duplicate source/sidecar IDs")
    selected = source[source.element_id.isin(sidecar.element_id)].sort_values("element_id")
    if set(selected.element_id) != set(sidecar.element_id):
        raise ValueError("Source is missing pilot IDs")
    validate_intervals(selected)
    if not selected.working_build.eq(config.GENOME_BUILD).all():
        raise ValueError("Human source is in the wrong assembly")
    records = extract_intervals(selected, genome, directory / "reference")
    coordinates = sidecar.set_index("element_id")
    for row in selected.itertuples(index=False):
        if records[row.element_id].upper() != row.sequence.upper():
            raise ValueError(f"Human reference differs from verified library DNA: {row.element_id}")
        for field in ("chrom", "start", "end", "strand"):
            if coordinates.loc[row.element_id, f"{config.GENOME_BUILD}_element_{field}"] != getattr(
                row, field
            ):
                raise ValueError(f"Sidecar/source coordinate mismatch for {row.element_id}")
    a, b = config.MOTIF_SPAN
    if a < pilot_fimo.FLANK_BP or b + 1 + pilot_fimo.FLANK_BP > config.ELEMENT_LENGTH:
        raise ValueError("Human comparison window does not fit inside the element")
    result = selected[["element_id"]].copy()
    for label, start, end in (
        ("core", a, b + 1),
        ("window", a - pilot_fimo.FLANK_BP, b + 1 + pilot_fimo.FLANK_BP),
        ("element", 0, config.ELEMENT_LENGTH),
    ):
        hits = score_records(
            {name: seq[start:end] for name, seq in records.items()}, directory / label
        )
        hits = hits.rename(columns={c: f"human_{label}_{c}" for c in hits if c != "element_id"})
        result = result.merge(hits, on="element_id", validate="one_to_one")
    return result


def mapping_requests(sidecar: pd.DataFrame, regional: pd.DataFrame, assembly: str) -> list[dict]:
    """Request evidence for both ordinary tracks and every accepted rescue flank."""
    requests = []
    for row in sidecar.to_dict("records"):
        for track in ("element", "motif"):
            if not row[f"{track if track == 'element' else 'motif_window'}_lifted"]:
                continue
            req = dict(element_id=row["element_id"], track=track)
            for label, prefix in (("source", config.GENOME_BUILD), ("target", assembly)):
                for key in ("chrom", "start", "end", "strand"):
                    value = row[f"{prefix}_{track}_{key}"]
                    req[f"{label}_{key}"] = int(value) if key in {"start", "end"} else value
            requests.append(req)
    for row in regional.loc[regional.valid_pair].to_dict("records"):
        flank = int(row["flank_bp"])
        for track, start, end in (
            ("left", row["source_start"] - flank, row["source_start"]),
            ("right", row["source_end"], row["source_end"] + flank),
        ):
            req = dict(
                element_id=row["element_id"],
                track=track,
                source_chrom=row["source_chrom"],
                source_start=int(start),
                source_end=int(end),
                source_strand="+",
            )
            for key in ("chrom", "start", "end", "strand"):
                value = row[f"{track}_{key}"]
                req[f"target_{key}"] = int(value) if key in {"start", "end"} else value
            requests.append(req)
    return requests


def alignment_summary(
    sidecar: pd.DataFrame, regional: pd.DataFrame, evidence: pd.DataFrame
) -> pd.DataFrame:
    """Keep alignment quality independent of source strength and mouse scores."""
    lookup = evidence.set_index(["element_id", "track"])
    rescued = regional.set_index("element_id")
    rows = []
    for source in sidecar.itertuples(index=False):
        row = dict(
            element_id=source.element_id,
            alignment_status="NO_MOTIF_CORRESPONDENCE",
            chain_consistent=False,
            motif_complete=False,
        )
        tracks = {}
        for track in ("element", "motif", "left", "right"):
            key = source.element_id, track
            if key not in lookup.index:
                continue
            hit = lookup.loc[key]
            tracks[track] = hit
            for field in (
                "chain_matches",
                "chain_ids",
                "aligned_bp",
                "source_unaligned_bp",
                "target_unaligned_bp",
                "source_left_trim",
                "source_right_trim",
                "block_count",
                "complete",
            ):
                row[f"{track}_{field}"] = hit[field]
        if "motif" in tracks and "element" in tracks:
            motif, element = tracks["motif"], tracks["element"]
            unique = motif.chain_matches == element.chain_matches == 1
            same = unique and motif.chain_ids == element.chain_ids
            row.update(chain_consistent=bool(same), motif_complete=bool(motif.complete))
            row["alignment_status"] = (
                "AMBIGUOUS_CHAIN"
                if not unique
                else "DISCORDANT_CHAINS"
                if not same
                else "INTACT_MOTIF"
                if motif.complete
                else "DISRUPTED_MOTIF"
            )
        elif source.element_id in rescued.index and rescued.loc[source.element_id, "valid_pair"]:
            left, right = tracks["left"], tracks["right"]
            unique = left.chain_matches == right.chain_matches == 1
            same = unique and left.chain_ids == right.chain_ids
            row.update(
                chain_consistent=bool(same),
                alignment_status=("REGIONAL_ANCHORS_ONLY" if same else "AMBIGUOUS_ANCHORS"),
            )
            row["left_inner_trim_bp"] = left.source_right_trim
            row["right_inner_trim_bp"] = right.source_left_trim
        elif source.call == "UNDETERMINED":
            row["alignment_status"] = "SINGLE_TRACK_ONLY"
        rows.append(row)
    return pd.DataFrame(rows)


def comparison(row: dict) -> tuple[str, str, str]:
    """Explicit eligibility gate; raw negatives never bypass uncertainty checks."""
    if row["human_core_status"] == "unscorable":
        return "SOURCE_UNSCORABLE", "human_core_unscorable", OrthologState.UNSCORABLE
    if row["human_core_status"] == "below_threshold":
        return "HUMAN_WEAK", "human_core_below_m_value", OrthologState.UNSCORABLE
    if row["human_core_status"] != "present":
        raise ValueError("Unknown human core status")
    if row["alignment_status"] != "INTACT_MOTIF":
        state = (
            OrthologState.FAILED_LIFT
            if row["alignment_status"] == "NO_MOTIF_CORRESPONDENCE"
            else OrthologState.UNSCORABLE
        )
        return "UNRESOLVED", row["alignment_status"].lower(), state
    if row["mouse_core_status"] == "present":
        return (
            "RETAINED_CANDIDATE",
            "strong_human_and_corresponding_mouse_site",
            OrthologState.RETAINED,
        )
    if (
        row["mouse_core_status"] == "below_threshold"
        and row["mouse_window_status"] == "below_threshold"
    ):
        return (
            "NONRETENTION_CANDIDATE",
            "strong_human_intact_correspondence_mouse_below_threshold",
            OrthologState.NO_MOTIF,
        )
    reason = (
        "nearby_mouse_motif_only"
        if row["mouse_window_status"] == "present"
        else "mouse_site_unscorable"
    )
    return "UNRESOLVED", reason, OrthologState.UNSCORABLE


def assemble(
    sidecar: pd.DataFrame,
    baseline: pd.DataFrame,
    human: pd.DataFrame,
    mouse_core: pd.DataFrame,
    quality: pd.DataFrame,
    regional: pd.DataFrame,
) -> pd.DataFrame:
    """Join all pilot IDs without dropping weak sources or overwriting observations."""
    result = sidecar.copy().rename(
        columns={
            "call": "mapping_call",
            "exclude_mouse_downstream": "mapping_exclude_mouse_downstream",
            "include_in_fasta": "mapping_include_in_fasta",
        }
    )
    if not mouse_core.element_id.is_unique or set(mouse_core.element_id) != set(
        sidecar.loc[sidecar.call.eq("MAPPED"), "element_id"]
    ):
        raise ValueError("Mouse core scores must cover exactly the MAPPED elements")
    if not regional.element_id.is_unique or set(regional.element_id) != set(
        sidecar.loc[sidecar.call.eq("ABSENT"), "element_id"]
    ):
        raise ValueError("Regional results must cover exactly the neither-track-lifted cohort")
    for frame in (baseline, human, quality):
        if not frame.element_id.is_unique or set(frame.element_id) != set(result.element_id):
            raise ValueError("Evidence must cover every pilot element exactly once")
    raw = baseline[["element_id", "motif_status", "score"]].rename(
        columns={"motif_status": "mouse_window_status", "score": "mouse_window_score"}
    )
    raw["mouse_window_status"] = raw.mouse_window_status.replace(
        {"absent": "below_threshold", "no_hit": "unscorable"}
    )
    core = mouse_core.rename(
        columns={c: f"mouse_core_{c}" for c in mouse_core if c != "element_id"}
    )
    for frame in (raw, human, core, quality):
        result = result.merge(frame, on="element_id", how="left", validate="one_to_one")
    regional_cols = [
        c
        for c in regional
        if c.startswith(("whole_", "center_", "gap_"))
        or c in {"element_id", "valid_pair", "rejection_reason", "search_space_sensitive"}
    ]
    result = result.merge(
        regional[regional_cols], on="element_id", how="left", validate="one_to_one"
    )
    calls = [comparison(row) for row in result.to_dict("records")]
    result[["comparison_status", "comparison_reason", "ortholog_state"]] = pd.DataFrame(
        calls, index=result.index
    )
    result["eligible_for_direct_comparison"] = result.comparison_status.isin(
        ["RETAINED_CANDIDATE", "NONRETENTION_CANDIDATE"]
    )
    result["exclude_mouse_downstream"] = ~result.eligible_for_direct_comparison
    result["candidate_nonretention"] = result.comparison_status.eq("NONRETENTION_CANDIDATE")
    assert len(result) == len(sidecar) and result.element_id.is_unique
    assert result.comparison_status.isin(COMPARISONS).all()
    assert result.loc[result.candidate_nonretention, "alignment_status"].eq("INTACT_MOTIF").all()
    assert result.loc[result.candidate_nonretention, "human_core_status"].eq("present").all()
    assert not result.loc[result.mapping_call.eq("ABSENT"), "candidate_nonretention"].any()
    return result.sort_values("element_id").reset_index(drop=True)


def render(result: pd.DataFrame, previous: pd.DataFrame | None) -> str:
    """Explain the exhaustive corrected counts and the limits of each comparison."""
    counts = result.comparison_status.value_counts().reindex(COMPARISONS, fill_value=0)
    alignment = (
        result.alignment_status.value_counts()
        .rename_axis("alignment_status")
        .reset_index(name="elements")
    )
    mapped = result[result.mapping_call.eq("MAPPED")]
    negatives = mapped[mapped.mouse_window_status.eq("below_threshold")]
    lines = [
        "# Corrected mouse pilot: source motif and alignment evidence",
        "",
        f"All **{len(result)}** elements retained. Intended motif: zero-based inclusive "
        f"{config.MOTIF_SPAN}; source-study FIMO positions 91–111 are one-based. "
        f"Shared matrix {config.JASPAR_ID}; threshold ≥{config.M_VALUE}.",
        "",
        "## Current comparison calls",
        "",
        rescue.markdown(counts.rename_axis("comparison").reset_index(name="elements")),
        "",
        "HUMAN_WEAK means a below-threshold human canonical site, not missing DNA. "
        "UNRESOLVED retains uncertain correspondence or scoring. Retention/nonretention "
        "are candidate sequence comparisons, not demonstrated REST binding or evolutionary loss. "
        "A human-positive/mouse-negative pair cannot establish the ancestral state.",
        "",
        "## Mapping and raw mouse scans",
        "",
        rescue.markdown(pd.crosstab(result.mapping_call, result.mouse_window_status).reset_index()),
        "",
        "The raw FIMO table's `absent` means only below threshold in the scanned region. "
        "Use this table's comparison_status, eligibility and evidence fields for contrasts; "
        "do not use raw ABSENT/absent labels as a loss variable.",
        "",
        "## Human starting sites",
        "",
        rescue.markdown(
            pd.DataFrame(
                {
                    mode: result[f"human_{mode}_status"].value_counts()
                    for mode in ("core", "window", "element")
                }
            )
            .fillna(0)
            .astype(int)
            .rename_axis("status")
            .reset_index()
        ),
        "",
        "The exact 21 bp canonical core controls eligibility. The separate 101 bp window "
        "and full 200 bp scans expose nearby motifs without substituting them for the "
        "intended site. All human sequences were freshly extracted and matched exactly "
        "to the verified source library. Five published plus/minus/weak/borderline "
        "reference sites independently check coordinates, motif bases and score scale.",
        "",
        "## Alignment evidence",
        "",
        rescue.markdown(alignment),
        "",
        f"Among {len(negatives)} MAPPED window negatives, "
        f"{int(negatives.human_core_status.eq('below_threshold').sum())} have weak human cores; "
        f"{int(negatives.comparison_status.eq('NONRETENTION_CANDIDATE').sum())} pass the direct "
        "comparison gates. The remaining rows retain their explicit reasons.",
        "",
        "INTACT_MOTIF requires a unique chain reproducing the motif mapping, the same "
        "chain for the element, all 21 human motif bases aligned and no extra target "
        "bases between them. Source bases missing from the alignment, target insertions, "
        "block counts, terminal trimming and chain IDs are retained. DISRUPTED_MOTIF "
        "does not mean technical failure: indel-mediated disruption can be real biology "
        "and needs review before entering a high-confidence contrast.",
        "",
        "## Flank rescue",
        "",
        "Regional rescue is recomputed from the corrected neither-track-lifted cohort. "
        "REGION_MOTIF_PRESENT / REGION_BELOW_THRESHOLD describe regional scans; neither "
        "certifies correspondence of the original motif. Accepted anchor pairs still "
        "receive chain-consistency and inner-edge-trimming diagnostics here. They remain "
        "unresolved at the motif level; NO_SYNTENY is never a loss. See pilot_mouse_rescue.md "
        "for both whole-gap and midpoint scans and all three flank sizes.",
        "",
        "## Downstream contract",
        "",
        "ortholog_state is an explicit conservative adapter to classify.py: retained and "
        "lifted_no_motif occur only for eligible direct comparisons; weak/unscorable human "
        "sites and unresolved mapped/regional evidence become lifted_unscorable; missing "
        "motif correspondence becomes failed_lift for strong human sites. Always retain "
        "human_core_status and comparison_reason to distinguish these exclusions. Weak "
        "human sites remain available for repression modeling; exclusion applies to the "
        "strong-human comparative contrast, not to the entire research dataset.",
        "",
        "## Prior run comparison",
        "",
    ]
    if previous is not None:
        pair = result.merge(
            previous[["element_id", "call", "motif_status"]], on="element_id", validate="one_to_one"
        )
        lines += [
            rescue.markdown(pd.crosstab(pair.call, pair.mapping_call).reset_index()),
            "",
            "The preserved baseline and its SHA-256 manifest are in "
            + str(SNAPSHOT.relative_to(config.ROOT))
            + ". "
            "comparison_to_previous.csv records every old/new raw call "
            "and new evidence-aware call.",
        ]
    else:
        lines += ["No prior snapshot supplied; old counts are not inferred."]
    lines += [
        "",
        "## Reproduction",
        "",
        "Run `make pilot-validate` to rebuild the pilot, both lift tracks, baseline FIMO, "
        "flank rescue and this evidence layer. No downloads are performed. Outputs are "
        "in data/processed/pilot_mouse_evidence/: calls.csv (all elements), human_scores.csv, "
        "mouse_core_scores.csv, alignment_blocks.csv, reference_check.csv and input hashes. "
        "`python -m src.pilot_evidence --check` detects stale inputs "
        "or altered published evidence. "
        "Raw and final calls are kept separately. Missing hits keep missing scores.",
        "",
    ]
    return "\n".join(lines)


def fingerprints(paths: list[Path]) -> dict[str, str]:
    """Hash source inputs and implementations for stale-output detection."""
    return {
        str(p.relative_to(config.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths
    }


def main() -> None:
    """Validate current raw outputs and publish the full evidence layer on success."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    manifest_path = OUTPUT / "manifest.json"
    if args.check:
        manifest = json.loads(manifest_path.read_text())
        for name, digest in (manifest["inputs"] | manifest["outputs"]).items():
            if hashlib.sha256((config.ROOT / name).read_bytes()).hexdigest() != digest:
                raise ValueError(
                    f"Evidence is stale or modified: {name}; rerun make pilot-validate"
                )
        for name, digest in manifest.get("tools", {}).items():
            if hashlib.sha256(Path(name).read_bytes()).hexdigest() != digest:
                raise ValueError(f"Executable changed since validation: {name}")
        print("Evidence inputs, executables and published outputs match their hashes.")
        return
    sidecar = pd.read_csv(pilot_orthologs.SIDECAR)
    baseline = pd.read_csv(pilot_fimo.SCORES)
    regional_path = Path(f"{RESCUE_PREFIX}.resolved.csv")
    regional = pd.read_csv(regional_path)
    source = pd.read_csv(config.HG38_ELEMENTS)
    if len(sidecar) != pilot_orthologs.EXPECTED_PILOT_COUNT or not sidecar.element_id.is_unique:
        raise ValueError("Expected the complete unique 390-element pilot")
    if not sidecar.min_match.eq(config.MINMATCH_CROSS_SPECIES).all():
        raise ValueError("Sidecar minMatch differs from config; rebuild pilot-orthologs")
    if set(regional.element_id) != set(sidecar.loc[sidecar.call.eq("ABSENT"), "element_id"]):
        raise ValueError("Rescue cohort differs from corrected sidecar")
    if (
        not sidecar.source_build.eq(config.GENOME_BUILD).all()
        or not sidecar.target_build.eq("mm39").all()
    ):
        raise ValueError("Unexpected pilot assemblies")
    # Refuse a stale sidecar even if IDs and high-level mapping counts happen to agree.
    for row in sidecar.itertuples(index=False):
        a, b = config.MOTIF_SPAN
        expected = (
            (row.hg38_element_start + a, row.hg38_element_start + b + 1)
            if row.hg38_element_strand == "+"
            else (row.hg38_element_end - b - 1, row.hg38_element_end - a)
        )
        if expected != (row.hg38_motif_start, row.hg38_motif_end):
            raise ValueError("Sidecar motif coordinates are stale; rebuild pilot-orthologs")
    sizes = pilot.parse_sizes((pilot.METADATA / "mm39.chrom.sizes").read_text())
    expected_windows = pilot_fimo.scoring_windows(sidecar, sizes).set_index("element_id")
    if set(baseline.element_id) != set(sidecar.element_id) or not baseline.element_id.is_unique:
        raise ValueError("Baseline IDs differ from the pilot")
    old = baseline.set_index("element_id")
    for row in sidecar.itertuples(index=False):
        if old.loc[row.element_id, "call"] != row.call:
            raise ValueError("Baseline mapping calls are stale")
        if row.call == "MAPPED":
            for col in ("window_chrom", "window_start", "window_end", "window_strand"):
                if old.loc[row.element_id, col] != expected_windows.loc[row.element_id, col]:
                    raise ValueError("Baseline scoring windows are stale")
    for path in (pilot_orthologs.SIDECAR, config.MOTIF_PATH, pilot_fimo.GENOME):
        if local_md5(path) not in pilot_fimo.REPORT.read_text():
            raise ValueError("Baseline report input hashes are stale")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    reference = check_reference(OUTPUT / "reference")
    human = human_scores(source, sidecar, config.HG38_TWOBIT, OUTPUT / "human")
    core_intervals = sidecar[sidecar.call.eq("MAPPED")][
        [
            "element_id",
            "mm39_motif_chrom",
            "mm39_motif_start",
            "mm39_motif_end",
            "mm39_motif_strand",
        ]
    ].rename(columns=lambda c: c.removeprefix("mm39_motif_"))
    core_intervals[["start", "end"]] = core_intervals[["start", "end"]].astype(int)
    core = score_records(
        extract_intervals(core_intervals, pilot_fimo.GENOME, OUTPUT / "mouse_core"),
        OUTPUT / "mouse_core",
    )
    evidence = audit_mappings(
        mapping_requests(sidecar, regional, "mm39"), config.CHAIN_PATHS["mm39"]
    )
    quality = alignment_summary(sidecar, regional, evidence)
    result = assemble(sidecar, baseline, human, core, quality, regional)
    previous_path = SNAPSHOT / "data/processed/pilot_mouse/fimo_scores.csv"
    previous = pd.read_csv(previous_path) if previous_path.is_file() else None
    report = render(result, previous)
    tables = {
        "calls.csv": result,
        "human_scores.csv": human,
        "mouse_core_scores.csv": core,
        "alignment_blocks.csv": evidence,
        "alignment_summary.csv": quality,
        "reference_check.csv": reference,
    }
    if previous is not None:
        tables["comparison_to_previous.csv"] = previous[
            ["element_id", "call", "motif_status", "score"]
        ].merge(
            result[
                [
                    "element_id",
                    "mapping_call",
                    "mouse_window_status",
                    "mouse_window_score",
                    "human_core_status",
                    "alignment_status",
                    "comparison_status",
                    "comparison_reason",
                ]
            ],
            on="element_id",
            validate="one_to_one",
        )
    inputs = [
        pilot_orthologs.SIDECAR,
        pilot_fimo.SCORES,
        pilot_fimo.REPORT,
        regional_path,
        config.HG38_ELEMENTS,
        config.HG38_TWOBIT,
        pilot_fimo.GENOME,
        config.MOTIF_PATH,
        config.CHAIN_PATHS["mm39"],
        REFERENCE,
        *sorted((config.ROOT / "src").glob("*.py")),
    ]
    manifest = dict(inputs=fingerprints(inputs))
    executables = [shutil.which(name) for name in ("liftOver", "twoBitInfo", "twoBitToFa", "fimo")]
    if not all(executables):
        raise FileNotFoundError("Required genomics executables are missing from PATH")
    manifest["tools"] = {
        name: hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in executables
    }
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".publish-", dir=OUTPUT) as directory:
        replacements = []
        for name, frame in tables.items():
            temp = Path(directory) / name
            frame.to_csv(temp, index=False)
            replacements.append((temp, OUTPUT / name))
        temp = Path(directory) / REPORT.name
        temp.write_text(report)
        replacements.append((temp, REPORT))
        manifest["outputs"] = {
            str(dest.relative_to(config.ROOT)): hashlib.sha256(src.read_bytes()).hexdigest()
            for src, dest in replacements
        }
        for src, dest in replacements:
            os.replace(src, dest)
        temp = Path(directory) / "manifest.json"
        temp.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        os.replace(temp, manifest_path)
    print(result.comparison_status.value_counts().to_string(), flush=True)
    print(f"Report: {REPORT}", flush=True)


if __name__ == "__main__":
    main()
