"""Known chain geometries expose defects that interval-level lift success hides."""

import pandas as pd
import pytest

from src.alignment_evidence import audit_mappings, chains
from src.pilot_evidence import alignment_summary


def request(start=90, end=111, target_start=90, target_end=111, strand="+", source_strand="+"):
    return dict(
        element_id="synthetic",
        track="motif",
        source_chrom="chrS",
        source_start=start,
        source_end=end,
        source_strand=source_strand,
        target_chrom="chrT",
        target_start=target_start,
        target_end=target_end,
        target_strand=strand,
    )


@pytest.mark.parametrize(
    "source_gap,target_gap,aligned,span",
    [
        (0, 0, 21, 21),
        (9, 0, 12, 12),
        (0, 24, 21, 45),
        (9, 24, 12, 36),
    ],
)
def test_indels_hidden_by_successful_interval_lifts(
    tmp_path, source_gap, target_gap, aligned, span
):
    path = tmp_path / "test.chain"
    path.write_text(
        f"chain 100 chrS 1000 + 0 {200 + source_gap} chrT 1000 + 0 {200 + target_gap} 1\n"
        f"100 {source_gap} {target_gap}\n100\n\n"
    )
    row = audit_mappings([request(target_end=90 + span)], path).iloc[0]
    assert row.aligned_bp == aligned
    assert row.source_unaligned_bp == 21 - aligned
    assert row.target_unaligned_bp == target_gap
    assert bool(row.complete) == (source_gap == target_gap == 0)


@pytest.mark.parametrize("source_strand,target_strand", [("+", "-"), ("-", "+")])
def test_reverse_target_and_reverse_oligo_orientation(tmp_path, source_strand, target_strand):
    path = tmp_path / "test.chain"
    path.write_text("chain 100 chrS 1000 + 0 200 chrT 1000 - 0 200 1\n200\n\n")
    row = audit_mappings(
        [
            request(
                target_start=889, target_end=910, strand=target_strand, source_strand=source_strand
            )
        ],
        path,
    ).iloc[0]
    assert row.complete
    assert row.aligned_bp == 21 and row.chain_matches == 1


def test_trimmed_anchor_keeps_inner_boundary_distance(tmp_path):
    path = tmp_path / "test.chain"
    path.write_text("chain 100 chrS 1000 + 0 100 chrT 1000 + 50 150 1\n100\n\n")
    row = audit_mappings([request(start=0, end=150, target_start=50, target_end=150)], path).iloc[0]
    assert row.source_right_trim == 50 and row.source_left_trim == 0
    assert not row.complete


def test_duplicate_chains_are_ambiguous_and_stale_mapping_fails(tmp_path):
    path = tmp_path / "test.chain"
    path.write_text(
        "".join(f"chain 100 chrS 1000 + 0 200 chrT 1000 + 0 200 {i}\n200\n\n" for i in (1, 2))
    )
    row = audit_mappings([request()], path).iloc[0]
    assert row.chain_matches == 2 and row.chain_ids == "1;2" and not row.complete
    assert pd.isna(row.aligned_bp)
    with pytest.raises(ValueError, match="No chain reproduces"):
        audit_mappings([request(target_end=112)], path)


def test_truncated_chain_cannot_certify_an_alignment(tmp_path):
    path = tmp_path / "bad.chain"
    path.write_text("chain 100 chrS 1000 + 0 200 chrT 1000 + 0 200 1\n100 0 0\n")
    with pytest.raises(ValueError, match="Truncated"):
        list(chains(path, {"chrS"}))


def test_different_parent_and_motif_chains_are_not_intact(tmp_path):
    path = tmp_path / "test.chain"
    path.write_text("chain 100 chrS 1000 + 0 200 chrT 1000 + 0 200 1\n200\n\n")
    evidence = audit_mappings(
        [
            request(),
            request(start=0, end=200, target_start=0, target_end=200) | dict(track="element"),
        ],
        path,
    )
    evidence.loc[evidence.track.eq("element"), "chain_ids"] = "different"
    sidecar = pd.DataFrame([dict(element_id="synthetic", call="MAPPED")])
    regional = pd.DataFrame(columns=["element_id", "valid_pair"])
    quality = alignment_summary(sidecar, regional, evidence)
    assert quality.loc[0, "alignment_status"] == "DISCORDANT_CHAINS"
    assert not quality.loc[0, "chain_consistent"]
