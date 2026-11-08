"""Protect the mapping-availability distinction and padded, oriented extraction."""

import random
import shutil
import struct

import pytest
from Bio import SeqIO
from Bio.Seq import Seq

from src import pilot_orthologs as orth


def interval(start=100, end=300, strand="+", chrom="chr1"):
    return dict(chrom=chrom, start=start, end=end, strand=strand, failure_reason="")


def tracks():
    """Each logical mapping outcome occurs once, including motif-only recovery."""
    source = {name: interval() for name in ("both", "element_only", "motif_only", "neither")}
    windows = orth.motif_windows(source)
    element = {
        "mapped": {name: interval() for name in ("both", "element_only")},
        "unmapped": {
            name: source[name] | {"failure_reason": "Partially deleted in new"}
            for name in ("motif_only", "neither")
        },
    }
    motif = {
        "mapped": {name: interval(191, 212) for name in ("both", "motif_only")},
        "unmapped": {
            name: windows[name] | {"failure_reason": "Deleted in new"}
            for name in ("element_only", "neither")
        },
    }
    return source, windows, element, motif


def test_window_coordinates_both_strands_and_width_validation():
    source = {"plus": interval(), "minus": interval(strand="-")}
    windows = orth.motif_windows(source)
    assert (windows["plus"]["start"], windows["plus"]["end"]) == (191, 212)
    assert (windows["minus"]["start"], windows["minus"]["end"]) == (188, 209)
    assert all(row["end"] - row["start"] == 21 for row in windows.values())
    source["plus"]["end"] = 301
    with pytest.raises(ValueError, match="200 bp"):
        orth.motif_windows(source)
    with pytest.raises(ValueError, match="21 bp"):
        orth.motif_windows({}, (90, 111))


def test_all_four_track_outcomes_remain_in_three_calls():
    audit = orth.classify_tracks(*tracks()).set_index("element_id")
    assert audit.call.value_counts().to_dict() == {"UNDETERMINED": 2, "MAPPED": 1, "ABSENT": 1}
    assert audit.loc["motif_only", "call_reason"] == "motif_only_element_unmapped"
    assert audit.loc["motif_only", "mm39_motif_width"] == 21
    assert audit.loc["element_only", "mm39_element_width"] == 200
    assert audit.loc[["element_only", "motif_only"], "exclude_mouse_downstream"].all()
    assert not audit.loc[["element_only", "motif_only"], "include_in_fasta"].any()
    assert audit.loc["neither", "mm39_element_width"] == ""


@pytest.mark.parametrize("change", [{"chrom": "chr2"}, {"strand": "-"}, {"start": 50, "end": 71}])
def test_both_tracks_must_describe_the_same_locus(change):
    source, windows, element, motif = tracks()
    motif["mapped"]["both"].update(change)
    with pytest.raises(ValueError, match="Discordant"):
        orth.classify_tracks(source, windows, element, motif)


def test_missing_ids_never_vanish():
    source, windows, element, motif = tracks()
    del motif["unmapped"]["neither"]
    with pytest.raises(AssertionError, match="partition"):
        orth.classify_tracks(source, windows, element, motif)


def test_padding_clipping_and_empty_nonmapped_fields():
    source, windows, element, motif = tracks()
    element["mapped"]["both"] = interval(0, 200)
    motif["mapped"]["both"] = interval(91, 112)
    audit = orth.classify_tracks(source, windows, element, motif)
    records, padded = orth.padded_records(audit, {"chr1": 210})
    assert records["both"]["start"] == 0 and records["both"]["end"] == 210
    row = padded.set_index("element_id").loc["both"]
    assert row.padding_left == 0 and row.padding_right == 10 and row.padding_clipped
    assert padded.loc[padded.call.ne("MAPPED"), "padded_start"].eq("").all()


@pytest.mark.integration
def test_real_extraction_orientation_length_flags_and_replay(tmp_path):
    if not shutil.which("twoBitToFa"):
        pytest.skip("Requires UCSC twoBitToFa")
    source, windows, element, motif = tracks()
    # A short, reverse-oriented mapped span still carries a mapped motif window.
    # It must be flagged rather than dropped from the MAPPED FASTA.
    element["mapped"]["both"] = interval(100, 115, "-")
    motif["mapped"]["both"] = interval(101, 112, "-")
    audit = orth.classify_tracks(source, windows, element, motif)
    rng = random.Random(2)
    dna = "".join(rng.choices("ACGT", k=500))
    name = b"chr1"
    genome = tmp_path / "synthetic.2bit"
    header = struct.pack("<IIII", 0x1A412743, 0, 1, 0)
    index = bytes([len(name)]) + name + struct.pack("<I", 16 + 1 + len(name) + 4)
    record = struct.pack("<IIIIII", len(dna), 0, 1, 0, len(dna), 0)
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytes(
        sum(code[base] << (6 - 2 * j) for j, base in enumerate(dna[i : i + 4]))
        for i in range(0, len(dna), 4)
    )
    genome.write_bytes(header + index + record + packed)
    candidate = tmp_path / "candidate.fa"
    result = orth.extract(audit, {"chr1": len(dna)}, genome, tmp_path, candidate)
    records = list(SeqIO.parse(candidate, "fasta"))
    assert len(records) == 1 and records[0].id == "both"
    assert str(records[0].seq) == str(Seq(dna[85:130]).reverse_complement())
    row = result.set_index("element_id").loc["both"]
    assert row.extracted_length == 45 and row.sequence_under_50bp
    assert row.width_outside_160_240
    before = candidate.read_bytes(), result.to_csv(index=False)
    repeated = orth.extract(audit, {"chr1": len(dna)}, genome, tmp_path, candidate)
    assert before == (candidate.read_bytes(), repeated.to_csv(index=False))
    unmapped = audit.loc[audit.call.ne("MAPPED")].reset_index(drop=True)
    result = orth.extract(unmapped, {"chr1": len(dna)}, genome, tmp_path, candidate)
    assert candidate.read_bytes() == b"" and result.extracted_length.eq("").all()
