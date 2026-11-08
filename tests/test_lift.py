"""Test the actual failure modes that could silently corrupt the assembly handoff.

All fixtures are synthetic and explicitly labeled; no network or production data
is required. The integration fixture exercises installed UCSC programs against a
small local chain and a masked .2bit written in the published binary format.
"""

import gzip
import random
import shutil
import struct
import zipfile

import pandas as pd
import pytest
from Bio.Seq import Seq

from src import fetch_lift_inputs, lift


def table(sequences, starts=None):
    """Create explicitly hg19 synthetic records with independent element IDs."""
    starts = starts if starts is not None else [i * 300 for i in range(len(sequences))]
    return pd.DataFrame(
        {
            "element_id": [f"synthetic_{i:02}" for i in range(len(sequences))],
            "chrom": "chr1",
            "start": starts,
            "end": [start + 200 for start in starts],
            "sequence": sequences,
            "source_build": "hg19",
            "project": "motif",
        }
    )


def write_repeats(path, intervals):
    """Write schema-correct rmsk rows; repeat scores/classes are irrelevant to overlap."""
    with gzip.open(path, "wt") as handle:
        for chrom, start, end in intervals:
            fields = [
                0,
                0,
                0,
                0,
                0,
                chrom,
                start,
                end,
                0,
                "+",
                "synthetic",
                "LINE",
                "synthetic",
                0,
                end - start,
                0,
                0,
            ]
            handle.write("\t".join(map(str, fields)) + "\n")


def test_width_and_source_build_rejected():
    frame = table(["A" * 200])
    frame.loc[0, "end"] = 201
    with pytest.raises(ValueError, match="synthetic_00,0,201,201"):
        lift.validate_input(frame)
    frame.loc[0, "end"] = 200
    frame.loc[0, "source_build"] = "hg38"
    with pytest.raises(ValueError, match="Expected source_build=hg19"):
        lift.validate_input(frame)
    frame.loc[0, "source_build"] = "hg19"
    frame["start"] = 0.5
    with pytest.raises(ValueError, match="refusing to round"):
        lift.validate_input(frame)


def test_aliases_are_explicit_not_fuzzy(tmp_path, monkeypatch):
    archive = tmp_path / "oligos.zip"
    source = pd.DataFrame(
        {
            "silencer": ["first", "second"],
            "chromosome": ["1", "1"],
            "start": [0, 300],
            "end": [200, 500],
            "project": ["motif", "noMotif"],
        }
    )
    monkeypatch.setattr(lift, "load_evidence", lambda: (source, None))
    monkeypatch.setattr(lift, "TABLE_S2", archive)
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr(
            lift.TABLE_S2_MEMBER,
            "oligo_name\tsequence\tlibrary\n"
            f"(first;explicit_variant)\t{'A' * 200}\tWG_REST\n"
            f"second\t{'C' * 200}\tWG_REST\n",
        )
    observed = lift.published_elements()
    assert observed.sequence_join.tolist() == ["explicit_alias", "exact_name"]
    source.loc[0, "silencer"] = "firs"
    with pytest.raises(ValueError, match="expected one explicit oligo alias, got 0"):
        lift.published_elements()


def test_adapters_are_independent_of_reference():
    prefix, suffix = "GATTACAG", "CTGATCGA"
    frame = table([prefix + "A" * 200 + suffix, prefix + "C" * 200 + suffix, "G" * 200])
    prepared = lift.prepare_cores(frame)
    assert prepared.sequence.tolist() == ["A" * 200, "C" * 200, "G" * 200]
    assert prepared.adapter_5p.tolist() == [prefix, prefix, ""]
    assert prepared.adapter_status.tolist() == ["shared_flanks_inferred"] * 2 + ["none_present"]
    singleton = lift.prepare_cores(frame.iloc[:1])
    assert singleton.iloc[0].adapter_status == "unresolved_long_oligo"
    explicit = lift.prepare_cores(frame.iloc[:1], prefix, suffix)
    assert explicit.iloc[0].sequence == "A" * 200
    # Homopolymer flanks have many possible split points: never guess one.
    ambiguous = lift.prepare_cores(table(["A" * 216, "A" * 8 + "C" * 190 + "A" * 18]))
    assert ambiguous.adapter_status.eq("unresolved_long_oligo").all()


@pytest.mark.parametrize(
    "reference,kind,differences",
    [
        ("A" * 199, "length_differs", ""),
        ("C" + "A" * 199, "one_base_difference", 1),
        ("CC" + "A" * 198, "multiple_base_differences", 2),
        ("C" * 200, "high_divergence_or_unrelated", 200),
        ("N" + "A" * 199, "ambiguous_bases", 1),
    ],
)
def test_mismatch_categories(reference, kind, differences):
    result = lift.verify_sequence("A" * 200, reference)
    assert not result["verified"]
    assert result["mismatch_type"] == kind
    assert result["differences"] == differences


def test_case_orientation_and_ambiguity():
    sequence = "A" * 93 + "C" * 107
    result = lift.verify_sequence(sequence.lower(), str(Seq(sequence).reverse_complement()))
    assert result["verified"] and result["best_orientation"] == "-"
    assert not lift.verify_sequence("N" * 200, "N" * 200)["verified"]


def test_unmapped_reasons_and_duplicate_detection(tmp_path):
    path = tmp_path / "unmapped.bed"
    path.write_text(
        "#Deleted in new\nchr1\t0\t200\tsynthetic_a\t0\t+\n"
        "#Partially deleted in new\nchr1\t300\t500\tsynthetic_b\t0\t-\n"
    )
    records = lift.parse_bed(path, unmapped=True)
    assert records["synthetic_b"]["failure_reason"] == "Partially deleted in new"
    with path.open("a") as handle:
        handle.write("#Deleted in new\nchr1\t0\t200\tsynthetic_a\t0\t+\n")
    with pytest.raises(ValueError, match="Multiple output records"):
        lift.parse_bed(path, unmapped=True)


def test_repeat_overlap_matches_brute_force(tmp_path):
    rng = random.Random(9)
    intervals = [("chr1", 50, 100), ("chr1", 60, 70), ("chr1", 500, 700)]
    intervals += [
        ("chr1", start, start + rng.randrange(1, 300))
        for start in [rng.randrange(0, 10000) for _ in range(100)]
    ]
    elements = table(["A" * 200] * 100, starts=[rng.randrange(0, 10000) for _ in range(100)])
    path = tmp_path / "rmsk.gz"
    write_repeats(path, intervals)
    observed, count = lift.repeat_overlap(elements, path)
    expected = [
        any(start < row.end and end > row.start for _, start, end in intervals)
        for row in elements.itertuples()
    ]
    assert observed.tolist() == expected and count == len(intervals)
    write_repeats(path, [("chr1", 200, 300)])
    observed, _ = lift.repeat_overlap(table(["A" * 200] * 3, [0, 100, 300]), path)
    assert observed.tolist() == [False, True, False]


def write_twobit(path, dna):
    """Tiny .2bit fixture with ALL bases soft-masked; -noMask must recover uppercase."""
    name = b"chr1"
    offset = 16 + 1 + len(name) + 4
    header = struct.pack("<IIII", 0x1A412743, 0, 1, 0)
    index = bytes([len(name)]) + name + struct.pack("<I", offset)
    record = struct.pack("<IIIIII", len(dna), 0, 1, 0, len(dna), 0)
    encoding = {"T": 0, "C": 1, "A": 2, "G": 3}
    padded = dna + "T" * ((-len(dna)) % 4)
    packed = bytes(
        sum(encoding[base] << (6 - 2 * j) for j, base in enumerate(padded[i : i + 4]))
        for i in range(0, len(padded), 4)
    )
    path.write_bytes(header + index + record + packed)


@pytest.mark.integration
def test_real_ucsc_all_buckets_orientation_and_replay(tmp_path):
    if not all(shutil.which(tool) for tool in ("liftOver", "twoBitToFa")):
        pytest.skip("Requires installed UCSC liftOver and twoBitToFa")
    rng = random.Random(17)
    dna = "".join(rng.choices("ACGT", k=7000))
    genome, chain, repeats = [
        tmp_path / name for name in ("synthetic.2bit", "test.chain", "rmsk.gz")
    ]
    write_twobit(genome, dna)
    chain.write_text("chain 1000 chr1 7000 + 0 6000 chr1 7000 + 10 6010 1\n6000\n\n")
    sequences = [dna[i * 300 + 10 : i * 300 + 210] for i in range(21)]
    sequences[0] = str(Seq(sequences[0]).reverse_complement())
    sequences[1] = ("A" if sequences[1][0] != "A" else "C") + sequences[1][1:]
    frame = table(sequences)
    frame["strand"] = "+"
    frame.loc[0, "strand"] = "-"
    write_repeats(repeats, [("chr1", 0, 1), ("chr1", 6000, 6200)])
    output, report = tmp_path / "verified.csv", tmp_path / "report.md"
    audit = lift.convert(frame, chain, genome, repeats, tmp_path / "audit", output, report)
    assert audit.bucket.value_counts().to_dict() == {
        "lifted+verified": 19,
        "lifted+mismatch": 1,
        "failed-to-lift": 1,
    }
    observed = pd.read_csv(output)
    assert observed.iloc[0].strand == "-"
    assert observed.iloc[0].sequence == sequences[0]
    assert observed.iloc[0].start == 10 and observed.iloc[0].hg19_start == 0
    assert "one_base_difference" in report.read_text()
    assert "LARGE" in report.read_text()
    paths = [output, report, *sorted((tmp_path / "audit").iterdir())]
    before = {path: path.read_bytes() for path in paths}
    lift.convert(frame.iloc[::-1], chain, genome, repeats, tmp_path / "audit", output, report)
    assert all(path.read_bytes() == value for path, value in before.items())
    # With half the source unmappable, the same defaults must stop before extraction.
    with pytest.raises(RuntimeError, match="LOW LIFT SUCCESS"):
        lift.convert(frame.iloc[[0, 20]], chain, genome, repeats, tmp_path / "low", output, report)
    assert not output.exists()


def test_resumable_download_provenance_and_checksum(tmp_path, monkeypatch):
    data = tmp_path / "data"
    data.mkdir()
    ledger = data / "README.md"
    ledger.write_text("existing rows\n")
    monkeypatch.setattr(fetch_lift_inputs.config, "ROOT", tmp_path)
    monkeypatch.setattr(fetch_lift_inputs.config, "DATA_DIR", data)
    path = data / "file"
    partial = data / "file.part"
    partial.write_bytes(b"test")
    digest = fetch_lift_inputs.local_md5(partial)
    # Complete partial files are validated without another request.
    fetch_lift_inputs.fetch("https://example.invalid/file", path, 4, digest)
    first = ledger.read_text()
    assert first.startswith("existing rows\n") and digest in first
    fetch_lift_inputs.fetch("https://example.invalid/file", path, 4, digest)
    assert ledger.read_text() == first
    path.write_bytes(b"oops")
    with pytest.raises(ValueError, match="MD5"):
        fetch_lift_inputs.fetch("https://example.invalid/file", path, 4, digest)
