"""Synthetic tests for the genome gate, strict BED writing, and seeded regression."""

import random
import shutil
import struct
import sys

import pandas as pd
import pytest
from Bio.Seq import Seq

from src import pilot


def elements(count=320):
    """Canonical chr8 records plus contrast rows that must never enter the pilot."""
    frame = pd.DataFrame(
        {
            "element_id": [f"synthetic_{i:03}" for i in range(count + 20)],
            "project": ["motif"] * count + ["noMotif"] * 20,
            "working_build": "hg38",
            "chrom": "chr8",
            "start": [i * 200 for i in range(count + 20)],
            "end": [(i + 1) * 200 for i in range(count + 20)],
            "strand": "+",
            "sequence": "A" * 200,
        }
    )
    return frame


def test_selection_and_widths(tmp_path):
    frame = elements()
    selected = pilot.choose_pilot(frame)
    assert len(selected) == 320 and selected.project.eq("motif").all()
    pd.testing.assert_frame_equal(selected, pilot.choose_pilot(frame.iloc[::-1]))
    path = tmp_path / "pilot.bed"
    pilot.write_bed(selected, path)
    assert all(len(line.split("\t")) == 6 for line in path.read_text().splitlines())
    frame.loc[0, "end"] = 201
    with pytest.raises(ValueError, match="synthetic_000"):
        pilot.choose_pilot(frame)
    frame.loc[0, "end"] = 200
    frame["working_build"] = "hg19"
    with pytest.raises(ValueError, match="must be hg38"):
        pilot.choose_pilot(frame)


@pytest.mark.parametrize("count", [149, 601])
def test_population_guard_suggests_without_switching(count):
    frame = elements(count)
    alternative = elements(350)
    alternative["chrom"] = "chr7"
    with pytest.raises(ValueError, match="Consider chr7"):
        pilot.choose_pilot(pd.concat([frame, alternative]))


def test_genome_failure_never_selects_or_publishes(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot, "BED", tmp_path / "pilot.bed")
    monkeypatch.setattr(pilot, "REPORT", tmp_path / "report.md")
    pilot.BED.write_text("stale result from previous run\n")

    def failure(assembly, executable):
        return dict(
            assembly=assembly,
            file_bytes=None,
            sequence_count=None,
            total_bases=None,
            expected_sequence_count=None,
            expected_total_bases=None,
            observed_md5="",
            expected_md5="",
            status="FAIL",
            errors=["synthetic missing genome"],
            sizes_url="synthetic",
        )

    monkeypatch.setattr(pilot, "verify_genome", failure)

    def forbidden(*args):
        pytest.fail("Selection must not run after a genome failure")

    monkeypatch.setattr(pilot, "choose_pilot", forbidden)
    monkeypatch.setattr(sys, "argv", ["pilot"])
    with pytest.raises(SystemExit) as exc:
        pilot.main()
    assert exc.value.code == 1
    assert not pilot.BED.exists()
    assert "NOT RUN" in pilot.REPORT.read_text()


@pytest.mark.integration
def test_actual_tools_replay_and_corrupted_genome(tmp_path, monkeypatch):
    """Run all six checks on small synthetic references, then detect silent DNA corruption.

    A payload mutation leaves every contig length intact. Requiring the publisher
    MD5 catches the error that an apparently successful twoBitInfo cannot catch.
    """
    info = pilot.info_executable()
    if not info or not shutil.which("twoBitToFa"):
        pytest.skip("Requires UCSC twoBitInfo and twoBitToFa")
    monkeypatch.setattr(pilot, "info_executable", lambda: info)
    genomes, metadata = tmp_path / "genomes", tmp_path / "metadata"
    genomes.mkdir()
    metadata.mkdir()
    monkeypatch.setattr(pilot, "GENOMES", genomes)
    monkeypatch.setattr(pilot, "METADATA", metadata)
    monkeypatch.setattr(pilot, "EVIDENCE", tmp_path / "evidence")
    monkeypatch.setattr(pilot, "REPORT", tmp_path / "report.md")
    monkeypatch.setattr(pilot, "BED", tmp_path / "pilot.bed")
    monkeypatch.setattr(pilot.config, "HG38_TWOBIT", genomes / "hg38.2bit")
    monkeypatch.setattr(pilot.config, "HG38_ELEMENTS", tmp_path / "elements.csv")
    # Provenance-helper behavior is covered separately; these are explicit fixtures,
    # never real downloaded genomes or claims about publisher data.
    monkeypatch.setattr(pilot, "download_small_file", lambda url, path, **kwargs: path)
    rng = random.Random(3)
    dna = "".join(rng.choices("ACGT", k=68000))
    name = b"chr8"
    header = struct.pack("<IIII", 0x1A412743, 0, 1, 0)
    index = bytes([len(name)]) + name + struct.pack("<I", 16 + 1 + len(name) + 4)
    record = struct.pack("<IIIIII", len(dna), 0, 1, 0, len(dna), 0)
    code = {"T": 0, "C": 1, "A": 2, "G": 3}
    packed = bytes(
        sum(code[base] << (6 - 2 * j) for j, base in enumerate(dna[i : i + 4]))
        for i in range(0, len(dna), 4)
    )
    for assembly in pilot.ASSEMBLIES:
        path = genomes / f"{assembly}.2bit"
        path.write_bytes(header + index + record + packed)
        (metadata / f"{assembly}.chrom.sizes").write_text(f"chr8\t{len(dna)}\n")
        (metadata / f"{assembly}.md5sum.txt").write_text(
            f"{pilot.local_md5(path)}  {assembly}.2bit\n"
        )
    frame = elements()
    frame["strand"] = ["-" if i % 2 else "+" for i in range(len(frame))]
    frame["sequence"] = [
        str(Seq(dna[r.start : r.end]).reverse_complement())
        if r.strand == "-"
        else dna[r.start : r.end]
        for r in frame.itertuples()
    ]
    frame.to_csv(pilot.config.HG38_ELEMENTS, index=False)
    monkeypatch.setattr(sys, "argv", ["pilot"])
    pilot.main()
    assert len(pilot.BED.read_text().splitlines()) == 320
    assert "10/10" in pilot.REPORT.read_text()
    artifacts = [pilot.BED, pilot.REPORT, *pilot.EVIDENCE.iterdir()]
    before = {path: path.read_bytes() for path in artifacts}
    pilot.main()
    assert all(path.read_bytes() == content for path, content in before.items())
    path = genomes / "mm39.2bit"
    damaged = bytearray(path.read_bytes())
    damaged[-1] ^= 1
    path.write_bytes(damaged)
    observed = pilot.verify_genome("mm39", info)
    assert observed["total_bases"] == observed["expected_total_bases"] == len(dna)
    assert observed["status"] == "FAIL"
    assert any("MD5 differs" in error for error in observed["errors"])
