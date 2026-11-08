"""Test step-specific gates, generic lift accounting, and sweep width boundaries."""

import shutil

import pytest

from src import pilot, pilot_mouse
from src.lift import lift_intervals


@pytest.mark.parametrize(
    "step,expected",
    [
        ("pilot", ("hg38",)),
        ("mouse", ("hg38", "mm39")),
        ("multi_species", ("hg38", "rheMac10", "mm39", "canFam6", "bosTau9", "monDom5")),
    ],
)
def test_step_verifies_only_its_dependencies(step, expected, monkeypatch):
    called = []

    def verify(assembly, executable):
        called.append(assembly)
        return {"assembly": assembly, "status": "FAIL" if assembly == "mm39" else "PASS"}

    monkeypatch.setattr(pilot, "verify_genome", verify)
    monkeypatch.setattr(pilot, "info_executable", lambda: "synthetic-twoBitInfo")
    results = pilot.verify_step(step)
    assert tuple(called) == expected
    assert len(results) == len(expected)
    # Scoping must never turn a required assembly failure into a pass.
    if "mm39" in expected:
        assert any(row["status"] == "FAIL" for row in results)
    else:
        assert all(row["status"] == "PASS" for row in results)


def test_width_boundary_and_empty_run():
    widths = [159, 160, 200, 240, 241]
    run = dict(
        min_match=0.2,
        mapped={str(i): {"start": 10, "end": 10 + width} for i, width in enumerate(widths)},
        unmapped={"failed": {}},
    )
    summary = pilot_mouse.summarize(run, 6)
    assert summary["median_width"] == 200
    assert summary["width_deviation_gt20pct"] == 2
    assert summary["width_within20pct"] == 3
    assert summary["lift_rate"] == 5 / 6
    run["unmapped"].update(run["mapped"])
    run["mapped"] = {}
    assert pilot_mouse.summarize(run, 6)["median_width"] is None
    with pytest.raises(AssertionError):
        pilot_mouse.summarize(run, 7)


@pytest.mark.integration
def test_generic_cross_species_lift_retains_failures(tmp_path):
    """A partial synthetic chain separates coverage thresholds without any genome fetch."""
    if not shutil.which("liftOver"):
        pytest.skip("Requires UCSC liftOver")
    bed = tmp_path / "input.bed"
    bed.write_text(
        "chr8\t0\t200\tsynthetic_partial\t0\t+\nchr8\t400\t600\tsynthetic_absent\t0\t-\n"
    )
    chain = tmp_path / "synthetic.chain"
    chain.write_text("chain 1000 chr8 1000 + 0 100 chr1 2000 + 100 200 1\n100\n\n")
    strict = lift_intervals(bed, "hg38", "mm39", chain, tmp_path / "strict", min_match=0.95)
    loose = lift_intervals(bed, "hg38", "mm39", chain, tmp_path / "loose", min_match=0.5)
    assert len(strict["unmapped"]) == 2 and not strict["mapped"]
    assert len(loose["mapped"]) == len(loose["unmapped"]) == 1
    mapped = loose["mapped"]["synthetic_partial"]
    assert mapped["end"] - mapped["start"] == 100 and mapped["chrom"] == "chr1"
    assert loose["unmapped"]["synthetic_absent"]["failure_reason"]
    assert "-minMatch=0.5" in loose["command"]
    assert "-multiple" not in loose["command"]
    original = lift_intervals(bed, "hg19", "hg38", chain, tmp_path / "default")
    assert not any(arg.startswith("-minMatch") for arg in original["command"])
    assert original["mapped_path"].name == "lifted.hg38.bed"
    assert original["unmapped_path"].name == "unmapped.hg19.bed"
    source = pilot_mouse.parse_bed(bed)
    audit = tmp_path / "audit.csv"
    pilot_mouse.write_audit(source, loose, audit)
    assert "synthetic_absent" in audit.read_text() and "failed-to-lift" in audit.read_text()


def test_failed_mouse_gate_cannot_run_liftover(tmp_path, monkeypatch):
    monkeypatch.setattr(pilot_mouse, "OUTPUT", tmp_path / "outputs")
    monkeypatch.setattr(pilot_mouse.pilot, "verify_step", lambda step: [{"status": "FAIL"}])
    monkeypatch.setattr(pilot_mouse, "render", lambda *args, **kwargs: "")
    monkeypatch.setattr("sys.argv", ["pilot_mouse"])

    def forbidden(*args, **kwargs):
        pytest.fail("Failed dependency gate must prevent liftOver")

    monkeypatch.setattr(pilot_mouse, "lift_intervals", forbidden)
    with pytest.raises(SystemExit) as error:
        pilot_mouse.main()
    assert error.value.code == 1
