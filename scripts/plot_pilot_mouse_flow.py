"""Recreate the pilot mapping/FIMO flow from current CSVs, without rerunning analysis.

Run from any directory: python scripts/plot_pilot_mouse_flow.py
"""

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as DrawingPath

ROOT = Path(__file__).resolve().parents[1]


def read_rows(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ids = [row["element_id"] for row in rows]
    if not rows or len(ids) != len(set(ids)):
        raise ValueError(f"Empty or duplicate element IDs: {path}")
    return {row["element_id"]: row for row in rows}


def ribbon(ax, x0, x1, top0, top1, width, color):
    bend = (x1 - x0) * 0.45
    vertices = [
        (x0, top0),
        (x0 + bend, top0),
        (x1 - bend, top1),
        (x1, top1),
        (x1, top1 + width),
        (x1 - bend, top1 + width),
        (x0 + bend, top0 + width),
        (x0, top0 + width),
        (x0, top0),
    ]
    codes = [
        DrawingPath.MOVETO,
        *([DrawingPath.CURVE4] * 3),
        DrawingPath.LINETO,
        *([DrawingPath.CURVE4] * 3),
        DrawingPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(DrawingPath(vertices, codes), facecolor=color, lw=0))


def main():
    mapped_rows = read_rows(ROOT / "data/processed/pilot_mouse/orthologs.mm39.csv")
    scores = read_rows(ROOT / "data/processed/pilot_mouse/fimo_scores.csv")
    if mapped_rows.keys() != scores.keys():
        raise ValueError("Mapping and scoring IDs differ")
    for element_id, row in mapped_rows.items():
        score = scores[element_id]
        if row["call"] != score["call"]:
            raise ValueError(f"Mapping/scoring calls disagree: {element_id}")
        allowed = {"present", "absent"} if row["call"] == "MAPPED" else {"unscored"}
        if score["motif_status"] not in allowed:
            raise ValueError(f"Unexpected score status: {element_id}")
    counts = Counter(row["call"] for row in mapped_rows.values())
    if set(counts) != {"MAPPED", "ABSENT", "UNDETERMINED"}:
        raise ValueError("Unexpected mapping categories")
    total = len(mapped_rows)
    mapped, absent, undetermined = (counts[k] for k in ("MAPPED", "ABSENT", "UNDETERMINED"))
    hits = Counter(row["motif_status"] for row in scores.values())
    present, negative = hits["present"], hits["absent"]
    if present + negative != mapped:
        raise ValueError("Scored total differs from mapped total")

    plt.rcParams.update(
        {"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none", "pdf.fonttype": 42}
    )
    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=0.035, right=0.98, top=0.83, bottom=0.21)
    fig.suptitle(
        "Mouse pilot: mapping and motif scans",
        x=0.055,
        y=0.96,
        ha="left",
        fontsize=18,
        fontweight="bold",
    )
    fig.text(0.055, 0.902, "Corrected motif coordinates · hg38 → mm39", color="#555555")
    ax.set(xlim=(0, 1), ylim=(total + 62, -42))
    ax.axis("off")
    purple, blue, amber, pink, green, lightblue = (
        "#352384",
        "#087BAF",
        "#E5A100",
        "#C96EAA",
        "#009D83",
        "#55B3DF",
    )
    ax.add_patch(Rectangle((0.025, 0), 0.065, total, color=purple, lw=0))
    ribbon(ax, 0.09, 0.44, 0, 0, mapped, blue)
    ribbon(ax, 0.09, 0.44, mapped, mapped + 18, absent, amber)
    ribbon(ax, 0.09, 0.44, mapped + absent, total + 34, undetermined, pink)
    ax.add_patch(Rectangle((0.44, 0), 0.30, present, color=green, lw=0))
    ax.add_patch(Rectangle((0.44, present + 3), 0.30, negative, color=lightblue, lw=0))
    ax.text(0.025, -13, f"chr8 RE1\nn = {total}", va="bottom")
    ax.text(0.34, -13, "LiftOver", ha="center", va="bottom")
    ax.text(0.59, -13, "FIMO\nof mapped", ha="center", va="bottom")
    ax.text(
        0.26,
        mapped / 2,
        f"Mapped\n{mapped} ({mapped / total:.1%})",
        ha="center",
        va="center",
        color="white",
    )
    ax.text(
        0.26,
        mapped + 9 + absent / 2,
        f"Neither track mapped\n{absent} ({absent / total:.1%})",
        ha="center",
        va="center",
        color="#252525",
    )
    ax.text(
        0.455,
        total + 34 + undetermined / 2,
        f"Undetermined\n{undetermined} ({undetermined / total:.1%})",
        va="center",
    )
    ax.text(0.755, present / 2, f"Motif present\n{present} ({present / mapped:.1%})", va="center")
    ax.text(
        0.59,
        present + 3 + negative / 2,
        f"Below motif threshold\n{negative} ({negative / mapped:.1%})",
        ha="center",
        va="center",
        color="#183E50",
    )
    fig.text(
        0.055,
        0.14,
        "Mapping percentages: all pilot elements. FIMO percentages: mapped elements only.",
        fontsize=10,
        color="#444444",
    )
    fig.text(
        0.055,
        0.10,
        "Raw window scans, not final retention/nonretention calls. "
        "Mapping failure does not prove loss.",
        fontsize=10,
        color="#444444",
    )
    fig.text(
        0.055,
        0.06,
        "Sources: orthologs.mm39.csv and fimo_scores.csv · MA0138.2 · score threshold ≥20.86",
        fontsize=9,
        color="#666666",
    )
    output = ROOT / "results/figures"
    output.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "svg", "pdf"):
        path = output / f"pilot_mouse_mapping_fimo_flow.{extension}"
        fig.savefig(path, dpi=300, facecolor="white")
        print(path.relative_to(ROOT))
    plt.close(fig)


if __name__ == "__main__":
    main()
