"""Recover base-level evidence behind interval liftOver results.

A successful BED lift is a bounding interval, not a guarantee that every source
base aligns. Chain blocks use source genomic order; negative target coordinates
are converted back to forward-genome coordinates before comparison. No species
names, biological loss rules, or sequence identity assumptions live here.
"""

import bisect
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd


def chains(path: Path, chromosomes: set[str]):
    """Stream selected source chromosomes with validated chain endpoints."""
    opener = gzip.open if path.suffix == ".gz" else open
    header, blocks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "chain":
                if header is not None:
                    raise ValueError("Unterminated chain before next header")
                if len(fields) != 13 or fields[4] != "+" or fields[9] not in {"+", "-"}:
                    raise ValueError("Malformed chain header or unsupported source strand")
                header = dict(
                    score=int(fields[1]),
                    chrom=fields[2],
                    source_size=int(fields[3]),
                    start=int(fields[5]),
                    end=int(fields[6]),
                    target_chrom=fields[7],
                    target_size=int(fields[8]),
                    strand=fields[9],
                    target_start=int(fields[10]),
                    target_end=int(fields[11]),
                    chain_id=fields[12],
                )
                source, target = header["start"], header["target_start"]
                blocks = []
            else:
                if header is None or len(fields) not in {1, 3}:
                    raise ValueError("Malformed chain block")
                size = int(fields[0])
                if size <= 0:
                    raise ValueError("Chain block must have positive size")
                if header["chrom"] in chromosomes:
                    blocks.append((source, source + size, target, target + size))
                source += size
                target += size
                if len(fields) == 3:
                    ds, dt = map(int, fields[1:])
                    if min(ds, dt) < 0:
                        raise ValueError("Negative chain gap")
                    source += ds
                    target += dt
                else:
                    if source != header["end"] or target != header["target_end"]:
                        raise ValueError("Chain blocks do not reach declared endpoints")
                    if blocks:
                        yield header, blocks
                    header = None
        if header is not None:
            raise ValueError("Truncated chain")


def intersection(
    header: dict, blocks: list, start: int, end: int, ends: list[int] | None = None
) -> dict | None:
    """Project only aligned bases; expose trimming and gaps instead of filling them."""
    if not 0 <= start < end:
        raise ValueError("Invalid source interval")
    ends = ends if ends is not None else [block[1] for block in blocks]
    pieces = []
    for index in range(bisect.bisect_right(ends, start), len(blocks)):
        a, b, c, _ = blocks[index]
        if a >= end:
            break
        lo, hi = max(a, start), min(b, end)
        if lo >= hi:
            continue
        x, y = c + lo - a, c + hi - a
        if header["strand"] == "-":
            x, y = header["target_size"] - y, header["target_size"] - x
        pieces.append((lo, hi, x, y))
    if not pieces:
        return None
    aligned = sum(b - a for a, b, _, _ in pieces)
    target_start = min(x for _, _, x, _ in pieces)
    target_end = max(y for _, _, _, y in pieces)
    return dict(
        chain_id=header["chain_id"],
        target_chrom=header["target_chrom"],
        target_start=target_start,
        target_end=target_end,
        strand=header["strand"],
        aligned_bp=aligned,
        source_unaligned_bp=end - start - aligned,
        target_unaligned_bp=target_end - target_start - aligned,
        source_left_trim=pieces[0][0] - start,
        source_right_trim=end - pieces[-1][1],
        block_count=len(pieces),
        complete=(aligned == end - start == target_end - target_start),
        blocks=json.dumps(pieces, separators=(",", ":")),
    )


def audit_mappings(requests: list[dict], chain: Path) -> pd.DataFrame:
    """Reproduce each published mapping and retain ambiguity instead of picking a chain.

    Requests contain element_id, track, source_{chrom,start,end,strand}, and
    target_{chrom,start,end,strand}. A missing match indicates stale inputs and
    fails; multiple matches remain explicitly ambiguous. Numeric block metrics
    are only published for a unique matching chain.
    """
    grouped = defaultdict(list)
    keys = set()
    matches = defaultdict(list)
    for req in requests:
        key = req["element_id"], req["track"]
        if key in keys:
            raise ValueError("Duplicate mapping-evidence request")
        keys.add(key)
        if req["source_strand"] not in {"+", "-"}:
            raise ValueError("Source orientation must be explicit")
        grouped[req["source_chrom"]].append(req)
    for header, blocks in chains(Path(chain), set(grouped)):
        ends = [block[1] for block in blocks]
        for req in grouped[header["chrom"]]:
            if req["source_end"] <= header["start"] or req["source_start"] >= header["end"]:
                continue
            strand = header["strand"]
            if req["source_strand"] == "-":
                strand = "+" if strand == "-" else "-"
            if req["target_chrom"] != header["target_chrom"] or req["target_strand"] != strand:
                continue
            row = intersection(header, blocks, req["source_start"], req["source_end"], ends)
            if row and (row["target_start"], row["target_end"]) == (
                req["target_start"],
                req["target_end"],
            ):
                matches[req["element_id"], req["track"]].append(row)
    rows = []
    metrics = (
        "aligned_bp",
        "source_unaligned_bp",
        "target_unaligned_bp",
        "source_left_trim",
        "source_right_trim",
        "block_count",
        "blocks",
    )
    for req in requests:
        hits = matches[req["element_id"], req["track"]]
        if not hits:
            raise ValueError(
                f"No chain reproduces published mapping: {req['element_id']} {req['track']}"
            )
        unique = len(hits) == 1
        row = req | dict(
            chain_matches=len(hits),
            chain_ids=";".join(sorted(h["chain_id"] for h in hits)),
            complete=unique and hits[0]["complete"],
        )
        row.update({key: hits[0][key] if unique else None for key in metrics})
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["element_id", "track"]).reset_index(drop=True)
