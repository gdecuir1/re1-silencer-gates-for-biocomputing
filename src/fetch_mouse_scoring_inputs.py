"""Acquire the raw inputs for mouse REST-motif scoring.

Nothing is downloaded on import. Run this module explicitly. It reuses the
existing small-file helper, the resumable ``fetch`` used for the human chain,
``pilot.fetch_metadata``, and ``scripts/fetch_genomes.sh`` for mm39 only.

This command does not write ``orthologs.mm39.csv`` or ``orthologs.mm39.fa``.
After these raw files exist, create them offline with ``make pilot-orthologs``.
"""

import argparse
import subprocess
from pathlib import Path

from src import config
from src.downloads import download_small_file
from src.fetch_lift_inputs import fetch, remote_size
from src.pilot import fetch_metadata

MATRIX_URL = "https://jaspar.elixir.no/api/v1/matrix/MA0138.2.meme"
CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm39.over.chain.gz"
MM39_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.2bit"
MM39_GENOME = config.RAW_DIR / "genomes" / "mm39.2bit"


def motif_ids(path: Path) -> set[str]:
    """Return MEME motif identifiers, using the same token as ``src.motifs``."""
    return {
        fields[1]
        for line in path.read_text().splitlines()
        if len(fields := line.split()) >= 2 and fields[0] == "MOTIF"
    }


def acquire_matrix() -> Path:
    """Download the JASPAR REST matrix and require motif id MA0138.2.

    The HTML matrix page is not a scoring file. Only the API MEME response is
    accepted. A response without that motif id is moved aside so a later run
    can replace it; the ledger row for the rejected bytes stays appended.
    """
    path = download_small_file(MATRIX_URL, config.MOTIF_PATH, fetch=True, max_bytes=1_000_000)
    if config.JASPAR_ID not in motif_ids(path):
        rejected = path.with_name(path.name + ".rejected")
        path.replace(rejected)
        raise ValueError(
            f"JASPAR response has no MOTIF {config.JASPAR_ID}; moved aside to {rejected.name}"
        )
    relative = path.relative_to(config.ROOT).as_posix()
    print(f"Matrix ready: {relative}", flush=True)
    return path


def acquire_chain() -> Path:
    """Keep a size-matched chain, or download it with the existing resumable fetch."""
    return fetch(CHAIN_URL, config.CHAIN_PATHS["mm39"], remote_size(CHAIN_URL))


def acquire_metadata() -> None:
    """Fetch the small hg38 and mm39 verification files the mouse gate reads offline."""
    errors = fetch_metadata(("hg38", "mm39"))
    if errors:
        raise RuntimeError("Genome metadata fetch failed: " + "; ".join(errors))
    print("Metadata ready: hg38 and mm39 chrom.sizes and md5sum.txt", flush=True)


def acquire_genome() -> Path:
    """Download mm39.2bit through the genome script, skipping a complete local file.

    The script owns resume, the publisher MD5, and provenance. Passing ``mm39``
    keeps the other target species out of this step. A present file is accepted
    when its byte size matches the published Content-Length; ``make pilot-orthologs``
    repeats the whole-file MD5 check before using it.
    """
    size = remote_size(MM39_URL)
    relative = MM39_GENOME.relative_to(config.ROOT).as_posix()
    if MM39_GENOME.is_file() and MM39_GENOME.stat().st_size == size:
        print(f"Cached and validated: {relative}, {size:,} bytes", flush=True)
        return MM39_GENOME
    script = config.ROOT / "scripts" / "fetch_genomes.sh"
    print(f"{relative}: {size:,} bytes; downloading with {script.name}", flush=True)
    completed = subprocess.run(
        ["bash", str(script), "--yes", "mm39"],
        cwd=config.ROOT,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"{script.name} exited {completed.returncode}")
    observed = MM39_GENOME.stat().st_size if MM39_GENOME.is_file() else 0
    if observed != size:
        raise RuntimeError(f"Incomplete {relative}: {observed} != {size} bytes")
    return MM39_GENOME


def main(argv: list[str] | None = None) -> None:
    """Acquire matrix, chain, verification metadata, then the mm39 genome."""
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    print(
        "Mouse scoring inputs. This does not write orthologs.mm39.csv or orthologs.mm39.fa.",
        flush=True,
    )
    acquire_matrix()
    acquire_chain()
    acquire_metadata()
    acquire_genome()
    print("Raw inputs are present. Next, offline: make pilot-orthologs", flush=True)


if __name__ == "__main__":
    main()
