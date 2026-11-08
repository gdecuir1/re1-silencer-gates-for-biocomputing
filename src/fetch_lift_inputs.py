"""Explicitly acquire the three reference inputs for the hg19 -> hg38 conversion.

Nothing is downloaded on import or by ``make lift``. Run this module explicitly;
``--genome`` additionally authorizes the full hg38 reference (about 835 MB).
Downloads resume into .part files, use HTTPS, and append provenance inside the
transfer function. A complete hg38 file must match the publisher's MD5. UCSC's
hg19 database directory has no md5sum.txt (HTTP 404 observed); the RepeatMasker
and chain digests therefore document local identity, not publisher verification.
"""

import argparse
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from src import config
from src.downloads import local_md5

BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"


def run_curl(arguments: list[str]) -> str:
    """Run curl without a shell and surface the complete diagnostic on failure."""
    result = subprocess.run(["curl", *arguments], capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"curl exited {result.returncode}: {result.stderr.strip()}")
    return result.stdout


def remote_size(url: str) -> int:
    """Require a published Content-Length before committing to a large transfer."""
    headers = run_curl(["-fIsSL", "--max-time", "60", url])
    sizes = [
        int(line.split(":", 1)[1])
        for line in headers.splitlines()
        if line.lower().startswith("content-length:")
    ]
    if not sizes or sizes[-1] <= 0:
        raise ValueError(f"No usable Content-Length for {url}")
    return sizes[-1]


def fetch(url: str, path: Path, size: int, expected_md5: str | None = None) -> Path:
    """Resume one download, validate it, then append one immutable ledger record.

    Existing complete files are accepted only with matching URL/digest provenance.
    Interrupted transfers retain .part bytes for resumption. Failed checksums stay
    in .part and cannot be mistaken for a usable reference. Never overwrite an
    existing final file. Repeated successful calls do not duplicate ledger rows.
    """
    path = path.resolve()
    if not url.startswith("https://") or not path.is_relative_to(config.DATA_DIR):
        raise ValueError("Require HTTPS and a destination under data/")
    relative = path.relative_to(config.ROOT).as_posix()
    ledger = config.DATA_DIR / "README.md"
    if path.exists():
        digest = local_md5(path)
        if path.stat().st_size != size or (expected_md5 and digest != expected_md5):
            raise ValueError(f"Existing reference fails size/MD5 validation: {relative}")
        if not any(
            f"`{relative}`" in row and url in row and digest in row
            for row in ledger.read_text().splitlines()
        ):
            raise ValueError(f"Existing reference has no matching provenance: {relative}")
        print(f"Cached and validated: {relative}", flush=True)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    # A fully transferred .part can be validated after an interruption just before
    # the rename. Otherwise continue from its exact current byte count.
    if not partial.exists() or partial.stat().st_size != size:
        run_curl(
            [
                "--fail",
                "--show-error",
                "--silent",
                "--location",
                "--continue-at",
                "-",
                "--retry",
                "3",
                "--connect-timeout",
                "30",
                "--max-time",
                "3600",
                "--output",
                str(partial),
                url,
            ]
        )
    if partial.stat().st_size != size:
        raise ValueError(f"Incomplete {relative}: {partial.stat().st_size} != {size} bytes")
    digest = local_md5(partial)
    if expected_md5 and digest != expected_md5:
        raise ValueError(f"MD5 mismatch for {relative}: observed {digest}; expected {expected_md5}")
    partial.replace(path)
    with ledger.open("a") as handle:
        handle.write(
            f"| `{relative}` | {url} | {datetime.now(UTC).isoformat(timespec='seconds')} "
            f"| {size} | `{digest}` |\n"
        )
    print(f"Downloaded: {relative}, {size:,} bytes, MD5 {digest}", flush=True)
    return path


def main() -> None:
    """Print the byte budget first; acquire only the explicitly requested inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", action="store_true", help="Also download the full hg38.2bit")
    args = parser.parse_args()
    specs = [
        (f"{BASE}/hg19/liftOver/hg19ToHg38.over.chain.gz", config.SOURCE_TO_WORKING_CHAIN),
        (f"{BASE}/hg19/database/rmsk.sql", config.HG19_REPEATS.parent / "hg19.rmsk.sql"),
        (f"{BASE}/hg19/database/rmsk.txt.gz", config.HG19_REPEATS),
    ]
    if args.genome:
        specs += [
            (f"{BASE}/hg38/bigZips/md5sum.txt", config.HG38_TWOBIT.parent / "hg38.md5sum.txt"),
            (f"{BASE}/hg38/bigZips/hg38.2bit", config.HG38_TWOBIT),
        ]
    sizes = [remote_size(url) for url, _ in specs]
    total = 0
    for (_url, path), size in zip(specs, sizes, strict=True):
        total += size
        print(f"{path.name}: {size:,} bytes; running total {total:,}", flush=True)
    print(f"Planned reference bytes (including caches): {total:,}", flush=True)
    for (url, path), size in zip(specs, sizes, strict=True):
        expected = None
        if path == config.HG38_TWOBIT:
            manifest = (path.parent / "hg38.md5sum.txt").read_text()
            matches = [
                line.split()[0]
                for line in manifest.splitlines()
                if len(line.split()) == 2 and line.split()[1].lstrip("*") == "hg38.2bit"
            ]
            if len(matches) != 1:
                raise ValueError("Require exactly one publisher MD5 for hg38.2bit")
            expected = matches[0]
        fetch(url, path, size, expected)


if __name__ == "__main__":
    main()
