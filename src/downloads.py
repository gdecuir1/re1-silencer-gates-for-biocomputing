"""Small, explicit downloads with an append-only provenance ledger.

This module never downloads anything at import time. It deliberately refuses files
above 100 MB, so it cannot accidentally acquire a complete reference genome. MD5
values are local identity records unless a publisher checksum is separately supplied.
"""

import hashlib
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from src import config


def local_md5(path: Path) -> str:
    """Compute an identity digest without loading the file into memory.

    Args:
        path: File whose bytes should be identified.

    Returns:
        Hexadecimal MD5 digest; this is not a security/authenticity guarantee.
    """
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "md5").hexdigest()


def download_small_file(
    url: str, destination: Path, *, fetch: bool = False, max_bytes: int = 10_000_000
) -> Path:
    """Read a verified local cache, or explicitly download and record a small file.

    Args:
        url: Exact HTTPS source URL, including any sequence-query coordinates.
        destination: Local file below the repository's data directory.
        fetch: Permit network access only when explicitly True.
        max_bytes: Maximum response size, at most 100,000,000 bytes.

    Returns:
        Local cache path with an existing matching provenance record.

    Raises:
        ValueError: The path/size is unsafe or the cache differs from its ledger.
        FileNotFoundError: A required cached file is absent and fetching is disabled.
        RuntimeError: curl fails. HTML responses with HTTP 200 still require format
            validation by the caller; a successful HTTP request is not a valid dataset.
    """
    destination = Path(destination).resolve()
    if not destination.is_relative_to(config.DATA_DIR.resolve()):
        raise ValueError("Downloads must remain below the repository data directory")
    if not url.startswith("https://") or not 0 < max_bytes <= 100_000_000:
        raise ValueError("Require HTTPS and a size limit no greater than 100 MB")
    ledger = config.DATA_DIR / "README.md"
    relative = destination.relative_to(config.ROOT).as_posix()
    if destination.exists():
        digest = local_md5(destination)
        # Caches must match both source URL and bytes, not merely a familiar filename.
        rows = ledger.read_text().splitlines()
        if not any(f"`{relative}`" in row and url in row and digest in row for row in rows):
            raise ValueError(f"Cache lacks matching provenance or has changed: {relative}")
        if destination.stat().st_size > max_bytes:
            raise ValueError(f"Cached file exceeds requested size limit: {relative}")
        return destination
    if not fetch:
        raise FileNotFoundError(f"Missing {relative}; rerun with --fetch for small evidence files")
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    if partial.exists():
        raise FileExistsError(f"Inspect an existing partial download before retrying: {partial}")
    try:
        # An argv list avoids shell interpretation of ';' in UCSC query URLs.
        result = subprocess.run(
            [
                "curl",
                "--fail",
                "--show-error",
                "--silent",
                "--location",
                "--max-time",
                "45",
                "--max-filesize",
                str(max_bytes),
                "--output",
                str(partial),
                url,
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            raise RuntimeError(f"Fetch failed for {url}: {result.stderr.strip()}")
        size = partial.stat().st_size
        if not 0 < size <= max_bytes:
            raise ValueError(f"Invalid download size: {size} bytes for {url}")
        digest = local_md5(partial)
        timestamp = datetime.now(UTC).isoformat(timespec="seconds")
        partial.replace(destination)
        # This write belongs inside the downloader so a caller cannot forget it.
        with ledger.open("a") as handle:
            handle.write(f"| `{relative}` | {url} | {timestamp} | {size} | `{digest}` |\n")
        return destination
    finally:
        # Only discard the temporary file this invocation created, never user data.
        if partial.exists():
            partial.unlink()
