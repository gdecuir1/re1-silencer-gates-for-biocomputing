"""Range-only extraction snapshots of verified element windows, never whole bigWigs.

Each cache contains only a chromosome's requested 200-base elements, in genomic
orientation. URL, HTTP identity, coordinates and array hashes guard replay.
"""

import argparse
import hashlib
import json
import os
import re
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import pyBigWig

from src import config
from src._utils import validate_intervals
from src.dataset import digest

CACHE = config.PROCESSED_DIR / "conservation_windows"


def track_identity(url):
    """Read one byte only after verifying a precise partial-content response."""
    request = Request(url, headers={"Range": "bytes=0-0", "Accept-Encoding": "identity"})
    with urlopen(request, timeout=30) as response:
        match = re.fullmatch(r"bytes 0-0/(\d+)", response.headers.get("Content-Range", ""))
        if response.status != 206 or not match:
            raise ValueError("Server ignored byte range; refusing full bigWig download")
        if len(response.read(1)) != 1:
            raise ValueError("Incomplete range probe")
        return dict(
            url=url,
            resolved_url=response.url,
            total_bytes=int(match[1]),
            etag=response.headers.get("ETag"),
            last_modified=response.headers.get("Last-Modified"),
            range_status=206,
        )


def verified_elements():
    """Admit only intervals identical to the successful sequence-verification audit."""
    elements = pd.read_csv(config.HG38_ELEMENTS)
    audit = pd.read_csv(config.PROCESSED_DIR / "liftover_hg19_to_hg38/element_audit.csv")
    validate_intervals(elements)
    if not elements.working_build.eq("hg38").all():
        raise ValueError("Conservation requires hg38 coordinates")
    verified = audit[audit.bucket.eq("lifted+verified")].set_index("element_id")
    if set(verified.index) != set(elements.element_id):
        raise ValueError("Verified element membership differs from the sequence audit")
    for col in ("project", "chrom", "start", "end", "strand", "sequence"):
        actual = elements.set_index("element_id")[col].sort_index()
        expected = verified[col].sort_index()
        if not actual.eq(expected).all():
            raise ValueError(f"Coordinate/sequence audit disagreement: {col}")
    return elements.sort_values(["chrom", "start", "end", "element_id"]).reset_index(drop=True)


def read_chunk(track, identity, group, cache=CACHE):
    """Read exact full element windows with one remote handle, or verify cached values."""
    chrom = group.chrom.iloc[0]
    key = hashlib.sha256(
        group[["element_id", "chrom", "start", "end", "strand"]].to_csv(index=False).encode()
    ).hexdigest()
    spec = dict(
        track=track,
        identity=identity,
        coordinates_sha256=key,
        shape=[len(group), config.ELEMENT_LENGTH],
        orientation="genomic_forward",
    )
    token = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:20]
    base = cache / f"{track}.{chrom}.{token}"
    path, receipt = Path(str(base) + ".npy"), Path(str(base) + ".json")
    if path.exists() and receipt.exists():
        saved = json.loads(receipt.read_text())
        if saved["spec"] != spec or saved["sha256"] != digest(path):
            raise ValueError(f"Conservation snapshot identity mismatch: {path}")
        if "original_receipt" in saved:
            original = config.ROOT / saved["original_receipt"]
            if digest(original) != saved["original_receipt_sha256"]:
                raise ValueError("Original transport receipt changed")
            if digest(original.with_suffix(".npy")) != saved["sha256"]:
                raise ValueError("Original transport values changed")
        values = np.load(path, allow_pickle=False)
        absent_chrom = saved["chromosome_absent"]
    else:
        if not pyBigWig.remote:
            raise RuntimeError("pyBigWig lacks remote HTTP support")
        values = np.full((len(group), config.ELEMENT_LENGTH), np.nan, dtype=np.float64)
        rows = list(group.itertuples(index=False))
        completed = 0
        for attempt in range(4):
            try:
                with pyBigWig.open(identity["url"]) as handle:
                    if not handle.isBigWig():
                        raise ValueError(f"Not a bigWig: {track}")
                    chromosomes = handle.chroms()
                    absent_chrom = chrom not in chromosomes
                    if not absent_chrom:
                        while completed < len(rows):
                            row = rows[completed]
                            if row.end > chromosomes[chrom]:
                                raise ValueError(
                                    f"Interval exceeds track chromosome: {row.element_id}"
                                )
                            queried = np.asarray(
                                handle.values(chrom, row.start, row.end), dtype=float
                            )
                            if queried.size not in (0, config.ELEMENT_LENGTH):
                                raise ValueError(f"Truncated bigWig response: {row.element_id}")
                            if queried.size:
                                values[completed] = queried
                            completed += 1
                            if completed % 200 == 0:
                                print(
                                    f"{track} {chrom}: {completed}/{len(group)} windows", flush=True
                                )
                break
            except (RuntimeError, OSError) as error:
                if attempt == 3:
                    raise RuntimeError(
                        f"Remote extraction failed for {track}/{chrom}; "
                        "failed requests are not no-data observations"
                    ) from error
                print(f"Retry {attempt + 1}: {track} {chrom} from window {completed}", flush=True)
                time.sleep(2**attempt)
        if np.isinf(values).any():
            raise ValueError("Infinite conservation score")
        if track == "phastCons" and ((values < 0) | (values > 1)).any():
            raise ValueError("phastCons outside probability range")
        cache.mkdir(parents=True, exist_ok=True)
        temp = Path(str(path) + ".tmp")
        with temp.open("wb") as stream:
            np.save(stream, values, allow_pickle=False)
        os.replace(temp, path)
        receipt.write_text(
            json.dumps(
                dict(spec=spec, sha256=digest(path), chromosome_absent=absent_chrom),
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    if values.shape != tuple(spec["shape"]) or np.isinf(values).any():
        raise ValueError("Invalid cached window array")
    return group.index.to_numpy(), values, absent_chrom, path, receipt


def reuse_transport_cache(elements, cache=CACHE):
    """Reuse complete snapshots only across HTTP/HTTPS for the exact same file identity.

    Original receipts remain intact and are referenced by SHA256 in copied receipts.
    Independently re-read overlapping snapshots must have identical array hashes.
    """
    identities = {name: track_identity(url) for name, url in config.BIGWIG_URLS.items()}
    reused, compared = 0, 0
    for receipt in sorted(cache.glob("*.json")):
        old = json.loads(receipt.read_text())
        spec = old["spec"]
        identity = identities[spec["track"]]
        before, after = urlparse(spec["identity"]["url"]), urlparse(identity["url"])
        if before.scheme == after.scheme:
            continue
        if {before.scheme, after.scheme} != {"http", "https"}:
            continue
        if before._replace(scheme="") != after._replace(scheme=""):
            continue
        if not identity.get("etag") or not identity.get("last_modified"):
            continue
        fields = ("etag", "last_modified", "total_bytes")
        if any(spec["identity"].get(k) != identity[k] for k in fields):
            raise ValueError("Same-file transport identities disagree; cannot reuse snapshot")
        old_path = receipt.with_suffix(".npy")
        if digest(old_path) != old["sha256"]:
            raise ValueError("Original snapshot hash mismatch")
        new_spec = dict(spec, identity=identity)
        chrom = receipt.name.split(".")[1]
        group = elements[elements.chrom.eq(chrom)]
        key = hashlib.sha256(
            group[["element_id", "chrom", "start", "end", "strand"]].to_csv(index=False).encode()
        ).hexdigest()
        if spec["coordinates_sha256"] != key:
            continue
        token = hashlib.sha256(json.dumps(new_spec, sort_keys=True).encode()).hexdigest()[:20]
        base = cache / f"{spec['track']}.{chrom}.{token}"
        path, target = Path(str(base) + ".npy"), Path(str(base) + ".json")
        if path.exists() and target.exists():
            if digest(path) != old["sha256"]:
                raise ValueError("Fresh HTTP values disagree with original HTTPS snapshot")
            compared += len(group)
            continue
        data = old_path.read_bytes()
        temp = Path(str(path) + ".alias.tmp")
        temp.write_bytes(data)
        os.replace(temp, path)
        record = dict(
            spec=new_spec,
            sha256=old["sha256"],
            chromosome_absent=old["chromosome_absent"],
            original_receipt=str(receipt.relative_to(config.ROOT)),
            original_receipt_sha256=digest(receipt),
        )
        target.write_text(json.dumps(record, sort_keys=True, indent=2) + "\n")
        reused += len(group)
    print(
        f"Transport cache: {compared} independently re-read windows identical; "
        f"{reused} complete windows reused with matching file identity",
        flush=True,
    )


def extract(elements, *, workers=4, cache=CACHE):
    """Fetch two tracks with bounded parallelism; requests stay at element scope."""
    identities = {name: track_identity(url) for name, url in config.BIGWIG_URLS.items()}
    arrays = {name: np.full((len(elements), config.ELEMENT_LENGTH), np.nan) for name in identities}
    absent = {name: np.zeros(len(elements), dtype=bool) for name in identities}
    files = []
    completion_messages = []
    # pyBigWig calls can retain the GIL during remote reads; processes bound the
    # concurrent connections without serializing all requests in Python threads.
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs = {
            pool.submit(read_chunk, name, identity, group, cache): (name, chrom)
            for name, identity in identities.items()
            for chrom, group in elements.groupby("chrom", sort=True)
        }
        for future in as_completed(jobs):
            name, chrom = jobs[future]
            try:
                idx, values, missing_chrom, path, receipt = future.result()
            except Exception:
                for job in jobs:
                    job.cancel()
                raise
            arrays[name][idx] = values
            absent[name][idx] = missing_chrom
            files.extend([path, receipt])
            saved = json.loads(receipt.read_text())
            if "original_receipt" in saved:
                original = config.ROOT / saved["original_receipt"]
                files.extend([original, original.with_suffix(".npy")])
            completion_messages.append(f"{name} {chrom}: {len(idx)} windows ready")
    for message in sorted(completion_messages):
        print(message, flush=True)
    # Do not combine observations across a remote source update.
    if identities != {name: track_identity(url) for name, url in config.BIGWIG_URLS.items()}:
        raise ValueError("Remote track changed during extraction")
    return arrays, absent, identities, sorted(files)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", action="store_true", help="Check HTTP identities only")
    parser.add_argument("--reuse-transport-cache", action="store_true")
    args = parser.parse_args()
    if args.probe:
        print(
            json.dumps(
                {name: track_identity(url) for name, url in config.BIGWIG_URLS.items()},
                indent=2,
                sort_keys=True,
            )
        )
    elif args.reuse_transport_cache:
        reuse_transport_cache(verified_elements())
    else:
        extract(verified_elements())


def load_cached(elements, manifest_path):
    """Strictly offline replay of an archived extraction: missing data never trigger HTTP."""
    manifest = json.loads(Path(manifest_path).read_text())
    identities = manifest["tracks"]
    if {k: v["url"] for k, v in identities.items()} != config.BIGWIG_URLS:
        raise ValueError("Configured tracks differ from the archived extraction")
    coordinate_path = str(config.HG38_ELEMENTS.relative_to(config.ROOT))
    if digest(config.HG38_ELEMENTS) != manifest["inputs"][coordinate_path]:
        raise ValueError("Verified coordinate catalog changed since extraction")
    files = [
        config.ROOT / p
        for p in manifest["inputs"]
        if p.startswith("data/processed/conservation_windows/")
    ]
    for path in files:
        if (
            not path.exists()
            or digest(path) != manifest["inputs"][str(path.relative_to(config.ROOT))]
        ):
            raise ValueError(
                f"Cached extraction missing or changed: {path}; remote reads prohibited"
            )
    arrays = {
        track: np.full((len(elements), config.ELEMENT_LENGTH), np.nan) for track in identities
    }
    absent = {track: np.zeros(len(elements), dtype=bool) for track in identities}
    filled = {track: np.zeros(len(elements), dtype=bool) for track in identities}
    for receipt in files:
        if receipt.suffix != ".json":
            continue
        saved = json.loads(receipt.read_text())
        spec = saved["spec"]
        track = spec["track"]
        if spec["identity"] != identities[track]:
            continue  # Original HTTPS receipts retained as provenance, not duplicate observations.
        chrom = receipt.name.split(".")[1]
        group = elements[elements.chrom.eq(chrom)]
        key = hashlib.sha256(
            group[["element_id", "chrom", "start", "end", "strand"]].to_csv(index=False).encode()
        ).hexdigest()
        if spec["coordinates_sha256"] != key or spec["orientation"] != "genomic_forward":
            raise ValueError("Cached window coordinates/orientation mismatch")
        path = receipt.with_suffix(".npy")
        if path not in files or digest(path) != saved["sha256"]:
            raise ValueError("Cached array is missing from verified extraction")
        values = np.load(path, allow_pickle=False)
        if values.shape != (len(group), config.ELEMENT_LENGTH) or np.isinf(values).any():
            raise ValueError("Invalid cached window shape or values")
        if filled[track][group.index].any():
            raise ValueError("Duplicate cached window observations")
        arrays[track][group.index] = values
        absent[track][group.index] = saved["chromosome_absent"]
        filled[track][group.index] = True
    if not all(mask.all() for mask in filled.values()):
        raise ValueError("Incomplete cached extraction; remote reads prohibited")
    return arrays, absent, identities, sorted(files)
