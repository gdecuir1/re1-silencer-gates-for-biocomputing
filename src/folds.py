"""CD-HIT-EST clustering and chromosome/cluster connected-component CV folds."""

import re
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from src import config
from src._utils import require_columns, run_command, validate_ids, validate_sequences, write_fasta


def parse_clusters(path: Path) -> pd.DataFrame:
    """Parse a CD-HIT .clstr file without truncating sequence identifiers.

    Args:
        path: CD-HIT-EST cluster membership file produced with '-d 0'.

    Returns:
        Unique element_id and cluster_id membership rows.
    """
    rows: list[tuple[str, str]] = []
    cluster: str | None = None
    seen: set[str] = set()
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        if line.startswith(">Cluster "):
            cluster = line.removeprefix(">Cluster ").strip()
            if not cluster or cluster in seen:
                raise ValueError("Invalid or duplicate CD-HIT cluster header")
            seen.add(cluster)
        else:
            match = re.fullmatch(r"\s*\d+\s+\d+nt, >(.+)\.\.\. (?:\*|at .+)\s*", line)
            if cluster is None or match is None:
                raise ValueError(f"Malformed CD-HIT membership line: {line}")
            rows.append((match.group(1), cluster))
    if not rows:
        raise ValueError("CD-HIT returned no cluster memberships")
    result = pd.DataFrame(rows, columns=["element_id", "cluster_id"])
    validate_ids(result["element_id"])
    return result


def cluster_sequences(
    elements: pd.DataFrame,
    *,
    identity: float = config.CDHIT_IDENTITY,
    coverage: float = config.CDHIT_COVERAGE,
    executable: str = "cd-hit-est",
    threads: int = 1,
) -> pd.DataFrame:
    """Cluster full DNA elements on both strands using CD-HIT-EST.

    Args:
        elements: Unique element_id and 200 bp sequence columns.
        identity: Global identity threshold, supported here from 0.8 through 1.0.
        coverage: Minimum alignment coverage of both sequences.
        executable: CD-HIT nucleotide executable.
        threads: Positive number of worker threads.

    Returns:
        One cluster membership per input element. CD-HIT is a greedy
        representative-based algorithm, not an exhaustive pairwise guarantee.
    """
    require_columns(elements, ["element_id", "sequence"])
    validate_ids(elements["element_id"])
    validate_sequences(elements["sequence"])
    if elements.empty:
        raise ValueError("Cannot cluster an empty element table")
    if not 0.8 <= identity <= 1 or not 0 < coverage <= 1:
        raise ValueError("Require 0.8 <= identity <= 1 and 0 < coverage <= 1")
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        raise ValueError("threads must be a positive integer")
    word_length = 10 if identity >= 0.9 else 8 if identity >= 0.88 else 7 if identity >= 0.85 else 5
    ordered = elements.sort_values("element_id")
    aliases = {f"seq_{index}": name for index, name in enumerate(ordered["element_id"])}
    records = dict(zip(aliases, ordered["sequence"].str.upper(), strict=True))
    with TemporaryDirectory(prefix="re1-cdhit-") as directory:
        temporary = Path(directory)
        fasta, output = temporary / "input.fa", temporary / "clusters.fa"
        write_fasta(records, fasta)
        run_command(
            [
                executable,
                "-i",
                str(fasta),
                "-o",
                str(output),
                "-c",
                str(identity),
                "-n",
                str(word_length),
                "-G",
                "1",
                "-aS",
                str(coverage),
                "-aL",
                str(coverage),
                "-r",
                "1",
                "-g",
                "1",
                "-d",
                "0",
                "-T",
                str(threads),
                "-M",
                "800",
            ]
        )
        clusters = parse_clusters(Path(f"{output}.clstr"))
    if set(clusters["element_id"]) != set(aliases):
        raise ValueError("CD-HIT must return every element exactly once")
    clusters["element_id"] = clusters["element_id"].map(aliases)
    return clusters


def assign_folds(elements: pd.DataFrame, *, n_splits: int = config.N_FOLDS) -> pd.DataFrame:
    """Assign whole connected chromosome/cluster components to balanced folds.

    Args:
        elements: One row per element with element_id, chrom, and cluster_id.
        n_splits: Number of nonempty folds, at least two.

    Returns:
        Assignments containing element_id, chrom, cluster_id, component_id, fold.
        Every chromosome and every cluster appears in exactly one fold. Connected
        components are greedily balanced by element count, with deterministic ties.

    Raises:
        ValueError: Too few independent components or incomplete/duplicate input.
    """
    require_columns(elements, ["element_id", "chrom", "cluster_id"])
    validate_ids(elements["element_id"])
    validate_ids(elements["chrom"], unique=False)
    validate_ids(elements["cluster_id"], unique=False)
    if not isinstance(n_splits, int) or isinstance(n_splits, bool) or n_splits < 2:
        raise ValueError("n_splits must be an integer >= 2")
    parent = {chrom: chrom for chrom in elements["chrom"].unique()}

    def find(chrom: str) -> str:
        """Return a component root, compressing its parent path.

        Args:
            chrom: Chromosome to locate.

        Returns:
            Root chromosome identifier.
        """
        while parent[chrom] != chrom:
            parent[chrom] = parent[parent[chrom]]
            chrom = parent[chrom]
        return chrom

    for _, group in elements.groupby("cluster_id", sort=True):
        chromosomes = sorted(group["chrom"].unique())
        for chrom in chromosomes[1:]:
            left, right = sorted([find(chromosomes[0]), find(chrom)])
            parent[right] = left
    result = elements[["element_id", "chrom", "cluster_id"]].copy()
    result["component_id"] = result["chrom"].map(find)
    sizes = result.groupby("component_id").size().to_dict()
    if len(sizes) < n_splits:
        raise ValueError(
            f"Only {len(sizes)} independent chromosome/cluster components for "
            f"{n_splits} folds; do not split components to force CV"
        )
    totals = [0] * n_splits
    allocation: dict[str, int] = {}
    for component in sorted(sizes, key=lambda key: (-sizes[key], key)):
        fold = min(range(n_splits), key=lambda index: (totals[index], index))
        allocation[component] = fold
        totals[fold] += sizes[component]
    result["fold"] = result["component_id"].map(allocation)
    for column in ("chrom", "cluster_id"):
        if result.groupby(column)["fold"].nunique().gt(1).any():
            raise RuntimeError(f"Internal leakage guard failure: {column} spans folds")
    return result.sort_values("element_id").reset_index(drop=True)


def build_folds(
    elements: pd.DataFrame,
    *,
    output: Path = config.FOLD_PATH,
    n_splits: int = config.N_FOLDS,
    identity: float = config.CDHIT_IDENTITY,
    coverage: float = config.CDHIT_COVERAGE,
    executable: str = "cd-hit-est",
) -> pd.DataFrame:
    """Cluster sequences, construct leakage-safe groups, and save a fold CSV.

    Args:
        elements: Unique element_id, chrom, and 200 bp sequence columns.
        output: CSV destination, default results/folds/assignments.csv.
        n_splits: Number of independent folds.
        identity: CD-HIT global identity threshold.
        coverage: CD-HIT alignment coverage of both sequences.
        executable: CD-HIT-EST executable.

    Returns:
        Assignments including the clustering thresholds for reproducibility.
    """
    require_columns(elements, ["element_id", "chrom", "sequence"])
    clusters = cluster_sequences(
        elements, identity=identity, coverage=coverage, executable=executable
    )
    joined = elements[["element_id", "chrom"]].merge(
        clusters, on="element_id", validate="one_to_one"
    )
    assignments = assign_folds(joined, n_splits=n_splits)
    assignments["cdhit_identity"] = identity
    assignments["cdhit_coverage"] = coverage
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    assignments.to_csv(output, index=False)
    return assignments
