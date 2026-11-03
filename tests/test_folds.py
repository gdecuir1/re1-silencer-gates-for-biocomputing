"""Exercise transitive chromosome/cluster constraints with synthetic memberships."""

from pathlib import Path

import pandas as pd
import pytest

from src import folds


def _memberships() -> pd.DataFrame:
    """Return a synthetic chromosome chain linked by two different clusters."""
    return pd.DataFrame(
        {
            "element_id": [f"synthetic_{i}" for i in range(7)],
            "chrom": ["chr1", "chr2", "chr2", "chr3", "chr4", "chr5", "chr6"],
            "cluster_id": ["a", "a", "b", "b", "c", "d", "e"],
        }
    )


def test_transitive_components_no_leakage() -> None:
    """Keep chromosomes 1, 2, and 3 together despite their different clusters."""
    result = folds.assign_folds(_memberships(), n_splits=3)
    assert result["fold"].nunique() == 3
    for column in ("chrom", "cluster_id", "element_id"):
        assert result.groupby(column)["fold"].nunique().eq(1).all()
    assert result[result["chrom"].isin(["chr1", "chr2", "chr3"])]["fold"].nunique() == 1
    for fold in result["fold"].unique():
        train, test = result[result["fold"] != fold], result[result["fold"] == fold]
        for column in ("element_id", "chrom", "cluster_id"):
            assert set(train[column]).isdisjoint(test[column])


def test_deterministic_under_row_permutation() -> None:
    """Fold assignments do not depend on row order or a random seed."""
    source = _memberships()
    pd.testing.assert_frame_equal(
        folds.assign_folds(source, n_splits=3),
        folds.assign_folds(source.sample(frac=1, random_state=7), n_splits=3),
    )


def test_insufficient_components() -> None:
    """Never split a connected component to satisfy a requested fold count."""
    with pytest.raises(ValueError, match="Only 4 independent"):
        folds.assign_folds(_memberships(), n_splits=5)


def test_duplicate_elements_rejected() -> None:
    """Require deduplication of cell-type rows before generating folds."""
    source = _memberships()
    with pytest.raises(ValueError, match="Duplicate element_id"):
        folds.assign_folds(pd.concat([source, source]), n_splits=2)


def test_missing_cluster_rejected() -> None:
    """Unclustered elements cannot silently become independent groups."""
    source = _memberships()
    source.loc[0, "cluster_id"] = None
    with pytest.raises(ValueError, match="nonempty"):
        folds.assign_folds(source, n_splits=2)


def test_parse_cluster_identifiers(tmp_path: Path) -> None:
    """Preserve long identifiers and reverse-strand CD-HIT memberships."""
    path = tmp_path / "synthetic.clstr"
    path.write_text(
        ">Cluster 0\n0\t200nt, >synthetic_identifier_longer_than_twenty... *\n"
        "1\t200nt, >synthetic_reverse... at -/99.00%\n>Cluster 1\n0\t200nt, >other... *\n"
    )
    result = folds.parse_clusters(path)
    assert result["cluster_id"].tolist() == ["0", "0", "1"]
    assert result.iloc[0]["element_id"] == "synthetic_identifier_longer_than_twenty"


def test_build_writes_csv(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Save fold definitions only to the explicitly supplied temporary destination."""
    memberships = _memberships()

    def synthetic_clusters(elements: pd.DataFrame, **kwargs: object) -> pd.DataFrame:
        """Supply known synthetic memberships without requiring CD-HIT."""
        return memberships[["element_id", "cluster_id"]]

    monkeypatch.setattr(folds, "cluster_sequences", synthetic_clusters)
    output = tmp_path / "folds" / "assignments.csv"
    result = folds.build_folds(memberships.assign(sequence="A" * 200), output=output, n_splits=3)
    pd.testing.assert_frame_equal(pd.read_csv(output), result)
