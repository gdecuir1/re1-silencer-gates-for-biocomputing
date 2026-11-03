"""Test retention precedence while keeping failed lifts distinguishable."""

from itertools import product

import pandas as pd
import pytest

from src import config
from src.classify import OrthologState, classify_element, classify_matrix


@pytest.mark.parametrize(
    "species,label",
    [
        ("monDom5", "therian_retention"),
        ("mm39", "eutherian_retention"),
        ("canFam6", "eutherian_retention"),
        ("bosTau9", "eutherian_retention"),
        ("rheMac10", "primate_retention"),
    ],
)
def test_single_retained_species(species: str, label: str) -> None:
    """Choose the correct observed ancestral depth for each species."""
    states = dict.fromkeys(config.TARGET_SPECIES, OrthologState.NO_MOTIF)
    states[species] = OrthologState.RETAINED
    call = classify_element(states)
    assert call.evolutionary_class == label
    assert call.retained_species == (species,)


def test_all_state_combinations() -> None:
    """Exhaustively check precedence and evidence partitioning across five species."""
    species = sorted(config.TARGET_SPECIES)
    for values in product(OrthologState, repeat=5):
        states = dict(zip(species, values, strict=True))
        call = classify_element(states)
        if states["monDom5"] == OrthologState.RETAINED:
            expected = "therian_retention"
        elif any(states[name] == OrthologState.RETAINED for name in ("mm39", "canFam6", "bosTau9")):
            expected = "eutherian_retention"
        elif states["rheMac10"] == OrthologState.RETAINED:
            expected = "primate_retention"
        elif all(state == OrthologState.NO_MOTIF for state in values):
            expected = "no_observed_outgroup_motif"
        else:
            expected = "unresolved"
        assert call.evolutionary_class == expected
        assert (
            set(
                call.retained_species
                + call.no_motif_species
                + call.failed_lift_species
                + call.unscorable_species
            )
            == config.TARGET_SPECIES
        )
        assert set(call.failed_lift_species) == {
            name for name, state in states.items() if state == OrthologState.FAILED_LIFT
        }


@pytest.mark.parametrize("value", [True, False, None, float("nan"), "absent"])
def test_ambiguous_states_rejected(value: object) -> None:
    """Do not collapse failed lifts into Boolean absence or missing values."""
    states = dict.fromkeys(config.TARGET_SPECIES, "retained")
    states["mm39"] = value
    with pytest.raises(ValueError, match="mm39"):
        classify_element(states)


def test_missing_species() -> None:
    """A missing species column is not an absence call."""
    with pytest.raises(ValueError, match="missing"):
        classify_element({})


def test_matrix_keeps_evidence_and_input() -> None:
    """Retain all states and flag incomplete evidence without mutating input."""
    matrix = pd.DataFrame(
        [{species: "retained" for species in config.TARGET_SPECIES}], index=["synthetic_element"]
    )
    matrix.loc["synthetic_element", "monDom5"] = "failed_lift"
    original = matrix.copy(deep=True)
    result = classify_matrix(matrix)
    assert result.loc["synthetic_element", "evolutionary_class"] == "eutherian_retention"
    assert not result.loc["synthetic_element", "evidence_complete"]
    assert result.loc["synthetic_element", "monDom5"] == "failed_lift"
    pd.testing.assert_frame_equal(matrix, original)
