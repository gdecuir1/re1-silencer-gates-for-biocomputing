"""Pure classification of explicit ortholog states; failed lift is not absence."""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum

import pandas as pd

from src import config
from src._utils import require_columns


class OrthologState(StrEnum):
    """Distinct outcomes that must be supplied for every target species."""

    RETAINED = "retained"
    NO_MOTIF = "lifted_no_motif"
    FAILED_LIFT = "failed_lift"
    UNSCORABLE = "lifted_unscorable"


@dataclass(frozen=True)
class Classification:
    """Observed retention class with separate absence and missing-evidence fields."""

    evolutionary_class: str
    retained_species: tuple[str, ...]
    no_motif_species: tuple[str, ...]
    failed_lift_species: tuple[str, ...]
    unscorable_species: tuple[str, ...]


def classify_element(states: Mapping[str, OrthologState | str]) -> Classification:
    """Apply the configured deepest-retention rule without inferring missing data.

    Args:
        states: Exactly the five target species mapped to explicit OrthologState
            values or their strings. Booleans/NaNs are deliberately rejected.

    Returns:
        Evidence class plus every species' outcome. With no retained motif,
        incomplete evidence is 'unresolved'; complete absence is
        'no_observed_outgroup_motif', not an assertion of human specificity.

    Raises:
        ValueError: A species is missing/extra or a state is invalid.
    """
    missing = config.TARGET_SPECIES - states.keys()
    extra = states.keys() - config.TARGET_SPECIES
    if missing or extra:
        raise ValueError(f"Species mismatch: missing={sorted(missing)}, extra={sorted(extra)}")
    normalized = {}
    for species, state in states.items():
        try:
            normalized[species] = OrthologState(state)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid ortholog state for {species}: {state!r}") from exc
    by_state = {
        state: tuple(sorted(species for species, value in normalized.items() if value == state))
        for state in OrthologState
    }
    retained = by_state[OrthologState.RETAINED]
    failed = by_state[OrthologState.FAILED_LIFT]
    unscorable = by_state[OrthologState.UNSCORABLE]
    label = "unresolved" if failed or unscorable else "no_observed_outgroup_motif"
    for candidate, species in config.EVOLUTION_RULES:
        if species.intersection(retained):
            label = candidate
            break
    return Classification(label, retained, by_state[OrthologState.NO_MOTIF], failed, unscorable)


def classify_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    """Classify an element-indexed species-state matrix without mutating it.

    Args:
        matrix: Unique, nonmissing element index and exactly five species columns.

    Returns:
        Copy of the input, retaining all states, plus evolutionary_class and
        evidence_complete columns. Missing evidence never becomes a negative call.
    """
    require_columns(matrix, config.TARGET_SPECIES)
    if set(matrix.columns) != config.TARGET_SPECIES:
        raise ValueError("Matrix must contain exactly the target species columns")
    if matrix.index.has_duplicates or matrix.index.hasnans:
        raise ValueError("Matrix requires a unique, nonmissing element index")
    calls = [classify_element(row.to_dict()) for _, row in matrix.iterrows()]
    result = matrix.copy()
    result["evolutionary_class"] = [call.evolutionary_class for call in calls]
    result["evidence_complete"] = [
        not (call.failed_lift_species or call.unscorable_species) for call in calls
    ]
    return result
