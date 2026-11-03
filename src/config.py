"""Repository paths, verified constants, and explicit user-supplied settings."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"
FOLDS_DIR = RESULTS_DIR / "folds"
NOTEBOOKS_DIR = ROOT / "notebooks"
MPRA_PATH = RAW_DIR / "mpra.tsv"
MOTIF_PATH = RAW_DIR / "MA0138.2.meme"
ELEMENT_FASTA = PROCESSED_DIR / "elements.fasta"
FOLD_PATH = FOLDS_DIR / "assignments.csv"
CHAIN_PATHS = {
    species: RAW_DIR / "chains" / f"hg38To{species[0].upper()}{species[1:]}.over.chain.gz"
    for species in ("rheMac10", "mm39", "canFam6", "bosTau9", "monDom5")
}
GENOME_PATHS = {species: RAW_DIR / "genomes" / f"{species}.fa" for species in CHAIN_PATHS}

ELEMENT_LENGTH = 200
MOTIF_SPAN = (91, 111)  # Zero-based, BOTH endpoints inclusive: 21 bases.
M_VALUE = 20.86
CELL_TYPES = ["GM12878", "HepG2", "K562", "SK-N-SH"]
JASPAR_ID = "MA0138.2"
GENOME_BUILD = "hg38"
TARGET_SPECIES = {"rheMac10", "mm39", "canFam6", "bosTau9", "monDom5"}
EXPECTED_ELEMENT_COUNT = 8436

# Fill these with actual source headers; None deliberately fails with a clear error.
# Canonical name -> source header. Coordinates, if supplied, must already be BED
# coordinates (zero-based, half-open). Strand describes the supplied sequence.
MPRA_COLUMNS: dict[str, str | None] = {"element_id": None, "sequence": None}
MPRA_REPRESSION_COLUMNS: dict[str, str | None] = dict.fromkeys(CELL_TYPES)
MPRA_LAYOUT = "wide"  # "long" uses cell_type and repression in MPRA_COLUMNS.
MPRA_SEPARATOR = "\t"

# Choose and record the actual hg38 track URLs in data/README.md before querying.
BIGWIG_URLS: dict[str, str | None] = {"phyloP": None, "phastCons": None}

# Operational defaults, NOT verified biological cutoffs. Record overrides.
CDHIT_IDENTITY = 0.95
CDHIT_COVERAGE = 0.95
N_FOLDS = 5

# Deepest observed retention first. Mouse, cow, and dog support the same
# human/outgroup ancestral depth; do not rank these branches against one another.
# Labels describe evidence, not proof of the element's evolutionary origin.
EVOLUTION_RULES = (
    ("therian_retention", frozenset({"monDom5"})),
    ("eutherian_retention", frozenset({"mm39", "canFam6", "bosTau9"})),
    ("primate_retention", frozenset({"rheMac10"})),
)
