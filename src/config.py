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
GENOME_BUILD = "hg38"  # Project analysis target; NOT the raw GEO native intervals.
# Empirically resolved using full Table S2 oligos at unchanged GEO coordinates:
# 47/47 stratified native oligos match hg19, versus 10/47 at hg38 coordinates.
# Scope: project='motif' and 'noMotif'. Do not assume variant-coordinate assemblies.
# Preserve source coordinates and explicitly lift/validate before hg38 workflows.
GEO_NATIVE_GENOME_BUILD = "hg19"
# Explicit names used by the one-time, sequence-verified assembly conversion.
SOURCE_BUILD = GEO_NATIVE_GENOME_BUILD
WORKING_BUILD = GENOME_BUILD
SOURCE_TO_WORKING_CHAIN = RAW_DIR / "hg19ToHg38.over.chain.gz"
HG38_TWOBIT = RAW_DIR / "genomes" / "hg38.2bit"
HG19_REPEATS = RAW_DIR / "repeats" / "hg19.rmsk.txt.gz"
HG38_ELEMENTS = PROCESSED_DIR / "elements_hg38.csv"
TARGET_SPECIES = {"rheMac10", "mm39", "canFam6", "bosTau9", "monDom5"}
EXPECTED_ELEMENT_COUNT = 8436
# Week 2 pilot: canonical elements only, using verified hg38 coordinates.
PILOT_CHROM = "chr8"
# chr8 sweep: 0.5 maps 227; 0.2/0.1 add mostly >20%-distorted spans (only 2 stable gains).
MINMATCH_CROSS_SPECIES = 0.5
MINMATCH_CROSS_SPECIES_RATIONALE = (
    "Recommend 0.5 for this mouse pilot: 227/390 map (58.21%), with 163 spans "
    "within 160–240 bp and 64 outside. At 0.95 only 63 map (56 within the band). "
    "Lowering 0.5 to 0.2 gains 21 IDs but loses one as duplicated (net +20), "
    "with only two additional in-band spans, "
    "while the out-of-band count rises from 64 to 82. Lowering 0.2 to 0.1 adds "
    "five mappings, all out of band (87 total); the in-band count remains 165. "
    "Thus 0.5 captures most of the width-stable coverage without the extra "
    "span distortion of 0.2/0.1. All 227 mappings are retained and flagged; "
    "this choice does not certify mapping accuracy or motif retention."
)

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
