PYTHON ?= python
CONDA ?= conda

.PHONY: env assembly-check lift pilot pilot-mouse pilot-orthologs pilot-fimo pilot-rescue pilot-validate scoring-validation dataset gate-metrics conservation validate-models test lint clean

# Environment and explicit input acquisition
# Data targets below are offline; fetch commands are separate:
# python -m src.assembly_check --fetch
# python -m src.fetch_lift_inputs [--genome]
# python -m src.fetch_mouse_scoring_inputs
# scripts/fetch_genomes.sh

env:
	$(CONDA) env create -f environment.yml

# Human coordinates: assembly resolution -> hg38 lift and sequence verification
# Source discovery is a recorded inventory, not an automated make target.

assembly-check:
	$(PYTHON) -m src.assembly_check

lift:
	$(PYTHON) -m src.lift

# Frozen mouse track: pilot -> element lift -> motif lift and extraction
# Requires the verified hg38 catalog from lift; never relifts or downloads implicitly.
# Canonical selection -> element lift -> motif lift, three-way calls and extraction.

pilot:
	$(PYTHON) -m src.pilot

pilot-mouse: pilot
	$(PYTHON) -m src.pilot_mouse

pilot-orthologs: pilot-mouse
	$(PYTHON) -m src.pilot_orthologs

# Frozen mouse track: FIMO -> flank rescue -> correspondence evidence
# FIMO and rescue consume existing inputs without rebuilding upstream stages.

pilot-fimo:
	$(PYTHON) -m src.pilot_fimo

pilot-rescue:
	$(PYTHON) -m src.rescue --assembly mm39 \
		--chain data/raw/chains/hg38ToMm39.over.chain.gz \
		--genome data/raw/genomes/mm39.2bit \
		--sidecar data/processed/pilot_mouse/orthologs.mm39.csv \
		--baseline-scores data/processed/pilot_mouse/fimo_scores.csv \
		--baseline-report results/tables/pilot_mouse_fimo.md \
		--output-prefix data/processed/pilot_mouse_rescue/rescue \
		--report results/tables/pilot_mouse_rescue.md

# Publish final evidence after dependencies, including under make -j.
pilot-validate: pilot-orthologs
	$(MAKE) pilot-fimo
	$(MAKE) pilot-rescue
	$(PYTHON) -m src.pilot_evidence

# Human activity: scoring/sign validation -> context-matched dataset
# The explicit continuous override is preserved; scoring/sign errors still fail.
# Dataset consumes cached validated scores; it does not silently rerun this audit.

scoring-validation:
	$(PYTHON) -m src.scoring_validation

dataset:
	$(PYTHON) -m src.dataset

# Reporter outcomes: matched ON/OFF -> operating regions -> cross-context ranks
# No enhancer pooling; separate scrambled and NegCtrl references.
gate-metrics:
	$(PYTHON) -m src.gate_metrics

# Conservation: cached features -> corrected sanity gate -> join -> endpoint previews
# Missing or changed extraction snapshots fail; this target never re-extracts.
conservation:
	$(PYTHON) -m src.conservation

# Held-out evaluation: chromosome/CD-HIT folds -> two fixed rank-linear models
# Cached inputs only; no tuning or model selection.
validate-models:
	$(PYTHON) -m src.heldout_validation

# Remaining production stages: cofactors, broader models, variant robustness,
# and species expansion. No placeholder target claims completion.

# Software checks

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

# Disposable cache cleanup (not an analysis stage)
clean:
	rm -rf .pytest_cache .ruff_cache
	find src tests notebooks -type d -name __pycache__ -prune -exec rm -rf {} +
