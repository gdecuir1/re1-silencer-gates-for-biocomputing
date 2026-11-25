PYTHON ?= python
CONDA ?= conda

.PHONY: env assembly-check lift pilot pilot-mouse pilot-orthologs pilot-fimo pilot-rescue pilot-validate scoring-validation dataset gate-metrics test lint clean

# Environment and explicit acquisition
# Data targets below are offline; fetch commands are separate:
# python -m src.assembly_check --fetch
# python -m src.fetch_lift_inputs [--genome]
# python -m src.fetch_mouse_scoring_inputs
# scripts/fetch_genomes.sh

env:
	$(CONDA) env create -f environment.yml

# Source assembly: cached oligos -> assembly evidence -> verified hg38 elements
# Source discovery is a recorded inventory, not an automated make target.

assembly-check:
	$(PYTHON) -m src.assembly_check

lift:
	$(PYTHON) -m src.lift

# Ortholog classification: complete and frozen at mouse-pilot scope
# Requires the verified hg38 catalog from lift; never relifts or downloads implicitly.
# Canonical selection -> element lift -> motif lift, three-way calls and extraction.

pilot:
	$(PYTHON) -m src.pilot

pilot-mouse: pilot
	$(PYTHON) -m src.pilot_mouse

pilot-orthologs: pilot-mouse
	$(PYTHON) -m src.pilot_orthologs

# Frozen classification: mapped-window FIMO -> flank regions -> correspondence audit
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

# Continuous analysis: canonical score/sign audit -> context-matched dataset
# The explicit continuous override is preserved; scoring/sign errors still fail.
# Dataset consumes cached validated scores; it does not silently rerun this audit.

scoring-validation:
	$(PYTHON) -m src.scoring_validation

dataset:
	$(PYTHON) -m src.dataset

# Reporter operating regions from matched scrambled/NegCtrl ON states; no enhancer pooling.
gate-metrics:
	$(PYTHON) -m src.gate_metrics

# Unimplemented production stages: conservation/cofactors -> grouped folds -> models
# -> variant robustness. No placeholder target claims completion.
# Conservation and fold helpers exist; species expansion also lacks a production target.

# Verification and disposable caches

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

clean:
	rm -rf .pytest_cache .ruff_cache
	find src tests notebooks -type d -name __pycache__ -prune -exec rm -rf {} +
