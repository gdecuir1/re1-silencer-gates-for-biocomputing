PYTHON ?= python
CONDA ?= conda

.PHONY: env test lint clean assembly-check lift pilot pilot-mouse pilot-orthologs pilot-fimo pilot-rescue pilot-validate

env:
	$(CONDA) env create -f environment.yml

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

# Offline replay of the assembly evidence. Fetching is a separate explicit command:
# python -m src.assembly_check --fetch (small evidence files only; no full genomes).
assembly-check:
	$(PYTHON) -m src.assembly_check

# Offline conversion. Inputs must already exist; this target never downloads genomes.
lift:
	$(PYTHON) -m src.lift

# Offline genome gate, canonical chr8 selection, and seeded sequence regression.
pilot:
	$(PYTHON) -m src.pilot

# Rebuild the hg38 pilot, then verify only hg38/mm39 and sweep mouse minMatch.
pilot-mouse: pilot
	$(PYTHON) -m src.pilot_mouse

# Track motif mapping separately, then extract padded jointly mapped mouse DNA.
pilot-orthologs: pilot-mouse
	$(PYTHON) -m src.pilot_orthologs

# Offline FIMO scan of the lifted mouse motif windows. Requires orthologs.mm39.csv.
pilot-fimo:
	$(PYTHON) -m src.pilot_fimo

# Separate outputs; requires the existing sidecar and successful MAPPED FIMO baseline.
pilot-rescue:
	$(PYTHON) -m src.rescue --assembly mm39 \
		--chain data/raw/chains/hg38ToMm39.over.chain.gz \
		--genome data/raw/genomes/mm39.2bit \
		--sidecar data/processed/pilot_mouse/orthologs.mm39.csv \
		--baseline-scores data/processed/pilot_mouse/fimo_scores.csv \
		--baseline-report results/tables/pilot_mouse_fimo.md \
		--output-prefix data/processed/pilot_mouse_rescue/rescue \
		--report results/tables/pilot_mouse_rescue.md

# Rebuild in dependency order even with make -j; publish the evidence layer last.
pilot-validate: pilot-orthologs
	$(MAKE) pilot-fimo
	$(MAKE) pilot-rescue
	$(PYTHON) -m src.pilot_evidence

# Only disposable Python/tool caches; never data, results, or environments.
clean:
	rm -rf .pytest_cache .ruff_cache
	find src tests notebooks -type d -name __pycache__ -prune -exec rm -rf {} +
