PYTHON ?= python
CONDA ?= conda

.PHONY: env test lint clean assembly-check

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

# Only disposable Python/tool caches; never data, results, or environments.
clean:
	rm -rf .pytest_cache .ruff_cache
	find src tests notebooks -type d -name __pycache__ -prune -exec rm -rf {} +
