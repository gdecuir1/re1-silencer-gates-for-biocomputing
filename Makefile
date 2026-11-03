PYTHON ?= python
CONDA ?= conda

.PHONY: env test lint clean

env:
	$(CONDA) env create -f environment.yml

test:
	$(PYTHON) -m pytest

lint:
	$(PYTHON) -m ruff check .

# Only disposable Python/tool caches; never data, results, or environments.
clean:
	rm -rf .pytest_cache .ruff_cache
	find src tests notebooks -type d -name __pycache__ -prune -exec rm -rf {} +
