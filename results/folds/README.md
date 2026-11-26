# Fixed chromosome/cluster folds

These assignments cover the union of complete-case above-threshold canonical scrambled-pair elements in the five passing contexts and the existing cross-context robustness target. Each element has one assignment across endpoints. Individual endpoints retain only their measured complete cases.

- `heldout_assignments.csv`: element ID, hg38 chromosome, CD-HIT cluster, connected component, zero-based fold, identity and coverage thresholds.
- `heldout_clusters.csv`: cluster sizes and chromosomes.
- `heldout_counts.csv`: element, cluster and component counts per fold.
- `heldout_manifest.json`: source and output SHA-256 hashes, parameters and software versions.
- `heldout_reproducibility.json`: two full runs with identical output hashes and unchanged input hashes.

Construction uses the existing `src/folds.py`: CD-HIT-EST at 0.95 global identity and 0.95 coverage of both full 200-base sequences, both strands, one thread and deterministic input ordering. Whole connected chromosome/cluster components are balanced across five folds; neither chromosomes nor clusters cross folds. Missing observations are not assigned artificial outcomes.

Run `make validate-models PYTHON=.venv/bin/python` with `cd-hit-est` on PATH. This uses local cached inputs; it does not extract conservation or download data. The report, all fold scores, model coefficients and out-of-fold predictions are linked from [heldout_validation.md](../tables/heldout_validation.md). The fold snapshot is versioned independently of ignored data and generated tables; the manifest identifies their exact bytes.
