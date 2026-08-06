# Virtual-knockout benchmark: GSE167595

## Purpose

Use the Ahr WT/KO mouse colonic stem-cell dataset as the first biological benchmark for `modules/27_virtual_knockout.R`. The virtual model must be fitted using WT cells only. Real KO cells are held out until every target, subset, gene-universe and network parameter is frozen.

## Required metadata

Standardize at least:

- `cell_id`;
- `sample_id` (five WT and five KO biological samples where recoverable from the processed object);
- `condition` with exactly `WT` and `KO`;
- a reviewed coherent stem/epithelial cell annotation used as `subset_column`;
- `annotation_confidence`.

## Suggested module configuration

```yaml
advanced_modules:
  virtual_knockout:
    enabled: true
    backend: sctenifoldknk
    target_genes: [Ahr]
    subset_column: cell_type_l1
    subset_values: [REVIEWED_STEM_CELL_LABEL]
    condition_column: condition
    baseline_conditions: [WT]
    analysis_unit: sample_stratified_consensus
    minimum_eligible_samples: 2
    minimum_cells_per_sample: 300
    maximum_cells_per_sample: 3000
    minimum_target_expressing_fraction: 0.01
    minimum_target_expressing_cells: 10
    subsampling_repeats: 3
    subsampling_fraction: 0.8
    minimum_genes: 500
    maximum_genes: 2000
    n_networks: 10
    cells_per_network: 300
    workers: 4
```

Values above are a starting benchmark configuration, not universal defaults. Freeze them before inspecting KO data and record any scientifically justified change before rerunning.

## Validation sequence

1. Run core QC/annotation with both conditions, but create the virtual-knockout input gate using WT only.
2. Confirm that every WT sample passes target-expression and cell-count audits or is explicitly `NOT_EVALUABLE`.
3. Confirm repeated runs use distinct outer cell subsets and that the run manifest records cell counts and seeds indirectly through deterministic run order/parameters.
4. Freeze `virtual_knockout_consensus_genes.csv` and pathway results.
5. Only then perform a real WT-versus-KO sample-aware pseudobulk comparison.
6. Compare virtual impact rankings with real-KO evidence using unsigned overlap, rank enrichment and pathway concordance. Do not compare them as signed predicted fold changes.
7. Record runtime, peak memory, failed runs, package versions and every independent figure QA status.

## Acceptance criteria

- no KO cell enters the virtual model;
- no target, cell subset, threshold or gene universe is changed after KO inspection;
- at least two eligible WT biological samples are required for `Replicate-consistent` wording;
- same-seed reruns preserve manifests and numerical checkpoints within declared tolerance;
- every requested virtual-knockout plot is an independent figure with Source Data;
- results remain labelled exploratory and unsigned.
