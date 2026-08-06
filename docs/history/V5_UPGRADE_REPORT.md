# V5 upgrade report — sample-aware virtual knockout

## Requirement implemented

V5 adds a complete virtual-knockout branch while preserving the V4 rule that every scientific display is an independent figure. The uploaded BHSC200121 code archive was inspected for reusable structure. It contained strong patterns for immutable inputs, timestamped branches, frozen parameters, Source Data, manifests and figure QA, but no dedicated virtual-knockout implementation. V5 therefore adds a new scientific module rather than relabelling unrelated project code.

## New execution module

`modules/27_virtual_knockout.R` implements an exploratory scTenifoldKnk workflow with:

- explicit target-gene declaration or target CSV;
- optional target-specific `subset_value` restriction;
- coherent cell-population and WT/control baseline gates;
- raw RNA count input and biological `sample_id` requirements;
- target-symbol, target-expression, sample-count and annotation-confidence audits;
- a frozen common gene universe with technical-gene exclusions;
- optional prior network with `regulators` and `targets` columns;
- per-biological-sample execution as the default;
- repeated **outer cell subsampling** before each backend call;
- optional balanced-pooled descriptive sensitivity analysis;
- backend-output schema validation;
- compatibility with historical/current WT/KO and X/Y network or manifold naming;
- sample-level repeat consensus and cross-sample recurrence;
- recurrent WT-network context, run similarity and optional pathway ORA;
- explicit `COMPLETED`, `FAILED`, `NOT_RUN` and `NOT_EVALUABLE` states.

## Why outer subsampling was added

The virtual-knockout backend controls parts of its own random process. V5 therefore draws an independent cell subset before every backend invocation using the predeclared `subsampling_fraction`. This makes the repeat analysis test sensitivity to observed-cell composition rather than assuming that changing only the external seed produces independent backend fits. The exact fraction is frozen and recorded in `virtual_knockout_parameters.yml`.

## Interpretation safeguards

The module treats scTenifoldKnk `distance`, `Z` and `FC` as measures derived from aligned network-state displacement. They are not converted into signed expression predictions. V5 forbids wording that a virtual knockout upregulated or downregulated a gene. It also forbids causal, viability, drug-response or therapeutic claims without real perturbation and orthogonal validation.

Evidence classes distinguish:

- `Replicate-consistent`;
- `Sample-restricted`;
- `Single-sample-descriptive`;
- `Pooled-descriptive`;
- `Not-supported`.

Cells are never treated as biological replicates.

## Independent figures added

For each evaluable target × cell subset, V5 can export separately:

1. unsigned regulatory-impact ranking;
2. biological-sample support heatmap;
3. target-expression/evaluability plot;
4. representative WT–KO manifold displacement plot;
5. recurrent target-centred WT-network context;
6. repeat-run similarity heatmap;
7. pathway over-representation plot.

Each figure uses the existing single-figure output contract: PDF, optional SVG, 600-dpi TIFF/PNG, 183-mm proof, thumbnail, plot RDS, parameter YAML, figure-specific Source Data and visual QA.

## Files added or materially changed

- `modules/27_virtual_knockout.R`;
- `recipes/27_virtual_knockout.md`;
- `templates/modules/27_virtual_knockout.yml`;
- `handbook/virtual_knockout_guide_zh.md`;
- `tests/VIRTUAL_KNOCKOUT_BENCHMARK_GSE167595.md`;
- `R/lib/advanced_figure_engine.R`;
- `R/11_figures.R`;
- `R/12_methods_legends.R`;
- `R/13_final_qa.R`;
- `run_all.R`;
- full/example config and JSON schema;
- dependency, static-validation and expected-output contracts.

## Validation boundary

The build environment does not contain `Rscript`. V5 has undergone static, configuration, route, file-contract, delimiter and archive checks, but has not completed a real scTenifoldKnk run here. The required next gate is local runtime testing, first on synthetic data and then on a public WT/KO benchmark such as GSE167595. Do not describe V5 as runtime validated until that gate passes.
