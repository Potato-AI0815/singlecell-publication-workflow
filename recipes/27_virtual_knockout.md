# Recipe 27 — exploratory virtual knockout

Use this recipe when the user asks for virtual knockout, in-silico knockout, GRN perturbation, target-gene perturbation prediction, or scTenifoldKnk.

## Required routing

Run `modules/27_virtual_knockout.R` through the normal advanced-module dispatcher. Do not write an ad hoc knockout script when this module is available.

## Required user/config inputs

1. One or more predeclared target genes in `target_genes` or `targets_path` (`target_gene,target_label,subset_value,priority,rationale`; `subset_value` is optional).
2. A biologically coherent cell population in `subset_values`.
3. Explicit WT/control `baseline_conditions` when multiple conditions are present.
4. Raw RNA counts retained in the annotated Seurat object.
5. Biological `sample_id` metadata.

## Scientific gates

- The target must exist in the selected gene universe.
- Each evaluable sample must pass minimum cell number, expressing-cell number and expressing-fraction gates.
- High/Medium annotation confidence is used by default.
- Sample-stratified repeated runs are the default; pooled analysis is descriptive only.
- Each repeat must draw an independent outer cell subset using `subsampling_fraction` before invoking the backend. This is a stability design and is recorded in the run manifest.
- At least two eligible biological samples are required before wording a result as replicate-consistent.
- Freeze targets, cell subset, baseline condition, gene universe and all network parameters before viewing results.

## Backend

Primary backend: `scTenifoldKnk`.

It constructs a WT scGRN, removes the target gene's outgoing edges, aligns the WT and KO network states, and reports differential-regulation distances. These distances are unsigned. Never describe them as predicted expression induction or repression.

Do not silently substitute:

- CellOracle, which requires a separate Python/GRN/TF-prior and trajectory-aware workflow;
- GEARS, which requires perturbational training data rather than WT-only observational scRNA-seq.

## Required outputs

- target-symbol audit;
- target-expression and target/subset evaluability audits;
- frozen gene universe;
- run manifest with completed/failed/unevaluable runs;
- all differential-regulation results;
- sample-level repeat consensus;
- cross-sample consensus and evidence class;
- manifold coordinates;
- recurrent WT-network edges;
- run-similarity table;
- optional pathway ORA;
- parameter YAML and module status.

## Independent figures only

For each target × cell subset, export separately when evaluable:

1. unsigned regulatory-impact ranking;
2. support-across-samples heatmap;
3. target-expression eligibility plot;
4. representative WT–KO manifold displacement plot;
5. recurrent target-centred WT-network context;
6. run-similarity heatmap;
7. pathway enrichment plot.

Never combine these into a composite panel.

## Interpretation contract

Use language such as:

> Virtual knockout of TARGET was associated with recurrent network-state displacement involving GENE SET in CELL TYPE across N of M evaluable samples.

Do not write:

> TARGET knockout upregulated/downregulated GENE.

Do not infer causality, viability, drug response or therapeutic benefit without real perturbation and orthogonal validation.
