# Architecture

## Core principle

The workflow separates four responsibilities:

1. **scientific result generation** — core stages and gated advanced modules;
2. **orchestration** — configuration, resume markers and failure propagation;
3. **presentation** — independent figure templates and Source Data;
4. **audit** — manifests, issue ledger and final QA.

## Execution sequence

```text
01 preflight
02 import/standardize
03 QC/doublets
04 normalize/reduce
05 integrate/cluster
06 annotate
07 markers
08 composition
09 pseudobulk DE
10 enrichment
advanced modules
11 independent figures
12 Methods/legends
13 final QA
```

## State model

Core stages use `.done` files only after successful completion. Advanced modules additionally write status rows and module-specific completion markers. A failed stage exits non-zero.

Resume mode reuses completed stages. Existing complete figure directories are marked `SKIPPED_EXISTING`; incomplete or unmanifested figure directories remain blocking errors.

## Data and statistical boundaries

- raw counts remain in the RNA assay;
- integration is used for alignment and visualization, not formal DE;
- formal comparisons operate on sample-level pseudobulk;
- missing biological combinations remain missing;
- advanced module subsets and parameters are explicit configuration inputs;
- all model-derived results retain interpretation disclaimers.
