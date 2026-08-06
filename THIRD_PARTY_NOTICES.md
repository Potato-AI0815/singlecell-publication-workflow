# Third-party notices

This repository's original orchestration, scientific gates, adapters, figure contracts and QA code are licensed separately from the tools they call.

No third-party package source code, model weight, proprietary database or public-dataset matrix is intentionally redistributed in this repository. Dependencies are installed from their upstream distribution channels.

## Major optional/runtime components

| Component | Use in this workflow | Upstream license status used for this release | Redistribution policy here |
|---|---|---|---|
| Seurat / SeuratObject | Core object model and analysis | MIT | External dependency; not vendored |
| CellChat | Communication inference | GPL-3.0 | External optional dependency; not vendored |
| hdWGCNA | Co-expression networks | GPL-3.0 | External optional dependency; not vendored |
| spacexr/RCTD | Spatial deconvolution | GPL-3.0 | External optional dependency; not vendored |
| CopyKAT | Expression-inferred CNA | Restrictive upstream copyright terms: academic/non-profit research use; commercial use/transfers require written approval; sale/redistribution prohibited without approval | External optional dependency; not redistributed; commercial users must obtain upstream permission |
| scTenifoldKnk | Virtual knockout backend | Upstream repository displays “All rights reserved” and no standard open-source license was confirmed | Manual user-installed optional backend; not auto-installed or redistributed |
| GeneNMF | Recurrent programs | MIT | External optional dependency; not vendored |
| Slingshot / tradeSeq / edgeR / SingleR / scDblFinder | Bioconductor methods | Read each installed package’s `License` field | External dependencies; not vendored |
| MSigDB-derived gene sets | Optional drug/signature resource | Separate resource terms apply | Generated/downloaded by the user; not bundled |

## Installed-environment audit

Run:

```bash
Rscript scripts/audit_installed_licenses.R installed_dependency_licenses.csv
```

The resulting table records the `License` and `URL` fields from the exact installed package versions. Review it before redistribution or paid deployment.

## Citations

Users must cite each upstream method actually used. The generated Methods and manifest files should be reviewed to determine the applicable citations.
