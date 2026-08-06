# SingleCell Publication Workflow

**Version:** `v0.1.0-alpha`  
**Status:** experimental public alpha  
**License:** GPL-3.0-or-later for this repository's original code; third-party tools remain under their own terms.

A configuration-driven, sample-aware and auditable workflow for scRNA-seq/snRNA-seq analysis and publication-ready **independent figures**.

The project does not invent replacements for Seurat, CellChat, CopyKAT, GeneNMF, hdWGCNA, Slingshot or other established methods. It provides the integration layer around them: input adaptation, experimental-design checks, scientific gates, reproducible execution, failure handling, Source Data, figure contracts and final QA.

> **Alpha warning** — This is research software. The core pipeline and selected tumour modules have completed targeted public-data runtime validation, but not every input type, tissue, platform or advanced module has been validated. All annotations and model-derived results require expert review.

[中文说明](README.zh-CN.md) · [Repository](https://github.com/Potato-AI0815/singlecell-publication-workflow) · [Validation matrix](VALIDATION_MATRIX.md) · [Known limitations](KNOWN_LIMITATIONS.md) · [Third-party notices](THIRD_PARTY_NOTICES.md)

## Maintainer

Maintained by **Potato-AI**.

Security reports: 939191004@qq.com

## Why this project exists

Single-cell analysis often fails at the integration layer rather than the algorithm layer: metadata are misread, cells are treated as biological replicates, paired designs are ignored, advanced analyses are forced despite inadequate data, figures cannot be traced to Source Data, and failed stages are mistaken for successful runs.

This workflow is designed to prevent those failures.

## Current capabilities

### Core workflow

- 10x matrix directory and 10x HDF5 import;
- Seurat RDS, SingleCellExperiment RDS and optional H5AD import;
- multi-sample metadata manifests;
- sample-aware QC and optional scDblFinder;
- LogNormalize or SCTransform v2;
- PCA, optional Seurat v5 RPCA integration, clustering and UMAP;
- marker-program annotation with competing-lineage evidence and uncertainty labels;
- optional SingleR reference support with graceful degradation;
- descriptive cluster markers and sample-level recovered-cell composition;
- raw-count `sample × cell type` pseudobulk edgeR analysis;
- enrichment from pseudobulk-ranked statistics;
- Methods, legends, manifests, session information and final QA.

### Gated advanced adapters

- expression-inferred CNA with CopyKAT or infercna;
- sample-wise CellChat and replicated condition contrasts;
- Slingshot and optional tradeSeq;
- recurrent sample-aware NMF meta-programs;
- sample-preserving hdWGCNA;
- spatial transcriptomics, label transfer and optional RCTD;
- exploratory drug-signature analysis;
- exploratory scTenifoldKnk virtual knockout.

Each advanced module reports exactly one of:

```text
COMPLETED
NOT_RUN
NOT_EVALUABLE
FAILED
```

`NOT_EVALUABLE` is a scientific safeguard, not a negative biological result.

## Independent-figure contract

Composite figures, montages and contact sheets are forbidden. Every scientific display receives its own directory containing:

```text
<figure_id>.pdf
<figure_id>.svg                    # when svglite is available
<figure_id>_600dpi.tiff
<figure_id>_600dpi.png
<figure_id>_proof183mm.png
<figure_id>_thumbnail.png
<figure_id>_plot.rds
<figure_id>_parameters.yml
<figure_id>_source_<table>.csv[.gz]
visual_QA.csv
VISUAL_QA_STATUS.txt
```

The plot RDS is a lightweight renderable proxy when the original ggplot environment would capture large Seurat objects. Scientific values remain in figure-specific Source Data.

## Runtime validation completed

### GSE160763 — mouse cortex/TBI

- 8 samples, 4 conditions;
- 40,666 raw cells and 38,215 QC-retained cells;
- SCTransform/RPCA core workflow completed;
- 18 independent figure directories;
- no composite figures;
- `n=2` differential analysis correctly returned `NOT_EVALUABLE`;
- final status: `PASS_WITH_WARNINGS`;
- dataset-specific validator: `PASS`.

### GSE132465 — paired colorectal cancer subset

Six paired patients, 12 samples, 20,074 input cells and 17,391 QC-retained cells:

- core tumour annotation: `PASS`;
- paired pseudobulk design `~ patient_id + condition`: `PASS`;
- sample-wise CopyKAT CNA: `PASS`;
- sample-wise CellChat: `PASS`;
- recurrent GeneNMF programs: `PASS`;
- hdWGCNA: partial validation through metacells, soft-power diagnostics and TOM construction; module detection was stopped before completion;
- drug-response analysis: not executed.

See [VALIDATION_MATRIX.md](VALIDATION_MATRIX.md) for exact scope.


> **Restricted optional backends:** CopyKAT and scTenifoldKnk are not installed by the default dependency installer and are not redistributed. Review upstream terms before enabling them.

## Quick start on Windows

```powershell
run_one_click.cmd `
  -Config "C:\path\to\config.yml" `
  -Mode fast `
  -RscriptPath "D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe"
```

Install dependencies when needed:

```powershell
run_one_click.cmd `
  -Config "C:\path\to\config.yml" `
  -Mode fast `
  -RscriptPath "C:\path\to\Rscript.exe" `
  -InstallPackages `
  -InstallGitHubPackages
```

CopyKAT upstream namespace repair is opt-in:

```powershell
run_one_click.cmd ... -RepairCopykat
```

The launcher propagates the R exit code and accepts success only when `14_qa/FINAL_STATUS.txt` contains `PASS` or `PASS_WITH_WARNINGS`.

## Direct R usage

```bash
Rscript install_dependencies.R
Rscript install_dependencies.R --install
Rscript install_dependencies.R --install-github
Rscript run_all.R /absolute/path/config.yml --mode=full
```

Start from `templates/config.example.yml` or `templates/config.full.yml` and validate it:

```bash
python tests/validate_config.py project/config.yml templates/config.schema.json
```

## Scientific safeguards

- biological samples or patients are the inferential units;
- cells are never substituted for biological replication;
- integrated values are not used for formal differential expression;
- missing sample–cell-type combinations remain `NA`, not zero;
- unsupported fine annotation remains `Ambiguous` or `Unresolved`;
- trajectory roots, gene sets, model subsets and thresholds are declared before outcome inspection;
- expression-inferred CNA is not DNA confirmation;
- CellChat probabilities are not direct physical interactions or causality;
- transcriptomic drug scores are not treatment recommendations;
- virtual-knockout distances are not expression direction or real perturbation evidence.

## Repository structure

```text
R/                         core stages and figure engines
modules/                   gated advanced adapters
recipes/                   task-specific execution contracts
templates/                 configuration schema and examples
resources/                 marker dictionaries and resource templates
validation_profiles/       public-data profiles and validators
docs/validation/           immutable validation evidence
handbook/                  Chinese technical handbooks
tests/                     static and smoke-test contracts
scripts/                   launch, repair and audit utilities
```

## What is not claimed

This alpha does not claim universal annotation accuracy, production-grade support for every platform, clinical validity, causal inference, or completed validation of trajectory, spatial, hdWGCNA, drug-response and virtual-knockout workflows.

## Citation and third-party methods

Cite this workflow using `CITATION.cff`, and cite every upstream method actually used in an analysis. This repository does not redistribute third-party R package source code or proprietary reference data. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Contributing

Bug reports that include the resolved config, stage log, package versions and `FINAL_STATUS.txt` are especially useful. See [CONTRIBUTING.md](CONTRIBUTING.md).

## Commercial use

Commercial use of the GPL-covered community version is permitted subject to the license. Paid analysis, deployment, support and training services can be offered around the project. Separate proprietary modules must not incorporate third-party code or data contrary to their licenses. See [COMMERCIAL_USE.md](COMMERCIAL_USE.md).
