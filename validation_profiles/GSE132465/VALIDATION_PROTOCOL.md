# GSE132465 validation protocol

## Purpose

Validate human tumour annotation, paired pseudobulk, expression-inferred CNA, sample-wise communication and recurrent malignant-cell programs.

## Scientific units

- `sample_id`: one tumour or normal library;
- `patient_id`: paired patient identifier;
- `condition`: `Normal` or `Tumor`;
- formal comparison: paired pseudobulk with patient blocking;
- tumour-program discovery subset: predeclared tumour epithelial cells.

## Run order

1. core workflow;
2. paired pseudobulk review;
3. CopyKAT, one sample at a time;
4. CellChat, one sample at a time;
5. GeneNMF on the predeclared tumour epithelial subset;
6. hdWGCNA only after NMF/core QA;
7. drug-signature analysis only with a documented directional signature resource.

## Data acquisition guidance

Prefer an already downloaded official GEO supplementary file. For large files, use resumable downloads and preserve the original archive. Record SHA-256. Do not commit expression matrices to GitHub.

## Interpretation boundaries

- CopyKAT outputs are expression-inferred CNA, not DNA confirmation.
- CellChat outputs are model-derived communication potential.
- NMF program names are assigned only after marker/pathway review.
- hdWGCNA soft-power fallback requires explicit review when the target scale-free threshold is not reached.
