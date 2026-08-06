# GSE160763 public-data validation protocol

## Dataset role

This dataset is the first core-workflow runtime benchmark. It tests flat GEO 10x preparation, eight-sample import, sample-level metadata joining, sample-aware QC, doublet handling, normalization, integration, clustering, provisional CNS annotation, composition, independent figures and QA.

It is not a valid benchmark for tumour CNA, spatial transcriptomics, drug-response prediction or gene virtual knockout.

## Design boundary

Each sequencing sample is a pool of cortex from three mice. There are two pooled sequencing replicates per group. The sequencing-level inferential unit is therefore the pool/library (`n=2` per condition), not six individual mice. Formal pseudobulk testing remains gated at three replicates and should return `NOT_EVALUABLE` in the second run.

## Run order

1. Extract this skill package under `C:\Users\YHN\Desktop\Qoder\GSE160763_analysis\skill\`.
2. Keep `GSE160763_RAW.tar` at `C:\Users\YHN\Desktop\Qoder\GSE160763_RAW.tar`.
3. Run the PowerShell wrapper in smoke mode.
4. Review `GSE160763_validation_report.md`, all proof PNG files and the issue ledger.
5. Fix runtime errors without changing biological thresholds based on desired results.
6. After smoke PASS/WARNING, run `-Mode sct` to test SCTransform/RPCA and the insufficient-replication gate.

## Expected sanity signals, not hard-coded results

- eight samples and 40,666 pre-QC barcodes;
- major cortex populations including neurons, astrocytes, oligodendrocytes, microglia and vascular cells;
- markedly reduced recovered microglia in PLX5622-treated samples is biologically plausible, but is not a substitute for formal inference;
- formal DE should not be declared with only two pooled libraries per group.

## Source anchors

- GEO accession: GSE160763
- PMID: 33452227
- DOI: 10.1523/JNEUROSCI.2469-20.2020
