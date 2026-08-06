# GSE132465 validation profile

This profile documents the six-pair colorectal-cancer validation used for `v0.1.0-alpha`.

## Recorded design

- source cohort: 23 tumours and 10 matched normal mucosa samples;
- validation subset: SMC01–SMC06, six paired patients, 12 samples;
- input cells: 20,074;
- QC-retained cells: 17,391;
- environment: Windows, R 4.5.3, 32 GB RAM.

## Validated stages

- core tumour annotation;
- paired pseudobulk with `~ patient_id + condition`;
- sample-wise CopyKAT;
- sample-wise CellChat;
- GeneNMF recurrent programs.

hdWGCNA was stopped after TOM construction. Drug-response analysis was not run.

The original local preparation script was not supplied with this release archive. `00_prepare_GSE132465_template.R` is therefore a conservative template, not a byte-for-byte copy of the locally validated script. The immutable run evidence is under `docs/validation/`.
