# Roadmap

## v0.1.x — alpha hardening

- complete hdWGCNA public-data benchmark;
- add per-sample CopyKAT checkpointing or explicit restart manifest;
- run scDblFinder and SingleR success-path benchmarks;
- improve dependency pinning with `renv.lock` from a validated machine;
- add golden schema checks for GSE160763 and GSE132465;
- add human PBMC benchmark.

## v0.2.0 — additional modalities

- Slingshot/tradeSeq differentiation benchmark;
- Visium/RCTD benchmark;
- WT-only virtual-knockout discovery followed by held-out real-KO comparison;
- Linux/Ubuntu runtime validation;
- H5AD import validation.

## v0.3.0 — usability

- interactive config generator;
- HTML result index linking figures, Source Data, Methods and QA;
- improved progress and runtime estimates;
- validated tissue profiles;
- project-level regression testing.

## v1.0.0 criteria

- every advertised module has a public-data benchmark;
- Windows and Linux validated;
- frozen dependency environment;
- installer and upgrade path tested;
- all proof images reviewed;
- licensing and redistribution audit completed;
- clear support and maintenance policy.
