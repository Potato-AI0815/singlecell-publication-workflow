# Build report — v0.1.0-alpha

## Source integration

The release was assembled from the v5.1 public-data package, the runtime-validated v5.1.2 SKILL contract, GSE160763 QA evidence, and the GSE132465 run/bug reports documenting v5.1.3 fixes.

## Runtime evidence incorporated

- GSE160763 core SCT/RPCA workflow: `PASS_WITH_WARNINGS`;
- GSE160763 profile validator: `PASS`;
- GSE132465 six-pair core workflow: `PASS_WITH_WARNINGS`;
- paired pseudobulk: validated;
- CopyKAT: validated after upstream namespace and Windows-core fixes;
- CellChat: validated after self-loop rendering fix;
- GeneNMF: validated;
- hdWGCNA: partial only;
- drug-response: not executed.

## Code hardening integrated

- Windows CopyKAT single-core safeguard;
- explicit CopyKAT namespace repair utility;
- runtime dependency load audit;
- resume-safe `SKIPPED_EXISTING` figure handling;
- current-run versus historical failure separation;
- human reference-label mapping;
- self-loop-safe network plotting;
- lightweight renderable plot proxies;
- CLI mode override and Windows one-click launcher.

## Validation performed in this build container

- Python/static source validation;
- YAML/JSON configuration validation;
- Python compilation;
- archive integrity and SHA-256 generation;
- forbidden composite-layout scan;
- repository document and required-file checks.

## Validation not performed in this build container

R is not available in this container. The reconstructed alpha archive has therefore not been re-executed here. Runtime claims refer only to the supplied local validation evidence. A clean regression run from the final release archive remains a mandatory pre-publication check.
