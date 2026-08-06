# Known limitations

## Product maturity

`v0.1.0-alpha` is research software. It has completed two targeted public-data validation programs but is not a universally validated analysis platform.

## Annotation

- Marker-program annotation is provisional and must be reviewed against tissue context, canonical markers and known biology.
- SingleR requires an available compatible reference; network or cache failure triggers a marker-only warning rather than a fabricated reference result.
- Tissue-specific fine annotation dictionaries are not comprehensive.
- Annotation agreement with an author label is not proof of biological truth.

## Statistics

- Formal between-condition testing requires biological replication and raw-count pseudobulk.
- Pooled animals remain one library-level replicate when pooling occurred before sequencing.
- `NOT_EVALUABLE` is expected when design gates fail.
- Multiple testing and small-sample power can produce no FDR-significant results even when descriptive directions exist.

## CopyKAT/CNA

- CopyKAT is expression-inferred CNA and can classify some histologically normal epithelial cells as aneuploid-like.
- CopyKAT is forced to one core on Windows because its `mclapply` path does not support multiple cores there.
- The tested upstream snapshot had a namespace packaging problem involving `sysdata.rda`; an explicit repair utility is provided, but upstream changes may alter this issue.
- CopyKAT currently lacks robust per-sample checkpointing in this wrapper; interruption can require re-running the module.
- DNA-based or orthogonal validation is required before genomic claims.

## CellChat

- CellChat results describe model-derived communication potential, not direct physical interaction, signaling flux or causality.
- Self-loop rendering is supported, but dense networks can remain visually crowded.
- Sample-wise execution can be slow, and condition contrasts require replicated samples.

## NMF and hdWGCNA

- GeneNMF meta-program recurrence depends on the selected population, rank grid, cell counts and sample representation.
- A program seen in only one sample must not be called recurrent.
- hdWGCNA can require substantial memory and runtime. The GSE132465 validation created a roughly 433 MB TOM for about 7,400 genes.
- No candidate soft power met an SFT threshold of 0.8 in the partial CRC run; fallback rules require review before module results are accepted.
- hdWGCNA module detection and final outputs have not yet completed public-data validation.

## Spatial, drug response and virtual knockout

- Spatial workflows are platform-specific and have not completed a dedicated public benchmark in this release.
- Drug signatures and model resources may have separate use and redistribution terms. Scores cannot be interpreted as treatment recommendations.
- scTenifoldKnk is an optional externally installed backend and is not redistributed. Its results are unsigned network-state displacement, not expression direction or real perturbation evidence.

## Figures

- Technical visual QA detects file and rendering problems but does not replace human aesthetic review.
- Plot RDS files may contain a lightweight grob-backed proxy rather than the original data-bearing ggplot. Full scientific values remain in Source Data.
- Composite figures are intentionally outside this project’s scope.

## Platforms and dependencies

- The recorded runtime validations used Windows and R 4.5.3.
- R/Bioconductor/GitHub package installation can change over time.
- Some optional packages may install but fail to load until transitive dependencies are installed.
- No guarantee is made that current upstream development branches remain API-compatible.
