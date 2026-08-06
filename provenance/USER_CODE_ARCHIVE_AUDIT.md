# User code archive engineering audit

The supplied archive contained 107 files; the external manifest listed 103 code entries, 87 unique SHA-256 values and 16 duplicate-content entries. Total indexed code size was 1592570 bytes.

## Reused engineering patterns

- multi-format PDF/SVG/TIFF/PNG export and 183-mm proof output;
- stable publication themes and panel tags;
- Seurat v5 layer-safe access and conservative `JoinLayers`;
- explicit Source Data pointers, parameter tables and export manifests;
- PASS/WARNING/FAIL QA tables;
- immutable timestamped output roots and hashes;
- separation of formal, historical, exploratory and presentation-only branches.

## Not generalized

BHSC200121-specific genes, thresholds, patient counts, MHC-I/APM conclusions and drift checkpoints remain only in the project profile. They are not defaults for unrelated datasets.

## Manifest role counts

| language   | code_role                      |   n |
|:-----------|:-------------------------------|----:|
| Python     | supporting_automation          |   2 |
| R          | QA_or_diagnostic               |  15 |
| R          | annotation_or_cell_selection   |  26 |
| R          | biological_downstream_analysis |   4 |
| R          | figure_or_table_generation     |  38 |
| R          | orchestration                  |   2 |
| R          | other_R_code                   |  14 |
| Shell      | supporting_automation          |   2 |

## Virtual-knockout review for V5

No file name, manifest role or script content in the supplied archive identified a dedicated scTenifoldKnk, CellOracle, GEARS or other virtual-knockout implementation. V5 therefore reused only the archive's engineering architecture—frozen parameters, immutable branches, Source Data, manifests, multi-format independent figures and QA—and added `modules/27_virtual_knockout.R` as a new exploratory scientific module. This prevents unrelated BHSC200121 code from being misrepresented as a validated knockout engine.
