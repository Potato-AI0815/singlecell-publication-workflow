# Package contents

Package: `singlecell-publication-workflow`
Version: `v0.1.0-alpha`

- Files: **137**
- Release type: experimental GitHub alpha
- Composite figures: forbidden by the execution contract

Use `PACKAGE_MANIFEST.csv` for per-file SHA-256 hashes (the manifest excludes its own self-hash).

## (root)

- `.gitattributes`
- `.gitignore`
- `BUILD_REPORT.md`
- `CHANGELOG.md`
- `CITATION.cff`
- `CODE_OF_CONDUCT.md`
- `COMMERCIAL_USE.md`
- `CONTENTS.md`
- `CONTRIBUTING.md`
- `COPYRIGHT`
- `DEPENDENCY_LICENSES.csv`
- `DISCLAIMER.md`
- `GITHUB_UPLOAD_GUIDE.md`
- `KNOWN_LIMITATIONS.md`
- `LICENSE`
- `MIGRATION_NOTES.md`
- `PACKAGE_MANIFEST.csv`
- `README.md`
- `README.zh-CN.md`
- `RELEASE_NOTES_v0.1.0-alpha.md`
- `ROADMAP.md`
- `SECURITY.md`
- `SKILL.md`
- `THIRD_PARTY_NOTICES.md`
- `VALIDATION_MATRIX.md`
- `VERSION`
- `freeze_environment.R`
- `install_dependencies.R`
- `run_all.R`
- `run_one_click.cmd`

## .github

- `.github/ISSUE_TEMPLATE/bug_report.yml`
- `.github/ISSUE_TEMPLATE/feature_request.yml`
- `.github/workflows/static-checks.yml`

## R

- `R/00_utils.R`
- `R/01_preflight.R`
- `R/02_import_standardize.R`
- `R/03_qc_doublets.R`
- `R/04_normalize_reduce.R`
- `R/05_integrate_cluster.R`
- `R/06_annotate.R`
- `R/07_markers.R`
- `R/08_composition.R`
- `R/09_pseudobulk_de.R`
- `R/10_enrichment.R`
- `R/11_figures.R`
- `R/12_methods_legends.R`
- `R/13_final_qa.R`
- `R/lib/advanced_figure_engine.R`
- `R/lib/figure_engine.R`

## docs

- `docs/ARCHITECTURE.md`
- `docs/RELEASE_CHECKLIST.md`
- `docs/history/README_V5.1_GSE160763.md`
- `docs/history/V4_UPGRADE_REPORT.md`
- `docs/history/V5.1_GSE160763_REPORT.md`
- `docs/history/V5_UPGRADE_REPORT.md`
- `docs/validation/GSE132465_RUN_AND_BUG_LOG.md`
- `docs/validation/GSE132465_TEST_FINDINGS.md`
- `docs/validation/GSE132465_TEST_REPORT_SOURCE.md`
- `docs/validation/GSE160763_final_QA_report.md`
- `docs/validation/GSE160763_validation_report.md`

## handbook

- `handbook/advanced_analysis_catalog_zh.md`
- `handbook/figure_engine_catalog_zh.md`
- `handbook/output_contract.md`
- `handbook/virtual_knockout_guide_zh.md`
- `handbook/visualization_cheatsheet_zh.md`

## modules

- `modules/20_cnv.R`
- `modules/21_communication_cellchat.R`
- `modules/22_trajectory_slingshot.R`
- `modules/23_nmf_programs.R`
- `modules/24_hdwgcna.R`
- `modules/25_spatial.R`
- `modules/26_drug_response.R`
- `modules/27_virtual_knockout.R`

## profiles

- `profiles/BHSC200121.md`

## provenance

- `provenance/CODE_PACKAGE_MANIFEST.csv`
- `provenance/USER_CODE_ARCHIVE_AUDIT.md`

## recipes

- `recipes/01_full_basic_analysis.md`
- `recipes/02_fast_first_pass.md`
- `recipes/03_reannotation_only.md`
- `recipes/04_compare_celltype.md`
- `recipes/05_figures_only.md`
- `recipes/06_audit_existing.md`
- `recipes/20_cnv.md`
- `recipes/21_communication.md`
- `recipes/22_trajectory.md`
- `recipes/23_nmf_programs.md`
- `recipes/24_hdwgcna.md`
- `recipes/25_spatial.md`
- `recipes/26_drug_response.md`
- `recipes/27_virtual_knockout.md`
- `recipes/30_publication_figure_engine.md`

## resources

- `resources/drug_signature_reversal_template.csv`
- `resources/drug_signature_score_template.csv`
- `resources/marker_dictionary_human_l1.csv`
- `resources/marker_dictionary_human_l2.csv`
- `resources/marker_dictionary_mouse_cortex_l1.csv`
- `resources/marker_dictionary_mouse_cortex_l2.csv`
- `resources/marker_dictionary_mouse_l1.csv`

## scripts

- `scripts/audit_installed_licenses.R`
- `scripts/run_one_click.ps1`

## templates

- `templates/config.example.yml`
- `templates/config.full.yml`
- `templates/config.schema.json`
- `templates/invocation_prompt.md`
- `templates/modules/20_cnv.yml`
- `templates/modules/21_communication.yml`
- `templates/modules/22_trajectory.yml`
- `templates/modules/23_nmf_programs.yml`
- `templates/modules/24_hdwgcna.yml`
- `templates/modules/25_spatial.yml`
- `templates/modules/26_drug_response.yml`
- `templates/modules/27_virtual_knockout.yml`
- `templates/modules/README.md`

## tests

- `tests/00_generate_smoke_input.R`
- `tests/02_validate_smoke_outputs.R`
- `tests/EXPECTED_OUTPUT_CONTRACT.csv`
- `tests/PUBLIC_DATA_VALIDATION_PLAN.md`
- `tests/VIRTUAL_KNOCKOUT_BENCHMARK_GSE167595.md`
- `tests/run_static_suite.sh`
- `tests/static_validate.py`
- `tests/validate_config.py`

## validation_profiles

- `validation_profiles/GSE132465/00_prepare_GSE132465_template.R`
- `validation_profiles/GSE132465/03_validate_GSE132465_outputs.R`
- `validation_profiles/GSE132465/README.md`
- `validation_profiles/GSE132465/VALIDATION_PROTOCOL.md`
- `validation_profiles/GSE132465/build_drug_signature_msigdbr.R`
- `validation_profiles/GSE132465/expected_checkpoints.csv`
- `validation_profiles/GSE160763/00_prepare_GSE160763.R`
- `validation_profiles/GSE160763/01_install_core_dependencies.R`
- `validation_profiles/GSE160763/02_run_validation.ps1`
- `validation_profiles/GSE160763/03_validate_GSE160763_outputs.R`
- `validation_profiles/GSE160763/LOCAL_MODEL_PROMPT.md`
- `validation_profiles/GSE160763/VALIDATION_PROTOCOL.md`
- `validation_profiles/GSE160763/config_01_core_smoke.yml`
- `validation_profiles/GSE160763/config_02_core_sct_gate_test.yml`
- `validation_profiles/GSE160763/expected_group_counts.csv`
- `validation_profiles/GSE160763/sample_metadata.csv`
