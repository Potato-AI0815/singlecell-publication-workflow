#!/usr/bin/env python3
from pathlib import Path
import sys, re, json

root = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1])
required = [
    'run_all.R','install_dependencies.R','run_one_click.cmd','scripts/run_one_click.ps1','SKILL.md','README.md','README.zh-CN.md','VERSION','LICENSE',
    'R/00_utils.R','R/01_preflight.R','R/02_import_standardize.R','R/03_qc_doublets.R',
    'R/04_normalize_reduce.R','R/05_integrate_cluster.R','R/06_annotate.R','R/07_markers.R',
    'R/08_composition.R','R/09_pseudobulk_de.R','R/10_enrichment.R','R/11_figures.R',
    'R/12_methods_legends.R','R/13_final_qa.R','R/lib/figure_engine.R','R/lib/advanced_figure_engine.R',
    'modules/20_cnv.R','modules/21_communication_cellchat.R','modules/22_trajectory_slingshot.R',
    'modules/23_nmf_programs.R','modules/24_hdwgcna.R','modules/25_spatial.R','modules/26_drug_response.R','modules/27_virtual_knockout.R',
    'templates/config.full.yml','templates/config.example.yml','templates/config.schema.json',
    'handbook/figure_engine_catalog_zh.md','handbook/advanced_analysis_catalog_zh.md','handbook/virtual_knockout_guide_zh.md',
    'recipes/23_nmf_programs.md','recipes/24_hdwgcna.md','recipes/25_spatial.md','recipes/26_drug_response.md','recipes/27_virtual_knockout.md',
    'recipes/30_publication_figure_engine.md',
    'templates/modules/README.md','templates/modules/20_cnv.yml','templates/modules/21_communication.yml',
    'templates/modules/22_trajectory.yml','templates/modules/23_nmf_programs.yml','templates/modules/24_hdwgcna.yml',
    'templates/modules/25_spatial.yml','templates/modules/26_drug_response.yml','templates/modules/27_virtual_knockout.yml',
    'tests/PUBLIC_DATA_VALIDATION_PLAN.md','tests/VIRTUAL_KNOCKOUT_BENCHMARK_GSE167595.md','tests/EXPECTED_OUTPUT_CONTRACT.csv','tests/run_static_suite.sh',
    'resources/drug_signature_reversal_template.csv','resources/drug_signature_score_template.csv',
    'VALIDATION_MATRIX.md','KNOWN_LIMITATIONS.md','THIRD_PARTY_NOTICES.md','DEPENDENCY_LICENSES.csv','CITATION.cff','CHANGELOG.md'
]
issues=[]
missing=[x for x in required if not (root/x).exists()]
if missing: issues.append('missing: '+', '.join(missing))

# Lightweight delimiter scan: quoted strings and comments are removed first.
for p in root.rglob('*.R'):
    t=p.read_text(errors='ignore')
    rel = p.relative_to(root)
    if 'validation_profiles' not in rel.parts and re.search(r'[A-Za-z]:/Users/|/home/[^/]+/|/data2/',t):
        issues.append(f'absolute project path: {rel}')
    if p.parent.name=='R' and re.match(r'\d\d_',p.name) and p.name!='00_utils.R' and 'stage_main <- function(ctx)' not in t:
        issues.append(f'missing stage_main: {p.relative_to(root)}')
    s=re.sub(r'"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'','',t)
    s='\n'.join(x.split('#')[0] for x in s.splitlines())
    for a,b in [('(',')'),('{','}'),('[',']')]:
        if s.count(a)!=s.count(b): issues.append(f'unbalanced {a}{b}: {p.relative_to(root)} {s.count(a)} {s.count(b)}')

core=(root/'R/lib/figure_engine.R').read_text(errors='ignore') if (root/'R/lib/figure_engine.R').exists() else ''
advanced=(root/'R/lib/advanced_figure_engine.R').read_text(errors='ignore') if (root/'R/lib/advanced_figure_engine.R').exists() else ''
core_fns=['plot_embedding_discrete','plot_embedding_continuous','plot_split_violin_box','plot_marker_dotplot',
          'plot_marker_heatmap','plot_trajectory_embedding','plot_pseudotime_trend','plot_volcano_single',
          'plot_enrichment_single','export_publication_figure','run_plot_qa']
adv_fns=['plot_advanced_heatmap','plot_rank_diagnostics','plot_program_prevalence_dotplot','plot_lollipop_ranking',
         'plot_cnv_heatmap','plot_cnv_scatter','plot_communication_bubble','plot_network_edges',
         'plot_soft_power_metric','plot_hub_genes','plot_spatial_discrete','plot_spatial_continuous',
         'plot_drug_score_heatmap','plot_drug_reversal_scatter','plot_metric_curve','plot_gene_test_ranking',
         'plot_group_score_distribution','plot_virtual_knockout_ranking','plot_virtual_knockout_sample_support',
         'plot_virtual_knockout_manifold','plot_virtual_knockout_network','plot_virtual_knockout_pathways',
         'plot_virtual_knockout_target_expression']
for fn in core_fns:
    if not re.search(rf'\b{re.escape(fn)}\s*<-\s*function\s*\(',core): issues.append(f'missing core figure function: {fn}')
for fn in adv_fns:
    if not re.search(rf'\b{re.escape(fn)}\s*<-\s*function\s*\(',advanced): issues.append(f'missing advanced figure function: {fn}')

run_all=(root/'run_all.R').read_text(errors='ignore') if (root/'run_all.R').exists() else ''
for token in ['advanced_figure_engine.R','20_cnv.R','21_communication_cellchat.R','22_trajectory_slingshot.R',
              '23_nmf_programs.R','24_hdwgcna.R','25_spatial.R','26_drug_response.R','27_virtual_knockout.R']:
    if token not in run_all: issues.append(f'run_all.R missing route: {token}')

# Result-directory aliases must remain aligned. The audit ledger intentionally
# aliases 13_logs while the numbered output contract stays unchanged.
utils_path = root/'R/00_utils.R'
if utils_path.exists():
    utils_text = utils_path.read_text(errors='ignore')
    keys_match = re.search(r'keys\s*<-\s*c\((.*?)\)', utils_text, flags=re.S)
    vals_match = re.search(r'vals\s*<-\s*c\((.*?)\)', utils_text, flags=re.S)
    if not keys_match or not vals_match:
        issues.append('result-directory key/value declarations are missing')
    else:
        dir_keys = re.findall(r'"([^"]+)"', keys_match.group(1))
        dir_vals = re.findall(r'"([^"]+)"', vals_match.group(1))
        if len(dir_keys) != len(dir_vals):
            issues.append(f'result-directory key/value length mismatch: {len(dir_keys)} keys vs {len(dir_vals)} values')
        else:
            dir_map = dict(zip(dir_keys, dir_vals))
            if dir_map.get('audit') != '13_logs': issues.append('audit directory must alias 13_logs')
            if dir_map.get('manifests') != '15_manifests': issues.append('manifests directory must remain 15_manifests')
            if dir_map.get('advanced') != '17_advanced_modules': issues.append('advanced directory must remain 17_advanced_modules')


module27=(root/'modules/27_virtual_knockout.R').read_text(errors='ignore') if (root/'modules/27_virtual_knockout.R').exists() else ''
for token in ['scTenifoldKnk::scTenifoldKnk','sample_stratified_consensus','virtual_knockout_consensus_genes.csv','virtual_knockout_target_expression_audit.csv','virtual_knockout_target_subset_evaluability.csv','subsampling_fraction','vk_prepare_backend_result','vk_get_tensor_network']:
    if token not in module27: issues.append(f'virtual-knockout module missing contract token: {token}')
if re.search(r'predicted\s+(up|down)[- ]?regulat', module27, flags=re.I):
    issues.append('virtual-knockout module contains forbidden directional-expression claim')

# No aggregate-layout generators may remain in executable R code.
all_r='\n'.join(p.read_text(errors='ignore') for p in root.rglob('*.R'))
if re.search(r'\bSeurat::JoinLayers\s*\(', all_r):
    issues.append('JoinLayers must be called from SeuratObject for Seurat 5 compatibility')
if not re.search(r'if \(nrow\(tab\) == 0L\) next', all_r):
    issues.append('empty figure source tables must not be registered in the manifest')
for pat,label in [(r'patchwork::wrap_plots\s*\(','patchwork::wrap_plots'),(r'cowplot::plot_grid\s*\(','cowplot::plot_grid'),
                  (r'ggpubr::ggarrange\s*\(','ggpubr::ggarrange'),(r'image_montage\s*\(','image_montage')]:
    if re.search(pat,all_r): issues.append(f'forbidden aggregate-layout generator: {label}')

cfg_path=root/'templates/config.full.yml'
if cfg_path.exists():
    cfg=cfg_path.read_text(errors='ignore')
    for key in ['cnv:','communication:','trajectory:','nmf:','hdWGCNA:','spatial:','drug_response:','virtual_knockout:',
                'export_individual_figures: true','composite_figures_forbidden: true','build_contact_sheet: false']:
        if key not in cfg: issues.append(f'config.full.yml missing: {key}')

schema_path=root/'templates/config.schema.json'
if schema_path.exists():
    try:
        schema=json.loads(schema_path.read_text())
        props=schema['properties']['advanced_modules']['properties']
        for key in ['cnv','communication','trajectory','nmf','hdWGCNA','spatial','drug_response','virtual_knockout']:
            if key not in props: issues.append(f'config schema missing advanced module: {key}')
        vkprops=props.get('virtual_knockout',{}).get('properties',{})
        if 'subsampling_fraction' not in vkprops: issues.append('config schema missing virtual_knockout.subsampling_fraction')
        fprops=schema['properties']['figures']['properties']
        if fprops.get('export_individual_figures',{}).get('const') is not True: issues.append('schema does not require independent figures')
        if fprops.get('composite_figures_forbidden',{}).get('const') is not True: issues.append('schema does not forbid aggregate layouts')
        if fprops.get('build_contact_sheet',{}).get('const') is not False: issues.append('schema does not disable contact sheet')
    except Exception as e: issues.append(f'invalid JSON schema: {e}')


# Alpha hardening contracts.
version=(root/'VERSION').read_text(errors='ignore').strip() if (root/'VERSION').exists() else ''
if version != '0.1.0-alpha': issues.append(f'unexpected VERSION: {version}')
checks = {
    'Windows CopyKAT core safeguard': (root/'modules/20_cnv.R', '.Platform$OS.type == "windows"'),
    'resume figure reuse': (root/'R/lib/figure_engine.R', 'SKIPPED_EXISTING'),
    'lightweight plot proxy': (root/'R/lib/figure_engine.R', 'make_renderable_plot_proxy'),
    'historical failure separation': (root/'R/13_final_qa.R', 'historical_FAIL'),
    'network self-loop handling': (root/'R/lib/advanced_figure_engine.R', 'n_self_loops'),
    'mode override': (root/'run_all.R', '--mode='),
}
for label,(path,token) in checks.items():
    if not path.exists() or token not in path.read_text(errors='ignore'):
        issues.append(f'missing alpha hardening: {label}')

print('STATIC_VALIDATION','PASS' if not issues else 'FAIL')
for x in issues: print('-',x)
sys.exit(1 if issues else 0)
