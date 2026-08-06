# GSE132465 测试运行 — 发现记录（持续更新）

数据集：GSE132465（23 例韩国 CRC 患者，23 肿瘤 + 10 配对正常，63,689 细胞）
测试子集：SMC01–SMC06 六对配对（12 样本，20,074 细胞，QC 后 17,391）
Skill：single-cell-publication-workflow v5.1.2
运行目录：results/GSE132465_CRC_6pair_20260805_120916

## 已确认问题

### B1【严重】copykat 上游打包 bug：`sysdata.rda` 放错位置
- 现象：每个样本 copykat 调用报 `object 'full.anno' not found`，CNA 模块 NOT_EVALUABLE。
- 根因：navinlabcode/copykat master（v1.2.5, 2026-06-12）将 `sysdata.rda`（含 full.anno、full.anno.mm10、DNA.hg20、cyclegenes）放在 `data/` 而非 `R/`，内部对象不进入命名空间。
- 修复：`mv data/sysdata.rda R/sysdata.rda` 后重装，命名空间对象恢复，copykat 可运行。
- SKILL 优化建议：install_dependencies.R 安装 copykat 后做 namespace 自检（exists('full.anno', asNamespace('copykat'))），失败时自动打补丁或给出明确指引；SKILL.md/文档注明该已知上游问题。

### B2【中】resume 模式下 stage 11 必然失败（图已存在）
- 现象：resume 新增高级模块后，stage 11 报 "The figure engine did not export any figure"。
- 根因：figure_engine 对非空图目录 `stop("Refusing to overwrite...")`（规则1保护），所有已存在图被拒 → 0 导出 → 11_figures.R:617 报错。
- 影响：与"逐模块 resume 测试"推荐工作流直接冲突。
- 修复建议：export_publication_figure 遇已存在非空图目录时返回 SKIPPED_EXISTING（保持 manifest 内重复检测为错误）；stage 11 允许 EXPORTED + SKIPPED_EXISTING 全部覆盖即通过。

### B3【中】map_reference_l1 缺少人类谱系规则（小鼠 CNS 偏向）
- 现象/推断：SingleR(HumanPrimaryCellAtlas) 的 "Epithelial cells" 无映射规则 → reference_unavailable；"Fibroblast" 被映射为小鼠专用标签 "Fibroblast/VLMC"、"Pericyte" 映射为 "Pericyte"（字典为 "Pericyte/smooth muscle"）→ 与字典标签不等 → agreement=conflict → 注释被降级 Ambiguous，而 DE 默认排除 Ambiguous 之外的规则会把关键基质群体从 pseudobulk 排除。
- 本次运行 SingleR 离线（celldex 下载失败）未实际触发，但代码路径明确。
- 修复建议：为 human_l1 字典补齐规则（Epithelial、Fibroblast、Pericyte/smooth muscle、Monocyte、Macrophage、Endothelial→已有、Smooth muscle 等），输出标签与 human_l1 字典一致；或按 species 分派映射表。

### B4【轻】药物响应资源仅有占位模板
- resources/drug_signature_reversal_template.csv 仅含 ExampleDrugA/GENE1 占位，signature_reversal 后端开箱无法做真实测试。
- 本次用 msigdbr v26.1.0 C2:CGP 的 _UP/_DN 成对集构造了 1,219 个方向性签名（269,867 行），含 CRC 相关（KANG_FLUOROURACIL_RESISTANCE、CALVET_IRINOTECAN_*、多个 CISPLATIN 集）。
- 优化建议：SKILL 补充"合规签名资源构建指引"（MSigDB/CMap 来源、方向性定义、provenance 要求），并将本脚本作为示例纳入。

### B5【轻】依赖链缺口（文档/检查层面）
- WGCNA（CRAN 二进制）安装成功但加载失败，因缺 Bioc impute/preprocessCore；install_dependencies 的 missing 检查能发现，但未提示"需先装 Bioc 依赖"。
- scDblFinder 因缺 bluster/KEGG.R 加载失败，QC 阶段优雅降级（WARNING，保留 singlets）——降级行为符合契约；但 install_dependencies 可加运行时加载自检。
- KEGG.R 在 R 4.5/Bioc 3.22 无可用版本（不影响 scDblFinder 加载）。

### B6【轻】GEO 大数据下载无指引
- NCBI FTP 单连接 ~2MB/min，个别 range 长时间停滞；建议 validation profile 文档注明"优先本地文件/镜像，分块下载 + 断点续传，FTP 协议备选"。

### B7【严重】copykat n.cores>1 在 Windows 直接报错
- 现象：`workers: 4` 时每样本 copykat 报 `'mc.cores' > 1 is not supported on Windows`（parallel::mclapply 在 Windows 的硬错误），全部样本跳过 → NOT_EVALUABLE。
- 已修复（modules/20_cnv.R）：`n.cores = if (.Platform$OS.type == "windows") 1L else (cfg$workers %||% 1)`。
- SKILL 优化建议：config 模板注释注明 Windows 单核限制。

### B8【中】issue ledger 历史 FAIL 条目永久阻塞已修复的运行
- 现象：resume 轮次修复缺陷后，final QA "No recorded blocking failure" 仍把修复前运行写入的历史 FAIL 条目计入 → FINAL_STATUS=FAIL（本测试中为 B2 修复前的两条 stage 11 失败记录）。
- 已修复（R/13_final_qa.R）：以本次运行 run_manifest start_time 为界，历史 FAIL 保留在台账与报告但不再阻塞；仅当前运行内的 FAIL 阻断；证据字段报告 current_FAIL / historical_FAIL 计数。

### B9【中】网络图模板在自分泌边上崩溃
- 现象：CellChat 轮次 Communication_network_Normal/Tumor 两图导出失败："end points must not be identical" → stage 13 FAIL。
- 根因：advanced_figure_engine.R plot_network_edges 用 geom_curve 画边，CellChat 天然含 source==target 的自分泌边 → 端点重合 → 图形设备报错。
- 已修复：自分泌边端点对称外移 0.16 渲染为小弧，parameters 记录 n_self_loops。

## 已实施的 SKILL 修复（本轮测试中）

- F1【B2】figure_engine.R：已存在的非空图目录返回 SKIPPED_EXISTING（不再 stop）；manifest 重复检测保留为错误；11_figures.R try_export_figure 记录 SKIPPED_EXISTING；stage 11 与 13_final_qa 接受 EXPORTED+SKIPPED_EXISTING。→ resume 逐模块工作流可用。
- F2【B7】modules/20_cnv.R：Windows 下 copykat 强制 n.cores=1。
- F3【B1】copykat 本地补丁：data/sysdata.rda → R/sysdata.rda 后重装（上游打包 bug，需在 SKILL 文档/安装脚本层面固化）。
- F4【B8】13_final_qa.R：历史 FAIL（早于本次运行 start_time）不再阻塞 FINAL_STATUS。
- F5【B3】06_annotate.R：map_reference_l1 增加人类规则（Epithelial/Fibroblast/Pericyte-smooth muscle/Monocyte/Macrophage/Megakaryocyte-platelet），20 个参考标签映射全部验证通过且无小鼠规则回归。
- F6【B9】advanced_figure_engine.R：plot_network_edges 自分泌边端点外移成弧，不再崩溃。
- 其他：install_dependencies.R 运行时可加载性审计；templates/modules/README.md 平台注意事项；SKILL.md v5.1.3 changelog；PACKAGE_MANIFEST/CONTENTS 更新；validation_profiles/GSE132465/ 完整档案（7 文件）。

## 正面验证（按阶段）

### 阶段 1 基础肿瘤注释 ✓
- 20,074 细胞导入，metadata 映射全部正确；count layer integer_like=TRUE。
- QC 保留 17,391/20,074（86.6%），12 样本全部 standard。
- RPCA 整合成功；自动选择 resolution 0.2 → 14 簇。
- 字典注释：T cell×4、B、Plasma、Monocyte、Mast、pDC-like、Fibroblast×2、Endothelial、Pericyte/smooth muscle、Epithelial；无 Unresolved。
- SingleR 离线优雅降级（WARNING + marker-only Medium），未中断流程。
- 与作者标签宏观一致（T 7,713 vs 7,874；上皮 2,605 vs 3,360；髓系 2,133 vs 2,365）。

### 阶段 2 配对 pseudobulk ✓
- 设计式 ~patient_id + condition 正确启用（metadata.paired: true）。
- 5 类 EVALUABLE（Epithelial/T/Fibroblast/Endothelial/B），5 类按 gate 规则 NOT_EVALUABLE（理由均正确：完整配对<3、单条件样本不足）。
- 上皮 Tumor vs Normal：上调 KRT23/KRT6A/KRT6B/BMP4/SAA1/REG1A；下调 PYY/GCG/CHGA/CA1/BEST4（分化丢失）——生物学合理。
- 显著基因数：T 1,645 / Fibroblast 1,321 / Endothelial 990 / Epithelial 872 / B 82。

### 阶段 3 CNA（copykat）✓（2026-08-06 完成）
- 12 样本全部通过 gate，17,391 细胞评分，COMPLETED。
- 肿瘤上皮 80.4% aneuploid-like（n=2,117）；同肿瘤内其他谱系 ≤8.8%（T 0.9%、浆 0.4%、内皮 0.6%）——恶性信号清晰。
- 正常黏膜上皮 37.9% aneuploid-like、部分正常谱系 20–40%：表达推断 CNA 在参照与候选同为二倍体池时的已知噪声，符合 SKILL "expression-inferred/需正交验证" 的解释边界，不是模块缺陷。
- 全细胞口径 Tumor 0.168 vs Normal 0.131；档案验证器 cnv 阶段 PASS（13 项）。
- 运行时长观察：Windows 单核 copykat 每样本约 30–90 分钟（12 样本一夜跑完）；文档已注明预期时长与"中断不得作为证据"。

### 阶段 4 CellChat ✓（2026-08-06 完成）
- 12 样本全部 COMPLETED（2–9 个支持组/样本；SMC05-T 仅 2 组但过 gate）。
- 6,247 条样本级互作（106 通路）；网络边/通路边/样本汇总表齐全。
- 条件对比按配对设计执行：paired=TRUE、n_pairs=6（患者均值配对 Wilcoxon）；2,335 对比行，FDR<0.05 为 0——n=6 对 + 2,335 重检验下符合统计预期，非缺陷；最小 p=0.036（T→T CLEC 下降、Epithelial→Monocyte APP 上升等方向合理）。
- 发现并修复 B9（自分泌边网络图崩溃）后 FINAL_STATUS=PASS_WITH_WARNINGS。

### 阶段 5 恶性细胞 NMF programs ✓（2026-08-06 完成）
- GeneNMF 后端成功（multiNMF + getMetaPrograms），识别 6 个 meta-programs。
- gate 行为正确：QC 后 SMC05-T 仅 2 个、SMC06-T 仅 25 个 TumorEpithelial 细胞（<100），被排除；SMC01–04-T（1,044/133/796/116）进入发现集——gate 使用的是注释对象内实际细胞数，非输入期计数。
- 程序生物学合理：MP1 增殖/细胞周期（CENPW/KIAA0101/BIRC5/UBE2C/MAD2L1）；MP3 热休克（HSPA6/HSPA1A/B/HSPH1/CRYAB）；MP4 分化（TFF1/TFF2/FABP1/KLK10/CEACAM1）；MP2 应激/炎症信号（LIF/EDN1/JUN/KLF6）；MP5 侵袭/炎症（MMP7/SAA1/COL17A1）；MP6 特异性低（recurrence 0.5）。
- recurrence 指标 0.5–1.0；细胞级评分回写全对象 + 按样本分组汇总齐全；契约产物（genes/metrics/similarity/weights/scored object）全部生成。

### 核心 QA
- FINAL_STATUS=PASS_WITH_WARNINGS；20+ 独立图（无 composite）；fgsea Hallmark 富集、组成分析、paired composition 图均导出。
