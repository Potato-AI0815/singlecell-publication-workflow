# GSE132465 测试运行与 BUG 修改日志

Skill：single-cell-publication-workflow（v5.1.2 → v5.1.3）
Skill 根目录：`C:\Users\YHN\Desktop\Qoder\single-cell-publication-skill-v5.1-GSE160763`
项目目录：`C:\Users\YHN\Desktop\Qoder\GSE132465_analysis`
结果目录：`results\GSE132465_CRC_6pair_20260805_120916`（FINAL_STATUS = PASS_WITH_WARNINGS）
运行环境：Windows，R 4.5.3，32GB RAM，i7-11800H，Rtools45
日志时间范围：2026-08-05 11:00 至 2026-08-06 14:23

---

## 一、运行日志（按轮次）

### 准备阶段（08-05 上午）
- 数据：GSE132465 = 23 例韩国 CRC 患者（23 肿瘤 + 10 配对正常黏膜，63,689 细胞，10x 3'）。两个 GEO 补充文件位于 `D:\qoderwork\`（本地已有，未重复下载）。
- 测试子集（用户选定）：SMC01–SMC06 六对配对，12 样本 / 20,074 细胞。
- 准备脚本：`input\prepare_GSE132465_6pair.R`（列选择式读取 123MB 密集矩阵，构建 raw-count Seurat 对象 + 样本清单 + `tumor_epithelial_flag` 恶性细胞预声明列）。QC 后 17,391 细胞。
- 依赖安装：核心包已备；补装 WGCNA（需 Bioc impute/preprocessCore 才能加载）、NMF、clustree、bluster；codeload tarball 本地安装 copykat、GeneNMF、CellChat；git clone 安装 hdWGCNA（dev 分支，需 GeneOverlap/tester）。

### Round 1：核心 full 运行（08-05 12:09–12:32，约 23 分钟）✓
阶段 01–13 全通。关键结果：
- 导入/QC：20,074 → 17,391 细胞（86.6%），12 样本全部 standard，integer_like=TRUE。
- RPCA 整合成功；自动选 resolution 0.2 → 14 簇。
- 注释：T/B/Plasma/Monocyte/Mast/pDC-like/Fibroblast/Endothelial/Pericyte/Epithelial，无 Unresolved；SingleR 离线（celldex 下载失败）优雅降级为 marker-only（Medium）。
- 配对 pseudobulk：设计式 `~patient_id + condition` 生效；5 类 EVALUABLE、5 类按 gate 正确 NOT_EVALUABLE；上皮 DE 生物学合理（KRT23/KRT6A/BMP4 上调，PYY/GCG/CHGA/CA1/BEST4 下调）。
- 20+ 独立图、fgsea 富集、组成分析齐全；FINAL_STATUS=PASS_WITH_WARNINGS。

### Round 2：CNA 第 1 次尝试（08-05 12:59–13:12）✗
copykat 每样本报 `object 'full.anno' not found` → 全部样本跳过 → NOT_EVALUABLE。
→ 定位为 copykat 上游打包 bug（B1），打补丁重装。

### Round 3：CNA 第 2 次尝试（08-05 13:09 起）✗
`workers: 4` 触发 `'mc.cores' > 1 is not supported on Windows` → NOT_EVALUABLE。
→ 模块内强制 Windows 单核（B7 修复）+ config workers:1。
同日还发现并修复 resume 图引擎阻塞 bug（B2）：11_figures 在图已存在时报 "did not export any figure"。

### Round 4：CNA 第 3 次尝试（08-05 21:12 重启）— 被关机中断 ✗
修复后重跑，夜间串行处理 12 样本；用户次日关机时仍在计算，进度丢失（copykat 无逐样本断点，全样本完成才写盘）。

### Round 5：CNA 第 4 次尝试（08-05 21:12 → 08-06 上午完成）✓
- 12 样本全部通过 gate，17,391 细胞评分，COMPLETED。
- 肿瘤上皮 80.4% aneuploid-like，同肿瘤微环境其他谱系 ≤8.8%；正常黏膜上皮 37.9%（表达推断 CNA 已知噪声，符合解释边界）。
- final QA 初次 FAIL：台账中 B2 修复前的两条历史 FAIL 被计入 → 修复 B8（历史 FAIL 不再阻塞）→ 重跑后 PASS_WITH_WARNINGS。
- 档案验证器 cnv 阶段 PASS（13 项）。

### Round 6：CellChat（08-06 约 11:4x–13:07，约 1.5 小时）✓
- 12 样本全部 COMPLETED（2–9 个支持组/样本）；6,247 条样本级互作、106 通路。
- 配对条件对比正确执行（paired=TRUE，n_pairs=6）；2,335 对比行，FDR<0.05 为 0（n=6 对 + 重检验下的统计预期）。
- 发现 B9：Communication_network_Normal/Tumor 两图因自分泌边（source==target）"end points must not be identical" 导出失败 → 修复后重跑收尾 → PASS_WITH_WARNINGS。

### Round 7：NMF programs（08-06 13:2x–13:30，约 8 分钟）✓
- GeneNMF 后端成功；gate 正确排除 SMC05-T（2 细胞）与 SMC06-T（25 细胞），SMC01–04-T 进入发现集。
- 6 个 meta-programs，生物学合理：MP1 增殖/细胞周期、MP3 热休克、MP4 分化、MP2 应激信号、MP5 侵袭/炎症、MP6 低特异性。
- 契约产物齐全（genes/metrics/similarity/weights/scored object）；PASS_WITH_WARNINGS。

### Round 8：hdWGCNA（08-06 13:30–14:23，49 分钟时应用户要求停止）◐ 部分完成
- 已完成：metacell 构建、soft-power 表（无 power 达到 SFT 0.8，truncated.R² 在 power 5 达 0.95）、TOM 计算（约 7,400 基因，433MB，13:45 写盘）。
- 停止时处于模块检测/聚类阶段；模块输出未写盘（无 ADV_hdWGCNA.done）。结果目录状态完好，重跑时该模块自动从头执行。

### 未执行
- 阶段 7 药物反应假设（signature_reversal；签名资源已构建：MSigDB C2:CGP 1,219 个方向性签名，269,867 行，含 5-FU/伊立替康/顺铂相关集）。

---

## 二、BUG 与修复日志

| 编号 | 严重度 | 问题 | 根因 | 修复 | 修改文件 |
|---|---|---|---|---|---|
| B1 | 严重 | copykat 每样本报 `object 'full.anno' not found` | 上游 navinlabcode/copykat v1.2.5 将 sysdata.rda 放在 data/ 而非 R/，内部对象不进命名空间 | 源码内移动文件后重装；install_dependencies.R 增加命名空间自检与指引 | 本地补丁 + install_dependencies.R |
| B2 | 中 | resume 轮次 stage 11 报 "did not export any figure" | 图引擎对已存在非空图目录 stop()，全部图被拒导出 | 已存在图返回 SKIPPED_EXISTING；stage 11/13 接受跳过状态；run 内重名仍报错 | R/lib/figure_engine.R、R/11_figures.R、R/13_final_qa.R |
| B3 | 中 | SingleR 人类参考映射缺规则，会把 Fibroblast/Pericyte 误判 conflict→Ambiguous | map_reference_l1 为小鼠 CNS 定制 | 追加 6 条人类规则并与 human_l1 字典对齐；20 标签映射测试全过 | R/06_annotate.R |
| B4 | 轻 | 药物签名资源仅占位模板 | 包内只有 ExampleDrug 模板 | 构建 msigdbr C2:CGP 方向性签名脚本并入档案 | validation_profiles/GSE132465/build_drug_signature_msigdbr.R |
| B5 | 轻 | WGCNA/scDblFinder 装上却无法加载 | 缺 Bioc 依赖（impute/preprocessCore/bluster）无提示 | install_dependencies.R 增加运行时可加载性审计 | install_dependencies.R |
| B6 | 轻 | GEO 大数据下载无指引 | 文档缺失 | 档案协议写入"本地优先/分块续传/ftp 备选" | validation_profiles/GSE132465/VALIDATION_PROTOCOL.md |
| B7 | 严重 | copykat workers>1 在 Windows 硬报错 | parallel::mclapply 不支持 Windows 多核 | 模块内按平台强制 n.cores=1 | modules/20_cnv.R |
| B8 | 中 | 历史 FAIL 台账条目永久阻塞已修复运行 | final QA 计入全部历史 FAIL | 以 run_manifest start_time 为界区分历史/当前 FAIL | R/13_final_qa.R |
| B9 | 中 | CellChat 网络图自分泌边崩溃 | geom_curve 端点重合 | 自分泌边端点外移渲染为小弧，记录 n_self_loops | R/lib/advanced_figure_engine.R |

其他改进（非 BUG）：SKILL.md 升级 v5.1.3 并记录全部变更；VERSION→5.1.3；CONTENTS.md 与 PACKAGE_MANIFEST.csv 更新；templates/modules/README.md 增加平台注意事项；新建 validation_profiles/GSE132465/ 完整档案（准备脚本、核心配置、模块 gate 片段、签名构建脚本、预期计数、验证脚本、协议文档）；档案验证器 core/cnv 阶段均 PASS。

---

## 三、当前状态与恢复指南

**已完成验证**：阶段 1 肿瘤注释 ✓、阶段 2 配对 pseudobulk ✓、阶段 3 CNA ✓、阶段 4 CellChat ✓、阶段 5 NMF ✓。
**部分完成**：阶段 6 hdWGCNA（TOM 已算，模块未出，被停止）。
**未开始**：阶段 7 药物反应。

恢复步骤（在 Skill 根目录执行）：

1. 重跑 hdWGCNA（config_core.yml 中 hdWGCNA 已 enabled）：
   删除 `13_logs/11.done 12.done 13.done`（若存在），然后
   `"D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe" run_all.R "C:/Users/YHN/Desktop/Qoder/GSE132465_analysis/config_core.yml" --mode=resume`
2. hdWGCNA 通过后，将 config_core.yml 的 `drug_response.enabled` 改为 true（signature_path 已指向构建好的签名），同样删除 11/12/13.done 后 resume。
3. 每轮结束后运行档案验证器：
   `Rscript validation_profiles\GSE132465\03_validate_GSE132465_outputs.R <结果目录> hdwgcna|drug|all`

注意：hdWGCNA 的 soft-power 表显示无 power 达到 0.8 scale-free 拟合（肿瘤转录组常见）；模块逻辑如何兜底选择 power 将在重跑时确认，必要时作为新发现记录。
