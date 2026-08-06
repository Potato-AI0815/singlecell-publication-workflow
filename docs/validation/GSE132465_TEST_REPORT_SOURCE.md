# single-cell-publication-workflow v5.1 → v5.1.3：GSE132465 测试运行与 SKILL 优化报告

日期：2026-08-05 至 2026-08-06
测试者：千问办公（自动执行）
Skill 根目录：C:\Users\YHN\Desktop\Qoder\single-cell-publication-skill-v5.1-GSE160763
项目目录：C:\Users\YHN\Desktop\Qoder\GSE132465_analysis
结果目录：results/GSE132465_CRC_6pair_20260805_120916

## 1. 测试目标与数据

按用户指定顺序验证：基础肿瘤注释 → 配对 pseudobulk → CNA → CellChat → 恶性细胞 NMF programs → 恶性细胞 hdWGCNA → 药物反应假设，随后进行 SKILL 优化。

数据集 GSE132465：23 例韩国结直肠癌患者（23 肿瘤 + 10 配对正常黏膜，63,689 细胞，10x 3'）。
测试规模（用户选定）：前 6 对配对患者 SMC01–SMC06，12 样本 / 20,074 细胞（QC 后 17,391）。
输入准备：validation_profiles/GSE132465/00_prepare_GSE132465.R（新增档案），从 GEO 原始 UMI 矩阵 + 作者注释表抽提；构建 raw-count Seurat 对象、样本清单、以及预声明的恶性细胞标记列 tumor_epithelial_flag（作者 Epithelial × Tumor）。

运行环境：Windows 10/11，R 4.5.3，32GB RAM，i7-11800H（8C/16T），Rtools45。

## 2. 各阶段验证结果

| 阶段 | 模块 | 状态 | 关键证据 |
|---|---|---|---|
| 1 | 基础肿瘤注释 | ✓ PASS | 14 簇；T/B/Plasma/Myeloid/Mast/Fibroblast/Endothelial/Pericyte/Epithelial；无 Unresolved；SingleR 离线优雅降级 |
| 2 | 配对 pseudobulk | ✓ PASS | ~patient_id+condition；5 类 EVALUABLE，5 类按 gate NOT_EVALUABLE；上皮 DE 生物学合理 |
| 3 | CNA（copykat） | ✓ PASS | 肿瘤上皮 80.4% aneuploid-like vs 微环境 ≤8.8%；档案验证器 13 项全过 |
| 4 | CellChat | （进行中/待填） | |
| 5 | NMF programs | （待填） | |
| 6 | hdWGCNA | （待填） | |
| 7 | 药物反应假设 | （待填） | |

## 3. 发现并修复的问题（B1–B8）

（详见 TEST_FINDINGS.md；此处摘要）

## 4. SKILL 优化清单（v5.1.3）

（修复文件、新增档案、文档更新）

## 5. 遗留事项与建议

（未测模块、资源建议）
