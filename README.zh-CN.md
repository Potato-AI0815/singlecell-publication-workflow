# SingleCell Publication Workflow v0.1.0-alpha

这是一个面向医学与生命科学研究者的、配置驱动、样本感知、可审计的单细胞分析工作流。

项目仓库：[Potato-AI0815/singlecell-publication-workflow](https://github.com/Potato-AI0815/singlecell-publication-workflow)  
维护者：**Potato-AI**  
安全问题请通过邮箱联系：939191004@qq.com

它的价值不是重新发明 Seurat、CellChat、CopyKAT、GeneNMF 或 hdWGCNA，而是把这些方法组织成一套受控执行系统：

```text
输入适配
→ 实验设计核对
→ QC/标准化/整合/聚类/初步注释
→ 样本级统计门控
→ 满足条件才运行高级模块
→ 独立投稿图与 Source Data
→ Methods、日志、manifest 和 Final QA
```

## Alpha 已验证范围

- GSE160763：小鼠脑皮层 8 个样本、40,666 个原始细胞；核心 SCT/RPCA 流程通过，QC 后 38,215 个细胞，18 个独立图，n=2 差异分析正确返回 `NOT_EVALUABLE`。
- GSE132465 六对 CRC 子集：基础肿瘤注释、配对 pseudobulk、CopyKAT、CellChat、GeneNMF 均完成真实运行验证。
- hdWGCNA：完成 metacell、soft power 与 TOM，模块检测阶段被人工停止，因此仅算部分验证。
- 轨迹、空间、药物敏感性和虚拟敲除：已有执行适配器，但尚未完成专用公共数据验证。


> **受限可选后端：** 默认依赖安装不会自动安装或再分发 CopyKAT 与 scTenifoldKnk。启用前必须审查其上游许可证与使用条件；CopyKAT 的修复安装只能通过用户显式传入 `-RepairCopykat` 触发。

## 一键运行

```powershell
run_one_click.cmd `
  -Config "C:\path\to\config.yml" `
  -Mode fast `
  -RscriptPath "D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe"
```

运行成功必须同时满足：

- R 退出码为 0；
- `14_qa/FINAL_STATUS.txt` 为 `PASS` 或 `PASS_WITH_WARNINGS`。

## 核心功能

- 10x、H5、Seurat RDS、SCE 和可选 H5AD 输入；
- 多样本 metadata 合并；
- 按样本 QC、可选 scDblFinder；
- LogNormalize 或 SCTransform v2；
- Seurat v5 RPCA、聚类和 UMAP；
- marker program＋竞争谱系＋可选 SingleR 的初步注释；
- cluster marker、样本级组成；
- `sample × cell type` pseudobulk edgeR；
- fgsea 富集；
- CNA、CellChat、轨迹、NMF、hdWGCNA、空间、药敏和虚拟敲除适配器；
- 每张图独立输出 PDF/SVG/TIFF/PNG、参数、Source Data 和视觉 QA；
- 完整 Methods、图注、软件版本、日志和 Final QA。

## 最重要的差异

这个工具不追求“无论如何都出结果”。数据前提不足时，模块应返回：

```text
NOT_EVALUABLE
```

它不会把细胞当作独立样本，也不会为了得到 P 值而降低预设重复数门槛。

## 当前边界

- 自动注释必须人工复核；
- 高级模型结果默认属于探索性证据；
- CNA 不是 DNA 证据；
- CellChat 不是因果或真实物理通讯；
- 药物分数不能用于临床用药；
- 虚拟敲除不能代替真实实验；
- 本版本仍是实验性 Alpha，不是经过所有平台与组织验证的生产级软件。

详细范围见 [VALIDATION_MATRIX.md](VALIDATION_MATRIX.md) 和 [KNOWN_LIMITATIONS.md](KNOWN_LIMITATIONS.md)。
