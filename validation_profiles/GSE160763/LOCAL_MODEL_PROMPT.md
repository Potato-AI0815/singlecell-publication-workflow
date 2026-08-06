# 给本地模型的固定提示词

你正在验证 `singlecell-publication-workflow-v0.1.0-alpha`，项目为小鼠大脑皮层 GSE160763。

严格执行：

1. 首轮只运行 `config_01_core_smoke.yml`，所有高级模块保持关闭。
2. 使用 `sample_id`/pooled cortex library 作为生物学重复；每组只有 2 个测序 pool，不得把每个 pool 中的 3 只小鼠拆成 6 个独立 scRNA-seq 重复。
3. 不得为了获得显著结果把 `minimum_replicates_per_condition` 从 3 降为 2。
4. `GSE160763_RAW.tar` 是只读输入；不得覆盖。
5. 先运行 `00_prepare_GSE160763.R`，确认 8 个样本、24 个输入文件和 40,666 个原始条码。
6. 使用小鼠皮层专用 L1/L2 marker dictionary；初步注释仍需人工核查。
7. 不生成综合大图，仅生成独立单图及各自 Source Data。
8. 遇到缺包、函数版本冲突、空图、元数据不匹配或 QA FAIL 时停止，不得跳过。
9. 最终报告必须列出：运行命令、result root、PASS/WARNING/FAIL、首个失败阶段、样本保留率、主要细胞类型、未解析比例、图形目录及需要修改的代码。

Windows 首轮命令：

```powershell
Set-ExecutionPolicy -Scope Process Bypass
& "<SKILL_ROOT>\validation_profiles\GSE160763\02_run_validation.ps1" -InstallPackages -Mode smoke
```
