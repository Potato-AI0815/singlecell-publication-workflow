# 投稿级单细胞可视化速查

## 必备图

1. QC：每样本细胞数、nFeature、nCount、线粒体比例，保留前后并列。
2. 聚类：UMAP by cluster、by sample、by condition；整合与未整合并列。
3. 注释：UMAP、canonical-marker dot plot、竞争谱系证据、置信度表。
4. Marker：每 cluster 的 top markers；热图优先使用 cluster/sample 平均，避免全部细胞热图。
5. 组成：样本级堆叠图和样本级散点/箱线图；不把细胞作为重复。
6. DE：pseudobulk volcano/MA；图中注明 target vs reference、样本数和模型。
7. 富集：来自 pseudobulk 排名；显示 NES、FDR 和 gene-set 来源。

## 排版

- 主图常用 183 mm 或 A4 横向工作画布；终稿文字至少 8 pt。
- 分类变量固定调色板；缺失值统一灰色；连续变量使用感知均匀或蓝–白–珊瑚发散色。
- 同一细胞类型在所有图中保持同色。
- 每幅图导出 PDF、SVG、600 dpi TIFF、PNG，并保存 plot RDS。
- 每个 panel 的数值必须能追溯到 Source Data。

## 禁止

- 彩虹色；
- 用 UMAP 距离作定量结论；
- 用星号替代模型和分母；
- 把缺失值画成零；
- 把 integrated/SCT residual 作为正式 pseudobulk 输入；
- 只展示支持预期结论的患者、cluster、方向或参数。

## V5 单图引擎补充

- cluster/细胞类型图可使用拆分不连通区域后的柔和轮廓；condition/sample 图默认不使用轮廓。
- 超大细胞数时只将细胞点栅格化，文字、轮廓、轨迹、图例和坐标箭头保持矢量。
- 连续变量图必须同时保存原始值与显示截断值，并在参数 YAML 中记录分位数和最终色标范围。
- split violin 只允许两个条件；更多条件改用普通 violin、箱线图或分面图。
- 轨迹箭头只能来自拟合模型的曲线坐标，禁止手工添加方向。
- 每个单图独立导出并带 visual QA；本工作流不生成综合大图、拼图或 contact sheet。
