# V5 独立单图引擎目录

## 强制规则

- 不生成综合大图、overview、montage 或 contact sheet。
- 不接受 patchwork/ggarrange 图对象。
- 每张图有独立目录、四种正式格式、proof PNG、plot RDS、参数、Source Data 和 QA。
- UMAP/tSNE 的轨迹箭头只能来自拟合模型。
- 颜色截断必须记录原始值与显示值。

## 核心单图

UMAP/tSNE 分类图、连续值图、split violin、DotPlot、marker heatmap、组成图、配对组成图、火山图、富集图、轨迹图和拟时序趋势图。

## 高级单图

CNV 热图/散点图、CellChat 气泡/网络/差异图、tradeSeq 诊断和基因排序、NMF rank/相似性/基因/活动/嵌入图、hdWGCNA power/连接度/hub/性状/树状/ME 图、空间分类/连续/去卷积/邻域图、药物 reversal/score/cluster/oncoPredict 图。


## 虚拟敲除单图模板

- `plot_virtual_knockout_ranking()`：无方向的共识网络影响排序。
- `plot_virtual_knockout_sample_support()`：跨样本与重复运行支持热图。
- `plot_virtual_knockout_target_expression()`：目标表达和样本可评价性。
- `plot_virtual_knockout_manifold()`：WT 与 KO 网络状态位移，箭头不代表表达方向。
- `plot_virtual_knockout_network()`：重复出现的目标中心 WT 网络背景。
- `plot_virtual_knockout_pathways()`：受影响基因的通路过度代表。

所有模板只生成单张图，禁止拼版。
