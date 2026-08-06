# Recipe: exploratory drug response

Select exactly one backend:

- `signature_reversal`: pseudobulk DE plus drug perturbation effects;
- `signature_score`: drug gene signatures scored in cells and aggregated by sample/cell type;
- `oncopredict`: declared training expression and response matrices plus pseudobulk expression.

Always audit resource provenance, direction and gene overlap. Export one ranking, contribution, heatmap, distribution or cluster plot per figure. State that outputs are hypothesis-generating and do not constitute treatment recommendations.
