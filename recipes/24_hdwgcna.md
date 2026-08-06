# Recipe: hdWGCNA

1. Declare a coherent cell population in `grouping_column` and `group_values`.
2. Require at least three biological samples and enough cells per sample.
3. Build metacells using both sample identity and the selected population.
4. Freeze the soft-power grid and network settings.
5. Run modules, eigengenes, connectivity, hub genes and sample summaries.
6. Run condition or trait analyses only at sample level and only when replicated.
7. Export every diagnostic and module display as an independent figure.
8. Describe modules as exploratory until preservation/external validation.
