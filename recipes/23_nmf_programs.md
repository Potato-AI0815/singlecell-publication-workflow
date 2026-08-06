# Recipe: recurrent NMF programs

1. Define a biologically coherent `subset_column` and `subset_values`.
2. Require multiple samples and adequate cells per sample.
3. Predeclare ranks and seed.
4. Enable `advanced_modules.nmf`.
5. Prefer GeneNMF; use fallback only when the preferred backend fails or is absent.
6. Review recurrence metrics and program genes before assigning biological names.
7. Export each rank diagnostic, similarity heatmap, gene display, activity summary and program embedding separately.
8. Report single-sample components as non-recurrent.
