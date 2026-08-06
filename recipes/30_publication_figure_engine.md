# Recipe: V5 independent publication figures

1. Require an accepted annotated object and/or generated Source Data.
2. Set:

```yaml
figures:
  export_individual_figures: true
  composite_figures_forbidden: true
  build_contact_sheet: false
```

3. Use only functions from `R/lib/figure_engine.R` and `R/lib/advanced_figure_engine.R`.
4. Export each scientific display with a unique `figure_id` and directory.
5. Reject patchwork, ggarrange, montage and collage objects.
6. Preserve raw and display-transformed values in Source Data.
7. Inspect every proof PNG directly.
8. Do not recompute annotations or statistics in `figures_only` mode.
