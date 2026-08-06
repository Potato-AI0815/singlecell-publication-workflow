# Migration notes to v0.1.0-alpha

This release changes the public version scheme from internal workflow revisions (`v5.x`) to semantic product versioning.

## Version mapping

- internal `v5.1.2` core plus GSE132465 runtime fixes → public `v0.1.0-alpha`.

## Important changes

- use `run_one_click.cmd` on Windows;
- `run_all.R` accepts `--mode=<mode>`;
- resume mode reuses complete existing figures as `SKIPPED_EXISTING`;
- plot RDS sidecars may be lightweight renderable proxies;
- CopyKAT is single-core on Windows;
- historical issue-ledger failures do not block a repaired resume run;
- human SingleR label mapping is aligned to the human L1 dictionary;
- CellChat network self-loops are rendered safely;
- optional backends with unclear redistribution rights are no longer auto-installed.

Existing result directories remain readable. Start a new timestamped result directory when changing scientific parameters or input data.
