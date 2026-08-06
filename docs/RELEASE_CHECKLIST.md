# GitHub release checklist

Before publishing `v0.1.0-alpha`:

- [ ] replace placeholder repository URL in `CITATION.cff`;
- [ ] confirm copyright holder/maintainer identity;
- [ ] review GPL-3.0-or-later selection;
- [ ] run installed-package license audit;
- [ ] confirm no third-party source code or restricted data are bundled;
- [ ] run static suite;
- [ ] run config schema validation;
- [ ] run a clean GSE160763 regression on the release archive;
- [ ] verify Windows launcher exit codes;
- [ ] inspect proof PNGs manually;
- [ ] confirm `VALIDATION_MATRIX.md` does not overstate incomplete modules;
- [ ] create Git tag `v0.1.0-alpha`;
- [ ] attach ZIP and SHA-256 to the GitHub release;
- [ ] use `RELEASE_NOTES_v0.1.0-alpha.md` as the release description.
