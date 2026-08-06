# GSE160763 validation-ready SKILL V5.1

This package is a runtime-hardened copy of V5 plus a project profile for GSE160763.

Key fixes before the first public-data run:

- fixed the preflight capability-length error after virtual knockout was added;
- added `modules/27_virtual_knockout.R` to the engine audit;
- made multi-10x import ignore incomplete/root directories;
- added one-row-per-sample metadata joining;
- made custom L1/L2 marker dictionaries operational for mouse and human data;
- added mouse cortex marker dictionaries and a two-stage validation protocol.

Start with `validation_profiles/GSE160763/VALIDATION_PROTOCOL.md`.
