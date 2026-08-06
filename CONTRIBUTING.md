# Contributing

Maintainer: Potato-AI  
Security reports: 939191004@qq.com

## Bug reports

Please include:

1. software version and commit;
2. operating system and R version;
3. sanitized resolved configuration;
4. `15_manifests/package_versions.csv` and `sessionInfo.txt`;
5. failing stage log from `13_logs/`;
6. `13_logs/issue_ledger.csv` or `14_qa/final_QA_report.md`;
7. minimal reproducible input when legally shareable;
8. whether the failure is reproducible from a new timestamped result directory.

Do not upload identifiable patient data or restricted datasets.

## Pull requests

- preserve sample/patient as the inferential unit;
- do not add cell-level formal group testing;
- add a scientific gate and output contract for new modules;
- return `NOT_EVALUABLE` rather than forcing unsupported analyses;
- export independent figures only;
- add Source Data, parameters and QA;
- update `VALIDATION_MATRIX.md` only after a real benchmark;
- record dependency license and citation requirements;
- run `python tests/static_validate.py .` and config validation.

## External contributions and future dual licensing

Contributors retain copyright unless a separate agreement is signed. The maintainer may introduce a contributor license agreement before accepting substantial external code if future dual licensing is planned.
