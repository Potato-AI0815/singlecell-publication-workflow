# GitHub upload guide

This package is prepared as **v0.1.0-alpha**. Before publishing, complete the release checklist.

## 1. Confirm public identity metadata

Edit at minimum:

- `CITATION.cff`: public author `Potato-AI`, no ORCID, and the repository URL;
- `COPYRIGHT`: legal copyright holder and year;
- `README.md` and `README.zh-CN.md`: repository URL and contact route;
- `SECURITY.md`: `939191004@qq.com` as the private security contact.

The public maintainer, author and brand identity for this release is **Potato-AI**.

## 2. Review licensing

The repository defaults to `GPL-3.0-or-later` for original project code. Third-party packages are external dependencies and remain under their upstream licenses. In particular, CopyKAT and scTenifoldKnk are not redistributed or installed by the default installer; review their terms before enabling them.

Read:

- `THIRD_PARTY_NOTICES.md`;
- `DEPENDENCY_LICENSES.csv`;
- `COMMERCIAL_USE.md`;
- `DISCLAIMER.md`.

## 3. Run local release checks

```powershell
python tests\static_validate.py .
python tests\validate_config.py templates\config.example.yml templates\config.schema.json
python tests\validate_config.py templates\config.full.yml templates\config.schema.json

& "D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe" install_dependencies.R --check
& "D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe" tests\00_generate_smoke_input.R C:\path\to\smoke_project
& "D:\Ruanjian\R-4.5.3\bin\x64\Rscript.exe" run_all.R C:\path\to\smoke_project\config.yml --mode=fast
```

Then inspect all proof PNG files and confirm `FINAL_STATUS.txt`.

## 4. Initialize and push

```bash
git init
git add .
git commit -m "Release v0.1.0-alpha"
git branch -M main
git remote add origin https://github.com/Potato-AI0815/singlecell-publication-workflow.git
git push -u origin main
```

## 5. Create the alpha release

Tag only after the uploaded repository passes GitHub Actions:

```bash
git tag -a v0.1.0-alpha -m "First public alpha release"
git push origin v0.1.0-alpha
```

Attach the packaged ZIP and SHA-256 file to the GitHub Release. Copy the text from `RELEASE_NOTES_v0.1.0-alpha.md`.

## 6. Public claims

Use the validation matrix exactly. Do not describe an implemented module as runtime validated unless `VALIDATION_MATRIX.md` marks it so.
