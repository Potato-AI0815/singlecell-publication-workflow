#!/usr/bin/env bash
set -euo pipefail
ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
python "$ROOT/tests/static_validate.py" "$ROOT"
python "$ROOT/tests/validate_config.py" "$ROOT/templates/config.full.yml"
python "$ROOT/tests/validate_config.py" "$ROOT/templates/config.example.yml"
python - "$ROOT" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
for rel in ("tests/static_validate.py", "tests/validate_config.py"):
    p = root / rel
    compile(p.read_text(encoding="utf-8"), str(p), "exec")
print("PYTHON_SYNTAX PASS")
PY
echo "STATIC_SUITE PASS"
