#!/usr/bin/env python3
from pathlib import Path
import sys,json,yaml
path=Path(sys.argv[1])
cfg=yaml.safe_load(path.read_text())
schema=json.loads((Path(__file__).resolve().parents[1]/"templates/config.schema.json").read_text())
try:
 import jsonschema
 jsonschema.validate(cfg,schema)
except ImportError:
 required=schema["required"]; missing=[x for x in required if x not in cfg]
 if missing: raise SystemExit("Missing top-level fields: "+", ".join(missing))
print("CONFIG_VALIDATION PASS")
