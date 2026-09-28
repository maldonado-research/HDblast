#!/usr/bin/env python3
"""Re-executes the workstream's controls.py unchanged except that the output goes to
verify/controls_rerun.json (the original CONTROLS.json is never touched)."""
import sys, runpy
from pathlib import Path
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
src = (ROOT / 'controls.py').read_text()
marker = "(HERE / 'CONTROLS.json')"
assert src.count(marker) == 1
src = src.replace(marker, "(Path(%r) / 'controls_rerun.json')" % str(HERE))
sys.path.insert(0, str(ROOT))
g = {'__name__': '__main__', '__file__': str(ROOT / 'controls.py')}
exec(compile(src, str(ROOT / 'controls.py'), 'exec'), g)
