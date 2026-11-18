"""Reproduce all archived interval and angular reports."""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

if __name__ == '__main__':
    flags = ['-O'] if sys.flags.optimize else []
    for script, report in json.loads((ROOT/'replays.json').read_text()):
        data = subprocess.check_output([sys.executable, '-I', '-B', *flags, script], cwd=ROOT)
        if data != (ROOT/report).read_bytes():
            raise ValueError(f'{script} differs from {report}')
        print(f'{script} reproduced {report}')
