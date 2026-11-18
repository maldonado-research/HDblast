"""Re-run each exact redshift calculation and check its archived output."""
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def replay():
    flags = ['-O'] if sys.flags.optimize else []
    for script, report in json.loads((ROOT / 'replays.json').read_text()):
        raw = subprocess.check_output([sys.executable, '-I', '-B', *flags, script], cwd=ROOT)
        if raw != (ROOT / report).read_bytes():
            raise ValueError(f'{script} disagrees with {report}')
        print(f'{script} reproduced {report}')


if __name__ == '__main__':
    replay()
