"""Replay the numerical reports using the published data in the capsule."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if __name__ == "__main__":
    flags = ["-O"] if sys.flags.optimize else []
    for executable, expected in json.loads((ROOT/"replays.json").read_text()):
        script = ROOT/executable
        actual = subprocess.check_output([sys.executable,"-I","-B",*flags,script.name],cwd=script.parent)
        if actual != (ROOT/expected).read_bytes():
            raise ValueError(f"{executable} differs from {expected}")
        print(f"{executable} reproduced {expected}")
