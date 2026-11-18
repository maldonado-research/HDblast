"""Replay each published arithmetic certificate using the current Python interpreter."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CERTIFICATES = (
    ("metric.py", "METRIC.json"),
    ("null_sky.py", "NULL_SKY.json"),
    ("angular_check.py", "ANGULAR.json"),
    ("curvature.py", "CURVATURE.json"),
    ("powers.py", "POWERS.json"),
    ("conformal_test.py", "CONFORMAL_TEST.json"),
    ("observability.py", "OBSERVABILITY.json"),
    ("second_order_freedom.py", "FREEDOM.json"),
)

if __name__ == "__main__":
    for script, report in CERTIFICATES:
        flags = ["-O"] if sys.flags.optimize else []
        data = subprocess.check_output([sys.executable, "-I", "-B", *flags, script], cwd=ROOT)
        if data != (ROOT / report).read_bytes():
            raise SystemExit(f"Output differs from {report}")
        print(f"Verified {report}")
