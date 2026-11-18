"""Verify the normalized affine-distance and finite-aperture response margins."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

# Publicly archived result in the 18 November 2026 joint-curvature edition.
SHA = "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42"
DEFAULT = Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118" / "METRIC.json"

def compute(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != SHA:
        raise ValueError("registered joint curvature report")
    target = json.loads(data)
    x0, x1 = F(49,32), F(99,64)
    a,b = map(F, target["windows"][1])
    if not a == x0 < x1 < b:
        raise ValueError("clock and strictly interior endpoint")
    low,high = map(F, target["uniform_null_abs_intervals"][1])
    # 0 <= t <= 1 and integral 2*t^2 dt = 2/3.
    plus = [low - 2*high/3, high - 2*low/3]
    if not F("0.051613") < plus[0] <= plus[1] < F("0.107100"):
        raise ValueError("ray margin")
    c = F(4,5)
    # Uniform cap has v1 uniform on [c,1], <v2^2>=(1-<v1^2>)/2.
    v1_squared = (1+c+c*c)/3
    angular = v1_squared - (1-v1_squared)/2
    if angular != F(18,25) or (1-c)/2 != F(1,10):
        raise ValueError("angular normalization")
    cap = [angular*v for v in plus]
    if not F("0.03716") < cap[0] <= cap[1] < F("0.07712"):
        raise ValueError("finite cap margin")
    return dict(input_sha256=SHA,clock=[x0,x1],remaining_clock_margin=b-x1,
                original_uniform_geometry=[low,high],normalized_axis_response=plus,
                cap_solid_angle_fraction=(1-c)/2,cap_quadrupole_mean=angular,
                normalized_cap_response=cap)

if __name__=="__main__":
    if len(sys.argv)>2: raise SystemExit("usage: optics.py [path to METRIC.json]")
    print(json.dumps(compute(Path(sys.argv[1]) if len(sys.argv)==2 else DEFAULT),
                     indent=2,default=str))
