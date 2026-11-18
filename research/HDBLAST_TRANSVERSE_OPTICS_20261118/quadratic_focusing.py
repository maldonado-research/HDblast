"""Exact arithmetic consequences of the written 5D optical-focusing proof.

This replay authenticates retained inputs and checks finite inequalities.
It does not prove the continuum Gauss/Codazzi or optical-ODE arguments.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import hashlib
import json

ROOT = Path(__file__).resolve().parent
PINS = {
    "METRIC.json": "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42",
    "SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
}

def require(condition, label):
    if not condition:
        raise ValueError(label)


def load(name):
    parent = ROOT.parent / "HDBLAST_JOINT_CURVATURE_20261118"
    path = parent / name if name == "METRIC.json" else parent / "inputs" / name
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != PINS[name]:
        raise ValueError("retained source hash: " + name)
    return json.loads(data)


def compute():
    metric, shell = load("METRIC.json"), load("SHELL_BOX.json")
    ell, x0, x1 = F(1, 10**6), F(49, 32), F(99, 64)
    wx0, wx1 = map(F, metric["windows"][1])
    require(wx0 == x0 < x1 < wx1, "quadratic_ledger_1")
    hl, hu = map(F, metric["uniform_null_abs_intervals"][1])
    c_exact = [hl - 2 * hu / 3, hu - 2 * hl / 3]
    cm, cp = F("0.051613"), F("0.107100")
    require(0 < cm < c_exact[0] < c_exact[1] < cp, "quadratic_ledger_2")
    s0, s1 = F(1, 8), F(7, 8)
    d = s1*s1*cm - s0*s0*cp
    require(d == F(2421937, 64000000) > F("0.03784"), "quadratic_ledger_3")
    # Integral ds/[s(1-s)] between s0,s1 = log 49.
    ratio = (s1 / (1-s1)) / (s0 / (1-s0))
    require(ratio == 49, "quadratic_ledger_4")
    exp4_partial = sum(F(4)**k / factorial(k) for k in range(8))
    require(exp4_partial > 49, "quadratic_ledger_5")  # Thus exp(4)>49 and log(49)<4.
    blo = F(shell["boxes"]["B_b"]["exact"][0])
    rhi = F(shell["boxes"]["rho_b"]["exact"][1])
    require(blo > 0 and rhi > 0, "quadratic_ledger_6")
    # Omega=ell*B/rho*(L/rho)^2, and L/rho>=ell*(x1-x0).
    omega_lower = ell**3 * (x1-x0)**2 * blo / rhi
    axis = omega_lower**2 * d**2 / 2
    cap_cut = F(4, 5)
    qmin = 2*cap_cut**2 - 1
    require(qmin == F(7, 25), "quadratic_ledger_7")
    mean_x2 = sum(cap_cut**k for k in range(3)) / 3
    mean_x4 = sum(cap_cut**k for k in range(5)) / 5
    mean_q2 = F(19,8)*mean_x4 - F(7,4)*mean_x2 + F(3,8)
    require(mean_q2 == F(10283, 18750), "quadratic_ledger_8")
    point_cap, mean_cap = qmin**2*axis, mean_q2*axis
    require(axis > F("2.9024e-50"), "quadratic_ledger_9")
    require(point_cap > F("2.2755e-51"), "quadratic_ledger_10")
    require(mean_cap > F("1.5917e-50"), "quadratic_ledger_11")
    # The normalized determinant deficit coefficient is three times chi.
    return {
        "status": "PASS_EXACT_RATIONAL_CONSEQUENCES",
        "scope": "conditional one-sided 5D screen Jacobi transfer along shell-tangent rays; tangent S2 caps only",
        "input_pins": PINS,
        "clock": [x0, x1],
        "prefix_affine_fractions": [s0, s1],
        "exact_uniform_prefix_C": c_exact,
        "rounded_prefix_C": [cm, cp],
        "normalized_strain_increment_lower": d,
        "reciprocal_weight_log_argument": ratio,
        "exp4_positive_partial_sum": exp4_partial,
        "Omega_lower": omega_lower,
        "cap_Q_abs_lower": qmin,
        "cap_Q_squared_mean": mean_q2,
        "radius_quadratic_deficit_lower": {
            "axis": axis, "each_cap_point": point_cap, "tangent_cap_average": mean_cap,
        },
        "volume_quadratic_deficit_lower": {
            "axis": 3*axis, "each_cap_point": 3*point_cap, "tangent_cap_average": 3*mean_cap,
        },
        "display_only": {
            "d": float(d), "Omega_lower": float(omega_lower),
            "axis_radius_coefficient_lower": float(axis),
            "cap_point_radius_coefficient_lower": float(point_cap),
            "cap_mean_radius_coefficient_lower": float(mean_cap),
        },
        "not_claimed": [
            "finite-amplitude numerical Taylor remainder",
            "positive-measure full bulk S3 sky aperture",
            "bulk crossing or shell reflection law",
            "experimental source-amplitude or affine calibration",
        ],
    }

if __name__ == "__main__":
    print(json.dumps(compute(), indent=2, default=str))
