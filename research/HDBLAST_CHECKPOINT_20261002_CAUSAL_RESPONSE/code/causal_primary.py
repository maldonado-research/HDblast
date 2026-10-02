#!/usr/bin/env python3
"""Frozen primary routes for a prescribed fixed-geometry linear response.

Importing this module performs no source, kernel or mode evaluation. Execute only
after the public registration is frozen. All results use y=a^2 deltaQ/epsilon.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import time
import traceback
import warnings

import mpmath
import numpy
import scipy
from scipy.integrate import IntegrationWarning, quad

ROOT = Path(__file__).resolve().parents[1]
PREF = 1.0 / (8.0 * math.pi**2)
GAMMA_E = float(numpy.euler_gamma)
DEADLINE = None


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def check_deadline():
    if DEADLINE is not None:
        require(time.monotonic() <= DEADLINE, "Primary wall-time budget exceeded")


def verify_registration(expected):
    registration = ROOT / "FULL_REGISTRATION.json"
    require(sha(registration) == expected, "Registration SHA256 does not match explicit pin")
    document = json.loads(registration.read_text())
    for filename, digest in document["files"].items():
        path = (ROOT / filename).resolve()
        require(path.is_relative_to(ROOT), "Registration path escapes package")
        require(sha(path) == digest, "Frozen input changed: " + filename)
    require("code/causal_primary.py" in document["files"], "Producer absent from registration")
    require("code/validate_causal.py" in document["files"], "Validator absent from registration")
    require("EXPERIMENT.json" in document["files"], "Experiment absent from registration")
    return document


def bump(u):
    if abs(u) >= 1.0:
        return 0.0
    return math.exp(1.0 - 1.0 / (1.0 - u * u))


def source(source_id, eta):
    u = eta + 4.0
    value = bump(u)
    if source_id == "positive_B":
        return value
    if source_id == "signed_uB":
        return u * value
    raise ValueError("Unknown frozen source ID")


def source_prime(source_id, eta):
    u = eta + 4.0
    if abs(u) >= 1.0:
        return 0.0
    value = bump(u)
    derivative = -2.0 * u * value / (1.0 - u * u)**2
    return derivative if source_id == "positive_B" else value + u * derivative


def integrate(function, left, right, config, **kwargs):
    check_deadline()
    if right <= left:
        return 0.0, 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        value, error = quad(function, left, right, epsabs=config["epsabs"],
                            epsrel=config["epsrel"], limit=config["limit"], **kwargs)
    check_deadline()
    require(math.isfinite(value) and math.isfinite(error), "Nonfinite quadrature")
    return float(value), float(error)


def logarithmic_memory(source_id, eta, config, input_factor_mutation=False):
    """Integrate f' log(eta-t); logarithmic endpoint handled by QAWSE."""
    if eta <= -5.0:
        return {"y": 0.0, "quadrature_error_estimate": 0.0}
    right = min(eta, -3.0)
    mass = math.sqrt(2.0) / (-eta)
    if input_factor_mutation:
        # Omitting a(t)^2 when converting physical delta x gives t^2 f(t).
        derivative = lambda t: t*t*source_prime(source_id, t) + 2*t*source(source_id, t)
        f0 = eta*eta*source(source_id, eta)
    else:
        derivative = lambda t: source_prime(source_id, t)
        f0 = source(source_id, eta)
    if right == eta:
        value, error = integrate(derivative, -5.0, right, config,
                                 weight="alg-logb", wvar=(0.0, 0.0))
    else:
        value, error = integrate(lambda t: derivative(t)*math.log(eta-t), -5.0, right, config)
    return {"y": -PREF*(value + f0*(math.log(mass) + GAMMA_E + 1.0)),
            "quadrature_error_estimate": PREF*error}


def source_difference(source_id, eta, tau):
    """Stable f(eta-tau)-f(eta), including signed-source cancellation."""
    if tau == 0.0:
        return 0.0
    u = eta + 4.0
    v = u - tau
    if abs(u) >= 1.0 or abs(v) >= 1.0:
        return source(source_id, eta-tau) - source(source_id, eta)
    gu = 1.0 - 1.0/(1.0-u*u)
    # (g(v)-g(u)) is evaluated as a rational difference, without subtracting g's.
    difference = tau*(2.0*u-tau)/((1.0-u*u)*(1.0-v*v))
    base = math.exp(gu)
    if abs(difference) < 0.5:
        delta_b = base*math.expm1(difference)
        return delta_b if source_id == "positive_B" else u*delta_b-tau*(base+delta_b)
    return source(source_id, eta-tau) - source(source_id, eta)


def finite_part_memory(source_id, eta, eta_i, config):
    if eta <= -5.0:
        return {"y": 0.0, "quadrature_error_estimate": 0.0}
    length = eta-eta_i
    f0 = source(source_id, eta)
    def integrand(tau):
        return -source_prime(source_id, eta) if tau == 0 else source_difference(source_id, eta, tau)/tau
    points = sorted({0.0, length, max(0.0, eta+3.0), eta+5.0})
    parts = [integrate(integrand, left, right, config) for left, right in zip(points, points[1:]) if right > left]
    value = math.fsum(pair[0] for pair in parts)
    error = math.fsum(pair[1] for pair in parts)
    mass = math.sqrt(2.0)/(-eta)
    return {"y": -PREF*(value + f0*(math.log(mass*length) + GAMMA_E + 1.0)),
            "quadrature_error_estimate": PREF*error}


def finite_k_interchanged(source_id, eta, cutoff, config):
    """Exact finite momentum integral after exchanging its two finite integrals.

    Direct 2*sin(K*tau)^2 avoids cancellation in 1-cos(2*K*tau).
    Every panel spans at most a half-period; no sampled oscillation test tunes it.
    """
    if eta <= -5.0:
        return {"K": cutoff, "y": 0.0, "quadrature_error_estimate": 0.0, "time_panels": 0}
    left, right = max(0.0, eta+3.0), eta+5.0
    count = max(1, math.ceil((right-left)*2.0*cutoff/math.pi))
    panel_config = dict(config, epsabs=config["epsabs"]/count)
    def integrand(tau):
        return 0.0 if tau == 0.0 else source(source_id, eta-tau)*2.0*math.sin(cutoff*tau)**2/tau
    pieces = []
    for index in range(count):
        l = left+(right-left)*index/count
        r = left+(right-left)*(index+1)/count
        pieces.append(integrate(integrand, l, r, panel_config))
    mass = math.sqrt(2.0)/(-eta)
    local = source(source_id, eta)*(math.asinh(cutoff/mass)-cutoff/math.hypot(cutoff, mass))
    return {"K": cutoff, "y": PREF*(local-math.fsum(p[0] for p in pieces)),
            "quadrature_error_estimate": PREF*math.fsum(p[1] for p in pieces), "time_panels": count}


def load_tail_module():
    path = ROOT/"theory"/"tail_bounds.py"
    spec = importlib.util.spec_from_file_location("registered_tail_bounds", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def controls(config):
    state = []
    for separation in (0.1, 0.4):
        occupation = lambda k: 0.01*bump((k-1.0)/0.25)
        integral, error = integrate(lambda k: occupation(k)*math.sin(2.0*k*separation), 0.75, 1.25, config)
        initial, initial_error = integrate(lambda k: k*occupation(k)*math.cos(2.0*k*separation), 0.75, 1.25, config)
        state.append({"separation": separation, "a2_extra_occupied_kernel": -integral/(2.0*math.pi**2),
                      "occupied_quadrature_error_estimate": error/(2.0*math.pi**2),
                      "a2_initial_beta_variance": initial/(2.0*math.pi**2),
                      "beta_quadrature_error_estimate": initial_error/(2.0*math.pi**2)})
    advanced, advanced_error = integrate(lambda t: source("positive_B", t)/(t+5.5), -5.0, -3.0, config)
    with mpmath.workdps(60):
        qx_memory = -(2*mpmath.euler+mpmath.log(2))/(16*mpmath.pi**2)
        qx_action = (mpmath.digamma(2)+mpmath.digamma(1)-mpmath.log(2)-1)/(16*mpmath.pi**2)
        analytic = {"memory_Qx": float(qx_memory), "common_action_Qx": float(qx_action),
                    "omit_plus_one_Qx": float(qx_memory+1/(8*mpmath.pi**2)),
                    "omit_euler_gamma_Qx": float(qx_memory+mpmath.euler/(8*mpmath.pi**2)),
                    "moving_reference_Qx": float(qx_action+1/(16*mpmath.pi**2)),
                    "absolute_difference": float(abs(qx_memory-qx_action))}
    return {"state_diagnostics": state, "analytic_stationary": analytic,
            "advanced_response_y_at_minus5p5": -PREF*advanced,
            "advanced_quadrature_error_estimate": PREF*advanced_error,
            "conformal_input_mutation_at_minus2p5": logarithmic_memory("positive_B", -2.5, config, True)}


def run(experiment, progress):
    tail_module = load_tail_module()
    config = experiment["primary_quadrature"]
    rows = []
    for source_id in experiment["source"]["ids"]:
        for eta in experiment["observations"]:
            check_deadline()
            active = {"source": source_id, "eta": eta}
            progress({"rows": rows, "active": active})
            continuum = logarithmic_memory(source_id, eta, config)
            finitepart = finite_part_memory(source_id, eta, experiment["geometry"]["eta_i"], config)
            tail = tail_module.derivative_tail_budget(source_id, eta)
            f0 = source(source_id, eta)
            finite = []
            for cutoff in experiment["cutoffs"]:
                entry = finite_k_interchanged(source_id, eta, cutoff, config)
                # The only certified numerical error component is this analytic UV bound.
                # Strictly before support, the full combined integrand is exactly zero.
                entry["analytic_tail_bound"] = tail_module.normalized_combined_tail_bound(source_id, eta, cutoff)
                finite.append(entry)
                progress({"rows": rows, "active": dict(active, continuum=continuum,
                          finite_part=finitepart, tail_derivative_budget=tail, finite_k=finite)})
            rows.append({"source": source_id, "eta": eta, "a": -1.0/eta,
                         "s_over_epsilon": f0, "continuum": continuum,
                         "finite_part": finitepart, "tail_derivative_budget": tail,
                         "finite_k": finite})
            progress({"rows": rows, "active": None})
    progress({"rows": rows, "active": {"stage": "state_and_analytic_controls"}})
    return {"schema_version": 1, "normalization": experiment["normalization"], "rows": rows,
            "controls": controls(config),
            "uncertainty_note": "QUADPACK errors are numerical estimates, not rigorous enclosures. The combined-tail bound is analytic; roundoff has a separately registered allowance."}


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--registration-sha256", required=True)
    parser.add_argument("--public-freeze-commit", required=True)
    args = parser.parse_args()
    verify_registration(args.registration_sha256)
    require(len(args.public_freeze_commit) == 40 and all(c in "0123456789abcdef" for c in args.public_freeze_commit),
            "An explicit full public freeze commit is required")
    output = args.output.resolve()
    require(not output.exists(), "Output path already exists; refusing overwrite")
    require(not output.is_relative_to(ROOT), "Fresh output must be outside frozen checkpoint")
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    provenance = {"started_at_utc": datetime.now(timezone.utc).isoformat(),
                  "registration_sha256": args.registration_sha256,
                  "public_freeze_commit": args.public_freeze_commit,
                  "producer_sha256": sha(__file__), "experiment_sha256": sha(ROOT/"EXPERIMENT.json"),
                  "command": sys.argv, "python": platform.python_version(),
                  "numpy": numpy.__version__, "scipy": scipy.__version__, "mpmath": mpmath.__version__}
    write_json(output/"started.json", provenance)
    try:
        experiment = json.loads((ROOT/"EXPERIMENT.json").read_text())
        for package, actual, expected in (("numpy", numpy.__version__, "2.2.6"),
                                          ("scipy", scipy.__version__, "1.15.3"),
                                          ("mpmath", mpmath.__version__, "1.3.0")):
            require(actual == expected, "Frozen dependency version mismatch: "+package)
        DEADLINE = start+experiment["execution"]["wall_time_budget_seconds"]
        def progress(data):
            write_json(output/"partial_results.json", dict(data, provenance=provenance,
                       elapsed_seconds=time.monotonic()-start, status="running"))
        result = run(experiment, progress)
        result["provenance"] = provenance
        result["elapsed_seconds"] = time.monotonic()-start
        require(result["elapsed_seconds"] <= experiment["execution"]["wall_time_budget_seconds"], "Primary wall-time budget exceeded")
        write_json(output/"results.json", result)
        write_json(output/"EXECUTION.json", {"status": "completed", "elapsed_seconds": result["elapsed_seconds"],
                                             "results_sha256": sha(output/"results.json")})
        print("Primary causal response written to", output)
    except BaseException as exc:
        write_json(output/"failure.json", {"exception": repr(exc), "traceback": traceback.format_exc(),
                                           "elapsed_seconds": time.monotonic()-start})
        raise


if __name__ == "__main__":
    main()
