#!/usr/bin/env python3
"""Registered primary matched scalar-variance, stress and current response.

No source, quadrature, mode or physical-response evaluation occurs on import.
The public prospective registration must be pinned before running this file.
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
import sympy
from scipy.integrate import IntegrationWarning, quad

from source_jet import source_jet_over_epsilon

ROOT = Path(__file__).resolve().parents[1]
PREF = 1/(8*math.pi**2)
DEADLINE = None


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, document):
    Path(path).write_text(json.dumps(document, indent=2, sort_keys=True, allow_nan=False)+"\n")


def check_deadline():
    if DEADLINE is not None:
        require(time.monotonic() <= DEADLINE, "Frozen primary wall-time budget exceeded")


def verify_registration(pin):
    path = ROOT/"FULL_REGISTRATION.json"
    require(sha(path) == pin, "Prospective registration SHA256 mismatch")
    registration = json.loads(path.read_text())
    for name, digest in registration["files"].items():
        target = (ROOT/name).resolve()
        require(target.is_relative_to(ROOT), "Registered path escapes checkpoint")
        require(sha(target) == digest, "Frozen file changed: "+name)
    required = {"code/stress_primary.py", "code/source_jet.py", "EXPERIMENT.json", "theory/stress_tail_bounds.py"}
    require(required <= set(registration["files"]), "Required primary input missing from prospective manifest")
    return registration


def integrate(function, left, right, config, **kwargs):
    check_deadline()
    if right <= left:
        return 0.0, 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        value, error = quad(function, left, right, epsabs=config["epsabs"],
                            epsrel=config["epsrel"], limit=config["limit"], **kwargs)
    check_deadline()
    require(math.isfinite(value) and math.isfinite(error), "Nonfinite primary quadrature")
    return float(value), float(error)


def log_history(source, eta, order, jet, config):
    """F[f^(order)] with its integrable logarithmic endpoint treated explicitly."""
    require(order in (1, 2, 3), "Unregistered logarithmic derivative order")
    if eta <= -5:
        return 0.0, 0.0
    right = min(eta, -3.0)
    derivative = lambda t: source_jet_over_epsilon(source, t)[order]
    if right == eta:
        value, error = integrate(derivative, -5.0, right, config, weight="alg-logb", wvar=(0.0, 0.0))
    else:
        value, error = integrate(lambda t: derivative(t)*math.log(eta-t), -5.0, right, config)
    M = math.sqrt(2.0)/(-eta)
    constant = math.log(M)+float(numpy.euler_gamma)+1
    # The compact source has f^(order-1)(eta_i)=0 exactly.
    return value+constant*jet[order-1], error


def continuum_q_jet(source, eta, jet, config):
    memories = [log_history(source, eta, n, jet, config) for n in (1, 2, 3)]
    L = -1/eta
    q = -PREF*memories[0][0]
    q_prime = -PREF*(memories[1][0]+L*jet[0])
    q_second = -PREF*(memories[2][0]+2*L*jet[1]+L*L*jet[0])
    return [q, q_prime, q_second], [PREF*entry[1] for entry in memories], memories


def local_coefficients(eta, cutoff=None):
    """Exact integrated fixed-reference coefficients, specialized to H=1,r=2."""
    L = a = -1/eta
    if cutoff is None:
        return {"v": 1.0, "Q0": 1/(12*math.pi**2), "Q0_prime": 0.0,
                "J5": 1/6, "J7": 1/30, "J9": 1/105,
                "density_contact_factor": 1.0, "density_contact_factor_prime": 0.0}
    M = math.sqrt(2.0)*a
    v = cutoff/math.hypot(cutoff, M)
    v_prime = -L*v*(1-v*v)
    baseline_polynomial = v*v/(2*(1+v))-v**3/48-v**5/16
    baseline_derivative = v*(2+v)/(2*(1+v)**2)-v*v/16-5*v**4/16
    return {"v": v, "v_prime": v_prime, "M": M,
            "A": math.asinh(cutoff/M)-v, "A_prime": -L*v**3,
            "A_second": L*L*v**3*(2-3*v*v),
            "Q0": baseline_polynomial/(2*math.pi**2),
            "Q0_prime": baseline_derivative*v_prime/(2*math.pi**2),
            "J5": v**3/6, "J7": (v**3/3-v**5/5)/4,
            "J9": (v**3/3-2*v**5/5+v**7/7)/8,
            "density_contact_factor": v**3,
            "density_contact_factor_prime": -3*L*v**3*(1-v*v)}


def finite_history(source, eta, cutoff, order, config):
    """The exact finite-band memory after exchanging two finite integrals."""
    if eta <= -5:
        return 0.0, 0.0, 0
    left, right = max(0.0, eta+3.0), eta+5.0
    count = max(1, math.ceil(2*cutoff*(right-left)/math.pi))
    panel_config = dict(config, epsabs=config["epsabs"]/count)
    def function(tau):
        if tau == 0:
            return 0.0
        return source_jet_over_epsilon(source, eta-tau)[order]*2*math.sin(cutoff*tau)**2/tau
    parts = []
    for panel in range(count):
        l = left+(right-left)*panel/count
        r = left+(right-left)*(panel+1)/count
        parts.append(integrate(function, l, r, panel_config))
    return math.fsum(p[0] for p in parts), math.fsum(p[1] for p in parts), count


def finite_q_jet(source, eta, cutoff, jet, coefficients, config):
    memories = [finite_history(source, eta, cutoff, n, config) for n in range(3)]
    A, Ap, App = (coefficients[key] for key in ("A", "A_prime", "A_second"))
    q = PREF*(jet[0]*A-memories[0][0])
    qp = PREF*(jet[1]*A+jet[0]*Ap-memories[1][0])
    qpp = PREF*(jet[2]*A+2*jet[1]*Ap+jet[0]*App-memories[2][0])
    return [q, qp, qpp], [PREF*m[1] for m in memories], memories


def closed_response(eta, jet, q_jet, q_errors, coefficients):
    """Use the direct subtraction reductions; neither stress is defined by Ward.

    Pressure uses the independently integrated direct fourth-order subtraction
    mismatch (J5,J7,J9), not a trace reconstruction.
    """
    a = L = -1/eta
    f, fp, fpp = jet[:3]
    q, qp, qpp = q_jet
    Q0 = coefficients["Q0"]
    J5, J7, J9 = (coefficients[name] for name in ("J5", "J7", "J9"))
    rho_local_base = (3*L*L*f-L*fp)/(96*math.pi**2)
    rho_contact = rho_local_base*coefficients["density_contact_factor"]
    pressure_contact = (70*L*L*f*J9-(30*L*L*f+10*L*fp)*J7
                        +(fpp-L*fp-9*L*L*f)*J5)/(48*math.pi**2)
    anomaly = ((fpp-12*L*L*f)*J5-(30*L*L*f+10*L*fp)*J7+70*L*L*f*J9)/(16*math.pi**2)
    rho_response = (3*L*L*q-L*qp)/2
    rho_mass_contact = a*a*Q0*f/2
    pressure_response = (qpp-3*L*qp-3*L*L*q)/6
    pressure_mass_contact = -a*a*Q0*f/6
    rho = rho_response+rho_mass_contact+rho_contact
    pressure = pressure_response+pressure_mass_contact+pressure_contact
    current_mass_contact = Q0*f/4
    values = {"q": q, "q_prime": qp, "q_second": qpp, "rho": rho,
              "p": pressure, "Q0": Q0, "anomaly": anomaly, "current": q+current_mass_contact}
    e0, e1, e2 = q_errors
    estimates = {"q": e0, "q_prime": e1, "q_second": e2,
                 "rho": (3*L*L*e0+L*e1)/2, "p": (e2+3*L*e1+3*L*L*e0)/6,
                 "current": e0, "Q0": 0.0, "anomaly": 0.0}
    variance_jet = [q/a**2, (qp-2*L*q)/a**2, (qpp-4*L*qp+2*L*L*q)/a**2]
    # An analytic derivative of the independently defined density is available
    # for separate algebraic checks. It is not obtained by rearranging Ward.
    rho_local_base_prime = (6*L**3*f+2*L*L*fp-L*fpp)/(96*math.pi**2)
    contact_prime = (rho_local_base_prime*coefficients["density_contact_factor"]
                     +rho_local_base*coefficients["density_contact_factor_prime"])
    rho_prime = (-3*L**3*q+3*L*L*qp-L*qpp/2
                 +a*a*((coefficients["Q0_prime"]/2-L*Q0)*f+Q0*fp/2)
                 +contact_prime-4*L*rho_contact)
    require(all(math.isfinite(value) for value in values.values()), "Nonfinite closed response")
    return {"values": values, "quadrature_error_estimates": estimates,
            "variance_jet_over_epsilon": variance_jet,
            "components": {"rho_variance_response": rho_response,
                           "rho_reduced_baseline_mass_contact": rho_mass_contact,
                           "rho_local_subtraction_contact": rho_contact,
                           "pressure_variance_response": pressure_response,
                           "pressure_reduced_baseline_mass_contact": pressure_mass_contact,
                           "pressure_direct_subtraction_contact": pressure_contact,
                           "current_mass_law_contact": current_mass_contact},
            "density_derivative": {"a4_delta_rho_prime_over_epsilon": rho_prime,
                                   "quadrature_error_estimate": 3*L**3*e0+3*L*L*e1+L*e2/2}}


def load_tail_module():
    path = ROOT/"theory/stress_tail_bounds.py"
    spec = importlib.util.spec_from_file_location("registered_stress_tail_bounds", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def primary_settings(experiment):
    """Translate the shared frozen registration, checking specialized geometry."""
    require(experiment["geometry"]["H"] == 1 and experiment["geometry"]["r"] == 2
            and experiment["geometry"]["xi"] == 0, "Producer is specialized to H=1,r=2,xi=0")
    require(experiment["initial_eta"] == -6 and experiment["final_eta"] == -1.5,
            "Unregistered time domain")
    require(experiment["source"]["ids"] == ["positive_B", "signed_uB"]
            and experiment["source"]["support"] == [-5, -3], "Unregistered source class")
    require(experiment["observation_eta"] == [-5.5, -4.5, -4, -3.5, -2.5, -1.5]
            and experiment["cutoffs"] == [64, 128, 256], "Unregistered response grid")
    return {"sources": experiment["source"]["ids"], "epsilon": experiment["source"]["epsilon"],
            "mass_law_b": experiment["source"]["mass_law_translation_b"],
            "observations": experiment["observation_eta"], "cutoffs": experiment["cutoffs"],
            "quadrature": experiment["primary"], "normalization": experiment["normalization"],
            "wall_time_budget_seconds": experiment["primary"]["producer_wall_seconds"]}


def run(settings, progress):
    tails = load_tail_module()
    config = settings["quadrature"]
    rows = []
    for source in settings["sources"]:
        for eta in settings["observations"]:
            check_deadline()
            progress({"rows": rows, "active": {"source": source, "eta": eta}})
            jet = source_jet_over_epsilon(source, eta)
            q, errors, memories = continuum_q_jet(source, eta, jet, config)
            coefficients = local_coefficients(eta)
            continuum = closed_response(eta, jet, q, errors, coefficients)
            continuum["local_coefficients"] = coefficients
            continuum["log_history_integrals_over_epsilon"] = [
                {"source_derivative_order": n+1, "value": item[0], "quadrature_error_estimate": item[1]}
                for n, item in enumerate(memories)]
            row = {"source": source, "eta": eta, "a": -1/eta,
                   "source_jet_over_epsilon": {"f"+str(n): value for n,value in enumerate(jet)},
                   "continuum": continuum, "finite_k": []}
            for cutoff in settings["cutoffs"]:
                coefficients = local_coefficients(eta, cutoff)
                q, errors, memories = finite_q_jet(source, eta, cutoff, jet, coefficients, config)
                finite = closed_response(eta, jet, q, errors, coefficients)
                finite.update(K=cutoff, local_coefficients=coefficients,
                              analytic_tail_bounds=tails.stress_tail_bounds(source, eta, cutoff),
                              finite_history_integrals_over_epsilon=[
                                  {"source_derivative_order": n, "value": item[0],
                                   "quadrature_error_estimate": item[1], "time_panels": item[2]}
                                  for n,item in enumerate(memories)])
                row["finite_k"].append(finite)
                progress({"rows": rows, "active": row})
            rows.append(row)
            progress({"rows": rows, "active": None})
    return {"schema_version": 1, "route": "logarithmic_variance_derivatives_and_direct_closed_stress",
            "normalization": settings["normalization"], "epsilon": settings["epsilon"], "mass_law_b": settings["mass_law_b"],
            "rows": rows, "uncertainty_scope": "Quadrature errors are numerical estimates and local coefficient arithmetic has floating-point uncertainty. Only separately recorded analytic removed-cutoff tails are interval enclosures.",
            "stress_definition": "Both density and pressure use direct matched subtraction reductions; neither is defined by the Ward identity or trace.",
            "scope": "Fixed-geometry linear homogeneous scalar-to-stress/current response; no metric kernels, coupled shell evolution, stability, particle yield or heating claim."}


def main():
    global DEADLINE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--registration-sha256", required=True)
    parser.add_argument("--public-freeze-commit", required=True)
    args = parser.parse_args()
    verify_registration(args.registration_sha256)
    require(len(args.public_freeze_commit) == 40 and all(c in "0123456789abcdef" for c in args.public_freeze_commit),
            "A full explicit public prospective freeze commit is required")
    settings = primary_settings(json.loads((ROOT/"EXPERIMENT.json").read_text()))
    versions = {"numpy": "2.2.6", "scipy": "1.15.3", "mpmath": "1.3.0", "sympy": "1.14.0"}
    for name, module in (("numpy", numpy), ("scipy", scipy), ("mpmath", mpmath), ("sympy", sympy)):
        require(module.__version__ == versions[name], "Frozen dependency mismatch: "+name)
    require(settings["mass_law_b"] == 1.0 and settings["epsilon"] == 1e-4, "Unregistered current/source amplitude")
    output = args.output.resolve()
    require(not output.exists(), "Refusing to overwrite a previous run")
    require(not output.is_relative_to(ROOT), "Fresh output must be outside frozen checkpoint")
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    DEADLINE = start+settings["wall_time_budget_seconds"]
    provenance = {"started_at_utc": datetime.now(timezone.utc).isoformat(), "command": sys.argv,
                  "registration_sha256": args.registration_sha256, "public_freeze_commit": args.public_freeze_commit,
                  "producer_sha256": sha(__file__), "source_jet_sha256": sha(ROOT/"code/source_jet.py"),
                  "experiment_sha256": sha(ROOT/"EXPERIMENT.json"), "python": platform.python_version(),
                  "dependencies": {name: module.__version__ for name,module in
                                   (("numpy", numpy), ("scipy", scipy), ("mpmath", mpmath), ("sympy", sympy))}}
    write_json(output/"started.json", provenance)
    try:
        def progress(data):
            write_json(output/"partial_results.json", dict(data, provenance=provenance,
                       elapsed_seconds=time.monotonic()-start, status="running"))
        result = run(settings, progress)
        check_deadline()
        result.update(provenance=provenance, elapsed_seconds=time.monotonic()-start, status="completed")
        write_json(output/"results.json", result)
        write_json(output/"EXECUTION.json", {"status": "completed", "elapsed_seconds": result["elapsed_seconds"],
                                             "results_sha256": sha(output/"results.json")})
        print("Primary matched stress result written to", output, flush=True)
    except BaseException as exc:
        write_json(output/"failure.json", {"exception": repr(exc), "traceback": traceback.format_exc(),
                                           "elapsed_seconds": time.monotonic()-start})
        raise


if __name__ == "__main__":
    main()
