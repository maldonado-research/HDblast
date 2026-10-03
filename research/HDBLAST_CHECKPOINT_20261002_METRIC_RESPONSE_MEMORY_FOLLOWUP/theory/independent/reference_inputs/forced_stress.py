#!/usr/bin/env python3
"""Prospectively frozen, independent forced-mode minimal-stress experiment."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

import numpy as np
from numpy.polynomial.legendre import leggauss


LD = np.longdouble
CD = np.clongdouble
EPSILON = LD("0.0001")
ETA_I, ETA_F = LD(-6), LD("-1.5")
OBSERVATIONS = (-5.5, -4.5, -4.0, -3.5, -2.5, -1.5)
CUTOFFS = (64, 128, 256)
SOURCES = ("positive_B", "signed_uB")
SETTINGS = (("coarse", 128, LD("0.5")), ("fine", 256, LD("0.25")))
QUANTITIES = ("q", "q_prime", "q_second", "rho", "p", "Q0", "anomaly", "current")
REFINEMENT_GATES = {"q": 2e-10, "q_prime": 2e-10, "q_second": 1e-8,
                    "rho": 1e-8, "p": 1e-8}
WRONSKIAN_GATE = 1e-10
WARD_ENDPOINT_GATE = 2e-6
WARD_REFINEMENT_GATE = 1e-6
WALL_BUDGET = 900.0
RSS_BUDGET_KIB = 256*1024


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest_provenance(directory, expected_sha):
    path = directory / "MANIFEST.json"
    require(sha(path) == expected_sha, "Independent manifest differs from public pin")
    manifest = json.loads(path.read_text())
    for entry in manifest["files"]:
        require(sha(directory / entry["path"]) == entry["sha256"],
                "Independent frozen input mismatch: " + entry["path"])
    return {"manifest_sha256": expected_sha, "files": manifest["files"]}


def add_poly(left, right):
    result = [0]*max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply_poly(left, right):
    result = [0]*(len(left)+len(right)-1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            result[i+j] += first*second
    return result


def bump_derivative_polynomials():
    # Exact integer algebra only: B^(n)=B P_n/(1-u^2)^(2n).
    polynomials = [[1]]
    for n in range(5):
        previous = polynomials[-1]
        derivative = [i*value for i, value in enumerate(previous)][1:] or [0]
        term1 = multiply_poly([1, 0, -2, 0, 1], derivative)
        term2 = multiply_poly([0, 4*n-2, 0, -4*n], previous)
        polynomials.append(add_poly(term1, term2))
    return tuple(tuple(p) for p in polynomials)


BUMP_POLYNOMIALS = bump_derivative_polynomials()


def source_jet(eta, source_id, order=5):
    require(source_id in SOURCES, "Unknown source")
    u = np.asarray(eta, dtype=LD)+LD(4)
    inside = np.abs(u) < 1
    selected = u[inside]
    denominator = 1-selected*selected
    bump = np.exp(1-1/denominator)
    derivatives = []
    for n in range(order+1):
        polynomial = np.zeros_like(selected)
        for coefficient in reversed(BUMP_POLYNOMIALS[n]):
            polynomial = polynomial*selected+coefficient
        out = np.zeros_like(u)
        out[inside] = bump*polynomial/denominator**(2*n)
        derivatives.append(out)
    if source_id == "signed_uB":
        derivatives = [u*derivatives[n]+(n*derivatives[n-1] if n else 0)
                       for n in range(order+1)]
    return np.stack(derivatives)


def oscillatory_propagator(k, distance):
    # Shared expm1 evaluation implements the exact propagator; no W projection.
    increment = np.expm1(CD(2j)*k*distance)
    phi = np.empty(np.broadcast_shapes(np.shape(k), np.shape(distance)), dtype=CD)
    nonzero = np.broadcast_to(k != 0, phi.shape)
    np.divide(increment, np.broadcast_to(CD(2j)*k, phi.shape), out=phi, where=nonzero)
    np.copyto(phi, np.broadcast_to(distance, phi.shape), where=~nonzero)
    return increment+1, phi


def momentum_rule(width):
    nodes, weights = leggauss(16)
    nodes, weights = nodes.astype(LD), weights.astype(LD)
    n_panels = int(LD(256)/width)
    left = np.arange(n_panels, dtype=LD)*width
    k = (left[:, None]+(nodes[None, :]+1)*width/2).ravel()
    wk = np.broadcast_to(weights[None, :]*width/2, (n_panels, 16)).ravel().copy()
    return k, wk


def resource_check(started):
    elapsed = time.monotonic()-started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    require(elapsed <= WALL_BUDGET, "Frozen 900-second wall budget exceeded")
    require(peak <= RSS_BUDGET_KIB, "Frozen 256-MiB peak-memory budget exceeded")
    return {"elapsed_seconds": elapsed, "peak_rss_kib": peak}


def direct_observables(eta, k, u, w, jet):
    """Two direct minimal stresses and complete fixed-r W2/W4 subtraction."""
    a = L = -1/eta
    mass2 = 2*a*a
    s, s1, s2 = EPSILON*jet[:3]
    omega = np.sqrt(k*k+mass2)
    omega1 = L*mass2/omega
    omega2 = 3*L*L*mass2/omega-L*L*mass2*mass2/omega**3
    U = -L*L/omega-omega2/(4*omega**2)+3*omega1**2/(8*omega**3)
    variation_W2 = s/(2*omega)
    variation_W2_1 = s1/(2*omega)-s*omega1/(2*omega**2)
    variation_W2_2 = (s2/(2*omega)-s1*omega1/omega**2
                      -s*omega2/(2*omega**2)+s*omega1**2/omega**3)
    variation_W4 = (-U*variation_W2/omega-variation_W2_2/(4*omega**2)
                    +omega2*variation_W2/(4*omega**3)
                    +3*omega1*variation_W2_1/(4*omega**3)
                    -3*omega1**2*variation_W2/(4*omega**4))
    b = -k*k/3-mass2
    c = 1-b/omega**2
    j4 = (L*variation_W2_1/omega**2-2*L*omega1*variation_W2/omega**3
          +omega1*variation_W2_1/(2*omega**3)
          -3*omega1**2*variation_W2/(4*omega**4))
    # Differentiate the COMPLETE inherited rho2+rho4 and p2+p4.
    # These unsimplified density terms retain the cancellation 2 U deltaW2/omega=s U/omega^2.
    sub_rho = (s/omega+2*U*variation_W2/omega-s*U/omega**2
               -L*L*variation_W2/omega**2+j4)/4
    sub_p = (c*variation_W2-s/omega+c*variation_W4
             +2*b*U*variation_W2/omega**3+s*U/omega**2
             -L*L*variation_W2/omega**2+j4)/4
    # Actual complex canonical modes and their forced variations.
    v = np.exp(-CD(1j)*k*eta)/np.sqrt(2*k)
    vp = -CD(1j)*k*v
    delta_v = v*u
    delta_vp = v*(w-CD(1j)*k*u)
    Dv = vp-L*v
    delta_Dv = delta_vp-L*delta_v
    delta_abs_v = 2*np.real(np.conj(v)*delta_v)
    delta_abs_Dv = 2*np.real(np.conj(Dv)*delta_Dv)
    baseline_abs_v = 1/(2*k)
    explicit_mass_contact = s*baseline_abs_v
    bare_rho = (delta_abs_Dv+(k*k+mass2)*delta_abs_v+explicit_mass_contact)/2
    bare_p = (delta_abs_Dv-(k*k/3+mass2)*delta_abs_v-explicit_mass_contact)/2
    # True u,w equations supply variance derivatives. No Wronskian substitution.
    A = u.real/k
    A1 = w.real/k
    A2 = np.real(CD(2j)*k*w-s)/k
    C = -s/(4*omega**3)
    C1 = -s1/(4*omega**3)+3*s*omega1/(4*omega**4)
    C2 = (-s2/(4*omega**3)+3*s1*omega1/(2*omega**4)
          +3*s*omega2/(4*omega**4)-3*s*omega1**2/omega**5)
    baseline_subtraction = 1/(2*omega)-U/(2*omega**2)
    # Direct local subtraction mismatch. It is a diagnostic, never a definition of p.
    anomaly = (sub_rho-3*sub_p-mass2*C-s*baseline_subtraction
               +(C2-2*L*C1-2*L*L*C)/2)
    pi = np.arccos(LD(-1))
    measure = k*k/(2*pi*pi)
    integrands = {
        "q": measure*(A-C)/EPSILON,
        "q_prime": measure*(A1-C1)/EPSILON,
        "q_second": measure*(A2-C2)/EPSILON,
        "rho": measure*(bare_rho-sub_rho)/EPSILON,
        "p": measure*(bare_p-sub_p)/EPSILON,
        "Q0": measure*(baseline_abs_v-baseline_subtraction)/(a*a),
        "anomaly": measure*anomaly/EPSILON,
    }
    integrands["current"] = integrands["q"]+jet[0]*integrands["Q0"]/4
    linear_W = delta_v*np.conj(vp)+v*np.conj(delta_vp)-delta_vp*np.conj(v)-vp*np.conj(delta_v)
    physical_W = np.max(np.abs(linear_W))/EPSILON
    canonical_energy = 2*np.real(np.conj(vp)*delta_vp)+k*k*delta_abs_v
    diagnostics = {"physical_wronskian_scaled": float(physical_W),
                   "canonical_energy_max_over_epsilon": float(np.max(np.abs(canonical_energy))/EPSILON)}
    archival = {"bare_rho": bare_rho/EPSILON, "sub_rho": sub_rho/EPSILON,
                "bare_p": bare_p/EPSILON, "sub_p": sub_p/EPSILON,
                "delta_W2": variation_W2, "delta_W4": variation_W4}
    return integrands, diagnostics, archival


def cutoff_sums(integrands, weights, width):
    weighted = np.stack([integrands[name] for name in QUANTITIES])*weights[None, :]
    values = np.empty((len(CUTOFFS), len(QUANTITIES)), dtype=LD)
    for i, cutoff in enumerate(CUTOFFS):
        count = int(LD(cutoff)/width)*16
        values[i] = np.sum(weighted[:, :count], axis=1, dtype=LD)
    require(np.all(np.isfinite(values)), "Nonfinite combined cutoff response")
    return values


def simpson_to(history, end_step, h):
    require(end_step % 2 == 0, "Observation must align with an even Simpson endpoint")
    return h/3*(history[0]+history[end_step]
                +4*np.sum(history[1:end_step:2], axis=0, dtype=LD)
                +2*np.sum(history[2:end_step:2], axis=0, dtype=LD))


def evolve(source_id, setting, output_dir, started, progress, active_arrays):
    name, steps_per_unit, width = setting
    h = LD(1)/steps_per_unit
    n_steps = int((ETA_F-ETA_I)*steps_per_unit)
    k, weights = momentum_rule(width)
    gx, gw = leggauss(8)
    local_time = (gx.astype(LD)+1)*h/2
    local_weight = gw.astype(LD)*h/2
    E, drift = oscillatory_propagator(k, h)
    E_local, drift_local = oscillatory_propagator(k[:, None], (h-local_time)[None, :])
    force_w = E_local*local_weight[None, :]
    force_u = drift_local*local_weight[None, :]
    del E_local, drift_local
    u, w = np.zeros(k.shape, dtype=CD), np.zeros(k.shape, dtype=CD)
    require(np.count_nonzero(u) == np.count_nonzero(w) == 0, "Nonzero initial forced modes")
    active_arrays.update({"k": k, "weights": weights, "u": u, "w": w})
    observation_steps = {int((LD(eta)-ETA_I)*steps_per_unit): eta for eta in OBSERVATIONS}
    history_values = np.empty((n_steps+1, len(CUTOFFS), len(QUANTITIES)), dtype=LD)
    history_ledger = np.empty((n_steps+1, len(CUTOFFS)), dtype=LD)
    history_jets = np.empty((n_steps+1, 6), dtype=LD)
    history_d_prime = np.empty(n_steps+1, dtype=LD)
    history_eta = ETA_I+np.arange(n_steps+1, dtype=LD)*h
    rows, archive = [], {"k": k, "momentum_weights": weights}
    progress["completed_observations"] = rows
    amplitude_peak = physical_peak = energy_peak = 0.0
    for step in range(n_steps+1):
        eta = history_eta[step]
        if step:
            local_source = EPSILON*source_jet(eta-h+local_time, source_id, 0)[0]
            forcing_u = np.einsum("ij,j->i", force_u, local_source, optimize=False)
            forcing_w = np.einsum("ij,j->i", force_w, local_source, optimize=False)
            u_new = u+drift*w-forcing_u
            w_new = E*w-forcing_w
            u, w = u_new, w_new
        progress.update({"step": step, "eta": float(eta)})
        active_arrays.update({"u": u, "w": w})
        residual = 2*u.real-w.imag/k
        active_arrays["wronskian_residual"] = residual
        amplitude_W = float(np.max(np.abs(residual))/EPSILON)
        progress["wronskian_scaled"] = amplitude_W if math.isfinite(amplitude_W) else None
        require(np.all(np.isfinite(u)) and np.all(np.isfinite(w)), "Nonfinite mode evolution")
        require(amplitude_W <= WRONSKIAN_GATE,
                f"Amplitude Wronskian failed at {source_id}/{name}, step {step}")
        amplitude_peak = max(amplitude_peak, amplitude_W)
        jet = source_jet(eta, source_id)
        integrands, diagnostics, archival = direct_observables(eta, k, u, w, jet)
        require(diagnostics["physical_wronskian_scaled"] <= WRONSKIAN_GATE,
                f"Physical Wronskian failed at {source_id}/{name}, step {step}")
        physical_peak = max(physical_peak, diagnostics["physical_wronskian_scaled"])
        energy_peak = max(energy_peak, diagnostics["canonical_energy_max_over_epsilon"])
        values = cutoff_sums(integrands, weights, width)
        history_values[step] = values
        history_jets[step] = jet
        a = L = -1/eta
        d_prime = EPSILON*(jet[1]-2*L*jet[0])/(a*a)
        history_d_prime[step] = d_prime
        q0_index, rho_index, p_index = (QUANTITIES.index(key) for key in ("Q0", "rho", "p"))
        history_ledger[step] = (a**4*values[:, q0_index]*d_prime/(2*EPSILON)
                                +L*values[:, rho_index]-3*L*values[:, p_index])
        active_arrays.update({"history_eta": history_eta[:step+1],
                              "history_values": history_values[:step+1],
                              "history_ledger": history_ledger[:step+1],
                              "source_jet": jet})
        if eta <= -5:
            perturbation_indices = [i for i, key in enumerate(QUANTITIES) if key != "Q0"]
            require(np.count_nonzero(u) == np.count_nonzero(w) == 0
                    and np.count_nonzero(values[:, perturbation_indices]) == 0,
                    "Strict source-free initial interval failed")
        if step % 32 == 0:
            resource_check(started)
        if step not in observation_steps:
            continue
        ledger = simpson_to(history_ledger, step, h)
        finite_k = []
        for cutoff_index, cutoff in enumerate(CUTOFFS):
            density = values[cutoff_index, rho_index]
            finite_k.append({"K": cutoff,
                             "values": {key: float(values[cutoff_index, j]) for j, key in enumerate(QUANTITIES)},
                             "ward": {"ledger": float(ledger[cutoff_index]),
                                      "direct_endpoint": float(density),
                                      "endpoint_difference": float(abs(ledger[cutoff_index]-density))}})
        row = {"source": source_id, "eta": float(eta), "a": float(a),
               "source_jet_over_epsilon": {f"f{i}": float(jet[i]) for i in range(6)},
               "wronskian_scaled": amplitude_W, **diagnostics, "finite_k": finite_k}
        rows.append(row)
        index = len(rows)-1
        archive[f"u_{index}"] = u.copy()
        archive[f"w_{index}"] = w.copy()
        archive[f"source_jet_{index}"] = jet.copy()
        for key in QUANTITIES:
            archive[f"integrand_{key}_{index}"] = integrands[key]
        for key, value in archival.items():
            archive[f"{key}_{index}"] = value
        # Retain completed snapshots for a possible exception later in this run.
        active_arrays.update({f"archived_{key}": value for key, value in archive.items()})
    require(len(rows) == len(OBSERVATIONS), "Missing forced-mode observation")
    archive.update({"history_eta": history_eta, "history_values": history_values,
                    "history_ledger_integrand": history_ledger,
                    "history_source_jet": history_jets, "history_d_prime": history_d_prime,
                    "history_quantity_names": np.array(QUANTITIES),
                    "history_cutoffs": np.array(CUTOFFS, dtype=np.int64)})
    archive_path = output_dir/f"stress_modes_{source_id}_{name}.npz"
    np.savez_compressed(archive_path, **archive)
    return {"source": source_id, "setting": name, "dt": float(h),
            "momentum_panel_width": float(width), "time_steps": n_steps,
            "momentum_nodes": len(k), "time_gauss_nodes": 8,
            "local_time_momentum_pairs": n_steps*8*len(k),
            "direct_stress_history_points": n_steps+1,
            "wronskian_max_scaled": amplitude_peak,
            "physical_wronskian_max_scaled": physical_peak,
            "canonical_energy_max_over_epsilon": energy_peak,
            "archive": {"path": archive_path.name, "sha256": sha(archive_path)},
            "rows": rows}


def combine_runs(runs):
    lookup = {(r["source"], r["setting"]): r for r in runs}
    maxima = {key: 0.0 for key in QUANTITIES}
    ward_max = ward_refinement_max = 0.0
    rows = []
    failures = []
    for source_id in SOURCES:
        for coarse, fine in zip(lookup[(source_id, "coarse")]["rows"], lookup[(source_id, "fine")]["rows"]):
            require(coarse["eta"] == fine["eta"], "Coarse/fine observation mismatch")
            row = {key: value for key, value in fine.items() if key != "finite_k"}
            row["finite_k"] = []
            for c, f in zip(coarse["finite_k"], fine["finite_k"]):
                differences = {key: abs(f["values"][key]-c["values"][key]) for key in QUANTITIES}
                for key, difference in differences.items():
                    maxima[key] = max(maxima[key], difference)
                    if key in REFINEMENT_GATES and difference > REFINEMENT_GATES[key]:
                        failures.append(f"{source_id}/{fine['eta']}/{f['K']}/{key} refinement")
                ward_difference = abs(f["ward"]["ledger"]-c["ward"]["ledger"])
                ward_max = max(ward_max, f["ward"]["endpoint_difference"])
                ward_refinement_max = max(ward_refinement_max, ward_difference)
                if f["ward"]["endpoint_difference"] > WARD_ENDPOINT_GATE:
                    failures.append(f"{source_id}/{fine['eta']}/{f['K']} Ward endpoint")
                if ward_difference > WARD_REFINEMENT_GATE:
                    failures.append(f"{source_id}/{fine['eta']}/{f['K']} Ward refinement")
                row["finite_k"].append({"K": f["K"], "values": f["values"],
                                         "coarse_values": c["values"],
                                         "refinement_differences": differences,
                                         "estimated_numerical_errors": {key: 2*value for key, value in differences.items()},
                                         "ward": {**f["ward"], "coarse_ledger": c["ward"]["ledger"],
                                                  "coarse_endpoint_difference": c["ward"]["endpoint_difference"],
                                                  "refinement_difference": ward_difference}})
            rows.append(row)
    return rows, maxima, ward_max, ward_refinement_max, failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-commit", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    require(len(args.freeze_commit) == 40
            and all(x in "0123456789abcdef" for x in args.freeze_commit), "Public freeze commit required")
    require(len(args.manifest_sha256) == 64, "Public manifest SHA-256 required")
    provenance = manifest_provenance(Path(__file__).resolve().parent, args.manifest_sha256)
    require(sys.version_info[:2] == (3, 12) and np.__version__ == "2.2.6", "Pinned Python/NumPy runtime required")
    require(np.finfo(LD).nmant >= 63, "Extended precision with at least 63 fraction bits required")
    checkpoint = Path(__file__).resolve().parents[1]
    require(not args.output_dir.resolve().is_relative_to(checkpoint),
            "Output must be external to the frozen checkpoint")
    require(not args.output_dir.exists(), "Refusing to replace existing output directory")
    args.output_dir.mkdir(parents=True)
    started = time.monotonic()
    report = {"schema_version": 1, "route": "independent_forced_modes_direct_minimal_stress",
              "freeze_commit": args.freeze_commit, "provenance": provenance,
              "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                          "platform": platform.platform(), "python_optimization": sys.flags.optimize,
                          "real_dtype": np.dtype(LD).name, "complex_dtype": np.dtype(CD).name,
                          "longdouble_nmant": np.finfo(LD).nmant,
                          "longdouble_eps": str(np.finfo(LD).eps)},
              "epsilon": float(EPSILON), "mass_law_b": 1,
              "normalizations": {"q": "a^2 deltaQ/epsilon", "q_prime": "(a^2 deltaQ)'/epsilon",
                                 "q_second": "(a^2 deltaQ)''/epsilon", "rho": "a^4 delta_rho/epsilon",
                                 "p": "a^4 delta_p/epsilon", "Q0": "physical Q0,K",
                                 "anomaly": "a^4 deltaA_K/epsilon", "current": "a^2 delta_j/epsilon"},
              "status": "running", "runs": [], "active_run": None}
    result_path = args.output_dir/"results.json"
    active_arrays = {}
    try:
        for source_id in SOURCES:
            for setting in SETTINGS:
                progress = {"source": source_id, "setting": setting[0], "step": 0, "eta": float(ETA_I)}
                report["active_run"] = progress
                result_path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
                run = evolve(source_id, setting, args.output_dir, started, progress, active_arrays)
                report["runs"].append(run)
                report["active_run"] = None
                active_arrays.clear()
                report["resources"] = resource_check(started)
                result_path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
                print(json.dumps({"source": source_id, "setting": setting[0],
                                  "status": "completed", "resources": report["resources"]}), flush=True)
        rows, maxima, ward_max, ward_refinement, failures = combine_runs(report["runs"])
        report.update({"rows": rows, "maximum_refinement_differences": maxima,
                       "maximum_ward_endpoint_difference": ward_max,
                       "maximum_ward_refinement_difference": ward_refinement,
                       "gate_failures": failures})
        require(not failures, "Registered internal gates failed: "+", ".join(failures))
        report["status"] = "passed_internal_gates"
        report["limits"] = [
            "Both stresses are direct physical minimal-mode bilinears with explicit mass contacts and complete W2/W4 subtraction.",
            "The Simpson Ward ledger is independently integrated from recorded direct stresses; it defines neither stress.",
            "Refinement and cross-route differences are numerical estimates, not certified integration error bounds.",
            "Stress and derivative UV-tail enclosures must be supplied separately; the variance tail alone is insufficient.",
            "One fixed geometry, unchanged BD state, linear response only; no shell evolution, stability, particle yield, or heating."]
    except Exception as exc:
        report.update({"status": "failed", "exception": {"type": type(exc).__name__,
                       "message": str(exc), "traceback": traceback.format_exc()}})
        if active_arrays:
            snapshot = args.output_dir/"failure_snapshot.npz"
            try:
                np.savez(snapshot, **active_arrays)
                report["failure_snapshot"] = {"path": snapshot.name, "sha256": sha(snapshot)}
            except Exception as snapshot_error:
                report["failure_snapshot_error"] = str(snapshot_error)
        raise
    finally:
        report["resources"] = {"elapsed_seconds": time.monotonic()-started,
                               "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        result_path.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": report["status"], "result": str(result_path),
                      "maximum_refinement_differences": maxima,
                      "maximum_ward_endpoint_difference": ward_max,
                      "maximum_ward_refinement_difference": ward_refinement}), flush=True)


if __name__ == "__main__":
    main()
