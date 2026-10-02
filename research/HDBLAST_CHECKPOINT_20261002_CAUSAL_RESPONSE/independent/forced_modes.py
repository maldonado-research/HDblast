#!/usr/bin/env python3
"""Independent, prospectively frozen forced canonical-mode calibration."""
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


EPSILON = 1e-4
ETA_I = -6.0
ETA_F = -1.5
OBSERVATIONS = (-5.5, -4.5, -4.0, -3.5, -2.5, -1.5)
CUTOFFS = (64, 128, 256)
SOURCES = ("positive_B", "signed_uB")
SETTINGS = (("coarse", 128, 0.5), ("fine", 256, 0.25))
CONVERGENCE_GATE = 2e-10
WRONSKIAN_GATE = 1e-10
WALL_BUDGET_SECONDS = 900.0
MEMORY_BUDGET_KIB = 128 * 1024


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_manifest(directory, expected_manifest_sha256):
    manifest_path = directory / "MANIFEST.json"
    require(sha256(manifest_path) == expected_manifest_sha256,
            "Independent manifest does not match the public prospective pin")
    manifest = json.loads(manifest_path.read_text())
    for entry in manifest["files"]:
        path = Path(entry["path"])
        if not path.is_absolute():
            path = directory / path
        require(sha256(path) == entry["sha256"], f"Source pin changed: {path}")
    return {"manifest_sha256": sha256(manifest_path), "files": manifest["files"]}


def source(eta, source_id):
    """Physical forcing s=a^2 delta x; support is exactly (-5,-3)."""
    u = np.asarray(eta, dtype=float) + 4.0
    result = np.zeros_like(u)
    inside = np.abs(u) < 1.0
    q = u[inside]
    result[inside] = EPSILON * np.exp(1.0 - 1.0 / (1.0 - q * q))
    if source_id == "signed_uB":
        result *= u
    elif source_id != "positive_B":
        raise ValueError(f"Unknown source ID: {source_id}")
    return result


def phi(k, distance):
    """Stable exact oscillatory drift coefficient, including its k=0 limit."""
    k = np.asarray(k)
    d = np.asarray(distance)
    z = 2j * k * d
    result = np.empty(np.broadcast_shapes(k.shape, d.shape), dtype=complex)
    nonzero = np.broadcast_to(k != 0, result.shape)
    numerator = np.expm1(z)
    denominator = np.broadcast_to(2j * k, result.shape)
    np.divide(numerator, denominator, out=result, where=nonzero)
    np.copyto(result, np.broadcast_to(d, result.shape), where=~nonzero)
    return result


def momentum_rule(panel_width):
    x, w = leggauss(16)
    n_panels = int(CUTOFFS[-1] / panel_width)
    left = np.arange(n_panels) * panel_width
    nodes = left[:, None] + (x[None, :] + 1.0) * (panel_width / 2.0)
    weights = np.broadcast_to(w[None, :] * (panel_width / 2.0), nodes.shape)
    return nodes.ravel(), weights.ravel(), n_panels


def check_resources(started):
    elapsed = time.monotonic() - started
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    require(elapsed <= WALL_BUDGET_SECONDS, "Frozen 15-minute wall budget exceeded")
    require(peak_kib <= MEMORY_BUDGET_KIB, "Frozen 128-MiB peak-memory budget exceeded")
    return {"elapsed_seconds": elapsed, "peak_rss_kib": peak_kib}


def evolve(source_id, setting, output_dir, started, progress, active_arrays):
    name, steps_per_unit, panel_width = setting
    h = 1.0 / steps_per_unit
    n_steps = int((ETA_F - ETA_I) * steps_per_unit)
    k, momentum_weights, n_panels = momentum_rule(panel_width)
    gx, gw = leggauss(8)
    local_t = (gx + 1.0) * h / 2.0
    local_w = gw * h / 2.0
    distance = h - local_t
    propagator = np.exp(2j * k * h)
    drift = phi(k, h)
    forcing_w = np.exp(2j * k[:, None] * distance[None, :]) * local_w[None, :]
    forcing_u = phi(k[:, None], distance[None, :]) * local_w[None, :]
    u = np.zeros(k.shape, dtype=complex)
    w = np.zeros(k.shape, dtype=complex)
    progress.update({"step": 0, "eta": ETA_I, "momentum_nodes": len(k),
                     "wronskian_scaled": 0.0})
    active_arrays.update({"k": k, "u": u, "w": w})
    require(np.count_nonzero(u) == 0 and np.count_nonzero(w) == 0,
            "Forced-mode initial data must be zero")
    observe_steps = {int((eta - ETA_I) * steps_per_unit): eta for eta in OBSERVATIONS}
    peak_wronskian = 0.0
    physical_peak = 0.0
    rows = []
    archive = {"k": k, "momentum_weights": momentum_weights}
    for step in range(1, n_steps + 1):
        eta_left = ETA_I + (step - 1) * h
        local_source = source(eta_left + local_t, source_id)
        # Independent local forced-ODE update; no full-history kernel is evaluated.
        source_u = np.einsum("ij,j->i", forcing_u, local_source, optimize=False)
        source_w = np.einsum("ij,j->i", forcing_w, local_source, optimize=False)
        u_new = u + drift * w - source_u
        w_new = propagator * w - source_w
        u, w = u_new, w_new
        progress.update({"step": step, "eta": ETA_I + step * h})
        active_arrays.update({"u": u, "w": w})
        residual = 2.0 * u.real - w.imag / k
        step_residual = float(np.max(np.abs(residual)) / EPSILON)
        progress["wronskian_scaled"] = step_residual if math.isfinite(step_residual) else None
        active_arrays["wronskian_residual"] = residual
        peak_wronskian = max(peak_wronskian, step_residual)
        require(np.all(np.isfinite(u)) and np.all(np.isfinite(w)), "Nonfinite mode state")
        require(step_residual <= WRONSKIAN_GATE,
                f"Linear amplitude Wronskian gate failed at {source_id}/{name} step {step}: {step_residual}")
        if step % 64 == 0:
            check_resources(started)
        if step not in observe_steps:
            continue
        eta = observe_steps[step]
        a = -1.0 / eta
        mass = a * math.sqrt(2.0)
        s0 = float(source(np.array(eta), source_id))
        # Reconstruct canonical modes to check the same invariant directly.
        v = np.exp(-1j * k * eta) / np.sqrt(2.0 * k)
        vp = -1j * k * v
        delta_v = v * u
        delta_vp = v * (w - 1j * k * u)
        physical_residual = (delta_v * vp.conj() + v * delta_vp.conj()
                             - delta_vp * v.conj() - vp * delta_v.conj())
        physical_scaled = float(np.max(np.abs(physical_residual)) / EPSILON)
        physical_peak = max(physical_peak, physical_scaled)
        require(physical_scaled <= WRONSKIAN_GATE, "Physical linear Wronskian gate failed")
        mode_term = 4.0 * k * u.real
        subtraction = s0 * k * k / (k * k + mass * mass) ** 1.5
        combined = mode_term + subtraction
        # Combine before summing. K endpoints contain complete GL panels.
        terms = (combined * momentum_weights).reshape(n_panels, 16).sum(axis=1)
        finite_k = []
        for cutoff in CUTOFFS:
            count = int(cutoff / panel_width)
            y = math.fsum(float(x) for x in terms[:count]) / (8.0 * math.pi**2 * EPSILON)
            require(math.isfinite(y), "Nonfinite finite-K response")
            if eta < -5.0:
                require(y == 0.0 and np.count_nonzero(u) == 0
                        and np.count_nonzero(w) == 0, "Strict pre-pulse causality failed")
            finite_k.append({"K": cutoff, "y": y})
        high = k >= 128.0
        cancel_stats = {
            "max_abs_mode_term_k_ge_128": float(np.max(np.abs(mode_term[high]))),
            "max_abs_subtraction_k_ge_128": float(np.max(np.abs(subtraction[high]))),
            "max_abs_combined_k_ge_128": float(np.max(np.abs(combined[high]))),
        }
        rows.append({"source": source_id, "eta": eta, "a": a,
                     "s_over_epsilon": s0 / EPSILON, "finite_k": finite_k,
                     "wronskian_scaled": step_residual,
                     "physical_wronskian_scaled": physical_scaled,
                     "tail_cancellation": cancel_stats})
        index = len(rows) - 1
        archive[f"u_{index}"] = u.copy()
        archive[f"w_{index}"] = w.copy()
        archive[f"subtraction_{index}"] = subtraction
        archive[f"combined_{index}"] = combined
    require(len(rows) == len(OBSERVATIONS), "Missing observation")
    archive_path = output_dir / f"modes_{source_id}_{name}.npz"
    np.savez_compressed(archive_path, **archive)
    return {"setting": name, "source": source_id, "dt": h,
            "momentum_panel_width": panel_width, "time_steps": n_steps,
            "momentum_nodes": len(k), "time_gauss_nodes": 8,
            "local_time_momentum_pairs": n_steps * 8 * len(k),
            "wronskian_max_scaled": peak_wronskian,
            "physical_wronskian_max_scaled": physical_peak,
            "archive": {"path": archive_path.name, "sha256": sha256(archive_path)},
            "rows": rows}


def combine_results(runs):
    lookup = {(run["source"], run["setting"]): run for run in runs}
    rows = []
    maximum_difference = 0.0
    for source_id in SOURCES:
        coarse = lookup[(source_id, "coarse")]["rows"]
        fine = lookup[(source_id, "fine")]["rows"]
        for c, f in zip(coarse, fine):
            require(c["eta"] == f["eta"], "Observation mismatch")
            row = {key: f[key] for key in ("source", "eta", "a", "s_over_epsilon",
                                          "wronskian_scaled", "physical_wronskian_scaled")}
            row["finite_k"] = []
            for ck, fk in zip(c["finite_k"], f["finite_k"]):
                require(ck["K"] == fk["K"], "Cutoff mismatch")
                difference = abs(fk["y"] - ck["y"])
                maximum_difference = max(maximum_difference, difference)
                row["finite_k"].append({"K": fk["K"], "y": fk["y"],
                                        "y_coarse": ck["y"],
                                        "refinement_difference": difference,
                                        "estimated_numerical_error": 2.0 * difference})
            rows.append(row)
    return rows, maximum_difference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-commit", required=True,
                        help="Root's public prospective freeze commit, before execution")
    parser.add_argument("--manifest-sha256", required=True,
                        help="Independent manifest SHA-256 pinned in the public registration")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "outputs")
    args = parser.parse_args()
    require(len(args.freeze_commit) in (40, 64)
            and all(c in "0123456789abcdef" for c in args.freeze_commit),
            "Provide the actual public freeze commit identifier")
    require(len(args.manifest_sha256) == 64
            and all(c in "0123456789abcdef" for c in args.manifest_sha256),
            "Provide the independently pinned manifest SHA-256")
    directory = Path(__file__).resolve().parent
    provenance = verify_manifest(directory, args.manifest_sha256)
    require(not args.output_dir.exists(), "Refusing to overwrite an existing output directory")
    args.output_dir.mkdir(parents=True)
    started = time.monotonic()
    report = {"schema_version": 1, "route": "independent_forced_canonical_modes",
              "freeze_commit": args.freeze_commit, "provenance": provenance,
              "runtime": {"python": sys.version, "numpy": np.__version__,
                          "platform": platform.platform()},
              "normalization": "y=a^2 deltaQ/epsilon", "epsilon": EPSILON,
              "status": "running", "runs": []}
    active_arrays = {}
    result_path = args.output_dir / "results.json"
    try:
        for source_id in SOURCES:
            for setting in SETTINGS:
                progress = {"source": source_id, "setting": setting[0], "status": "running"}
                report["active_run"] = progress
                result_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
                run = evolve(source_id, setting, args.output_dir, started, progress, active_arrays)
                report["runs"].append(run)
                report["active_run"] = None
                active_arrays.clear()
                report["resources"] = check_resources(started)
                result_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
                print(json.dumps({"source": source_id, "setting": setting[0],
                                  "state": "completed", "resources": report["resources"]}), flush=True)
        rows, maximum_difference = combine_results(report["runs"])
        report["rows"] = rows
        report["maximum_refinement_difference"] = maximum_difference
        require(maximum_difference <= CONVERGENCE_GATE,
                "Frozen finite-K refinement gate failed; preserve this diagnostic")
        report["status"] = "passed_internal_gates"
        report["limitations"] = [
            "Refinement differences are numerical estimates, not certified integration bounds.",
            "Independent primary-route finite-K agreement and actual omitted-tail bounds are separate gates.",
            "Fixed geometry scalar variance only: no stress, stability, heating, or particle yield."]
    except Exception as exc:
        report["status"] = "failed"
        report["exception"] = {"type": type(exc).__name__, "message": str(exc),
                               "traceback": traceback.format_exc()}
        if active_arrays:
            snapshot_path = args.output_dir / "failure_snapshot.npz"
            try:
                # No copies and no compression: preserve the current state even
                # when the failed gate concerns the memory or wall-time budget.
                np.savez(snapshot_path, **active_arrays)
                report["failure_snapshot"] = {"path": snapshot_path.name,
                                               "sha256": sha256(snapshot_path)}
            except Exception as snapshot_error:
                report["failure_snapshot_error"] = str(snapshot_error)
        raise
    finally:
        report["resources"] = {"elapsed_seconds": time.monotonic() - started,
                               "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        result_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": report["status"], "result": str(result_path),
                      "maximum_refinement_difference": maximum_difference}), flush=True)


if __name__ == "__main__":
    main()
