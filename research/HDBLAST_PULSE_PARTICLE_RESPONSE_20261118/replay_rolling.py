"""Verify the matched rolling B1 event and the finite-band mode comparison."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

CAPSULE = Path(__file__).resolve().parent
REPO = CAPSULE.parents[1]
CLASSICAL = CAPSULE / "classical"
RUNS = CLASSICAL / "runs"
GRIDS = ("coarse", "fine")
BRANCHES = ("quiet", "0.001", "-0.001")
CUTS = ("1.95", "2.02", "2.1", "2.15")


def read(path):
    return json.loads(path.read_bytes())


def require(ok, label):
    if not ok:
        raise ValueError(label)


def close(actual, expected, label, rel=3e-6, absolute=1e-11):
    require(math.isfinite(actual) and math.isfinite(expected)
            and math.isclose(actual, expected, rel_tol=rel, abs_tol=absolute),
            f"{label}: got {actual}, expected {expected}")


def invoke(script, *args, stdout=None):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    opts = ["-O"] if sys.flags.optimize else []
    subprocess.run([sys.executable, "-I", "-B", *opts, script, *args],
                   stdout=stdout, check=True, env=env)


def compare_history(reference, new):
    import numpy as np
    with np.load(reference / "rolling_timeseries.npz", allow_pickle=False) as stored:
        with np.load(new / "rolling_timeseries.npz", allow_pickle=False) as generated:
            expected = stored["rec"]
            cols = list(stored["cols"])
            actual = generated["rec"][:, [list(generated["cols"]).index(k) for k in cols]]
    require(expected.shape == actual.shape, "Changed history shape")
    for j, key in enumerate(cols):
        label = str(new) + "/" + str(key)
        x, y = actual[:, j], expected[:, j]
        require(np.array_equal(np.isnan(x), np.isnan(y)), label+" missing-value pattern")
        use = np.isfinite(y)
        require(np.all(np.isfinite(x[use])), label+" has nonfinite samples")
        if key == "Hmax_z":
            # The location of a discrete maximum can change cells when two sites tie.
            # Its norm and near-shell norm are compared in the corresponding columns.
            require(np.all((x[use] >= -16.01) & (x[use] <= 0.01)), label)
        else:
            relative = 5e-3 if key in ("Hmax","Hmax_near","Hshell","Mmax","Mmax_near","Mshell") else 5e-6
            absolute = 3e-9 if relative == 5e-3 else 4e-10
            np.testing.assert_allclose(x[use], y[use], rtol=relative, atol=absolute, err_msg=label)
    actual = read(new / "rolling_summary.json")
    expected = read(reference / "rolling_summary.json")
    require(actual["stop_reason"] == "T_final" and len(actual["matter"]["events"]) == 1,
            "Incomplete RK4 branch or wrong number of events")
    event = actual["matter"]["events"][0]
    old = expected["matter"]["events"][0]
    for key in ("T_star","s_star","A_star","q","v_star","D_over_q","G_min_screen",
                "number_factor","n_hat_at_crossing"):
        close(event[key], old[key], key, rel=2e-5, absolute=1e-9)
    require(event["screen_pass"] is old["screen_pass"], "Local screen changed")
    if reference.name != "quiet":
        initial = read(new/"initial.json")
        ref_initial = read(reference/"initial.json")
        for key in ("delta_phi_compensator","matched_delta_y_left","min_abs_D","min_abs_Q",
                    "H_relative_delta_max","M_relative_delta_max"):
            close(initial[key], ref_initial[key], key, rel=5e-3, absolute=3e-9)


def bulk_grid(work, grid):
    import numpy as np
    driver = CLASSICAL / "run_rolling.py"
    invoke(driver, REPO, work, "--grid", grid, "--branch=capture")
    with np.load(RUNS/grid/"checkpoint.npz", allow_pickle=False) as stored:
        with np.load(work/grid/"checkpoint.npz", allow_pickle=False) as new:
            np.testing.assert_allclose(new["hist"], stored["hist"], rtol=5e-7, atol=3e-10)
    for branch in BRANCHES:
        invoke(driver, REPO, work, "--grid", grid, "--branch="+branch)
        compare_history(RUNS/grid/branch, work/grid/branch)
        print(f"Reproduced rolling RK4 shell history, event and constraints for {grid}/{branch}", flush=True)


def bulk(work):
    with ThreadPoolExecutor(max_workers=2) as pool:
        for result in [pool.submit(bulk_grid, work, grid) for grid in GRIDS]:
            result.result()




def check_receipts(label, grid, params, rolling):
    start, spacing, rtol, order = params
    reports = {}
    for branch in BRANCHES:
        d = read(CLASSICAL/"modes"/f"{label}-{grid}-{branch}.json")
        event = rolling[grid][branch]
        require(d["grid"] == grid and d["branch"] == branch and d["basis_order"] == order
                and (d["start"], d["k_step"], d["rtol"]) == (start, spacing, rtol),
                "Mixed mesh, initial state, WKB order or quadrature settings")
        require(tuple(d["ends"]) == CUTS and 0 < d["start_adiabaticity"] < 0.012
                and d["max_lna_drift"] < 5e-8, "Rolling trajectory control")
        close(d["M5_over_H0"], 6.49e-6**(-1/3), "Gravitational scale")
        # Analytic integral of k² exp[-pi(k/a*)²/q] / (2 pi²), from zero to 8.
        u = math.pi / (d["event_a_star"]**2 * event["q"])
        leading = (math.sqrt(math.pi)*math.erf(8*math.sqrt(u))/(4*u**1.5)
                   - 4*math.exp(-64*u)/u)/(2*math.pi**2)
        close(d["event_leading_comoving_N_8"], leading, "Same-band Gaussian integral", rel=3e-7)
        close(d["event_comoving_N_8"], leading*event["derivative_factor"],
              "Same-band derivative-corrected integral", rel=3e-7)
        for cut, row in d["ends"].items():
            require(row["band_E_upper_H0"] < d["M5_over_H0"] and row["wronskian_max"] < 6e-8,
                    "Frequency screen or normalization failed")
            close(row["N8_phys_density"]*row["a"]**3, row["comoving_N_8"], "Comoving definition")
            close(row["N8_half_grid"], row["comoving_N_8"], "Momentum quadrature", rel=8e-6)
        reports[branch] = d

    quiet, plus, minus = (reports[k] for k in BRANCHES)
    require(minus["event_comoving_N_8"] < quiet["event_comoving_N_8"] < plus["event_comoving_N_8"],
            "Corrected event ordering")
    require(plus["event_leading_comoving_N_8"] < quiet["event_leading_comoving_N_8"]
            < minus["event_leading_comoving_N_8"], "Leading event ordering")
    for cut in CUTS:
        q, p, m = [r["ends"][cut]["comoving_N_8"] for r in (quiet, plus, minus)]
        require(-.0022 < p/q-1 < -.0019 and .0096 < m/q-1 < .0102, "Direct free-mode contrasts")


def quick():
    reports = read(CLASSICAL/"rolling_results.json")
    leading = read(CLASSICAL/"leading_coupled.json")
    baseline = (.8, .25, 2e-9, 2)
    for grid in GRIDS:
        for branch in BRANCHES:
            d = reports[grid][branch]
            src = read(RUNS/grid/branch/"rolling_summary.json")
            require(src["stop_reason"] == "T_final" and len(src["matter"]["events"]) == 1,
                    "Wrong event count or incomplete branch")
            event = src["matter"]["events"][0]
            require(src["matter"]["G"] == 100 and src["matter"]["phistar"] == .5
                    and src["matter"]["b"] == 6.49e-6 and src["matter"]["y"] == 1.
                    and src["matter"]["mu0"] == 0. and src["matter"]["apply_D"] is True,
                    "Original B1 couplings")
            require(1.210 < d["T_star"] < 1.255 and d["passes_screen"] is (branch != "0.001"),
                    "Arrival causal window or local validity")
            for key in ("T_star","s_star","q","v_star","D_over_q"):
                close(event[key], d[key], key)
            close(100*abs(d["v_star"]), d["q"], "Original event q")
            close(d["derivative_factor"], 1 + d["D_over_q"], "Corrected event factor")
            limit = 2e-7 if grid == "coarse" else 5e-12
            require(d["pre_arrival_max_abs_v_difference"] < limit, "Pre-arrival velocity leakage")
            if branch != "quiet":
                cons = d["construction"]
                require(abs(cons["matched_delta_y_left"]) < 1e-13 and cons["amplitude"] == float(branch)
                        and cons["H_relative_delta_max"] < 2e-5 and cons["M_relative_delta_max"] < 9e-5,
                        "Constraint-matched initial data")
            modified = read(RUNS/grid/("leading-"+branch)/"rolling_summary.json")
            require(modified["matter"]["apply_D"] is False
                    and modified["matter"]["events"][0]["number_factor"] == 1,
                    "Leading-only coupled model setting")
        check_receipts("base", grid, baseline, reports)
        for t in ("2.1", "2.15"):
            q, p, m = (leading[grid][branch][t]["R"] for branch in BRANCHES)
            require(p < q < m, "Leading-only coupled radiation contrast")
    for label, params in (("basis", (.8,.25,2e-9,0)), ("start", (.5,.25,2e-9,2)),
                          ("strict", (.8,.125,2e-10,2))):
        check_receipts(label, "fine", params, reports)
    print("Verified 15 rolling mode reports, same-band event integrals and matched controls", flush=True)


def full(work):
    import numpy as np
    expected_summary = read(CLASSICAL/"rolling_results.json")
    summary_file = work/"summary.json"
    invoke(CLASSICAL/"summarize_rolling.py", RUNS, summary_file)
    summary = read(summary_file)
    leading = read(CLASSICAL/"leading_coupled.json")
    for grid in GRIDS:
        for branch in BRANCHES:
            for key in ("T_star","s_star","v_star","derivative_factor","n_target","R_at_s2_1",
                        "Qhat_at_s2_1","proper_nhat_at_s2_1"):
                close(summary[grid][branch][key], expected_summary[grid][branch][key], key)
            npz = RUNS/grid/("leading-"+branch)/"rolling_timeseries.npz"
            with np.load(npz, allow_pickle=False) as archive:
                cols = dict(zip(archive["cols"], archive["rec"].T))
            for cut, stored in leading[grid][branch].items():
                for k,value in stored.items():
                    close(float(np.interp(float(cut),cols["H0tau"],cols[k])), value, k)
            out = work/f"{grid}-{branch}.json"
            invoke(CLASSICAL/"rolling_modes.py", RUNS, out, "--grid", grid, "--branch="+branch)
            actual = read(out)
            stored = read(CLASSICAL/"modes"/f"base-{grid}-{branch}.json")
            for k in ("event_a_star","event_comoving_N_8","event_leading_comoving_N_8","start_adiabaticity"):
                close(actual[k],stored[k],k)
            for t in CUTS:
                row = actual["ends"][t]
                for k in ("comoving_N_8","N8_half_grid","N8_phys_density","N8_rho_quasiparticle",
                          "band_E_upper_H0","a","mass_H0","endpoint_adiabaticity"):
                    close(row[k],stored["ends"][t][k],k)
                require(row["wronskian_max"] < 6e-8, "Replayed mode normalization")
            print(f"Reproduced same-band rolling free modes for {grid}/{branch}", flush=True)
    audit = CAPSULE/"rolling_audit"
    out = work/"energy.json"
    with out.open("wb") as stream:
        invoke(audit/"energy_budget.py", RUNS, stdout=stream)
    actual, original = read(out), read(audit/"ENERGY.json")
    for grid in GRIDS:
        for branch in BRANCHES:
            case = actual[grid]["branches"][branch]
            require(max(map(abs, case["checks"].values())) < .002, "Independent matter budget")
            for k in ("q","s_star","D_over_q","number_factor","n_hat_at_crossing"):
                close(case["event"][k],original[grid]["branches"][branch]["event"][k],k)
        for branch in BRANCHES[1:]:
            for k in ("recorded_R_pct","recorded_X_pct","transport_R_pct","transport_R_no_D_pct",
                      "exact_fixed_wall_R_no_D_pct"):
                close(actual[grid]["contrasts"][branch]["cuts"]["2.15"][k],
                      original[grid]["contrasts"][branch]["cuts"]["2.15"][k],k)
    print("Reproduced rolling timing, leading-only coupled radiation and independent transport budget", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="rerun free modes and radiation transport (NumPy, SciPy)")
    parser.add_argument("--bulk", action="store_true", help="rerun the six original-prescription 5D bulk evolutions from authenticated sources (NumPy, SciPy)")
    args = parser.parse_args()
    quick()
    if args.full or args.bulk:
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            if args.full:
                full(work)
            if args.bulk:
                bulk(work)


if __name__ == "__main__":
    main()
