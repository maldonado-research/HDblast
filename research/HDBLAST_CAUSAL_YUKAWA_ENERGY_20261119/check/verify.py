"""Cross-check the causal budget against two independent calculations."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE = Path(__file__).resolve().parent


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def build_singleton(native, solver):
    q, weights, initial, background = solver.read_native(native, "fine")
    j = 11
    return solver.run_budget(q[j:j+1], np.ones(1),
                              {key: tuple(row[j:j+1] for row in value) for key, value in initial.items()},
                              background, 60., order=24, steps=4096)


def verify(spectra, finite, singleton, budgets, refinement):
    comparisons = {}
    for section, window in (("grids", "wide"), ("narrow_grids", "narrow")):
        for grid, cuts in spectra[section].items():
            for cutoff, run in cuts.items():
                key = f"{grid}_{window}_{cutoff}"
                if key not in budgets["runs"]:
                    continue
                for state, independent in run["states"].items():
                    primary = budgets["runs"][key][state]["bare"]["checkpoints"]["cut"]["energy"]["daughter"]
                    relative = (primary-independent["energy"])/independent["energy"]
                    comparisons[f"{key}_{state}"] = dict(primary=primary, independent=independent["energy"],
                                                          relative_difference=relative)
                    if abs(relative) > 1e-7:
                        raise ValueError(("Independent daughter-energy mismatch", key, state, relative))
    if len(comparisons) != 30:
        raise ValueError(("Native spectral coverage", len(comparisons)))
    finite_checks = {}
    for state, schemes in finite["runs"].items():
        for scheme, runs in schemes.items():
            coarse, fine = runs
            ratio = (fine["proper_time_steps"]/coarse["proper_time_steps"])**2
            extrapolate = lambda section, key: (ratio*fine[section][key]-coarse[section][key])/(ratio-1)
            reference = singleton[state][scheme]["checkpoints"]["cut"]
            endpoint_difference = {key: extrapolate("order_y2_endpoint", key)-
                                  reference["energy"]["daughter" if key=="bath" else key]
                                  for key in ("parent", "bath")}
            work_reference = dict(parent_background=reference["work"]["mass"]+reference["work"]["scale_parent"],
                                  bath_redshift=reference["work"]["scale_daughter"],
                                  switch=reference["work"]["switch"],
                                  counterterm=reference["work"]["scale_counterterm"]+
                                              reference["work"]["switch_counterterm"])
            work_difference = {key: extrapolate("order_y2_work", key)-value
                               for key, value in work_reference.items()}
            spectral_difference = extrapolate("order_y2_endpoint", "bath")-fine["spectrum_bath_energy"]
            error_ratio = coarse["series_minus_spectrum"]/fine["series_minus_spectrum"]
            balance_extrapolated = (ratio*fine["energy_balance_residual"]-coarse["energy_balance_residual"])/(ratio-1)
            weak_y_error = max(abs(value) for value in fine["extrapolation_minus_series"].values())
            test = dict(endpoint_difference=endpoint_difference, work_difference=work_difference,
                        daughter_difference_from_direct_spectrum=spectral_difference,
                        step_halving_spectral_error_ratio=error_ratio,
                        extrapolated_energy_balance_residual=balance_extrapolated,
                        finite_y_extrapolation_max_difference=weak_y_error,
                        max_symplectic_residual=max(run["max_symplectic_residual"] for run in runs),
                        max_order_y2_commutator_residual=max(run["order_y2_commutator_residual"] for run in runs))
            if (max(map(abs, endpoint_difference.values())) > 1e-7 or
                    max(map(abs, work_difference.values())) > 1e-7 or
                    abs(spectral_difference) > 1e-8 or abs(error_ratio-4) > .02 or
                    abs(balance_extrapolated) > 1e-7 or weak_y_error > 2e-7 or
                    test["max_symplectic_residual"] > 1e-10 or
                    test["max_order_y2_commutator_residual"] > 1e-10):
                raise ValueError(("Finite bath check", state, scheme, test))
            finite_checks[f"{state}_{scheme}"] = test
    refinement_difference = {state: (refinement["states"][state]["energy"]-
                                     spectra["grids"]["fine"]["120"]["states"][state]["energy"])
                             for state in refinement["states"]}
    if max(map(abs, refinement_difference.values())) > 1e-8:
        raise ValueError(("Independent spectrum refinement", refinement_difference))
    return dict(passed=True, native_energy_comparisons=len(comparisons),
                max_native_energy_relative_difference=max(abs(row["relative_difference"]) for row in comparisons.values()),
                native_spectral_comparisons=comparisons, finite_bath_checks=finite_checks,
                spectrum_refinement_energy_difference=refinement_difference,
                interpretation="Finite-y exactness refers to the quadratic oscillator-bath realization; only its order-y² kernel is identified with the Yukawa pair channel.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native", required=True, type=Path)
    ap.add_argument("--solver", type=Path, default=HERE.parent/"solver/budget.py")
    ap.add_argument("--budgets", type=Path, default=HERE.parent/"solver/BUDGETS.json")
    ap.add_argument("--output", type=Path, default=HERE/"CHECKS.json")
    ap.add_argument("--full", action="store_true", help="Recompute both independent methods and the singleton solver")
    args = ap.parse_args()
    if args.full:
        subprocess.run([sys.executable, str(HERE/"native_spectrum.py"), "--native", str(args.native),
                        "--output", str(HERE/"SPECTRUM.json")], check=True)
        subprocess.run([sys.executable, str(HERE/"finite_bath.py"), "--native", str(args.native),
                        "--output", str(HERE/"FINITE_BATH.json")], check=True)
    if args.full or not (HERE/"SINGLETON_SOLVER.json").exists():
        singleton = build_singleton(args.native, module(args.solver, "causal_budget"))
        (HERE/"SINGLETON_SOLVER.json").write_text(json.dumps(singleton, indent=2, allow_nan=False)+"\n")
    if args.full or not (HERE/"SPECTRUM_REFINED.json").exists():
        reference = module(HERE/"native_spectrum.py", "direct_spectral_reference")
        refined = reference.spectrum(*reference.read_native(args.native), cutoff=120., order=128, samples=8192)
        (HERE/"SPECTRUM_REFINED.json").write_text(json.dumps(refined, indent=2, allow_nan=False)+"\n")
    read = lambda path: json.loads(path.read_text())
    report = verify(read(HERE/"SPECTRUM.json"), read(HERE/"FINITE_BATH.json"),
                    read(HERE/"SINGLETON_SOLVER.json"), read(args.budgets), read(HERE/"SPECTRUM_REFINED.json"))
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(f"Passed {report['native_energy_comparisons']} native energy comparisons and six exact finite-bath controls; "
          f"max native relative difference {report['max_native_energy_relative_difference']:.3g}.")


if __name__ == "__main__":
    main()
