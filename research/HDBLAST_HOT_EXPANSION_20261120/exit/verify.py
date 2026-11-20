"""Check the hot interval and its energy source from the recorded samples."""
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def crossing(left, right, left_value, right_value):
    weight = -left_value / (right_value-left_value)
    return {k: left[k]+weight*(right[k]-left[k]) for k in left}


def close(actual, expected, tolerance, label):
    if abs(actual-expected) > tolerance:
        raise ValueError(f"{label}: got {actual:.12g}, expected {expected:.12g} ± {tolerance:.2g}")


def verify(data):
    run = data["runs"]["primary"]
    ps = data["parameters"]["primary"]
    columns = data["trace_columns"]
    rows = [dict(zip(columns, vals)) for vals in run["trace"]]
    if data["solver_sha256"] != hashlib.sha256((HERE/"growth_forward.cpp").read_bytes()).hexdigest():
        raise ValueError("Solver hash mismatch")
    # Independently reconstruct the rest-frame EOS random energy carried by
    # each laboratory particle, then check the conservative gas ledger.
    max_intrinsic_error = 0.
    for r in rows:
        theta = r["theta"]
        internal = 1.5*theta + 2.25*theta*theta/(1+math.hypot(1,1.5*theta))
        close(r["intrinsic_random"], r["rest"]*internal, 2e-10, "thermal EOS")
        gamma = math.hypot(1,r["u"])
        reference = r["rest"]*((1+internal+theta)*gamma-theta/gamma)
        close(r["gas"], reference, 2e-10, "gas energy")
        max_intrinsic_error = max(max_intrinsic_error, abs(r["intrinsic_random"]-r["rest"]*internal))
    # Locate crossings directly in exported states, independently of the native
    # every-time-step event counters. Linear interpolation leaves sampling error.
    def margin(x):
        return x["intrinsic_random"]-.9*x["e"]
    exit_field = ps["exit_field"] + ps["exit_width"]
    i = next(j for j,r in enumerate(rows) if r["psi"]>=exit_field and r["H"]>0 and margin(r)>=0)
    hot = crossing(rows[i-1], rows[i], margin(rows[i-1]), margin(rows[i]))
    j = next(j for j in range(i+1,len(rows)) if rows[j]["H"]<=0)
    turnaround = crossing(rows[j-1], rows[j], rows[j-1]["H"], rows[j]["H"])
    if any(r["V"] != 0 or r["rad"] != 0 or margin(r)<0 for r in rows[i:j]):
        raise ValueError("Sampled hot interval was interrupted")
    target = 2*run["first_ninety_intrinsic_radius"]
    k = next(k for k in range(i+1,j) if rows[k]["a"]>=target)
    doubled = crossing(rows[k-1], rows[k], rows[k-1]["a"]-target, rows[k]["a"]-target)
    bulk_term = doubled["mass"]/doubled["a"]**4
    positive_friedmann = doubled["e"]*(2+doubled["e"])+bulk_term
    growth = math.log(turnaround["a"]/hot["a"])
    close(growth, run["ninety_intrinsic_expanding_efolds"], .008, "hot growth")
    close(doubled["intrinsic_random"]/doubled["e"], run["doubling_intrinsic_fraction"], 2e-5, "doubling share")
    close(doubled["t"], run["first_hot_doubling_tau"], 1e-3, "doubling time")
    close(run["peak_gas_total_a3"]+run["peak_gas_pressure_work"],
          run["peak_daughter_work_a3"]+run["peak_decay_energy_a3"], 2e-7, "particle ledger")
    close(run["end_comoving_sheet_energy"]+run["end_charge_work_a4"],
          run["comoving_wave_to_sheet"]+run["comoving_scalar_to_sheet"], 2e-6, "Maxwell ledger")
    if run["ninety_intrinsic_intervals"]!=1 or growth <= math.log(2):
        raise ValueError("No unbroken radius doubling")
    return {"solver_sha256":data["solver_sha256"], "maximum_thermal_reconstruction_error":max_intrinsic_error,
            "sample_interpolated_hot_efolds":growth,"sample_interpolated_hot_time":hot["t"],
            "sample_interpolated_turnaround_time":turnaround["t"],"sample_interpolated_doubling_time":doubled["t"],
            "sample_interpolated_doubling_share":doubled["intrinsic_random"]/doubled["e"],
            "doubling_bulk_gravity_fraction":bulk_term/positive_friedmann,
            "doubling_total_wall_pressure_over_density":doubled["p_total"]/doubled["e"],
            "work_to_direct_decay_energy":run["peak_daughter_work_a3"]/run["peak_decay_energy_a3"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, nargs="?", default=HERE/"EXIT.json")
    print(json.dumps(verify(json.loads(parser.parse_args().data.read_text())), indent=2))
