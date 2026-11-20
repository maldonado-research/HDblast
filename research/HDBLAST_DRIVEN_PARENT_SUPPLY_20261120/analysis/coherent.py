"""Homogeneous coherent action, local pressure, and complete-cycle dust tests."""
import argparse
from bisect import bisect_left
import json
import math
from pathlib import Path


# Four-point Gauss-Legendre per saved time interval.
QUAD = [
    (-0.8611363115940526, 0.34785484513745385),
    (-0.33998104358485626, 0.6521451548625461),
    (0.33998104358485626, 0.6521451548625461),
    (0.8611363115940526, 0.34785484513745385),
]


def sample_scalar(left, right, t, key, slope):
    dt = right["t"] - left["t"]
    x = (t - left["t"]) / dt
    y0, y1 = left[key], right[key]
    d0, d1 = slope(left) * dt, slope(right) * dt
    b = 3 * (y1 - y0) - 2 * d0 - d1
    c = 2 * (y0 - y1) + d0 + d1
    return y0 + x * (d0 + x * (b + x * c)), (d0 + x * (2 * b + x * 3 * c)) / dt


class FreeEpoch:
    """Free, constant-mass subinterval of the numerical trajectory."""

    def __init__(self, rows, start, end, mass):
        self.rows = rows
        self.times = [r["t"] for r in rows]
        self.start = start
        self.end = end
        self.mass = mass
        self.knots = [start, *(t for t in self.times if start < t < end), end]
        if start < self.times[0] or end > self.times[-1] or end <= start:
            raise ValueError("Free epoch outside trace")

    def state(self, t):
        i = max(1, min(len(self.rows) - 1, bisect_left(self.times, t)))
        l, r = self.rows[i-1:i+1]
        a, adot = sample_scalar(l, r, t, "a", lambda x: x["a"] * x["H"])
        chi, vel = sample_scalar(l, r, t, "chi", lambda x: x["chi_dot"])
        e = (vel ** 2 + (self.mass * chi) ** 2) / 2
        p = (vel ** 2 - (self.mass * chi) ** 2) / 2
        action = a ** 3 * e / self.mass
        corr = 1.5 * adot / (a * self.mass) * a ** 3 * chi * vel
        return dict(a=a, H=adot/a, chi=chi, vel=vel, e=e, p=p, I=action, J=action+corr)

    def turning_points(self):
        # All resolved turning points are separated by many saved steps.
        roots = []
        for l, r in zip(self.knots, self.knots[1:]):
            left, right = l, r
            vl, vr = self.state(left)["vel"], self.state(right)["vel"]
            if vl == 0:
                roots.append(left)
            if vl * vr < 0:
                for _ in range(46):
                    m = (left + right) / 2.
                    if self.state(m)["vel"] * vl <= 0:
                        right = m
                    else:
                        left = m
                roots.append((left+right)/2)
        return roots

    def integrated(self, start, end):
        # Measures are proper time, proper time times volume, and its
        # expansion-work weight. Radius and field use cubic Hermite data.
        keys = ("I", "p", "e", "volume_p", "volume_e", "work", "work_weight", "H")
        sums = dict.fromkeys(keys, 0.)
        knots = [start, *(t for t in self.knots if start < t < end), end]
        for l, r in zip(knots, knots[1:]):
            mid, half = (l + r) / 2, (r - l) / 2
            for x, weight in QUAD:
                u = self.state(mid + x * half)
                v = u["a"] ** 3
                vals = (u["I"], u["p"], u["e"], v*u["p"], v*u["e"],
                        v*u["H"]*u["p"], v*u["H"]*u["e"], u["H"])
                for k, val in zip(keys, vals):
                    sums[k] += val * weight * half
        return sums


def whole_sphere():
    # The source-only model has zero charge. Both examples are conversions
    # conditional on the local coupling of a putative charged daughter model.
    mu = 1.
    target = 8 * 1e-4 / (2 * mu * .95)
    per_action = 2 * math.pi ** 2 * target
    examples = {}
    for label, charge, gauge in (("previous", 1e8, 13200.04), ("two_stage", 2e8, 26.04)):
        s0 = 3 * charge ** 2 * mu ** 2 / (4 * math.pi * gauge)
        examples[label] = dict(Q=charge, K=gauge, S0_at_alpha_one=s0, N_target_at_alpha_one=per_action*s0)
    return dict(target_action=target, N_target_over_S0=per_action,
                illustrative_local_couplings=examples)


def epoch_report(epoch, target):
    rows = [epoch.state(t) for t in epoch.knots]
    actions = [r["I"] / target for r in rows]
    corrected = [r["J"] / target for r in rows]
    full = epoch.integrated(epoch.start, epoch.end)
    extrema = epoch.turning_points()
    periods = (len(extrema)-1)//2 if extrema else 0
    result = dict(
        start=epoch.start, end=epoch.end, mass=epoch.mass, sampled_points=len(rows),
        mass_phase_in_cycles=epoch.mass * (epoch.end-epoch.start)/(2*math.pi),
        action_over_target_min=min(actions), action_over_target_max=max(actions),
        action_over_target_mean=full["I"] / (epoch.end-epoch.start) / target,
        action_over_target_start=actions[0], action_over_target_end=actions[-1],
        corrected_action_over_target_min=min(corrected), corrected_action_over_target_max=max(corrected),
        max_abs_H_over_M=max(abs(r["H"]) for r in rows) / epoch.mass,
        total_volume_weighted_pressure_over_energy=full["volume_p"]/full["volume_e"],
        turning_points=extrema, complete_periods=periods,
    )
    if periods:
        start, end = extrema[0], extrema[2*periods]
        integrals = epoch.integrated(start, end)
        is_, ie = epoch.state(start)["I"], epoch.state(end)["I"]
        cycles = dict(start=start, end=end, action_over_target_start=is_/target,
                      action_over_target_end=ie/target, fractional_action_change=ie/is_-1,
                      proper_time_pressure_over_energy=integrals["p"]/integrals["e"],
                      comoving_pressure_over_energy=integrals["volume_p"]/integrals["volume_e"],
                      expansion_weighted_pressure_over_energy=integrals["work"]/integrals["work_weight"],
                      action_identity_residual_over_target=(ie-is_+3*integrals["work"]/epoch.mass)/target)
        result["all_complete_cycles"] = cycles
    return result


def late_native(data, name):
    r = data["runs"][name]
    p = data["parameters"][name]
    target = 8e-4 / (2 * p["daughter_mass"] * p["initial_mass_ratio"])
    intervals = [dict(zip(data["half_cycle_columns"], v)) for v in r["half_cycles"]]
    trace = [dict(zip(data["trace_columns"], v)) for v in r["trace"]]
    turn = r["maximum_tau"]
    left = next(s for s in reversed(trace) if s["t"] < turn)
    right = next(s for s in trace if s["t"] >= turn)
    # On each side of the first maximum a is monotone. Constant mass, no
    # decay and |p|<=e give |ln I(turn)/I(row)|<=3|ln a(turn)/a(row)|.
    brackets = []
    for row in (left, right):
        dilation = (r["maximum_radius"] / row["a"]) ** 3
        ratio = row["a"] ** 3 * row["rho_chi"] / row["M"] / target
        brackets.append((ratio / dilation, ratio * dilation))
    lo, hi = max(pair[0] for pair in brackets), min(pair[1] for pair in brackets)
    if lo > hi:
        raise ArithmeticError("Turnaround action brackets do not overlap")
    return dict(
        daughter_mass=p["daughter_mass"], drive_amplitude=p["drive_amplitude"],
        drive_width=p["drive_width"], opening_pressure_over_energy=r["opening_pressure"]/r["opening_density"],
        opening_mass_WKB=r["opening_mass_WKB"], opening_H_over_M=r["opening_H_over_M"],
        opening_action_over_target=r["opening_target_ratio"],
        mass_complete_action_over_target=r["mass_complete_target_ratio"],
        complete_half_cycles=len(intervals), complete_periods=len(intervals)//2,
        min_half_cycle_end_action_over_target=min(z["action"] for z in intervals)/target,
        max_half_cycle_end_action_over_target=max(z["action"] for z in intervals)/target,
        max_abs_half_cycle_time_weighted_pressure_ratio=max(abs(z["mean_w_time"]) for z in intervals),
        max_abs_half_cycle_volume_weighted_pressure_ratio=max(abs(z["mean_w_volume"]) for z in intervals),
        turn_radius=r["maximum_radius"], turn_time=turn, times_bracketing_turn=[left["t"],right["t"]],
        turnaround_action_ratio_lower=lo, turnaround_action_ratio_upper=hi,
    )


def analyze(data):
    ref = whole_sphere()
    target = ref["target_action"]
    r = data["runs"]["primary"]
    tr = [dict(zip(data["trace_columns"], row)) for row in r["trace"]]
    mass = r["drive_end_mass"]
    before = FreeEpoch(tr, r["drive_end_time"], r["before_mass_time"], mass)
    opening = dict(pressure_over_energy=r["opening_pressure"]/r["opening_density"],
                   H_over_M=r["opening_H_over_M"], mass_sweep_over_M2=r["opening_mass_WKB"],
                   coherent_target_ratio=r["opening_action"]/target,
                   N_coherent_over_S0=2*math.pi**2*r["opening_action"])
    opening["instantaneous_mixing_over_frequency"] = (
        3* r["opening_H"]/r["opening_mass"] + r["opening_mass_WKB"])/2
    end = tr[-1]
    comoving_e = end["a"]**3*end["rho_chi"]
    ledger = dict(t=end["t"], comoving_energy=comoving_e, work_from_moving_minimum=end["drive_work"],
                  work_from_mass=end["mass_work"], pressure_work_spent=end["chi_pressure_work"],
                  residual=comoving_e + end["chi_pressure_work"]-end["mass_work"]-end["drive_work"])
    result = dict(whole_sphere=ref, opening=opening, end_of_primary_run=ledger,
                  source_off_before_gate=epoch_report(before, target))
    if "stable" in data["runs"]:
        stable = data["runs"]["stable"]
        sr = [dict(zip(data["trace_columns"], row)) for row in stable["trace"]]
        upper_mass = max(s["M"] for s in sr)
        plateau = next(s for s in sr if s["t"] > r["opening_time"] and s["M"] == upper_mass)
        last = max(s["t"] for s in sr if s["t"] <= stable["maximum_tau"])
        after = FreeEpoch(sr, plateau["t"], last, upper_mass)
        result["source_off_after_gate_expanding"] = epoch_report(after,target)
        # A half-spaced companion checks the virial and action integrals.
        coarse = sr[::2]
        if coarse[-1] is not sr[-1]:
            coarse.append(sr[-1])
        check = epoch_report(FreeEpoch(coarse, plateau["t"],last,upper_mass),target)
        result["after_gate_coarsened"] = {key:check[key] for key in
            ("complete_periods","action_over_target_mean","turning_points","all_complete_cycles")
            if key in check}
    result["native_postgate_cycles_and_turnaround"] = {name: late_native(data, name)
        for name in ("stable", "heavy_stable") if name in data["runs"]}
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1]/"coupled/PARENT.json")
    p.add_argument("--output", type=Path)
    cli = p.parse_args()
    result = json.dumps(analyze(json.loads(cli.source.read_bytes())),indent=2,sort_keys=True,allow_nan=False)+"\n"
    if cli.output:
        cli.output.write_text(result)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
