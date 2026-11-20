"""Check scalar decay, charge-limited heat, and finite-time scales on a coupled trajectory."""
import argparse
from bisect import bisect_left
import json
import math
from pathlib import Path


def interpolate(rows, key, value, fields):
    values = [r[key] for r in rows]
    index = bisect_left(values, value)
    if index >= len(rows) or value < values[0]:
        return None
    if values[index] == value:
        return {c: rows[index][c] for c in fields}
    left, right = rows[index - 1], rows[index]
    z = (value - left[key]) / (right[key] - left[key])
    return {c: left[c] + z * (right[c] - left[c]) for c in fields}


def characteristic_scales(p, r):
    mu, kappa = p["daughter_mass"], p["kappa"]
    charge = abs(p["charge_to_mass"])
    R = p["final_mass_ratio"]
    mass, excess = 2 * mu * R, 2 * mu * (R - 1)
    width = kappa * R * (1 - 1 / R**2) ** 1.5
    d = dict(mu=mu, kappa=kappa, absolute_Q=charge, summed_yukawa_alpha=kappa/mu,
             maximum_width=width, plateau_width_over_M=width/mass,
             plateau_width_over_threshold_excess=width/excess, plateau_excess=excess)
    if "opening_mass_WKB" in r:
        d.update(opening_time=r["opening_time"], opening_abs_H_over_M=r["opening_H_over_M"],
                 opening_abs_Mdot_over_M2=r["opening_mass_WKB"])
    if "exit_complete_field_density" in r and charge > 0:
        # The late cap of the base gauge function is 0.04, and this
        # preparation has completed both monotone stages on this slice.
        K = .04 + p["final_gauge_increment"]
        E = math.sqrt(2 * r["exit_complete_field_density"] * r["exit_complete_electric_fraction"] / K)
        if E > 0:
            force = charge * E
            d["electric_at_exit"] = {
                "K": K, "abs_X": E, "electric_share": r["exit_complete_electric_fraction"],
                "abs_eE_over_m2": force/mu,
                "flat_constant_field_exponent": math.pi*mu/force,
                "inverse_force_in_wall_time": 1/force,
                "inverse_sqrt_abs_eE_in_wall_time": 1/math.sqrt(mu*force),
            }
            if "opening_time" in r:
                d["electric_at_exit"]["proper_time_from_exit_to_opening"] = r["opening_time"] - r["exit_complete_tau"]
    return d


def budget_at(r, p, rows, point):
    kappa, charge, mu = p["kappa"], abs(p["charge_to_mass"]), p["daughter_mass"]
    exposure = r[point+"_exposure"]
    actual_work = r[point+"_electric_work"]
    injection = r[point+"_birth_kinetic"]
    d = dict(time=r[point+"_time"], radius=r[point+"_radius"],
             rest=r[point+"_rest"],
             birth_energy=r[point+"_birth_energy"],
             birth_kinetic=injection, field_exposure=exposure,
             electric_work=actual_work, speed_limited_work=charge*exposure,
             spent_pressure_work=r[point+"_gas_pressure_work"],
             integrated_excess_mass_phase=r[point+"_formation_phase"])
    if kappa > 0:
        d.update(A=injection/kappa, B=exposure/kappa)
    if exposure > 0 and charge > 0:
        d.update(absolute_work_over_speed_limit=abs(actual_work)/(charge*exposure),
                 min_Q_alpha_to_reproduce_actual_electric_work=kappa/mu*abs(actual_work)/exposure)
    if rows and kappa > 0:
        row = next((s for s in rows if s["t"] >= r[point+"_time"] and s["H"] >= 0
                    and s["intrinsic"] >= .9 * s["e"]), None)
        if row and row["field_exposure"] > 0:
            # Keep primitives and cumulative budgets on the same actual
            # saved sample. Independently interpolated primitives can break
            # the small gas-energy inequality near a rapidly varying event.
            A_row, B_row = row["birth_kinetic"]/kappa, row["field_exposure"]/kappa
            required = row["intrinsic"] * row["a"]**3 + row["gas_pressure_work"]
            required90 = .9*row["e"]*row["a"]**3 + row["gas_pressure_work"]
            d["saved_hot_state_budget"] = dict(
                time=row["t"], radius=row["a"], actual_internal_fraction=row["intrinsic"]/row["e"],
                A=A_row, B=B_row, required_intrinsic_plus_pressure_work=required,
                upper_from_actual_work=row["birth_kinetic"]+row["electric_work"],
                min_Q_at_same_alpha_heat=max(0., required-row["birth_kinetic"])/row["field_exposure"],
                min_alpha_at_same_Q_heat=required/(mu*(A_row+charge*B_row)),
                min_Q_alpha_for_measured_heat_if_alpha_le_one=max(0., required/mu-A_row)/B_row,
                min_Q_alpha_for_90pct_heat_if_alpha_le_one=max(0., required90/mu-A_row)/B_row,
            )
    return d


def trace_audit(rows, p, r):
    mu, charge = p["daughter_mass"], abs(p["charge_to_mass"])
    pair_balance = birth_balance = gas_balance = scalar_balance = current_balance = 0.
    heat_excess = speed_excess = 0.
    for s in rows:
        volume = s["a"]**3
        pair_balance = max(pair_balance, abs(s["D"]-2*mu*s["pairs"]))
        birth_balance = max(birth_balance, abs(s["D"]+s["birth_kinetic"]-s["decay_energy"]))
        gas_balance = max(gas_balance, abs(volume*s["gas"] + s["gas_pressure_work"]
                                      - s["electric_work"] - s["decay_energy"]))
        scalar_balance = max(scalar_balance, abs(volume*s["rho_chi"] + s["chi_pressure_work"]
                                     + s["decay_energy"] - s["drive_work"] - s["mass_work"]))
        heat_excess = max(heat_excess, volume*s["intrinsic"]+s["gas_pressure_work"] -
                          s["electric_work"]-s["birth_kinetic"])
        speed_excess = max(speed_excess, abs(s["electric_work"])-charge*s["field_exposure"])
        current_balance = max(current_balance, abs(volume*s["current"] -
                                                p["charge_to_mass"]*s["D"]*s["velocity"]))
    result = dict(pair_number_error=pair_balance, birth_energy_error=birth_balance,
                  gas_comoving_energy_error=gas_balance,
                  scalar_comoving_energy_error=scalar_balance, current_identity_error=current_balance,
                  largest_positive_heat_excess=heat_excess, largest_positive_speed_limit_excess=speed_excess)
    if "opening_time" not in r:
        return result
    deadline = r.get("first_hot_time", r.get("maximum_tau",rows[-1]["t"]))
    result["phase_reference_time"] = deadline
    reference = interpolate(rows, "t", deadline, ("D","formation_phase"))
    if not reference or reference["D"]<=0:
        return result
    result["formed_rest_at_reference"] = reference["D"]
    result["reference_excess_mass_phase"] = reference["formation_phase"]
    for mark, name in ((1.,"one_radian"),(2*math.pi,"one_cycle")):
        if reference["formation_phase"] >= mark:
            sample=interpolate(rows,"formation_phase",mark,("t","D"))
            result[name] = dict(time=sample["t"],fraction_of_reference_rest_already_born=sample["D"]/reference["D"])
        else:
            result[name]=dict(time=None,fraction_of_reference_rest_already_born=1.)
    selected=[s for s in rows if r["opening_time"] <= s["t"] <= deadline]
    if selected and charge:
        result["max_sampled_abs_eE_over_m2_open_to_reference"]=max(charge*abs(s["X"])/mu for s in selected)
        result["max_sampled_total_damping_over_two_M"]=max((3*s["H"]+s["Gamma"])/(2*s["M"]) for s in selected)
    if selected and "first_hot_exposure" in r and r["first_hot_exposure"] > 0:
        active = [s for s in selected if s["Gamma"] > 0 and s["decay_power"] > 0]
        if active:
            # A measured sample maximum is a lower bound on the supremum.
            # At alpha_g<=1 and all fixed masses within 5%, it yields an
            # intentionally permissive upper bound on Q*sqrt(alpha_sum).
            max_ratio = max(abs(s["chi"]) / math.sqrt(s["K"]) for s in active)
            if max_ratio > 0:
                upper = .05 / max_ratio
                lower = (p["kappa"]/mu*abs(r["first_hot_electric_work"])/r["first_hot_exposure"])
                result["fixed_mass_three_flavor_test"] = dict(
                    fractional_shift_allowance=.05, maximum_local_gauge_alpha=1.,
                    maximum_summed_yukawa_alpha=1., positive_birth_samples=len(active),
                    sampled_max_abs_chi_over_sqrt_K=max_ratio, upper_Q_alpha=upper,
                    lower_Q_alpha_for_recorded_work=lower, lower_over_upper=lower/upper)
    return result


def audit(data):
    out = {}
    for name, r in data["runs"].items():
        p=data["parameters"][name]
        rows=[dict(zip(data["trace_columns"], v)) for v in r.get("trace", [])]
        d=characteristic_scales(p,r)
        for point in ("first_hot","doubling"):
            if point+"_exposure" in r:
                d[point]=budget_at(r,p,rows,point)
        if "first_hot_aH" in r:
            f = .9
            d["general_geometric_upper_factor"]=(1+r["first_hot_aH"]**2)**(1/(4*f-2))
        if rows:
            d["trace_audit"]=trace_audit(rows,p,r)
        out[name]=d
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1]/"coupled/JOINED.json")
    parser.add_argument("--output", type=Path)
    args=parser.parse_args()
    result=json.dumps(audit(json.loads(args.source.read_bytes())),sort_keys=True,indent=2,allow_nan=False)+"\n"
    if args.output:
        args.output.write_text(result)
    else:
        print(result,end="")


if __name__=="__main__":
    main()
