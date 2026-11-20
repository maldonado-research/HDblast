"""Verify field stress, energy exchange, coherent supply and its post-ramp bound."""
from pathlib import Path
import argparse
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent


def stable_bound(run, columns, p):
    """Positive expansion, fixed mass and zero drive imply I >= I_left (a_left/a_right)^3."""
    target=.0008/(2*p["daughter_mass"]*p["initial_mass_ratio"])
    rows=[dict(zip(columns,x)) for x in run["trace"]]
    left_a=run["mass_complete_radius"]
    left_ratio=run["mass_complete_target_ratio"]
    minimum=left_ratio
    for r in rows:
        if not (run["mass_complete_time"]<r["t"]<run["maximum_tau"]):
            continue
        if abs(r["M"]-2*p["daughter_mass"]*p["final_mass_ratio"])>1e-12 or r["v"]!=0 or r["H"]<=0:
            raise ValueError("Bound requires freely expanding constant-mass oscillator")
        if r["a"]<left_a:
            raise ValueError("Unexpected radius decrease")
        minimum=min(minimum,left_ratio*(left_a/r["a"])**3)
        left_a=r["a"]
        left_ratio=r["a"]**3*r["rho_chi"]/(r["M"]*target)
    return min(minimum,left_ratio*(left_a/run["maximum_radius"])**3)


def verify(data):
    source=HERE/"production_forward.cpp"
    if hashlib.sha256(source.read_bytes()).hexdigest()!=data["solver_sha256"]:
        raise ValueError("Solver hash mismatch")
    results={}
    for key in ("primary","stable","heavy_stable"):
        p=data["parameters"][key]
        run=data["runs"][key]
        rows=[dict(zip(data["trace_columns"],r)) for r in run["trace"]]
        initial=rows[0]
        if any(initial[k]!=0 for k in ("psi","pi","chi","chi_dot","v","rho_chi","e","rad","B","X")):
            raise ValueError("Nonempty initial wall excitations")
        if initial["mass"]!=4. or initial["a"]!=2. or initial["H"]!=0.:
            raise ValueError("Incorrect initial geometry")
        target=.0008/(2*p["daughter_mass"]*p["initial_mass_ratio"])
        for r in rows:
            kinetic=.5*r["chi_dot"]**2
            potential=.5*(r["M"]*(r["chi"]-r["v"]))**2
            if max(abs(kinetic+potential-r["rho_chi"]),abs(kinetic-potential-r["p_chi"]))>1e-11:
                raise ValueError("Incorrect oscillator stress")
        err=max(abs(run[f"{e}_radius"]**3*run[f"{e}_density"] +run[f"{e}_pressure_work"]
                    -run[f"{e}_drive_work"]-run[f"{e}_mass_work"])
                for e in ("drive_end","opening"))
        if err>5e-9 or run["driver_balance_error"]>5e-8 or run["oscillator_balance_error"]>2e-11:
            raise ValueError("Scalar energy ledger failed")
        if abs(run["opening_action"]/target-run["opening_target_ratio"])>2e-10:
            raise ValueError("Incorrect abundance normalization")
        results[key]={"opening_ratio":run["opening_target_ratio"],"event_ledger_error":err}
        if key!="primary":
            cycles=run["half_cycles"]
            values=[c[3]/target for c in cycles]
            bound=stable_bound(run,data["trace_columns"],p)
            if bound<=1 or max(abs(c[5]) for c in cycles)>1e-7:
                raise ValueError("Late supply or volume-weighted virial check failed")
            results[key].update(complete_half_cycles=len(cycles),turning_point_ratio_range=[min(values),max(values)],
                                lower_bound_all_postramp_expansion=bound)
    for key in ("drive_end_action","opening_action","opening_radius","opening_time"):
        if data["runs"]["flip_drive"][key]!=data["runs"]["primary"][key]:
            raise ValueError("Sign-symmetry check failed")
    if data["runs"]["no_drive"]["opening_action"]!=0. or "bounce_radius" in data["runs"]["no_pulse"]:
        raise ValueError("Null controls failed")
    return results


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data",nargs="?",type=Path,default=HERE/"PARENT.json")
    print(json.dumps(verify(json.loads(parser.parse_args().data.read_text())),indent=2))
