"""Compare the simulated stress and energy with both counterstreaming and wall-rest thermal gases."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import kve

def fluid(theta, stream_u, enthalpy):
    gamma=float(np.hypot(1.,stream_u))
    # Per wall-frame particle, with the number density Gamma*n_rest.
    energy=enthalpy*gamma-theta/gamma
    sideways=theta/gamma
    axial=sideways+enthalpy*stream_u**2/gamma
    trace=2*sideways+axial
    return dict(mean_gamma=float(energy),transverse_pressure_each_over_trace=float(sideways/trace),
                longitudinal_pressure_over_trace=float(axial/trace))

def compare(observed,model):
    transverse=model["transverse_pressure_each_over_trace"]
    return dict(mean_energy_over_boosted_prediction=observed["mean_gamma"]/model["mean_gamma"],
        least_principal_stress_fraction_over_prediction=observed["stress_eigen_min"]/transverse,
        middle_principal_stress_fraction_over_prediction=observed["stress_eigen_mid"]/transverse,
        largest_principal_stress_fraction_minus_prediction_points=100*(observed["stress_eigen_max"]-
                                                                  model["longitudinal_pressure_over_trace"]),
        ratio_of_two_smaller_stresses=observed["stress_eigen_min"]/observed["stress_eigen_mid"])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history",type=Path,required=True)
    p.add_argument("--calibration",type=Path,required=True);p.add_argument("--axial-budget",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.history/"markers.json").read_text())
    cal=json.loads(args.calibration.read_text());budget=json.loads(args.axial_budget.read_text())
    result={}
    for label in ("dedicated_half","shared_half"):
        marker=meta["runs"][label]["markers"]["0.01"]
        with np.load(args.history/(label+".npz")) as h:
            u=float(np.interp(marker["time"],h["t"],h["u"]))
        theta=marker["theta"]
        # Exact relativistic thermal gas and the algebraic Taub–Mathews closure of the wall code.
        exact=fluid(theta,u,kve(3,1/theta)/kve(2,1/theta))
        approximate=fluid(theta,u,2.5*theta+np.hypot(1.,1.5*theta))
        rest=fluid(theta,0.,kve(3,1/theta)/kve(2,1/theta))
        selected=cal[label]["best_fine"]
        resolutions={r['resolution']:compare(r["result"],exact) for r in
                     [selected,*cal[label]["independent_resolution"]]}
        bound=budget[label]["0.01"]
        if not np.isclose(rest["mean_gamma"],bound["thermal_mean_gamma"],rtol=1e-13):
            raise ValueError("Incompatible thermal equations of state")
        result[label]=dict(momentum_drift_u=u,bulk_lorentz_gamma=float(np.hypot(1.,u)),temperature=theta,
            source_amplitude=selected["amplitude"],k_physical=selected["k_physical"],
            wall_rest_maxwell_juttner=rest,boosted_maxwell_juttner=exact,
            boosted_taub_mathews=approximate,comparisons_by_resolution=resolutions,
            fine_comparison_to_archive_closure=compare(selected["result"],approximate),
            seven_mode_floor_over_boosted_thermal_mean=bound["energy_over_thermal_lower"]["7"]*
                    rest["mean_gamma"]/exact["mean_gamma"])
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
if __name__=="__main__":main()
