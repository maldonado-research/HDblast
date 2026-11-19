"""Finite-band spectator trajectories immediately across actual first passage.

This imports the frozen dynamic Hamiltonian solver without modifying it.
Finite-stretch differences are not excess power caused by minimum passage.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("frozen_spectator", HERE.parent/"dynamic/spectator.py")
dynamic = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dynamic)

HALF_WIDTH = .005


class JoinedGeometry:
    def __init__(self, before, after):
        self.before, self.after = before, after
        self.crossing = before.t[-1]

    def sol(self, t):
        return (self.before if t <= self.crossing else self.after).sol(t)


def geometry_through_passage(module, branch, amplitude, degree, nodes,
                             tolerance=2e-12, max_step=.003):
    geometry = module.Evolution(branch, degree, nodes, monitor_nodes=6*degree+1)
    initial = geometry.initial(amplitude, 0.)
    sign = branch.fitted.sign
    target = sign*(branch.fitted.tmin-branch.top)
    def crossing(t, state):
        q = geometry.monitor @ state.reshape(6, geometry.n)[4]
        return target-max(sign*q)
    crossing.terminal, crossing.direction = True, -1
    control = dict(method="DOP853", rtol=tolerance, atol=2e-15,
                   max_step=max_step, dense_output=True)
    before = solve_ivp(geometry.equation, (0., .63), initial, events=crossing, **control)
    if not before.success or not len(before.t_events[0]):
        raise ArithmeticError("Actual first passage was not found inside the controlled interval")
    xstar = float(before.t_events[0][0])
    after = solve_ivp(geometry.equation, (xstar, xstar+HALF_WIDTH), before.y[:, -1], **control)
    if not after.success:
        raise ArithmeticError(after.message)
    solution = JoinedGeometry(before, after)
    chi = 0. if sign*amplitude>0 else np.pi/2
    at_worldline = geometry.matrices(np.array([np.cos(chi)]))[0][0]
    dense_y = np.cos(np.linspace(0., np.pi, 1537))
    dense_basis = geometry.matrices(dense_y)[0]
    low, high = branch.fitted.theta[[0, -1]]-branch.top
    samples = {}
    constraints = []
    margins = []
    for t in np.linspace(0., xstar+HALF_WIDTH, 49):
        state = solution.sol(float(t))
        record = geometry.record(float(t), state)
        constraints.append(record["constraints"])
        q = dense_basis @ state.reshape(6, geometry.n)[4]
        margins.append(float(min(q.min()-low, high-q.max())))
    for label, t in (("before", xstar-HALF_WIDTH), ("crossing",xstar), ("after",xstar+HALF_WIDTH)):
        state = solution.sol(t)
        qcoeff = state.reshape(6, geometry.n)[4]
        q = dense_basis @ qcoeff
        signed_worldline = sign*(at_worldline @ qcoeff)
        record = geometry.record(t, state)
        samples[label] = dict(x=t, signed_worldline_distance_to_minimum=float(target-signed_worldline),
                              global_signed_distance_to_minimum=float(target-max(sign*q)),
                              worldline_field=float(at_worldline @ qcoeff),
                              volume=record["volume"], H_volume=record["H_volume"],
                              minimum_spatial_volume_ratio=record["minimum_spatial_volume_ratio"])
    max_constraints = {key:max(row[key] for row in constraints) for key in constraints[0]}
    if min(margins) <= 1e-4:
        raise ArithmeticError("First-passage neighbourhood left the fitted-domain guard")
    return geometry, solution, dict(crossing_x=xstar, first_worldline_chi=chi,
        geometry_degree=degree, geometry_nodes=nodes, geometry_tolerance=tolerance,
        geometry_max_step=max_step, sample_times={key:row["x"] for key,row in samples.items()},
        geometry_samples=samples, minimum_sampled_fitted_domain_margin=min(margins),
        maximum_sampled_constraints=max_constraints)


def samples_for(hamiltonian, times, tolerance=2e-11, max_step=.006):
    values, ode = dynamic.propagate(hamiltonian, list(times.values()), tolerance=tolerance, max_step=max_step)
    return {label:values[f"{t:g}"] for label,t in times.items()}, ode.nfev


def compare_samples(primary, other):
    return {label:dynamic.compare_rows(primary[label],other[label]) for label in primary}


def stretch(samples):
    first,last = samples["before"],samples["after"]
    return dict(duration=last["x"]-first["x"],
        low_modes={key:dict(occupation_difference=last["low_modes"][key]["occupation"]-
                                                   first["low_modes"][key]["occupation"],
                           off_diagonal_contribution_difference=
                           last["low_modes"][key]["from_other_initial_modes"]-
                           first["low_modes"][key]["from_other_initial_modes"])
                   for key in first["low_modes"]},
        off_diagonal_beta_squared_sum_difference=last["off_diagonal_beta_squared_sum"]-
                                                first["off_diagonal_beta_squared_sum"])


def calculate(gravity):
    module = dynamic.intrinsic.gravity_module(gravity)
    key,branch = next((key,branch) for key,branch in module.constraint.load_branches(gravity/"previous").items()
                      if branch.fitted.pm<0)
    result = dict(branch=key, sample_half_width=HALF_WIDTH, primary_band_maximum=36,
        angular_ell=2, radial_parity="even n-ell", azimuthal_blocks="one fixed m block",
        convention="frozen dynamic solver: instantaneous original-field Hamiltonian in Einstein-frame geodesic proper time",
        interpretation="finite-stretch mode trajectories through first passage; no attribution of excess power to crossing, no UV sum or backreaction",
        kinetic_weight="B(Theta)=1/f(Theta); same-geometry Z=1 control",
        crossing_method="terminal first zero of minimum offset minus spatial maximum in the roll direction; pole-inclusive event grid and independent dense spatial check",
        runs={})
    for amplitude in (.001,-.001):
        geometry,solution,entry = geometry_through_passage(module,branch,amplitude,64,160)
        refined_geometry,refined_solution,refined_entry = geometry_through_passage(
            module,branch,amplitude,96,256,tolerance=2e-13,max_step=.0015)
        time_geometry,time_solution,time_entry = geometry_through_passage(
            module,branch,amplitude,64,160,tolerance=2e-13,max_step=.0015)
        times=entry["sample_times"]
        entry["refined_geometry"] = refined_entry
        entry["time_refined_geometry"] = time_entry
        entry["kinetic_controls"] = {}
        for kinetic in ("matter","geometric"):
            ham=dynamic.FieldHamiltonian(geometry,solution,36,192,kinetic)
            primary,nfev=samples_for(ham,times)
            run=entry["kinetic_controls"][kinetic]=dict(labels=ham.labels.tolist(),samples=primary,
                ode_evaluations=nfev,finite_stretch_difference=stretch(primary),band_refinements=[])
            for cutoff in (18,24,48):
                other,_=samples_for(dynamic.FieldHamiltonian(geometry,solution,cutoff,192,kinetic),times)
                run["band_refinements"].append(dict(maximum_n=cutoff,
                    comparisons=compare_samples(primary,other),
                    off_diagonal_beta_squared_sums={label:row["off_diagonal_beta_squared_sum"]
                                                   for label,row in other.items()}))
            other,_=samples_for(ham,times,tolerance=2e-12,max_step=.003)
            run["field_time_refinement"]=compare_samples(primary,other)
            other,_=samples_for(dynamic.FieldHamiltonian(geometry,solution,36,272,kinetic),times)
            run["quadrature_refinement"]=compare_samples(primary,other)
            other,_=samples_for(dynamic.FieldHamiltonian(refined_geometry,refined_solution,36,272,kinetic),times)
            run["spatial_geometry_refinement"]=compare_samples(primary,other)
            other,_=samples_for(dynamic.FieldHamiltonian(time_geometry,time_solution,36,192,kinetic),times)
            run["geometry_time_refinement"]=compare_samples(primary,other)
        result["runs"][f"{amplitude:g}"]=entry
        print(key,f"{amplitude:g}","crossing",entry["crossing_x"],"passage modes complete",flush=True)
    result["checks"]=verify(result)
    return result


def verify(result):
    event_errors=[]
    canonical=[]
    refinements=[]
    high_bands=[]
    constraints=[]
    for entry in result["runs"].values():
        for name in ("refined_geometry","time_refined_geometry"):
            event_errors.append(abs(entry[name]["crossing_x"]-entry["crossing_x"]))
        for samples in (entry["geometry_samples"],entry["refined_geometry"]["geometry_samples"]):
            if not (samples["before"]["global_signed_distance_to_minimum"]>0 and
                    samples["after"]["global_signed_distance_to_minimum"]<0):
                raise ArithmeticError("Samples do not bracket actual first passage")
            event_errors.append(abs(samples["crossing"]["global_signed_distance_to_minimum"]))
        for row in (entry,entry["refined_geometry"],entry["time_refined_geometry"]):
            constraints.extend(row["maximum_sampled_constraints"][key] for key in
                 ("max_abs_hamiltonian_over_6h2","max_abs_physical_momentum_over_h2"))
        for run in entry["kinetic_controls"].values():
            for sample in run["samples"].values():
                canonical.extend(sample["commutators"].values())
                canonical.extend([sample["bogoliubov_column_norm_error"],sample["bogoliubov_column_pair_error"]])
            for key in ("field_time_refinement","quadrature_refinement","spatial_geometry_refinement","geometry_time_refinement"):
                refinements.append(dynamic.comparison_maximum(run[key]))
            for band in run["band_refinements"]:
                if band["maximum_n"]>=24:
                    high_bands.append(dynamic.comparison_maximum(band["comparisons"]))
    checks=dict(maximum_event_time_or_field_distance_error=max(event_errors),
                maximum_canonical_error=max(canonical),
                maximum_time_quadrature_or_geometry_low_occupation_absolute_difference=max(refinements),
                maximum_n24_or_higher_low_occupation_absolute_difference=max(high_bands),
                maximum_sampled_normalized_physical_constraint=max(constraints))
    limits=dict(maximum_event_time_or_field_distance_error=2e-9,
                maximum_canonical_error=2e-9,
                maximum_time_quadrature_or_geometry_low_occupation_absolute_difference=2e-8,
                maximum_n24_or_higher_low_occupation_absolute_difference=2e-8,
                maximum_sampled_normalized_physical_constraint=1e-6)
    checks["required_limits"]=limits
    checks["passed"]=all(checks[name]<limit for name,limit in limits.items())
    if not checks["passed"]:
        raise ArithmeticError(f"Through-passage controls failed: {checks}")
    return checks


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gravity",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    data=calculate(args.gravity.resolve())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
