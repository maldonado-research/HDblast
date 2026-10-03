#!/usr/bin/env python3
"""Independent homogeneous metric response from forced canonical modes.

No physical producer can run without both prospective registration and the
independent code/input manifest pinned by that registration.  --preflight only
inspects runtime, integer algebra, and the scheduled array dimensions.
"""
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

LD, CD = np.longdouble, np.clongdouble
EPSILON = LD("0.0001")
ETA_I, ETA_F = LD(-6), LD("-1.5")
OBSERVATIONS = (-5.5, -4.5, -4.0, -3.5, -2.5, -1.5)
CUTOFFS = (64, 128, 256)
SOURCES = ("positive_B", "signed_uB")
SETTINGS = (("coarse", 128, LD("0.5")), ("fine", 256, LD("0.25")))
QUANTITIES = ("q", "q_prime", "q_second", "rho", "p", "Q0", "rho0", "p0", "current")
WALL_BUDGET, RSS_BUDGET_KIB = 900.0, 256*1024
WRONSKIAN_GATE = 1e-10
REFINEMENT_GATES = {"q": 2e-10, "q_prime": 2e-10, "q_second": 1e-8,
                    "rho": 1e-8, "p": 1e-8, "current": 2e-10}
WARD_ENDPOINT_GATE, WARD_REFINEMENT_GATE = 2e-6, 1e-6


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def runtime_preflight():
    require(sys.version_info[:2] == (3, 12), "Pinned Python 3.12 required")
    require(np.__version__ == "2.2.6", "Pinned NumPy 2.2.6 required")
    require(sys.flags.optimize == 0, "Unoptimized pinned runtime required")
    require(np.finfo(LD).nmant >= 63, "Extended real/complex precision required")
    require(np.empty(1, dtype=CD).real.dtype == np.dtype(LD), "Complex precision mismatch")
    return {"python": platform.python_version(), "numpy": np.__version__,
            "platform": platform.platform(), "python_optimization": sys.flags.optimize,
            "real_dtype": np.dtype(LD).name, "complex_dtype": np.dtype(CD).name,
            "longdouble_nmant": np.finfo(LD).nmant, "longdouble_eps": str(np.finfo(LD).eps)}


def frozen_configuration():
    return {"epsilon": "0.0001", "eta_initial": "-6", "eta_final": "-1.5",
            "observations": list(OBSERVATIONS), "cutoffs": list(CUTOFFS), "sources": list(SOURCES),
            "settings": [{"name": n, "steps_per_unit": dt, "momentum_panel_width": str(w)}
                         for n, dt, w in SETTINGS],
            "momentum_gauss_nodes": 16, "local_time_gauss_nodes": 8,
            "mass_squared": 2, "xi": 0, "a0": "-1/eta", "H": 1,
            "state": "incoming_BD", "phi": "fixed", "mass_law_b": 1,
            "wronskian_projection": False, "quantities": list(QUANTITIES),
            "wall_budget_seconds": WALL_BUDGET, "rss_budget_kib": RSS_BUDGET_KIB}


def provenance_gate(args):
    require(len(args.freeze_commit) == 40 and
            all(c in "0123456789abcdef" for c in args.freeze_commit), "Public freeze commit required")
    for label, digest in (("registration", args.registration_sha256), ("manifest", args.manifest_sha256)):
        require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest),
                label + " SHA-256 required")
    require(sha(args.registration) == args.registration_sha256, "Registration SHA mismatch")
    registration = json.loads(args.registration.read_text())
    require(registration.get("schema_version") == 1,"Unknown registration schema")
    require(registration["independent_manifest_sha256"] == args.manifest_sha256,
            "Registration does not pin the supplied independent manifest")
    require(registration["frozen_configuration"]["independent_configuration"] == frozen_configuration(),
            "Registered independent settings do not match this producer")
    gates = registration["frozen_gates"]["independent"]
    require(gates == {"wronskian": WRONSKIAN_GATE, "refinement": REFINEMENT_GATES,
                     "ward_endpoint": WARD_ENDPOINT_GATE, "ward_refinement": WARD_REFINEMENT_GATE},
            "Registered independent gates do not match this producer")
    directory = Path(__file__).resolve().parent
    manifest_path = directory / "MANIFEST.json"
    require(sha(manifest_path) == args.manifest_sha256, "Independent manifest differs from registration")
    manifest = json.loads(manifest_path.read_text())
    listed = set()
    for entry in manifest["files"]:
        relative = Path(entry["path"])
        require(not relative.is_absolute() and ".." not in relative.parts,
                "Unsafe independent manifest path")
        path = directory / relative
        require(path.resolve().is_relative_to(directory), "Independent manifest escapes code directory")
        require(sha(path) == entry["sha256"], "Independent frozen input mismatch: " + str(relative))
        listed.add(str(relative))
    require({"forced_metric.py", "metric_wkb.py", "stable_baselines.py"}.issubset(listed),
            "Independent executable/module pin missing")
    checkpoint = args.registration.resolve().parent
    require(directory.resolve().is_relative_to(checkpoint), "Independent inputs must be in registered checkpoint")
    files = registration["files"]
    experiment_path = checkpoint/"EXPERIMENT.json"
    require(files.get("EXPERIMENT.json") == sha(experiment_path), "Experiment bytes are not registered")
    require(registration["frozen_configuration"] == json.loads(experiment_path.read_text()),
            "Registration configuration differs from registered experiment")
    require(registration["frozen_gates"] == registration["frozen_configuration"]["gates"],
            "Registration gate objects differ")
    for relative_string,digest in files.items():
        relative = Path(relative_string)
        require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe registration file path")
        path = checkpoint/relative
        require(path.resolve().is_relative_to(checkpoint), "Registered path escapes checkpoint")
        require(sha(path) == digest, "Registered frozen file mismatch: "+str(relative))
    require(files[str(directory.relative_to(checkpoint)/"MANIFEST.json")] == args.manifest_sha256,
            "Registration file table omits/mismatches independent manifest")
    for entry in manifest["files"]:
        relative = str(directory.relative_to(checkpoint)/entry["path"])
        require(files.get(relative) == entry["sha256"], "Independent manifest differs from full registration")
    return {"registration_sha256": args.registration_sha256,
            "manifest_sha256": args.manifest_sha256, "files": manifest["files"]}


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
    # Exact integer algebra: B^(n)=B P_n/(1-u^2)^(2n).
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


def forcing_jet(eta, jet, order=3):
    L = -LD(1)/np.asarray(eta, dtype=LD)
    # d^j L = j! L^(j+1), d^j L^2=(j+1)! L^(j+2).
    return np.stack([sum(math.comb(n, j)*(4*math.factorial(j+1)*L**(j+2)*jet[n-j]
                                           -2*math.factorial(j)*L**(j+1)*jet[n-j+1])
                         for j in range(n+1))-jet[n+2] for n in range(order+1)])


def oscillatory_propagator(k, distance):
    increment = np.expm1(CD(2j)*k*distance)
    phi = np.empty(np.broadcast_shapes(np.shape(k), np.shape(distance)), dtype=CD)
    nonzero = np.broadcast_to(k != 0, phi.shape)
    np.divide(increment, np.broadcast_to(CD(2j)*k, phi.shape), out=phi, where=nonzero)
    np.copyto(phi, np.broadcast_to(distance, phi.shape), where=~nonzero)
    return increment+1, phi


def momentum_rule(width):
    nodes, weights = leggauss(16)
    nodes, weights = nodes.astype(LD), weights.astype(LD)
    panels = int(LD(256)/width)
    left = np.arange(panels, dtype=LD)*width
    k = (left[:, None]+(nodes[None, :]+1)*width/2).ravel()
    wk = np.broadcast_to(weights[None, :]*width/2, (panels, 16)).ravel().copy()
    return k, wk


def resource_check(started):
    elapsed = time.monotonic()-started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    require(elapsed <= WALL_BUDGET, "Frozen 900-second wall budget exceeded")
    require(peak <= RSS_BUDGET_KIB, "Frozen 256-MiB peak-memory budget exceeded")
    return {"elapsed_seconds": elapsed, "peak_rss_kib": peak}


def frequency_jets(L, k, hjet):
    """Directional time derivatives from omega^2=k^2+M; no differencing."""
    M = 2*L*L
    mass = [math.factorial(n+1)*M*L**n for n in range(5)]
    delta_mass = [2*sum(math.comb(n,j)*mass[j]*hjet[n-j] for j in range(n+1)) for n in range(5)]
    omega = [np.sqrt(k*k+M)]
    delta = [delta_mass[0]/(2*omega[0])]
    for n in range(1,5):
        numerator = mass[n]-sum(math.comb(n,j)*omega[j]*omega[n-j] for j in range(1,n))
        omega.append(numerator/(2*omega[0]))
        variation = delta_mass[n]-sum(math.comb(n,j)*(delta[j]*omega[n-j]+omega[j]*delta[n-j])
                                    for j in range(1,n))-2*delta[0]*omega[n]
        delta.append(variation/(2*omega[0]))
    return omega, delta


def direct_observables(eta, k, u, w, jet, retain=False):
    from metric_wkb import evaluate_wkb
    from stable_baselines import baseline_differences
    a = L = -LD(1)/eta
    M = 2*a*a
    local_contacts = retain or np.count_nonzero(jet) != 0
    if local_contacts:
        omega, delta_omega = frequency_jets(L, k, jet)
        C = (2*L*L,4*L**3,12*L**4)
        delta_C = (jet[2]+2*L*jet[1],jet[3]+2*L*jet[2]+2*L*L*jet[1],
                   jet[4]+2*L*jet[3]+4*L*L*jet[2]+4*L**3*jet[1])
        sub = evaluate_wkb(k*k,L,M,omega,delta_omega,C,delta_C,jet[1],2*M*jet[0])
    else:
        # Outside the exact compact support every directional local
        # subtraction and prefactor contact vanishes identically. Baseline
        # operators are still integrated directly, with stable exact algebra.
        omega = (np.sqrt(k*k+M),)
    # Actual complex incoming-BD modes and unconstrained forced variations.
    v = np.exp(-CD(1j)*k*eta)/np.sqrt(2*k)
    vp = -CD(1j)*k*v
    delta_v = v*u
    delta_vp = v*(w-CD(1j)*k*u)
    Dv = vp-L*v
    delta_D_mode = delta_vp-L*delta_v
    delta_D_operator = -EPSILON*jet[1]*v
    delta_D = delta_D_mode+delta_D_operator
    delta_abs_v = 2*np.real(np.conj(v)*delta_v)
    delta_abs_D_mode = 2*np.real(np.conj(Dv)*delta_D_mode)
    delta_abs_D_operator = 2*np.real(np.conj(Dv)*delta_D_operator)
    base_abs_v = np.real(np.conj(v)*v)
    base_abs_D = np.real(np.conj(Dv)*Dv)
    delta_mass_contact = 2*EPSILON*M*jet[0]*base_abs_v
    bare_R0 = (base_abs_D+(k*k+M)*base_abs_v)/2
    bare_P0 = (base_abs_D-(k*k/3+M)*base_abs_v)/2
    bare_dR_mode = (delta_abs_D_mode+(k*k+M)*delta_abs_v)/(2*EPSILON)
    bare_dP_mode = (delta_abs_D_mode-(k*k/3+M)*delta_abs_v)/(2*EPSILON)
    bare_dR_operator = delta_abs_D_operator/(2*EPSILON)
    bare_dP_operator = bare_dR_operator
    bare_dR_mass = delta_mass_contact/(2*EPSILON)
    bare_dP_mass = -bare_dR_mass
    bare_dR = bare_dR_mode+bare_dR_operator+bare_dR_mass
    bare_dP = bare_dP_mode+bare_dP_operator+bare_dP_mass
    bare_A = base_abs_v
    A = delta_abs_v/EPSILON
    # Stable exact amplitude equations give the derivatives of this same
    # incoming-BD bilinear, whose reference modulus is 1/(2k).
    A1 = w.real/k/EPSILON
    g = forcing_jet(eta, jet, 0)[0]
    A2 = np.real(CD(2j)*k*w-EPSILON*g)/k/EPSILON
    # Exact rationalized bare-minus-subtraction baselines remove the leading
    # ultraviolet cancellation before quadrature; literal pieces remain raw.
    Qcanonical0,R0,P0 = baseline_differences(k,omega[0],a)
    if local_contacts:
        # S has variance grades0/2; R,P have stress grades0/2/4.
        S, S1, S2 = sub["S"], sub["S_prime"], sub["S_second"]
        dS, dS1, dS2 = sub["deltaS"], sub["deltaS_prime"], sub["deltaS_second"]
        q = A-dS-2*jet[0]*Qcanonical0
        q1 = A1-dS1-2*jet[1]*Qcanonical0+2*jet[0]*S1
        q2 = A2-dS2-2*jet[2]*Qcanonical0+4*jet[1]*S1+2*jet[0]*S2
        scaled_rho = bare_dR-sub["deltaR"]-4*jet[0]*R0
        scaled_p = bare_dP-sub["deltaP"]-4*jet[0]*P0
    else:
        q,q1,q2 = A,A1,A2
        scaled_rho,scaled_p = bare_dR,bare_dP
    pi = np.arccos(LD(-1))
    measure = k*k/(2*pi*pi)
    integrands = {"q": measure*q, "q_prime": measure*q1, "q_second": measure*q2,
                  "rho": measure*scaled_rho, "p": measure*scaled_p,
                  "Q0": measure*Qcanonical0/(a*a), "rho0": measure*R0/a**4,
                  "p0": measure*P0/a**4, "current": measure*q}
    linear_W = delta_v*np.conj(vp)+v*np.conj(delta_vp)-delta_vp*np.conj(v)-vp*np.conj(delta_v)
    canonical_energy = 2*np.real(np.conj(vp)*delta_vp)+k*k*delta_abs_v
    diagnostics = {"physical_wronskian_scaled": float(np.max(np.abs(linear_W))/EPSILON),
                   "canonical_energy_variation_max_over_epsilon": float(np.max(np.abs(canonical_energy))/EPSILON)}
    archive = {}
    if retain:
        archive = {"v": v, "v_prime": vp, "delta_v": delta_v, "delta_v_prime": delta_vp,
                   "D_v": Dv, "delta_D_v": delta_D, "delta_D_v_mode": delta_D_mode,
                   "delta_D_v_operator": delta_D_operator, "bare_variance0": bare_A,
                   "bare_R0": bare_R0, "bare_P0": bare_P0,
                   "bare_delta_R_mode": bare_dR_mode, "bare_delta_P_mode": bare_dP_mode,
                   "bare_delta_R_operator": bare_dR_operator, "bare_delta_P_operator": bare_dP_operator,
                   "bare_delta_R_mass": bare_dR_mass, "bare_delta_P_mass": bare_dP_mass,
                   "bare_delta_R_scale": -4*jet[0]*bare_R0,
                   "bare_delta_P_scale": -4*jet[0]*bare_P0,
                   "bare_delta_Q_scale": -2*jet[0]*bare_A,
                   "sub_delta_R_scale": -4*jet[0]*sub["R"],
                   "sub_delta_P_scale": -4*jet[0]*sub["P"],
                   "sub_delta_Q_scale": -2*jet[0]*S,
                   "bare_delta_R": bare_dR, "bare_delta_P": bare_dP,
                   "bare_delta_A": A, "bare_delta_A_prime": A1, "bare_delta_A_second": A2,
                   "baseline_Q_combined": Qcanonical0,"baseline_R_combined": R0,"baseline_P_combined": P0,
                   "baseline_Q_literal": bare_A-S,
                   "baseline_R_literal": bare_R0-sub["R"],"baseline_P_literal": bare_P0-sub["P"],
                   "canonical_energy_variation": canonical_energy/EPSILON,
                   "linear_wronskian": linear_W/EPSILON,
                   "source_jet": jet.copy(), "forcing_jet": forcing_jet(eta,jet)}
        archive.update({"C_"+str(n): np.asarray(value,dtype=LD) for n,value in enumerate(C)})
        archive.update({"delta_C_"+str(n): np.asarray(value,dtype=LD) for n,value in enumerate(delta_C)})
        archive.update({"omega_"+str(n): arr for n,arr in enumerate(omega)})
        archive.update({"delta_omega_"+str(n): arr for n,arr in enumerate(delta_omega)})
        archive.update({"sub_"+key: value for key,value in sub.items()})
    return integrands, diagnostics, archive


def cutoff_sums(integrands, weights, width):
    weighted = np.stack([integrands[name] for name in QUANTITIES])*weights[None,:]
    result = np.empty((len(CUTOFFS),len(QUANTITIES)), dtype=LD)
    for i, cutoff in enumerate(CUTOFFS):
        result[i] = np.sum(weighted[:,:int(LD(cutoff)/width)*16], axis=1, dtype=LD)
    require(np.all(np.isfinite(result)), "Nonfinite combined cutoff responses/baselines")
    return result


def simpson_to(history, end_step, dt):
    require(end_step%2 == 0, "Observation must align with an even Simpson endpoint")
    return dt/3*(history[0]+history[end_step]
                 +4*np.sum(history[1:end_step:2], axis=0, dtype=LD)
                 +2*np.sum(history[2:end_step:2], axis=0, dtype=LD))


def evolve(source_id, setting, output_dir, started, progress, active_arrays):
    name, steps_per_unit, width = setting
    dt = LD(1)/steps_per_unit
    steps = int((ETA_F-ETA_I)*steps_per_unit)
    k, weights = momentum_rule(width)
    gx, gw = leggauss(8)
    local_time = (gx.astype(LD)+1)*dt/2
    local_weight = gw.astype(LD)*dt/2
    E, drift = oscillatory_propagator(k, dt)
    E_local, drift_local = oscillatory_propagator(k[:,None], (dt-local_time)[None,:])
    force_w, force_u = E_local*local_weight[None,:], drift_local*local_weight[None,:]
    del E_local, drift_local
    u, w = np.zeros(k.shape,dtype=CD), np.zeros(k.shape,dtype=CD)
    observation_steps = {int((LD(eta)-ETA_I)*steps_per_unit): eta for eta in OBSERVATIONS}
    history_eta = ETA_I+np.arange(steps+1,dtype=LD)*dt
    history_values = np.empty((steps+1,len(CUTOFFS),len(QUANTITIES)),dtype=LD)
    history_ledger = np.empty((steps+1,len(CUTOFFS)),dtype=LD)
    history_baseline_contact = np.empty_like(history_ledger)
    history_hjet = np.empty((steps+1,6),dtype=LD)
    history_gjet = np.empty((steps+1,4),dtype=LD)
    history_geometry = np.empty((steps+1,10),dtype=LD)
    archive = {"k": k, "momentum_weights": weights}
    rows = []
    progress["completed_observations"] = rows
    amplitude_peak = physical_peak = energy_peak = 0.0
    qr, rr, pp, br, bp = (QUANTITIES.index(key) for key in ("q","rho","p","rho0","p0"))
    for step,eta in enumerate(history_eta):
        if step:
            previous = eta-dt
            if eta <= -5 or previous >= -3:
                # The interval is wholly outside the exact compact support.
                u,w = u+drift*w,E*w
            else:
                local_h = source_jet(previous+local_time,source_id,2)
                local_g = EPSILON*forcing_jet(previous+local_time,local_h,0)[0]
                forcing_u = np.einsum("ij,j->i",force_u,local_g,optimize=False)
                forcing_w = np.einsum("ij,j->i",force_w,local_g,optimize=False)
                u,w = u+drift*w-forcing_u,E*w-forcing_w
        progress.update({"step": step,"eta": float(eta)})
        residual = 2*u.real-w.imag/k
        amplitude_W = float(np.max(np.abs(residual))/EPSILON)
        require(np.all(np.isfinite(u)) and np.all(np.isfinite(w)), "Nonfinite forced evolution")
        require(amplitude_W <= WRONSKIAN_GATE, "Amplitude Wronskian failed")
        amplitude_peak = max(amplitude_peak,amplitude_W)
        jet = source_jet(eta,source_id)
        retain = step in observation_steps
        integrands,diagnostics,contact_archive = direct_observables(eta,k,u,w,jet,retain)
        require(diagnostics["physical_wronskian_scaled"] <= WRONSKIAN_GATE, "Physical Wronskian failed")
        physical_peak = max(physical_peak,diagnostics["physical_wronskian_scaled"])
        energy_peak = max(energy_peak,diagnostics["canonical_energy_variation_max_over_epsilon"])
        values = cutoff_sums(integrands,weights,width)
        a = L = -LD(1)/eta
        history_values[step] = values
        history_hjet[step] = jet
        history_gjet[step] = forcing_jet(eta,jet)
        history_geometry[step] = [a,L,L*L,2*L**3,6*L**4,2*a*a,4*a*a*L,
                                  12*a*a*L*L,48*a*a*L**3,240*a*a*L**4]
        baseline_contact = 3*jet[1]*a**4*(values[:,br]+values[:,bp])
        history_baseline_contact[step] = baseline_contact
        history_ledger[step] = L*values[:,rr]-3*L*values[:,pp]-baseline_contact
        active_arrays.update({"k": k,"momentum_weights": weights,"u": u,"w": w,
                              "history_eta": history_eta[:step+1],"history_values": history_values[:step+1],
                              "history_ledger_integrand": history_ledger[:step+1],
                              "history_source_jet": history_hjet[:step+1],
                              "history_forcing_jet": history_gjet[:step+1]})
        if eta <= -5:
            response_indices = [QUANTITIES.index(key) for key in ("q","q_prime","q_second","rho","p","current")]
            require(np.count_nonzero(u) == np.count_nonzero(w) == 0 and
                    np.count_nonzero(values[:,response_indices]) == 0, "Strict source-free initial interval failed")
        if step%32 == 0:
            resource_check(started)
        if not retain:
            continue
        ledger = simpson_to(history_ledger,step,dt)
        finite_k = [{"K": cutoff,"values": {key: float(values[i,j]) for j,key in enumerate(QUANTITIES)},
                     "ward": {"ledger": float(ledger[i]),"direct_endpoint": float(values[i,rr]),
                              "endpoint_difference": float(abs(ledger[i]-values[i,rr]))}}
                    for i,cutoff in enumerate(CUTOFFS)]
        rows.append({"source": source_id,"eta": float(eta),"a": float(a),
                     "source_jet_over_epsilon": {"h"+str(j): float(jet[j]) for j in range(6)},
                     "forcing_jet_over_epsilon": {"g"+str(j): float(history_gjet[step,j]) for j in range(4)},
                     "wronskian_scaled": amplitude_W,**diagnostics,"finite_k": finite_k})
        index = len(rows)-1
        archive["u_"+str(index)],archive["w_"+str(index)] = u.copy(),w.copy()
        archive.update({key+"_"+str(index): value for key,value in contact_archive.items()})
        archive.update({"integrand_"+key+"_"+str(index): value for key,value in integrands.items()})
        active_arrays.update({"archived_"+key: value for key,value in archive.items()})
    require(len(rows) == len(OBSERVATIONS),"Missing observation")
    archive.update({"history_eta": history_eta,"history_values": history_values,
                    "history_ledger_integrand": history_ledger,"history_baseline_contact": history_baseline_contact,
                    "history_source_jet": history_hjet,"history_forcing_jet": history_gjet,
                    "history_geometry": history_geometry,"history_quantity_names": np.array(QUANTITIES),
                    "history_geometry_names": np.array(("a0","L","L_prime","L_second","L_third",
                                                         "M","M_prime","M_second","M_third","M_fourth")),
                    "history_cutoffs": np.array(CUTOFFS,dtype=np.int64),
                    "observation_eta": np.array(OBSERVATIONS,dtype=LD)})
    path = output_dir/f"metric_modes_{source_id}_{name}.npz"
    np.savez_compressed(path,**archive)
    resource_check(started)
    return {"source": source_id,"setting": name,"dt": float(dt),"momentum_panel_width": float(width),
            "time_steps": steps,"momentum_nodes": len(k),"time_gauss_nodes": 8,
            "local_time_momentum_pairs": steps*8*len(k),"direct_stress_history_points": steps+1,
            "nonzero_support_local_time_momentum_pairs": 2*steps_per_unit*8*len(k),
            "wronskian_max_scaled": amplitude_peak,"physical_wronskian_max_scaled": physical_peak,
            "canonical_energy_variation_max_over_epsilon": energy_peak,
            "archive": {"path": path.name,"sha256": sha(path)},"rows": rows}


def combine_runs(runs):
    lookup = {(r["source"],r["setting"]): r for r in runs}
    maxima = {key: 0.0 for key in QUANTITIES}
    ward_max = ward_refinement_max = 0.0
    rows,failures = [],[]
    for source in SOURCES:
        for coarse,fine in zip(lookup[(source,"coarse")]["rows"],lookup[(source,"fine")]["rows"]):
            require(coarse["eta"] == fine["eta"],"Coarse/fine observation mismatch")
            row = {key:value for key,value in fine.items() if key != "finite_k"}
            row["finite_k"] = []
            for c,f in zip(coarse["finite_k"],fine["finite_k"]):
                diff = {key: abs(f["values"][key]-c["values"][key]) for key in QUANTITIES}
                for key,value in diff.items():
                    maxima[key] = max(maxima[key],value)
                    if key in REFINEMENT_GATES and value > REFINEMENT_GATES[key]:
                        failures.append(f"{source}/{fine['eta']}/{f['K']}/{key} refinement")
                ward_diff = abs(f["ward"]["ledger"]-c["ward"]["ledger"])
                ward_max = max(ward_max,f["ward"]["endpoint_difference"])
                ward_refinement_max = max(ward_refinement_max,ward_diff)
                if f["ward"]["endpoint_difference"] > WARD_ENDPOINT_GATE:
                    failures.append(f"{source}/{fine['eta']}/{f['K']} Ward endpoint")
                if ward_diff > WARD_REFINEMENT_GATE:
                    failures.append(f"{source}/{fine['eta']}/{f['K']} Ward refinement")
                row["finite_k"].append({"K": f["K"],"values": f["values"],"coarse_values": c["values"],
                                         "refinement_differences": diff,
                                         "estimated_numerical_errors": {key:2*value for key,value in diff.items()},
                                         "ward": {**f["ward"],"coarse_ledger": c["ward"]["ledger"],
                                                  "coarse_endpoint_difference": c["ward"]["endpoint_difference"],
                                                  "refinement_difference": ward_diff}})
            rows.append(row)
    return rows,maxima,ward_max,ward_refinement_max,failures


def preflight():
    runtime = runtime_preflight()
    require(BUMP_POLYNOMIALS[1] == (0,-2),"Integer bump derivative recurrence failed")
    inventory = [{"source": s,"setting": n,"time_steps": int((ETA_F-ETA_I)*dt),
                  "momentum_nodes": int(LD(256)/width)*16,
                  "local_time_momentum_pairs": int((ETA_F-ETA_I)*dt)*int(LD(256)/width)*16*8}
                 for s in SOURCES for n,dt,width in SETTINGS]
    return {"status": "pure_input_preflight_passed","physical_evaluations": 0,
            "runtime": runtime,"configuration": frozen_configuration(),"scheduled_runs": inventory,
            "limits": ["Timing is not measured before public freeze.",
                       "All physical source values, contacts, modes and integrals remain unexecuted."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight",action="store_true")
    parser.add_argument("--freeze-commit")
    parser.add_argument("--registration",type=Path)
    parser.add_argument("--registration-sha256")
    parser.add_argument("--manifest-sha256")
    parser.add_argument("--output-dir",type=Path)
    args = parser.parse_args()
    if args.preflight:
        require(not any((args.freeze_commit,args.registration,args.registration_sha256,
                         args.manifest_sha256,args.output_dir)),"Preflight cannot produce physical output")
        print(json.dumps(preflight(),indent=2,allow_nan=False))
        return
    require(all((args.freeze_commit,args.registration,args.registration_sha256,
                 args.manifest_sha256,args.output_dir)),"All prospective freeze/provenance/output arguments required")
    runtime = runtime_preflight()
    provenance = provenance_gate(args)
    checkpoint = args.registration.resolve().parent
    require(not args.output_dir.resolve().is_relative_to(checkpoint),"Output must be external to frozen checkpoint")
    require(not args.output_dir.exists(),"Refusing to replace output directory")
    args.output_dir.mkdir(parents=True)
    started = time.monotonic()
    report = {"schema_version": 1,"route": "independent_metric_forced_modes_direct_minimal_stress",
              "freeze_commit": args.freeze_commit,"provenance": provenance,"runtime": runtime,
              "configuration": frozen_configuration(),"epsilon": float(EPSILON),"mass_law_b": 1,
              "normalizations": {"q": "a0^2 deltaQ/epsilon","q_prime": "(a0^2 deltaQ)'/epsilon",
                                 "q_second": "(a0^2 deltaQ)''/epsilon","rho": "a0^4 delta_rho/epsilon",
                                 "p": "a0^4 delta_p/epsilon","Q0": "physical Q0,K",
                                 "rho0": "physical rho0,K","p0": "physical p0,K",
                                 "current": "a0^2 delta_j/epsilon=q, fixed phi, b=1, x_phi=2"},
              "status": "running","runs": [],"active_run": None}
    result_path = args.output_dir/"results.json"
    active_arrays = {}
    try:
        for source in SOURCES:
            for setting in SETTINGS:
                progress = {"source": source,"setting": setting[0],"step": 0,"eta": float(ETA_I)}
                report["active_run"] = progress
                result_path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
                report["runs"].append(evolve(source,setting,args.output_dir,started,progress,active_arrays))
                report["active_run"] = None
                active_arrays.clear()
                report["resources"] = resource_check(started)
                result_path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
                print(json.dumps({"source": source,"setting": setting[0],"status": "completed",
                                  "resources": report["resources"]}),flush=True)
        rows,maxima,ward_max,ward_refinement,failures = combine_runs(report["runs"])
        report.update({"rows": rows,"maximum_refinement_differences": maxima,
                       "maximum_ward_endpoint_difference": ward_max,
                       "maximum_ward_refinement_difference": ward_refinement,"gate_failures": failures})
        require(not failures,"Registered internal gates failed: "+", ".join(failures))
        report["status"] = "passed_internal_gates"
        report["limits"] = ["Both stresses are direct physical minimal-mode bilinears with operator, mass and scale contacts.",
                            "The complete metric direction of W2/W4 and J2/J4 is evaluated independently.",
                            "The Ward ledger uses direct stresses and actual finite-K baseline rho0+p0; it defines neither stress.",
                            "Refinement and cross-route estimates are not certified integration error bounds.",
                            "Canonical energy variation is a diagnostic only and establishes no particle yield or heating.",
                            "One fixed background and incoming BD state, fixed physical mass, linear homogeneous metric response only."]
    except Exception as exc:
        report.update({"status": "failed","exception": {"type": type(exc).__name__,"message": str(exc),
                                                            "traceback": traceback.format_exc()}})
        if active_arrays:
            snapshot = args.output_dir/"failure_snapshot.npz"
            try:
                np.savez(snapshot,**active_arrays)
                report["failure_snapshot"] = {"path": snapshot.name,"sha256": sha(snapshot)}
            except Exception as snapshot_error:
                report["failure_snapshot_error"] = str(snapshot_error)
        raise
    finally:
        report["resources"] = {"elapsed_seconds": time.monotonic()-started,
                               "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        result_path.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status": report["status"],"result": str(result_path),
                      "maximum_refinement_differences": maxima,
                      "maximum_ward_endpoint_difference": ward_max,
                      "maximum_ward_refinement_difference": ward_refinement}),flush=True)


if __name__ == "__main__":
    main()
