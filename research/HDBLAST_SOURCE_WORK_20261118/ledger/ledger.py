"""Proper-time energy and pressure of a finite band on the simulated B1 shell."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicSpline


ARMS = {"quiet": "classical/runs/{grid}/quiet",
        "near_plus": "classical/runs/{grid}/0.001",
        "near_minus": "classical/runs/{grid}/-0.001",
        "far_plus": "far/runs/{grid}/0.001"}
START = 0.80
END = 2.10
G = 100.
STAR = 0.5


def history(source, grid, arm):
    checkpoint = source / f"classical/runs/{grid}/checkpoint.npz"
    branch = source / ARMS[arm].format(grid=grid)
    with np.load(checkpoint, allow_pickle=False) as cp:
        hist = cp["hist"]
    with np.load(branch/"rolling_timeseries.npz", allow_pickle=False) as data:
        cols = dict(zip(data["cols"], data["rec"].T, strict=True))
    event = json.loads((source/f"classical/runs/{grid}/quiet/rolling_summary.json").read_bytes())["matter"]["events"][0]
    split = cols["H0tau"][0]
    early = hist[:, 0] < split-1e-10
    times = np.r_[hist[early, 0], cols["H0tau"]]
    phi = CubicSpline(times, np.r_[hist[early, 1], cols["phi_b"]])
    hubble = CubicSpline(times, np.r_[hist[early, 2], cols["H_over_H0"]])
    primitive = hubble.antiderivative()
    offset = cols["A_b"][0] - event["A_star"] - primitive(split)
    return phi, hubble, lambda t: np.exp(offset + primitive(t))


def analyzer(k, phi, hubble, scale, s):
    """Second-adiabatic-order analyzer and its slope, defined as in the input study."""
    def omega(t):
        return np.sqrt((k/scale(t))**2 + (G*(phi(t)-STAR))**2
                       - 2.25*hubble(t)**2 - 1.5*hubble(t, 1))

    def frequency(t):
        ds = .001
        coeff = omega(t)
        before, after = omega(t-ds), omega(t+ds)
        slope = (after-before)/(2*ds)
        curvature = (after-2*coeff+before)/ds**2
        return coeff - curvature/(4*coeff**2) + 3*slope**2/(8*coeff**3)

    w = frequency(s)
    derivative = (frequency(s+.001) - frequency(s-.001))/.002
    x = 1/np.sqrt(2*w)
    return x, (-1j*w - .5*derivative/w)*x


def physical_terms(C, K, k2, m2):
    return (K+(k2+m2)*C)/2, (K-(m2+k2/3)*C)/2


def work_integrals(t, measure, dm2, h, variance, pressure):
    return {"mass_work":float(simpson((measure@variance)*dm2/2, x=t)),
            "expansion_work":float(simpson(-3*h*(measure@pressure), x=t))}


def radiation_era_control():
    # a(s)=sqrt(1+2H s), a''(eta)=0, m=m0/a.  The rescaled mode a*chi
    # has constant conformal frequency and is known for all proper times.
    k = np.array([2., 6.])
    measure = np.array([1.25, .75])
    m0, H = 4., .6
    t = np.linspace(0., 1.5, 1025)
    a = np.sqrt(1+2*H*t)
    h = H/a**2
    m2 = (m0/a)**2
    om = np.hypot(k, m0)[:,None]
    q = np.exp(-1j*om*(a-1)/H)/(a*np.sqrt(2*om))
    p = a**3*(-h-1j*om/a)*q
    C = a**3*abs(q)**2
    K = abs(p)**2/a**3
    energy, press = physical_terms(C, K, (k[:,None]/a)**2, m2)
    result = work_integrals(t, measure, -2*h*m2, h, C, press)
    gap = 1/a[-1]-1
    expected_mass = float(measure@(m0**2*gap/(2*om[:,0])))
    expected_expansion = float(measure@(k**2*gap/(2*om[:,0]) + H**2*(a[-1]**-3-1)/(4*om[:,0])))
    delta = float(measure@(energy[:,-1]-energy[:,0]))
    errors = [abs(result["mass_work"]-expected_mass),
              abs(result["expansion_work"]-expected_expansion),
              abs(delta-expected_mass-expected_expansion)]
    if max(errors)>3e-11:
        raise ValueError(f"Exact radiation-era control failed: {errors}")
    return {"mass_work":expected_mass, "expansion_work":expected_expansion,
            "total_change":delta, "max_absolute_error":max(errors)}


def calculate(source, grid, arm, order, steps):
    phi, hubble, scale = history(source, grid, arm)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    k = 4*(1 + nodes)
    measure = 4*weights*k**2/(2*np.pi**2)
    t = np.linspace(START, END, steps+1)
    x, xp = analyzer(k, phi, hubble, scale, START)
    # q=chi and p=a^3 dchi/ds evolve by the physical field Hamiltonian.
    a0 = scale(START)
    q = x/a0**1.5
    p = a0**1.5*(xp-1.5*hubble(START)*x)

    def equation(s, state):
        a = scale(s)
        w2 = (k/a)**2 + (G*(phi(s)-STAR))**2
        q, p = np.split(state, 2)
        return np.r_[p/a**3, -a**3*w2*q]

    sol = solve_ivp(equation, (START, END), np.r_[q, p], method="DOP853",
                    rtol=6e-11, atol=2e-13, max_step=.0025, t_eval=t)
    if not sol.success:
        raise ValueError(sol.message)
    q, p = np.split(sol.y, 2)
    a, h = scale(t), hubble(t)
    m2 = (G*(phi(t)-STAR))**2
    dm2 = 2*G**2*(phi(t)-STAR)*phi(t, 1)
    k2 = (k[:, None]/a)**2
    w2 = k2 + m2
    w = np.sqrt(w2)
    C = a**3*abs(q)**2
    K = abs(p)**2/a**3
    energy, pressure = physical_terms(C, K, k2, m2)
    zero = w/2
    C_vac = .5/w
    C_occ = (energy-zero)/w2
    p_vac = zero*k2/(3*w2)
    p_occ = (energy-zero)*k2/(3*w2)
    norm_error = np.max(abs(1j*(q.conj()*p-p.conj()*q)-1))
    if norm_error>2e-8:
        raise ValueError(f"Mode norm {norm_error}")

    pieces = {name:work_integrals(t, measure, dm2, h, var, pres) for name, var, pres in (
        ("full", C, pressure), ("instantaneous_vacuum", C_vac, p_vac),
        ("diagonal_excess", C_occ, p_occ),
        ("coherence", C-C_vac-C_occ, pressure-p_vac-p_occ))}
    energies = {}
    for name, term in (("full", energy),("instantaneous_vacuum",zero),
                       ("instantaneous_excess", energy-zero)):
        e0, ef = float(measure@term[:,0]), float(measure@term[:,-1])
        energies[name] = {"initial": e0, "final": ef, "change": ef-e0}

    # The finite-time number is a separate observable. Match both analyzers used here.
    xe, xpe = analyzer(k, phi, hubble, scale, END)
    X = a[-1]**1.5*q[:,-1]
    Xp = p[:,-1]/a[-1]**1.5 + 1.5*h[-1]*X
    adiabatic_occupancy = abs(xe*Xp-xpe*X)**2
    number_adiabatic = float(measure@adiabatic_occupancy)
    energy_adiabatic = float(measure@(adiabatic_occupancy*w[:,-1]))
    close_full = energies["full"]["change"]-sum(pieces["full"].values())
    close_zero = energies["instantaneous_vacuum"]["change"]-sum(pieces["instantaneous_vacuum"].values())
    if max(abs(close_full),abs(close_zero))>2e-5:
        raise ValueError(f"Ward balance failed: {close_full}, {close_zero}")
    arm_report = {"energies":energies, "work":pieces, "ward_residual":close_full, "vacuum_residual":close_zero,
            "max_mode_norm_error":float(norm_error), "final_scale_factor":float(a[-1]),
            "number_second_order":number_adiabatic, "energy_second_order_occupation":energy_adiabatic,
            "number_instantaneous":float(measure@(energy[:,-1]/w[:,-1]-.5)),
            "pressure_comoving_final":float(measure@pressure[:,-1])}
    trace = {"time":t, "dm2":dm2, "C":measure@C, "h":h, "pressure":measure@pressure,
             "energy":measure@energy, "C_excess":measure@(C-C_vac),
             "pressure_excess":measure@(pressure-p_vac)}
    return arm_report, trace


def paired_attribution(quiet, pulse):
    """Exact symmetric difference of two products on their common proper clock."""
    t = quiet["time"]
    def split(x, y, factor):
        mean_x = (quiet[x]+pulse[x])/2
        mean_y = (quiet[y]+pulse[y])/2
        drive = factor*mean_y*(pulse[x]-quiet[x])
        response = factor*mean_x*(pulse[y]-quiet[y])
        return drive, response

    mass_rates = split("dm2","C",.5)
    pressure_rates = split("h","pressure",-3)
    excess_mass_rates = split("dm2","C_excess",.5)
    excess_pressure_rates = split("h","pressure_excess",-3)
    def total(rate):
        return float(simpson(rate, x=t))
    result = {"mass_drive":total(mass_rates[0]), "mass_variance":total(mass_rates[1]),
              "hubble_rate":total(pressure_rates[0]), "pressure":total(pressure_rates[1]),
              "excess_mass_drive":total(excess_mass_rates[0]),
              "excess_mass_variance":total(excess_mass_rates[1]),
              "excess_hubble_rate":total(excess_pressure_rates[0]),
              "excess_pressure":total(excess_pressure_rates[1])}
    # At the index cuts the partial Simpson integrations always have an even interval count.
    trajectory = []
    for i in range(0, len(t), (len(t)-1)//32):
        energy = float(pulse["energy"][i] - quiet["energy"][i])
        m = float(simpson((mass_rates[0]+mass_rates[1])[:i+1], x=t[:i+1])) if i else 0.
        p = float(simpson((pressure_rates[0]+pressure_rates[1])[:i+1], x=t[:i+1])) if i else 0.
        trajectory.append({"s":float(t[i]), "energy":energy, "mass_work":m, "expansion_work":p})
    return {"integrals":result, "proper_clock_samples":trajectory}


def differences(result, quiet):
    delta = {}
    for key in ("energies", "work"):
        delta[key]={name:{f:val-quiet[key][name][f] for f,val in group.items()}
                    for name,group in result[key].items()}
    for field in ("number_second_order", "energy_second_order_occupation",
                  "number_instantaneous", "pressure_comoving_final", "ward_residual"):
        delta[field]=result[field]-quiet[field]
    return delta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path, help="Previous release's HDBLAST_PULSE_PARTICLE_RESPONSE_20261118 directory")
    p.add_argument("--grid", choices=["fine","coarse"], default="fine")
    p.add_argument("--order", type=int, default=64)
    p.add_argument("--steps", type=int, default=4096)
    args=p.parse_args()
    if args.steps%64 or args.steps<256 or args.order<8:
        raise ValueError("The time-grid size must be a multiple of 64 and at least 256; at least eight momentum nodes are required")
    pins=json.loads(Path(__file__).with_name("INPUTS.json").read_bytes())
    for name, digest in pins.items():
        if args.grid in Path(name).parts and hashlib.sha256((args.source/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f"Unexpected input: {name}")
    pairs={arm:calculate(args.source,args.grid,arm,args.order,args.steps) for arm in ARMS}
    raw={arm:pair[0] for arm,pair in pairs.items()}
    delta={arm:differences(raw[arm],raw["quiet"]) for arm in ARMS if arm!="quiet"}
    for arm in delta:
        attribution = paired_attribution(pairs["quiet"][1],pairs[arm][1])
        components = attribution["integrals"]
        for key, expected in (("mass_work", components["mass_drive"]+components["mass_variance"]),
                              ("expansion_work", components["hubble_rate"]+components["pressure"])):
            if abs(delta[arm]["work"]["full"][key]-expected)>1e-9:
                raise ValueError(f"Paired-source identity failed for {arm} {key}")
        delta[arm]["paired_attribution"] = attribution
    print(json.dumps({"grid":args.grid, "quadrature_order":args.order, "time_intervals":args.steps,
                      "start":START, "end":END, "band":[0,8], "radiation_era_control":radiation_era_control(), "arms":raw, "pulse_minus_quiet":delta},indent=2,allow_nan=False))


if __name__=="__main__":
    main()
