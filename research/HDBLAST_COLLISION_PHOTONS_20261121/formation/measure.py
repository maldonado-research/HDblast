"""Event spacing, external acceleration, and energy processed by accepted pairs."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np


BANDS = ((1.e-5, 1.e-4), (1.e-4, 1.e-3), (1.e-3, 1.e-2), (.01, .1))


def neighbor_gaps(times, pairs, birth, end):
    """Return time since the predecessor and until the successor for each leg."""
    tags = pairs.ravel()
    wall_time = np.repeat(times, 2)
    order = np.lexsort((wall_time, tags))
    sorted_times, sorted_tags = wall_time[order], tags[order]
    repeated = sorted_tags[1:] == sorted_tags[:-1]
    before = birth[sorted_tags].copy()
    after = np.full(len(order), end)
    before[1:] = np.where(repeated, sorted_times[:-1], before[1:])
    after[:-1] = np.where(repeated, sorted_times[1:], after[:-1])
    gaps = []
    for delta in (sorted_times - before, after - sorted_times):
        if np.any(delta < 0.):
            raise ValueError("Event outside charged track")
        scattered = np.empty_like(delta)
        scattered[order] = delta
        gaps.append(scattered.reshape(-1, 2))
    return gaps


def local_chi(p, gamma, electric, magnetic, force_over_mass):
    length = np.linalg.norm(p, axis=-1, keepdims=True)
    direction = np.divide(p, length, out=np.zeros_like(p), where=length>0.)
    parallel = np.sum(direction * electric[:, None], axis=-1)
    perpendicular = electric[:, None] - direction * parallel[..., None]
    rest_field = gamma[..., None]*perpendicular + np.cross(p, magnetic[:, None])
    return force_over_mass * np.sqrt(parallel**2 + np.sum(rest_field**2, axis=-1))


def force_intervals(wall, wall_e, wall_b, times, charge, tolerance):
    # Endpoint maximum dominates both linearly interpolated field norms summed
    # throughout a step, including on an arbitrarily short partial interval.
    strength = np.linalg.norm(wall_e, axis=1) + np.linalg.norm(wall_b, axis=1)
    slope = abs(charge)*np.maximum(strength[:-1], strength[1:])
    cumulative = np.r_[0., np.cumsum(np.diff(wall)*slope)]
    if np.any(slope <= 0.):
        raise ValueError("This inverse assumes a nonzero field envelope on each cell")
    at_event = np.interp(times, wall, cumulative)
    start = np.interp(at_event-tolerance, cumulative, wall)
    stop = np.interp(at_event+tolerance, cumulative, wall)
    return times-start, stop-times


def portion(floor, lo, hi, weights):
    # At leading soft order the energy (not the photon count) is uniform in xi.
    available = np.clip((hi-floor)/(hi-lo), 0., 1.)
    return float(np.average(available, weights=weights))


def quantiles(values, weights):
    order = np.argsort(values)
    total = np.cumsum(weights[order])
    return [float(values[order[np.searchsorted(total, q*total[-1])]]) for q in (.5, .9, .99)]


def measure(data, coefficient, hard, interpolate):
    s, cosine, t = data["s"], data["mu"], data["t"]
    mass, charge = float(data["daughter_mass_L"]), float(data["charge_to_mass"])
    kinetic = (s-4.)/s
    if np.any(kinetic <= 0.) or np.any(kinetic >= 1.):
        raise ValueError("Expected massive, moving, accepted incoming particles")
    beta = np.sqrt(kinetic)
    alpha = 3.*charge**2*mass**2/(4.*np.pi*float(data["action"])*float(data["K"]))
    energy_in, energy_out = data["E_in"], data["E_out"]
    processed = data["weight"]*data["a"]*energy_in.sum(axis=1)/data["gas"]
    soft = processed*coefficient(s, cosine, alpha)*kinetic/2.
    if np.any(soft < 0.) or soft.sum() <= 0.:
        raise ValueError("A positive soft emission energy is required")
    pairs = np.column_stack((data["i"], data["j"]))
    previous, following = neighbor_gaps(t, pairs, data["birth"], data["time"])
    wall = data["wall_t"]
    electric, magnetic = [np.column_stack([np.interp(t, wall, data[key][:,k]) for k in range(3)])
                          for key in ("wall_E", "wall_B")]
    force = [local_chi(data[key], en, electric, magnetic, charge/mass)
             for key, en in (("p_in",energy_in),("p_out",energy_out))]
    # With xi = 2 omega*/(m sqrt(s) beta²), max one-radian proper
    # formation time per mass unit is (1+beta)/(xi beta²).
    factor = (1.+beta)/kinetic
    with np.errstate(divide="ignore"):
        time_factor = factor/mass * np.maximum((energy_in/previous).max(axis=1),
                                               (energy_out/following).max(axis=1))
    past_safe, future_safe = force_intervals(wall,data["wall_E"], data["wall_B"],t,charge,.1)
    with np.errstate(divide="ignore"):
        envelope_factor = factor/mass*np.maximum(energy_in.max(axis=1)/past_safe,
                                                 energy_out.max(axis=1)/future_safe)
    bend_factor = factor*np.maximum(force[0].max(axis=1), force[1].max(axis=1))
    bands = []
    for lo, hi in BANDS:
        bands.append(dict(xi_lo=lo, xi_hi=hi,
                          one_radian_time=portion(time_factor,lo,hi,soft),
                          ten_radian_time=portion(10.*time_factor,lo,hi,soft),
                          local_drift_below_one_tenth=portion(10.*bend_factor,lo,hi,soft),
                          integrated_force_envelope_and_ten_radians=portion(
                              np.maximum(10.*time_factor,envelope_factor),lo,hi,soft),
                          both_ten_radians_and_local_drift=portion(
                              10.*np.maximum(time_factor,bend_factor),lo,hi,soft)))
    # For a hard real photon the external line has virtuality +/-2 p.k.
    # At emission before scattering its exact energy-mismatch time is at most
    # E_in/(p.k). For emission after scattering, E_pair/(p_out.k) bounds that
    # time for every possible daughter momentum after hard recoil.
    quarter = .25
    incoming_clock = (factor/(mass*quarter))[:,None]*energy_in
    outgoing_clock = factor/(mass*quarter)*energy_in.sum(axis=1)
    time_in = np.all(10.*incoming_clock <= previous,axis=1)
    time_out = 10.*outgoing_clock <= following.min(axis=1)
    field_in = incoming_clock.max(axis=1) <= past_safe
    field_out = outgoing_clock <= future_safe
    raw = np.column_stack((hard["above_quarter_isr_energy_over_alpha"],
                           hard["above_quarter_fsr_energy_over_alpha"]))
    energy = alpha*processed[:,None]*interpolate(np.asarray(hard["log_s"]),raw[:,:,None],s)[:,:,0]
    total = energy.sum(axis=1)
    both_times = time_in & time_out
    both_fields = field_in & field_out
    robust = dict(energy_over_gas=float(total.sum()),
                  ten_radian_spacing_fraction=float(total[both_times].sum()/total.sum()),
                  field_margin_fraction=float(total[both_fields].sum()/total.sum()),
                  joint_both_sides_fraction=float(total[both_times & both_fields].sum()/total.sum()),
                  isr_joint_incoming_fraction=float(energy[time_in & field_in,0].sum()/energy[:,0].sum()),
                  fsr_joint_outgoing_surrogate_fraction=float(energy[time_out & field_out,1].sum()/energy[:,1].sum()))
    return dict(label=str(data["label"]), events=len(t), soft_weight_sum=float(soft.sum()),
                above_quarter_hard_band=robust, processed_pair_energy_over_gas=float(processed.sum()),
                average_shell_fraction_for_ten_percent=float(.1/processed.sum()),
                soft_weighted_xi_clock_quantiles=quantiles(time_factor,soft),
                soft_weighted_xi_local_bend_quantiles=quantiles(bend_factor,soft),
                band_energy_fractions=bands)


def checks():
    times = np.array([.1,.2,.4])
    pairs = np.array([[0,1],[0,3],[2,1]])
    gaps = neighbor_gaps(times,pairs,np.array([0.,0.,.3,.15]),.5)
    np.testing.assert_allclose(gaps[0], [[.1,.1],[.1,.05],[.1,.3]])
    np.testing.assert_allclose(gaps[1], [[.1,.3],[.3,.3],[.1,.1]])
    dt = force_intervals(np.array([0.,.5,1.]),np.array([[3.,0.,0.]]*3),
                         np.zeros((3,3)),np.array([.25,.75]),2.,.6)
    np.testing.assert_allclose(dt,[ [.1,.1],[.1,.1] ])
    weights = np.array([1.,2.,1.])
    np.testing.assert_allclose(portion(np.array([0.,.5,1.]),0.,1.,weights),.5)
    p = np.array([[[0.,0.,12.],[0.,0.,-12.]]])
    e, b = np.array([[0.,0.,3.]]), np.array([[0.,2.,0.]])
    np.testing.assert_allclose(local_chi(p,np.sqrt(145.)*np.ones((1,2)),e,b,1.),np.sqrt(585.))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--events", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--hard", required=True, type=Path)
    args = parser.parse_args()
    checks()
    coefficient = runpy.run_path(str(Path(__file__).resolve().parents[1]/"soft"/"soft_photons.py"))["coefficient"]
    hard_data = json.loads(args.hard.read_bytes())
    expected = {r["tape"]: r["bands"][-1]["fraction_of_endpoint_gas"] for r in hard_data["collision_histories"]}
    interpolate = runpy.run_path(str(Path(__file__).resolve().parents[1]/"hard"/"emission.py"))["interpolation"]
    rows = []
    for path in sorted(args.events.glob("*.npz")):
        with np.load(path, allow_pickle=False) as snapshot:
            row = measure(snapshot,coefficient,hard_data["energy_grid"],interpolate)
            if not np.isclose(row["above_quarter_hard_band"]["energy_over_gas"],expected[path.name],rtol=3.e-12,atol=0.):
                raise ArithmeticError("Inconsistent hard photon normalization")
            rows.append(row)
    if not rows:
        raise ValueError("No collision tapes")
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps({"runs":rows},indent=2,sort_keys=True,allow_nan=False)+"\n")
