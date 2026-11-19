"""First-contact geometry, scalar rest energy and regional stress-energy balance."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import roots_legendre


HALF = np.pi / 2
TIMES = (1e-6, 1e-5, 1e-4, .001, .01)
CAP_EDGE = np.pi / 3


def quadrature(lo, hi, nodes):
    point, weight = roots_legendre(nodes)
    return lo + (hi - lo) * (1 + point) / 2, (hi - lo) * weight / 2


def load(gravity):
    spec = importlib.util.spec_from_file_location('gravity_evolve', gravity / 'evolution/evolve.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def angular_grid(evolution, chi):
    chi = np.atleast_1d(chi)
    y = np.cos(chi)
    return chi, y, np.sin(chi)**2, *evolution.matrices(y)


def local(evolution, t, z, grid):
    chi, y, w, U, Uy, Uyy = grid
    state = z.reshape(6, evolution.n)
    m, d, b, e, q, pi = (U @ state.T).T
    my, dy, qy = (Uy @ state[[0, 1, 4]].T).T
    qyy = Uyy @ state[4]
    a, H, _ = evolution.background(t)
    ell = m + w * d
    inv = np.exp(-2*ell)/a**2
    s = -np.sqrt(inv) * np.sin(chi) * qy
    J = np.exp(3*m + w*d)
    u, _ = evolution.branch.potential(q)
    alpha, beta = H + b + w*e, H + b
    ly = my + w*dy - 2*y*d
    lap = inv*(w*qyy - 3*y*qy + w*(2*my - ly)*qy)
    return dict(q=q, pi=pi, s=s, u=u, J=J, alpha=alpha, beta=beta, m=m,
                w=w, qchichi=w*qyy-y*qy, lap=lap)


def roots_on_half(function, points=257):
    angles = np.linspace(0., HALF, points)
    values = function(angles)
    brackets = np.flatnonzero(values[:-1]*values[1:] < 0.)
    return [brentq(lambda a: float(function([a])[0]), angles[i], angles[i+1],
                   xtol=1e-14, rtol=1e-14) for i in brackets]


def volume_integral(evolution, z, intervals, nodes=32):
    coeff = z.reshape(6, evolution.n)
    value = 0.
    for lo, hi in intervals:
        if hi <= lo:
            continue
        chi, weights = quadrature(lo, hi, nodes)
        basis = evolution.matrices(np.cos(chi))[0]
        m, d = (basis @ coeff[:2].T).T
        w = np.sin(chi)**2
        value += weights @ (w*np.exp(3*m+w*d))
    return value


def partition(points, condition):
    return [(lo, hi) for lo, hi in zip(points, points[1:]) if condition((lo+hi)/2)]


def field(evolution, z, angles):
    return evolution.matrices(np.cos(angles))[0] @ z.reshape(6, evolution.n)[4]


def crossed(evolution, z, qmin):
    difference = lambda chi: field(evolution, z, chi)-qmin
    roots = roots_on_half(difference)
    intervals = partition([0., *roots, HALF], lambda chi: difference([chi])[0] > 0.)
    total = volume_integral(evolution, z, [(0., HALF)], 80)
    return dict(fraction=volume_integral(evolution, z, intervals)/total,
                fronts_chi_in_northern_hemisphere=roots)


def extrema(metric, intervals):
    candidates = []
    for lo, hi in intervals:
        chi = np.linspace(lo, hi, 193)
        values = metric(chi)
        candidates.extend([values[0], values[-1]])
        for sign in (1., -1.):
            scaled = sign*values
            for i in np.flatnonzero((scaled[1:-1] < scaled[:-2]) &
                                   (scaled[1:-1] < scaled[2:]))+1:
                opt = minimize_scalar(lambda x: sign * float(metric([x])[0]),
                                      method='bounded', bounds=(chi[i-1],chi[i+1]),
                                      options={'xatol': 1e-13})
                candidates.append(sign*opt.fun)
    return [float(min(candidates)), float(max(candidates))]


def well(evolution, t, z):
    """The archived fiducial well is the phi_m +/- 0.1 interval."""
    branch = evolution.branch
    lo, hi = sorted(np.array(branch.fitted.barriers)-branch.top)
    qs = lambda chi: field(evolution, z, chi)
    boundaries = sorted([0., HALF, *roots_on_half(lambda chi: qs(chi)-lo),
                         *roots_on_half(lambda chi: qs(chi)-hi)])
    intervals = partition(boundaries, lambda chi: lo < qs([chi])[0] < hi)
    total = volume_integral(evolution, z, [(0., HALF)], 80)
    fraction = volume_integral(evolution, z, intervals)/total
    u_min = float(branch.potential(branch.fitted.tmin-branch.top)[0])
    capture = float(min(branch.potential(np.array([lo, hi]))[0]) - u_min)
    def measures(chi):
        v = local(evolution, t, z, angular_grid(evolution, chi))
        X = v['pi']**2 - v['s']**2
        return X, (v['u']-u_min + X/6)/capture, v['s']**2/v['pi']**2
    return dict(proper_volume_fraction=fraction, intervals_chi=intervals,
                capture_window_energy_above_min=capture,
                pi_squared_minus_spatial_gradient_squared=extrema(lambda c: measures(c)[0], intervals),
                rest_energy_above_min_over_capture_energy=extrema(lambda c: measures(c)[1], intervals),
                spatial_gradient_squared_over_pi_squared=extrema(lambda c: measures(c)[2], intervals))


def balance(evolution, sol, end, temporal_nodes):
    """Integrate each side of dE/dx = -(radial+angular work rate) - outward flux."""
    branch = evolution.branch
    a0 = evolution.background(0)[0]
    regions = {}
    for name, (lo, hi), orientation in [('two_polar_caps', (0., CAP_EDGE), 1.),
                                         ('equatorial_belt', (CAP_EDGE, HALF), -1.)]:
        chi, wc = quadrature(lo, hi, 64)
        region_grid = angular_grid(evolution, chi)
        edge_grid = angular_grid(evolution, [CAP_EDGE])
        def rates(t):
            z = sol.sol(t)
            v = local(evolution, t, z, region_grid)
            a = evolution.background(t)[0]
            volume = (4/np.pi)*(a/a0)**3 * wc*v['w']*v['J']/branch.u_top
            kin = v['pi']**2/6
            grad = v['s']**2/6
            rho = kin + grad + v['u']
            pr = kin + grad - v['u']
            pt = kin - grad - v['u']
            work = volume @ (v['alpha']*pr + 2*v['beta']*pt)
            edge = local(evolution, t, z, edge_grid)
            flux_density = -edge['pi'][0]*edge['s'][0]/3
            area = (4/np.pi)*a*a/a0**3 * edge['w'][0]*np.exp(2*edge['m'][0])
            return np.array([volume@rho, work, orientation*area*flux_density/branch.u_top])
        begin_energy, final_energy = rates(0.)[0], rates(end)[0]
        ts, wt = quadrature(0., end, temporal_nodes)
        _, work, flux = wt @ np.array([rates(t) for t in ts])
        delta = final_energy-begin_energy
        regions[name] = dict(initial_energy=begin_energy, final_energy=final_energy,
                             energy_change=delta, integrated_pressure_work=work,
                             integrated_outward_energy_flux=flux,
                             balance_residual=delta+work+flux)
    return regions


def one_run(module, branch, amplitude, degree, step, tolerance, temporal_nodes):
    evolution = module.Evolution(branch, degree, 4*degree, 257)
    qmin = branch.fitted.tmin - branch.top
    site = 0. if amplitude>0 else HALF
    site_grid = angular_grid(evolution, [site])
    def at_min(t, z):
        return (site_grid[3] @ z.reshape(6, evolution.n)[4])[0]-qmin
    at_min.terminal = True
    at_min.direction = 1
    opts = dict(method='DOP853', rtol=tolerance, atol=2e-15, max_step=step, dense_output=True)
    sol = solve_ivp(evolution.equation, (0., .70), evolution.initial(amplitude, 0.),
                    events=at_min, **opts)
    if not sol.success or not len(sol.t_events[0]):
        raise ValueError('Failed to reach a minimum')
    t = sol.t_events[0][0]
    z = sol.sol(t)
    post = solve_ivp(evolution.equation, (t, t+TIMES[-1]), z, **opts)
    if not post.success:
        raise ValueError(post.message)
    at = local(evolution, t, z, site_grid)
    pi = at['pi'][0]
    qcc = at['qchichi'][0]
    total = volume_integral(evolution, z, [(0., HALF)], 80)
    k = np.sqrt(2*pi / -qcc)
    dimensionality = 3 if site==0. else 1
    power = dimensionality/2
    coefficient = (at['J'][0]/total)*k**dimensionality/dimensionality
    points = []
    for dt in TIMES:
        result = crossed(evolution, post.sol(t+dt), qmin)
        if len(result['fronts_chi_in_northern_hemisphere']) != 1:
            raise ArithmeticError('Unexpected early level-set topology')
        points.append(dict(delta_x=dt, **result,
                           measured_over_leading=result['fraction']/(coefficient*dt**power)))
    first = well(evolution, t, z)
    after = well(evolution, t+TIMES[-1], post.sol(t+TIMES[-1]))
    fraction = np.array([p['fraction'] for p in points])
    exponents = np.diff(np.log(fraction))/np.diff(np.log(TIMES))
    # Check for another position ahead of the claimed first-touch site at every
    # solver step and additionally sample local Einstein constraints throughout.
    site_lag = max(float(max(evolution.monitor@v.reshape(6,evolution.n)[4]) -
                         (site_grid[3]@v.reshape(6,evolution.n)[4])[0]) for v in sol.y.T)
    audits = [evolution.record(x, sol.sol(x))['constraints'] for x in np.linspace(0, t, 17)]
    audits += [evolution.record(x, post.sol(x))['constraints'] for x in np.linspace(t, t+TIMES[-1], 5)]
    constraints = {name:max(r[name] for r in audits) for name in audits[0]}
    return dict(amplitude=amplitude, degree=degree, projection_nodes=4*degree, crossing_x=t,
                site='poles' if site==0. else 'equator', leading_volume_power=power,
                leading_volume_coefficient=coefficient, q_chichi=qcc, normal_velocity=pi,
                dimensionless_angular_front_coefficient=k,
                event_laplacian=at['lap'][0], local_theta=(at['alpha']+2*at['beta'])[0],
                local_kinetic_over_capture_energy=(pi*pi/6)/first['capture_window_energy_above_min'],
                largest_sampled_field_excess_over_first_site=site_lag,
                fractions_past_minimum=points, adjacent_log_slopes=exponents.tolist(),
                in_well_at_first_contact=first, in_well_0p01_later=after,
                energy_balance_to_first_contact=balance(evolution, sol, t, temporal_nodes),
                maximum_sampled_constraints=constraints)


def calculate(gravity):
    module = load(gravity)
    branch = module.constraint.load_branches(gravity/'previous')['phi_m=-0.5']
    controls = []
    runs = []
    for A in (.001, -.001, .0005, -.0005):
        primary = one_run(module, branch, A, 96, .003, 8e-12, 64)
        other = one_run(module, branch, A, 64, .004, 2e-11, 48)
        first = primary['in_well_at_first_contact']
        max_balance = max(abs(r['balance_residual']) for r in primary['energy_balance_to_first_contact'].values())
        control = dict(amplitude=A, other_degree=64, primary_degree=96,
                       crossing_x_difference=primary['crossing_x']-other['crossing_x'],
                       leading_coefficient_relative_difference=primary['leading_volume_coefficient']/
                           other['leading_volume_coefficient']-1,
                       max_crossed_volume_relative_difference=max(abs(a['fraction']/b['fraction']-1)
                         for a,b in zip(primary['fractions_past_minimum'],other['fractions_past_minimum'])),
                       max_in_well_volume_absolute_difference=max(abs(primary[w]['proper_volume_fraction']-
                           other[w]['proper_volume_fraction']) for w in ('in_well_at_first_contact','in_well_0p01_later')),
                       min_rest_energy_ratio_difference=first['rest_energy_above_min_over_capture_energy'][0]-
                           other['in_well_at_first_contact']['rest_energy_above_min_over_capture_energy'][0],
                       integrated_cap_flux_difference=primary['energy_balance_to_first_contact']['two_polar_caps']
                           ['integrated_outward_energy_flux']-other['energy_balance_to_first_contact']['two_polar_caps']
                           ['integrated_outward_energy_flux'], max_abs_regional_work_flux_residual=max_balance)
        if (max_balance>2e-9 or abs(control['crossing_x_difference'])>1e-7 or
                abs(control['leading_coefficient_relative_difference'])>2e-6 or
                control['max_crossed_volume_relative_difference']>2e-5 or
                control['max_in_well_volume_absolute_difference']>5e-7 or
                abs(control['min_rest_energy_ratio_difference'])>2e-6 or
                abs(control['integrated_cap_flux_difference'])>2e-8 or
                primary['largest_sampled_field_excess_over_first_site']>1e-10 or
                abs(primary['adjacent_log_slopes'][0]-primary['leading_volume_power'])>.001 or
                abs(primary['fractions_past_minimum'][0]['measured_over_leading']-1)>.001):
            raise ArithmeticError(('Patch verification failed', control, primary['adjacent_log_slopes']))
        runs.append(primary)
        controls.append(control)
    return dict(branch='phi_m=-0.5', initial_profile='A U2, raw offset zero',
                q_min=branch.fitted.tmin-branch.top,
                well_bounds_q=sorted((np.array(branch.fitted.barriers)-branch.top).tolist()),
                primary_max_step=.003, primary_relative_tolerance=8e-12,
                comparison_max_step=.004, comparison_relative_tolerance=2e-11,
                runs=runs, refinement=controls)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gravity', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    data = calculate(args.gravity)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
