"""Local expansion law at the conserved-photon double root of the round constraint."""
from math import expm1, sqrt

from scipy.integrate import solve_ivp


def _integrate(branch, initial_q, force, ui, alpha, time_scale, log_radius_scale, find_maximum):
    """Evolve the exact field and Raychaudhuri equations in scaled variables."""
    r_initial = ui-alpha
    scaled_field = force*time_scale**2
    scaled_rate = time_scale**2/log_radius_scale

    def derivatives(s, state):
        displacement, velocity, log_radius, expansion = state
        q = initial_q+scaled_field*displacement
        normalized_force = float(branch.potential(q)[1])/force
        n = log_radius_scale*log_radius
        # h' = a_i^-2 exp(-2N) - 2 r_i exp(-4N) - v^2/2.
        # Combine the constants first so the O(alpha) acceleration remains visible.
        curvature_radiation = (ui+r_initial)*expm1(-2*n) - 2*r_initial*expm1(-4*n)
        dh = alpha + curvature_radiation - .5*(force*time_scale*velocity)**2
        return (velocity, normalized_force-3*log_radius_scale*expansion*velocity,
                expansion, scaled_rate*dh)

    def maximum(s, state):
        return state[3]

    maximum.direction = -1
    maximum.terminal = True
    result = solve_ivp(derivatives, (0., 2. if find_maximum else 1.), [0.]*4,
                       method='DOP853', rtol=2e-12, atol=2e-14, max_step=.025,
                       events=maximum if find_maximum else None)
    if not result.success:
        raise ArithmeticError(result.message)
    if find_maximum and (len(result.t_events[0]) != 1 or result.t[-1] <= 0.):
        raise ArithmeticError('Missing expanding turnaround')
    return float(result.t[-1]), float(result.y[2,-1])


def critical_check(branch, crosscheck_time=None):
    """Check the quartic critical maximum and the smaller subcritical expansion."""
    q = branch.fitted.sign*branch.fitted.dq
    ui, force = map(float, branch.potential(q))
    if ui <= 0 or force == 0:
        raise ValueError('Local test requires positive energy and a moving-field release force')
    critical = []
    for x in (.01, .003, .001, .0003):
        scale = force**2*x**4/24
        _, n = _integrate(branch, q, force, ui, 0., x, scale, False)
        # log(a/a_i) and a/a_i-1 have the same quartic leading coefficient.
        critical.append(dict(x=x, log_radius=scale*n, quartic_coefficient_ratio=-n))
    subcritical = []
    for gap in (1e-8, 1e-10, 1e-12, 1e-14, 1e-16):
        delta = sqrt(gap)
        alpha = 2*ui*delta/(1+delta)
        time_scale = sqrt(6*alpha)/abs(force)
        n_scale = 1.5*(alpha/force)**2
        s, n = _integrate(branch, q, force, ui, alpha, time_scale, n_scale, True)
        subcritical.append(dict(gap=gap, time_of_maximum=s*time_scale,
                                log_radius_gain=n*n_scale, time_ratio=s, gain_ratio=n))
    for ratios, tolerance in [(tuple(r['quartic_coefficient_ratio'] for r in critical), 5e-5),
                              (tuple(r['time_ratio'] for r in subcritical), 1e-3),
                              (tuple(r['gain_ratio'] for r in subcritical), 1e-3)]:
        errors = [abs(r-1) for r in ratios]
        if not all(earlier > later for earlier, later in zip(errors, errors[1:])) or errors[-1] > tolerance:
            raise ArithmeticError('Approach to critical asymptote failed')
    if crosscheck_time is not None and abs(crosscheck_time-subcritical[0]['time_of_maximum']) > 1e-7:
        raise ArithmeticError('Unscaled and scaled turnaround times disagree')
    return dict(force_at_release=force, critical=critical, subcritical=subcritical)
