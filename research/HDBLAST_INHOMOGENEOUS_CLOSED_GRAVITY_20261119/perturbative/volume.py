"""Second-order proper-volume response of a time-symmetric S^3 Einstein–scalar slice."""
import argparse
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import quad, solve_ivp
from scipy.special import hyp2f1


TIMES = (0., .03, .06, .09, .12)
K2 = 8


def load_potentials(previous):
    """Read only the archived potential and its stationary-point parameters."""
    spec = importlib.util.spec_from_file_location("prior_zonal", previous / "nonlinear/nonlinear.py")
    prior = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(prior)
    hilltop, model, refs, tables = prior.prior_model(previous / "prior")
    for pm in (.5, -.5):
        key = f"phi_m={pm}"
        potential = hilltop.BG.Potential(model, tables, pm, refs[key])
        top, h, mu2, _ = hilltop.field_top(potential)
        yield key, potential, top, h, mu2


def concavity_on_segments(potential, top, radius):
    """Bound each spline segment's curvature with exact rational arithmetic.

    Scale each segment to t in [0,1]. A polynomial c(t) obeys
    c(t) <= c[0] + sum(max(0,c[j])) for j >= 1.
    """
    def upper(poly):
        return poly.coef[0] + sum(max(0, c) for c in poly.coef[1:])

    def diff(poly):
        # Polynomial.deriv rescales its coefficient array in floating point.
        return Polynomial([j*c for j, c in enumerate(poly.coef) if j])

    theta = potential.spline.x
    segments = np.flatnonzero((theta[1:] >= top-radius) & (theta[:-1] <= top+radius))
    curvature_bounds = []
    v0, b, pm = (Fraction(float(v)) for v in (potential.v0, potential.b, potential.pm))
    for i in segments:
        width = Fraction(float(theta[i+1])) - Fraction(float(theta[i]))
        phi, f = (Polynomial([Fraction(float(c))*width**j for j, c in enumerate(cs[::-1])])
                  for cs in potential.spline.c[:, i, :].T)
        if v0 <= 0 or upper(-f) >= 0:
            raise ArithmeticError("Uncertified positive potential denominator")
        numerator = (1-phi/pm)**2 * (1+b*phi)
        fp, dp = diff(f), diff(numerator)
        curvature_numerator = (diff(dp)*f**2 - 4*dp*f*fp
                               - 2*numerator*f*diff(fp) + 6*numerator*fp**2)
        bound = upper(curvature_numerator)
        if not isinstance(bound, Fraction) or bound >= 0:
            raise ArithmeticError("Segment concavity failed")
        curvature_bounds.append(bound / (v0*width**2*upper(f)**4))
    return dict(symmetric_field_radius=radius, spline_segments_tested=len(segments),
                maximum_segment_upper_bound_on_u_second_derivative=float(max(curvature_bounds)))


def transfer(x, h, mu2):
    """Even linear n=2 mode and its analytic proper-time derivative."""
    sh, ch = math.sinh(h*x), math.cosh(h*x)
    disc = math.sqrt(9 + 4*mu2/h**2)
    p = 1 + (3 + disc)/4
    q = 1 + (3 - disc)/4
    z = -sh**2
    f = hyp2f1(p, q, .5, z)
    fz = 2*p*q*hyp2f1(p+1, q+1, 1.5, z)
    return ch**2*f, 2*h*sh*ch*(f - ch**2*fz)


def force(T, Tp, mu2):
    return -mu2*T**2/6 - Tp**2/3


def integrate_response(h, mu2, times):
    B0 = (mu2/h**2 - K2)/12

    def rhs(x, state):
        T, Tp, B, Bp = state
        H = h*math.tanh(h*x)
        inv_a2 = (h/math.cosh(h*x))**2
        return Tp, (mu2 - K2*inv_a2)*T - 3*H*Tp, Bp, force(T, Tp, mu2) - 2*H*Bp

    ode = solve_ivp(rhs, (0., times[-1]), (1., 0., B0, 0.), method="DOP853",
                    rtol=4e-13, atol=2e-14, max_step=.015, t_eval=times)
    if not ode.success:
        raise ArithmeticError(ode.message)
    rows = []
    for x, state in zip(ode.t, ode.y.T):
        T, Tp = transfer(x, h, mu2)
        H = h*math.tanh(h*x)
        inv_a2 = (h/math.cosh(h*x))**2

        def source(s):
            return math.cosh(h*s)**2*force(*transfer(s, h, mu2), mu2)

        def kernel(s):
            # sinh(h(x-s))/(cosh(hx)cosh(hs)) is tanh(hx)-tanh(hs).
            return math.sinh(h*(x-s))/(h*math.cosh(h*x)*math.cosh(h*s))

        B = B0 + quad(lambda s: source(s)*kernel(s), 0., x, epsabs=1e-13, epsrel=1e-13)[0]
        Bp = quad(source, 0., x, epsabs=1e-13, epsrel=1e-13)[0] / math.cosh(h*x)**2
        # The integrated Hamiltonian constraint gives a first integral.
        rho2 = (Tp**2 + (K2*inv_a2 - mu2)*T**2)/6
        first_integral_error = 2*H*Bp - 2*B*inv_a2 - rho2
        difference = np.abs(state - [T, Tp, B, Bp])
        rows.append(dict(x=x, T=T, Tp=Tp, radius_fractional_coefficient=B,
                         volume_fractional_coefficient=3*B, H_volume_difference_coefficient=Bp,
                         fractional_growth_since_initial_coefficient=B-B0,
                         acceleration_over_radius_difference_coefficient=force(T, Tp, mu2),
                         energy_density_shift_coefficient=rho2,
                         averaged_Hamiltonian_absolute_error=abs(first_integral_error),
                         absolute_ODE_minus_quadrature_error=dict(zip(("T", "Tp", "B", "Bp"), difference))))
        if abs(first_integral_error) > 2e-10 or max(difference) > 2e-10:
            raise ArithmeticError("Second-order volume cross-check failed")
    return rows


def calculate(previous):
    branches = {}
    for key, pot, top, h, mu2 in load_potentials(previous):
        cap = .04 if pot.pm > 0 else .08
        branches[key] = dict(h=h, mu_squared=mu2, times=integrate_response(h, mu2, TIMES),
                             potential_concavity=concavity_on_segments(pot, top, cap))
    return dict(amplitude_convention="q(0,chi)=A U_2(cos chi), q_x(0,chi)=0",
                coefficient_convention="X(A)=X(0)+A^2 C_X+O(|A|^3), at fixed proper time",
                branches=branches)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.previous), indent=2, allow_nan=False) + "\n")
