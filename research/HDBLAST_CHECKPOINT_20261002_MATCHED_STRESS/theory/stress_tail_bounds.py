"""Deferred constructive UV bounds for registered fixed-geometry stress.

Import performs no pulse/source/response evaluation. Call only after public
freeze. H=1, r=2, b=1; q=a^2 deltaQ/epsilon, rho,p use a^4/epsilon,
and the current uses a^2/epsilon. Only the omitted k>K band is bounded here.
"""
from fractions import Fraction
from functools import lru_cache
from math import inf, nextafter
from pathlib import Path
import importlib.util


INTERVAL_DECIMAL_PRECISION = 60
OBSERVATIONS = tuple(Fraction(v) for v in ("-5.5", "-4.5", "-4", "-3.5", "-2.5", "-1.5"))


@lru_cache(maxsize=1)
def _pulse_library():
    path = Path(__file__).with_name("pulse_derivative_budgets.py")
    if not path.is_file():
        path = Path(__file__).parent/"pulse-norms"/"pulse_derivative_budgets.py"
    spec = importlib.util.spec_from_file_location("registered_stress_pulse_budgets", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rational(iv, value):
    return iv.mpf(value.numerator)/iv.mpf(value.denominator)


def _upper(interval):
    if interval == 0:
        return 0.0
    return nextafter(float(interval.b), inf)


def _pair(interval):
    if interval == 0:
        return [0.0, 0.0]
    return [nextafter(float(interval.a), -inf), nextafter(float(interval.b), inf)]


def stress_tail_bounds(source, eta, cutoff):
    """Return rigorous constructive omitted-band envelopes, in stated units.

    q, q_prime, q_second: errors in q=a^2 deltaQ/epsilon and its conformal
    derivatives. rho,p: errors in a^4 delta(rho,p)/epsilon. current: error in
    a^2 delta j/epsilon=q+Q0*f/4 for b=1,r=2. Q0_difference is the physical
    reference variance Q0-Q0K, not divided by the source amplitude.

    Bounds use exact root isolation/total variation and directed interval
    arithmetic. They do not certify time/momentum quadrature or mode error.
    """
    time = Fraction(str(eta))
    if time not in OBSERVATIONS:
        raise ValueError("Time is outside the frozen six-observation set")
    k_fraction = Fraction(str(cutoff))
    if k_fraction <= 0:
        raise ValueError("The cutoff must be positive")
    from mpmath.ctx_iv import MPIntervalContext
    iv = MPIntervalContext()
    iv.dps = INTERVAL_DECIMAL_PRECISION
    pulse = _pulse_library().pulse_interval_data(iv, source, time)
    f, f1, f2 = pulse["jets"][:3]
    norms = pulse["budgets"]
    K = _rational(iv, k_fraction)
    a = -1/_rational(iv, time)
    L = a
    M2 = 2*a*a
    radial = iv.sqrt(K*K+M2)
    v = K/radial
    one_minus_v = M2/(radial*(radial+K))
    one_minus_v3 = one_minus_v*(1+v+v*v)
    second_contact_factor = one_minus_v*(3*v**4+3*v**3+v*v+v+1)
    C = 1/(8*iv.pi**2)
    basic = [(3*abs(pulse["jets"][j])*M2/2+norms[j]/4)/(16*iv.pi**2*K*K)
             for j in range(3)]
    q0 = basic[0]
    q1 = basic[1]+C*L*abs(f)*one_minus_v3
    q2 = basic[2]+C*(2*L*abs(f1)*one_minus_v3+L*L*abs(f)*second_contact_factor)
    # Exact positive reference tail, written in a stable factored form.
    reference_tail = (one_minus_v**2*(v+2)*(3*v**3+3*v*v+10*v+4)
                      /(96*iv.pi**2*(v+1)))
    density_contact = L*(3*L*f-f1)/(96*iv.pi**2)
    density_contact_tail = density_contact*one_minus_v3
    # r^((n-3)/2)*(J_n(infinity)-J_n(K/a)), with r=2.
    moment5 = one_minus_v3/3
    moment7 = one_minus_v**2*(3*v**3+6*v*v+4*v+2)/15
    moment9 = one_minus_v**3*(15*v**4+45*v**3+48*v*v+24*v+8)/105
    anomaly_tail = ((f2-12*L*L*f)*moment5/16
                    -(30*L*L*f+10*L*f1)*moment7/32
                    +70*L*L*f*moment9/64)/(2*iv.pi**2)
    pressure_contact_tail = (anomaly_tail+density_contact_tail)/3
    rho = (3*L*L*q0+L*q1)/2+a*a*abs(f)*reference_tail/2+abs(density_contact_tail)
    pressure = ((q2+3*L*q1+3*L*L*q0)/6+a*a*abs(f)*reference_tail/6
                +abs(pressure_contact_tail))
    current = q0+abs(f)*reference_tail/4
    bounds = {"q":q0, "q_prime":q1, "q_second":q2,
              "rho":rho, "p":pressure, "current":current,
              "Q0_difference":reference_tail,
              "anomaly":abs(anomaly_tail)}
    result = {name:_upper(value) for name,value in bounds.items()}
    result.update({
        "source":source, "eta":str(eta), "K":str(cutoff),
        "normalizations":{
            "q":"a^2 deltaQ/epsilon; q_prime and q_second are conformal derivatives",
            "rho":"a^4 delta_rho/epsilon", "p":"a^4 delta_p/epsilon",
            "current":"a^2 delta_j/epsilon; b=1,r=2,delta_x=br delta_phi",
            "Q0_difference":"physical reference variance Q0-Q0K",
            "anomaly":"a^4 |delta_A_infinity-delta_A_K|/epsilon"},
        "interval_dps":INTERVAL_DECIMAL_PRECISION,
        "source_jet_intervals_over_epsilon":{f"f{j}":_pair(value) for j,value in enumerate(pulse["jets"])},
        "derivative_budget_intervals":{f"N{j}":_pair(value) for j,value in enumerate(norms)},
        "basic_variance_functional_tail_upper":{f"f{j}":_upper(value) for j,value in enumerate(basic)},
        "density_local_tail_interval":_pair(density_contact_tail),
        "pressure_local_tail_interval":_pair(pressure_contact_tail),
        "anomaly_tail_interval":_pair(anomaly_tail),
        "pulse_details":pulse.get("details",{}),
        "uncertainty_scope":"Constructive analytic omitted-momentum-band bound only; quadrature, mode integration, derivative implementation and floating-point evaluation of responses need separate evidence."
    })
    return result
