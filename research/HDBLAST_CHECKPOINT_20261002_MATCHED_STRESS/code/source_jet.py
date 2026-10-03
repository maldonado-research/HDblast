"""Exact analytic source derivatives; importing performs no source sampling."""
from __future__ import annotations

import math


# B^(n)=B*P_n(u)/(1-u^2)^(2n), with exactly derived integer coefficients.
# P_(n+1)=(1-u²)² P_n' + 4n*u*(1-u²)*P_n - 2u*P_n.
POLYNOMIALS = (
    (1,),
    (-2, 0),
    (6, 0, 0, 0, -2),
    (-24, 0, -12, 0, 40, 0, -12, 0),
    (120, 0, 180, 0, -528, 0, 232, 0, 24, 0, -12),
    (-720, 0, -2160, 0, 6120, 0, -2400, 0, -2112, 0, 1360, 0, -120, 0),
)


def horner(coefficients, x):
    result = 0.0
    for coefficient in coefficients:
        result = result*x+coefficient
    return result


def source_jet_over_epsilon(source, eta):
    """Return f and five conformal-time derivatives of f=s/epsilon."""
    if source not in ("positive_B", "signed_uB"):
        raise ValueError("Source is not in the prospective registration")
    u = eta+4.0
    if abs(u) >= 1.0:
        return [0.0]*6
    denominator = (1.0-u)*(1.0+u)
    bump = math.exp(-u*u/denominator)
    if bump == 0.0:
        # This also avoids an underflowed zero multiplied by a large reciprocal.
        return [0.0]*6
    derivatives = [bump*horner(coefficients, u)/denominator**(2*n)
                   for n, coefficients in enumerate(POLYNOMIALS)]
    if source == "positive_B":
        return derivatives
    return [u*derivatives[n]+(n*derivatives[n-1] if n else 0.0) for n in range(6)]
