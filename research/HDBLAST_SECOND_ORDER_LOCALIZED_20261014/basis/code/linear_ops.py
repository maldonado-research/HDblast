"""First-order (linearized) constraint operators P_H and P_M at the
time-symmetric background, built from the same eps-jet engine."""
import sympy as sp
from fields import constraints

from engine import J3

def PH_PM(h=None, psi=0, pi=None, varpi=0):
    z = sp.zeros(4)
    J3.FIRST_ONLY = True
    try:
        Hm, M = constraints(h if h is not None else z, psi,
                            pi if pi is not None else z, varpi)
    finally:
        J3.FIRST_ONLY = False
    return Hm.b, [M[i].b for i in range(4)]
