"""Intrinsic shell scalar Hessian; second parameter derivative convention.

Inputs may be exact rationals or outward intervals. Background shell scalar
must be tangentially constant. eta_mat=v Y includes the embedding pullback.
"""
from fractions import Fraction as Q

Y2=(Q(1,6),Q(1,2),Q(1))
GRAD2=(Q(4,3),Q(2),Q(-4))

def coefficients(v,vt,rho,c):
    """I'' coefficients in (1,Z,F), not the epsilon-squared coefficient."""
    return tuple(2*(-vt**2*y+v**2*g/c**2)/rho**2 for y,g in zip(Y2,GRAD2))

def mean(v,vt,rho,c):
    return (-vt**2+8*v**2/c**2)/(3*rho**2)

def pole(vt,rho):
    return -2*vt**2/rho**2

def angular_equator(v,rho,c):
    return 8*v**2/(rho**2*c**2)

def nodal_ratio(v,vt,c):
    """Use only when v and vt are certified nonzero; ratio for Y^2/S."""
    return 4*v**2/(c**2*vt**2+4*v**2)

def material_variation(chi,p,zeta):
    return chi+p*zeta
