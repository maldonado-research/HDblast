"""Initial five means for the explicit 2026-11-01 B-mean preparation.

Pure algebra adapter; no old initial-data oracle or evolution is imported.
Arguments can be exact numbers or interval objects supporting arithmetic.
All A, B, QM are second parameter derivatives, not quadratic coefficients.
The prime/y suffix always means outward background proper radial derivative.
"""
from fractions import Fraction

ORDER = ("K0", "L0", "N0", "J0", "w0")
GRAM = (1, 2, 1, 3, 1)
ELL = Fraction(1, 1000000)


def radial_metric(f, fy, xi, xiy, xiyy, a, ay):
    """Return A,A_y,B,B_y for h=f gamma0+Lie_(xi d_y) gamma0."""
    return f + 2*xiy, fy + 2*xiyy, f + 2*a*xi, fy + 2*ay*xi + 2*a*xiy


def initial_means(*, rho, a, u, ut, uy, qm, A, Ay, B, By, ell=ELL):
    """Return U0, d_x U0, d_r U0 and physical/gauge intermediate means.

    x=t/ell, r=z/ell, dy=-rho dz. Frame: e_y=d_y, e_0=rho^-1 d_t,
    e_a=(rho cosh t)^-1 E_a. Angular probability average gives <Y^2>=1/6.
    N0 is the covariant timelike frame metric component, with no sign flip.
    """
    zero = u*0
    us = u*ut/rho
    h_gauge = Ay/2 - 3*By/2 + 4*a*A - 3*a*B
    U0 = (A, zero, -u*u/3, B, zero)
    V0 = tuple(ell*v for v in (
        2*rho*qm/(3*a) + 32*u*ut/3,
        rho*(h_gauge - 13*u*uy/3 - a*u*u/3),
        -2*rho*qm/(3*a) + 4*u*ut/3,
        -4*u*ut/3,
        zero,
    ))
    U0_r = tuple(-ell*rho*v for v in (Ay, zero, -2*u*uy/3, By, zero))
    return {
        "order": ORDER,
        "U0": U0, "V0": V0, "U0_r": U0_r,
        "K2_rr": qm/(3*a) + 6*us, "K2_angular": -us,
        "dot_N2": rho*rho*(qm/(3*a) - us),
        "dot_beta_y2": rho*rho*(h_gauge - 37*u*uy/3 - 97*a*u*u/3),
    }
