"""Pure second harmonic-gauge adapter; see CONDITIONED_GAUGE.md.

Separate from the frozen original mean_initial.py. Only arithmetic on supplied
actual initial profiles is performed; no old initial oracle is imported.
"""
from fractions import Fraction

ORDER = ("K0", "L0", "N0", "J0", "w0")
GRAM = (1, 2, 1, 3, 1)
ELL = Fraction(1, 1000000)


def conditioned_initial_means(*, rho, a, u, ut, uy, qm, f, fy,
                              xi, xiy, p, py, ell=ELL):
    """Return second derivative means, scaled time velocity and radial slope.

    Generator initially Z^y=-xi, Z^t=0, Z^y_t=0, Z^t_t=a xi. Its subsequent
    homogeneous monopole gauge wave has Z^y=0, Z^t_y=0 on the shell.
    x=t/ell, r=z/ell, dy=-rho dz. Inputs are exact numbers or intervals.
    """
    zero=u*0
    U0=(f,zero,-u*u/3,f,-p*xi)
    V0=tuple(ell*v for v in (
        2*rho*qm/(3*a)+32*u*ut/3,
        rho*(-fy+a*f-13*u*uy/3-a*u*u/3),
        -2*rho*qm/(3*a)+4*u*ut/3,
        -4*u*ut/3,
        zero,
    ))
    U0_r=tuple(-ell*rho*v for v in (fy,zero,-2*u*uy/3,fy,-py*xi-p*xiy))
    return {"order":ORDER,"U0":U0,"V0":V0,"U0_r":U0_r}
