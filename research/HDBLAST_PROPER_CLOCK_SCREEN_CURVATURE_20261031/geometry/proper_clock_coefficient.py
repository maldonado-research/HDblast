"""Proper-time Gaussian-label screen coefficient, exact/interval compatible.

All suffix derivatives are ordinary physical t,y derivatives. ``xi_t`` and
``xi_y`` are coordinate position variations multiplying Y, at fixed proper
time and initial labels, not vector components in an orthonormal frame.
No family existence or numerical enclosure is supplied by these functions.
"""
from fractions import Fraction as Q


def radial_quadratic_F(m, b):
    rho = b['rho']; R = rho*b['c']; theta = b['theta']
    eta = R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2-m['B']/2
    return ((-m['A']+4*m['M']-m['L'])*m['My']-2*m['M']*m['Ly'])/4 \
        +m['C']*(m['A']-m['L'])/(2*R) \
        +(m['E']/R-m['Mt']/(2*rho))*eta


def time_quadratic_F(m, b):
    rho = b['rho']; R = rho*b['c']; theta = b['theta']
    return ((m['D']+4*m['M']-m['L'])*m['Mt']-2*m['M']*m['Lt'])/(4*rho) \
        -m['E']*(m['D']+m['L'])/(2*R) \
        +(m['C']/R-m['My']/2)*(R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2+m['B']/2)


def first_weyl_invariants(m, b):
    """At=g+e0, Ay=ey-g, K=e0/ey mixed Hessian coefficient."""
    rho = b['rho']; R = rho*b['c']; theta = b['theta']; a = b['a']
    At = -(m['Mtt']+theta*m['Mt'])/(2*rho**2) \
        +m['Et']/(rho*R)-(m['D']+m['L'])/(2*R**2)
    Ay = -(m['Myy']+a*m['My'])/2+m['Cy']/R \
        -theta*m['Mt']/(2*rho**2)+theta*m['E']/(rho*R) \
        +(m['L']-m['A'])/(2*R**2)
    K = -(m['Myt']+theta*m['My'])/(2*rho) \
        +(m['Ct']+theta*m['C'])/(2*rho*R)+m['Ey']/(2*R)-m['B']/(2*R**2)
    return {'At': At, 'Ay': Ay, 'K': K}


def coefficients(m, b, ux, ur, cosh_chi, sinh_chi, xi_t, xi_y):
    """Return the full second-derivative coefficient of Tplus-Tminus."""
    rho = b['rho']; R = rho*b['c']; c = cosh_chi; s = sinh_chi
    Ft = time_quadratic_F(m,b); Fy = radial_quadratic_F(m,b)
    w = first_weyl_invariants(m,b)
    baseline = 4*(c*Ft+s*Fy)/R
    tilt = -4*(rho*xi_t*(c*w['At']+s*w['K'])
               +xi_y*(s*w['Ay']+c*w['K']))/R
    linear = Q(3,2)*(c*ux-s*ur)/(b['ell']*rho**2*b['c'])
    total = linear+baseline+tilt
    return {'Ft':Ft,'Fy':Fy,**w,'linear_TT':linear,
            'boosted_quadratic':baseline,'geodesic_tilt':tilt,
            'plus':total,'minus':-total}
