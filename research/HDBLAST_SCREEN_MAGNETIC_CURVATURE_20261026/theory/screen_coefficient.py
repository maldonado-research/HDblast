"""Registered-coordinate screen H'' coefficient, exact arithmetic compatible.

Inputs m,b have the actual_first_jets.actual_jet dictionaries' conventions.
The function returns the coefficient of Tplus-Tminus; ux is the equal TT
second-response x derivative. No interval comparisons or sign claims here.
"""
from fractions import Fraction as Q


def quadratic_F(m,b):
    rho=b['rho']; R=rho*b['c']; theta=b['theta']
    return ((m['D']+4*m['M']-m['L'])*m['Mt']-2*m['M']*m['Lt'])/(4*rho) \
        -m['E']*(m['D']+m['L'])/(2*R) \
        +(m['C']/R-m['My']/2)*(R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2+m['B']/2)


def coefficients(m,b,ux):
    R=b['rho']*b['c']
    F=quadratic_F(m,b)
    quadratic=4*F/R
    linear=Q(3,2)*ux/(b['ell']*b['rho']**2*b['c'])
    return {'F':F,'quadratic':quadratic,'linear':linear,'plus':linear+quadratic,'minus':-(linear+quadratic)}
