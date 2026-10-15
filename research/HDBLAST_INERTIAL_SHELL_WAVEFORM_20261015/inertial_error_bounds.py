"""Rational consequences of the independently proved characteristic estimate.

Coefficient bounds and analytic solution/support hypotheses are inputs proved
or stated separately. This evaluator alone is not a PDE existence proof.
"""
from fractions import Fraction as F


def constants():
    ell, m, M, B, gamma, K, a0 = F(1, 10**6), F(1, 500), F(102000), F(250), F(250), F(3*10**7), F(35, 4)
    T = 3 * ell
    U, r, s = 1+T, M*T*T/2, B*T*T/(2*m)
    determinant = (1-r)*(1-s)-T*r/m
    if not (r < 1 and s < 1 and determinant > 0):
        raise ValueError('positive majorant conditions failed')
    A = U*(1-s+T/m)/determinant
    V = U/(m*determinant)
    R0 = F(3, 8)*(K+a0*T)*ell**2+F(3, 4)*M*ell**2
    Itail = R0+M*V*T*T/4+M*M*A*T**3/6
    J1 = F(3, 4)*M*ell**2+T*R0+M*V*T**3/12+M*M*A*T**4/24
    Ev = (J1/m+a0*V*T**4/4)/(1-s)
    global_error = ell*(M*ell+Itail+(1+gamma*T)*Ev)
    tail_error = ell*(Itail+(1+gamma*T)*Ev)
    eta, tail, zeta = F('1.1e-7'), F('1.7e-10'), F('1.100002e-7')
    if not (global_error < eta and tail_error < tail and eta**2+tail**2 < zeta**2):
        raise ValueError('published error rounding failed')
    return dict(ell=ell,T=T,m=m,M=M,B=B,Gamma=gamma,Kz=K,a0=a0,U=U,r=r,s=s,
                determinant=determinant,A=A,V=V,R0=R0,incoming_tail=Itail,
                history_integral=J1,boundary_error_velocity=Ev,
                global_error=global_error,tail_error=tail_error,
                published_global=eta,published_tail=tail,published_norm_error=zeta)


def cell_error(left, right, q_upper, bounds=None):
    """Uniform H-error on one closed cell; Q_upper encloses its complete Q range."""
    left, right, q_upper = F(left), F(right), F(q_upper)
    if not (0 <= left <= right <= 3 and 0 <= q_upper <= 1):
        raise ValueError('invalid continuous error cell')
    if right <= 1:
        return F(0)
    b = constants() if bounds is None else bounds
    t = b['ell']*right
    incoming = b['M']*t*q_upper/2+b['R0']+b['M']*b['V']*t*t/4+b['M']**2*b['A']*t**3/6
    adaptive = b['ell']*(incoming+(1+b['Gamma']*b['T'])*b['boundary_error_velocity'])
    return min(adaptive, b['global_error'])
