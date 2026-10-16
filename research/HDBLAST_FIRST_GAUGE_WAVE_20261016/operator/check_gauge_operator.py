#!/usr/bin/env python3
"""Exact all-point operator identities and rational size certificate; stdlib only.

This proves finite algebraic identities, not a solution of a PDE. Matrix caps
are implications of the explicitly reported whole-collar background box.
No assert statements, floating-point arithmetic, or sampled point tests.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys


class P:
    """Sparse rational Laurent polynomials in named independent symbols."""
    def __init__(self, data=0):
        if isinstance(data, P):
            self.d = dict(data.d)
        elif isinstance(data, dict):
            self.d = {m: Q(c) for m, c in data.items() if c}
        else:
            self.d = {(): Q(data)} if data else {}

    @staticmethod
    def var(name):
        return P({((name, 1),): 1})

    def __add__(self, other):
        out = dict(self.d)
        for m, c in P(other).d.items():
            out[m] = out.get(m, Q(0)) + c
        return P(out)

    __radd__ = __add__

    def __neg__(self):
        return P({m: -c for m, c in self.d.items()})

    def __sub__(self, other):
        return self + (-P(other))

    def __rsub__(self, other):
        return P(other) - self

    def __mul__(self, other):
        out = {}
        for m, c in self.d.items():
            for n, d in P(other).d.items():
                powers = dict(m)
                for name, e in n:
                    powers[name] = powers.get(name, 0) + e
                key = tuple(sorted((name, e) for name, e in powers.items() if e))
                out[key] = out.get(key, Q(0)) + c*d
        return P(out)

    __rmul__ = __mul__

    def __pow__(self, n):
        if n < 0:
            if len(self.d) != 1:
                raise ValueError('negative power requires a Laurent monomial')
            (m, c), = self.d.items()
            return P({tuple((name, e*n) for name, e in m): c**n})
        out = P(1)
        for _ in range(n):
            out *= self
        return out

    def derivative(self, name):
        out = P()
        for m, c in self.d.items():
            powers = dict(m)
            e = powers.get(name, 0)
            if e:
                powers[name] = e-1
                key = tuple(sorted((v, k) for v, k in powers.items() if k))
                out += P({key: c*e})
        return out

    def sub(self, mapping):
        out = P()
        for m, c in self.d.items():
            term = P(c)
            for name, e in m:
                term *= P(mapping.get(name, P.var(name)))**e
            out += term
        return out


def vars_(names):
    return tuple(P.var(s) for s in names.split())


def require(condition, message):
    if not condition:
        raise ValueError(message)


MUTATIONS = {
    'source_sign': 'longitudinal_radial_defect',
    'radial_drift': 'inward_radial_drift',
    'tangential_drift': 'inward_tangential_drift',
    'angular_curvature': 'desitter_angular_component',
    'missing_tangential_coupling': 'warped_tangential_component',
    'reflection_drift': 'reflection_B_parity',
    'omit_master_time_factor': 'reconstruction_source',
    'contraction_underbound': 'contractive_cap',
}


def build(injection=None):
    identities = []
    controls = []

    def zero(name, lhs, rhs, error=None):
        residual = P(lhs)-rhs
        if injection and MUTATIONS[injection] == name:
            require(error is not None, 'missing negative-control witness')
            residual += error
        require(not residual.d, 'IDENTITY_FAILED:'+name)
        identities.append(name)
        if error is not None:
            require(bool((P(lhs)-rhs+error).d), 'ineffective control:'+name)
            controls.append(name)

    rho, a, ay, u, uy = vars_('rho a ay u uy')
    # DC(h)_y from div(h)_y - one-half trace derivative, longitudinal h.
    div_long = -4*uy + 4*a*(-4*u) - a*(8*u)
    zero('longitudinal_radial_defect', div_long-2*uy,
         -6*(uy+4*a*u), 12*(uy+4*a*u))
    ua = P.var('ua')
    zero('longitudinal_tangential_defect', 2*ua-2*ua, 0)

    # Warp tensor-divergence calculation, with all barred jets independent.
    X, Xy, Xyy, boxX, divV, divVy = vars_('X Xy Xyy boxX divV divVy')
    hyy, trtan = 2*Xy, 8*a*X+2*divV
    divh_y = 2*Xyy + rho**-2*boxX+divVy+4*a*hyy-a*trtan
    halftrace_y = Xyy+4*ay*X+4*a*Xy+divVy
    zero('warped_radial_component', divh_y-halftrace_y,
         Xyy+4*a*Xy+rho**-2*boxX-2*a*divV-(4*ay+8*a*a)*X)
    Xa, Xya, Vay, Vayy, boxVa, gradaDiv, Va = vars_(
        'Xa Xya Vay Vayy boxVa gradaDiv Va')
    hya = Xa+rho*rho*Vay
    partial_y_hya = Xya+2*a*rho*rho*Vay+rho*rho*Vayy
    divbar_tan = 2*a*Xa+boxVa+gradaDiv+3*Va
    divh_a = partial_y_hya-a*hya+divbar_tan+5*a*hya
    zero('warped_tangential_component', rho**-2*(divh_a-(Xya+4*a*Xa+gradaDiv)),
         Vayy+6*a*Vay+rho**-2*(boxVa+3*Va+2*a*Xa), 2*a*rho**-2*Xa)

    # Four-dimensional tensor-divergence calculation on de Sitter.
    H, c2, lam, b, bt, btt, d, dt, dtt = vars_(
        'H c2 lam b bt btt d dt dtt')
    S = 1-H*H
    # Use c2^{-1} before applying c2^{-1}=1-H^2 to the final identity.
    htt, hti, trspace = -2*bt, c2*dt-b, 6*H*b-2*lam*d
    div_t = 2*btt-lam*c2**-1*hti-3*H*htt-H*trspace
    halftrace_t = btt+3*S*b+3*H*bt-lam*dt
    Ct = -(div_t-halftrace_t)
    # c2 appears only as c2^{-1}; multiply by c2 and set c2*S=1
    # equivalently substitute inverse powers directly in a small normalizer.
    def desitter_normal(poly):
        out = P()
        for mon, coefficient in P(poly).d.items():
            powers = dict(mon)
            e = powers.pop('c2', 0)
            require(e <= 0, 'uncancelled positive c2 power')
            out += P({tuple(sorted(powers.items())): coefficient})*S**(-e)
        return out
    time_expected = -btt-3*H*bt+(3+3*H*H-lam*S)*b-2*lam*H*d
    zero('desitter_time_component', desitter_normal(Ct), time_expected)
    dt_hti = c2*dtt+2*H*c2*dt-bt
    # D^j Hess_ji Y=(2-lambda) D_i Y; four connection terms are retained.
    div_i = -dt_hti+H*hti+2*H*b+2*(2-lam)*d-4*H*hti
    halftrace_i = bt+3*H*b-lam*d
    angular_expected = -dtt-5*H*dt+(4-lam)*S*d+2*H*S*b
    zero('desitter_angular_component', desitter_normal(c2**-1*(div_i-halftrace_i)),
         angular_expected, 2*S*d)

    # Independent scalar-gradient and Killing checks on the component operator.
    chi, chit, chitt, chittt = vars_('chi chit chitt chittt')
    rates = {'H': S, 'chi': chit, 'chit': chitt, 'chitt': chittt}
    def Dt(poly):
        return sum((P(poly).derivative(name)*rate for name, rate in rates.items()), P())
    gradb, gradd = -chit, S*chi
    gradient_jets = {'b':gradb,'bt':Dt(gradb),'btt':Dt(Dt(gradb)),
                     'd':gradd,'dt':Dt(gradd),'dtt':Dt(Dt(gradd))}
    scalar = -chitt-3*H*chit-lam*S*chi+6*chi
    zero('desitter_scalar_gradient_time',time_expected.sub(gradient_jets),-Dt(scalar))
    zero('desitter_scalar_gradient_angular',angular_expected.sub(gradient_jets),S*scalar)
    boost = {'lam':3,'b':1,'bt':0,'btt':0,'d':H,'dt':S,'dtt':-2*H*S}
    zero('desitter_boost_time',time_expected.sub(boost),0)
    zero('desitter_boost_angular',angular_expected.sub(boost),0)

    fz, fzz = vars_('fz fzz')
    fyy = rho**-2*fzz+a*rho**-1*fz
    fy = -rho**-1*fz
    zero('inward_radial_drift',rho*rho*(fyy+4*a*fy),fzz-3*a*rho*fz,a*rho*fz)
    zero('inward_tangential_drift',rho*rho*(fyy+6*a*fy),fzz-5*a*rho*fz,a*rho*fz)
    uz = P.var('uz')
    zero('inward_source',-6*rho*rho*(-rho**-1*uz+4*a*u),6*rho*uz-24*a*rho*rho*u)

    # Obtain RHS from the original y amplitude equations, not from matrices.
    f, ft, ftt, bz, bzz, dz, dzz = vars_('f ft ftt bz bzz dz dzz')
    first_y_equation = (fzz-3*a*rho*fz-ftt-3*H*ft-lam*S*f
         -2*a*rho*rho*(bt+3*H*b-lam*d)-rho*rho*(4*ay+8*a*a)*f
         -6*rho*rho*(-rho**-1*uz+4*a*u))
    second_y_equation = (bzz-5*a*rho*bz-btt-3*H*bt+
         (3+3*H*H-lam*S)*b-2*lam*H*d-2*a*ft)
    third_y_equation = (dzz-5*a*rho*dz-dtt-5*H*dt+(4-lam)*S*d+2*H*S*b+2*a*S*f)
    A = [[-3*H,-2*a*rho*rho,0],[-2*a,-3*H,0],[0,0,-5*H]]
    B = [[-3*a*rho,0,0],[0,-5*a*rho,0],[0,0,-5*a*rho]]
    D = [[-lam*S-rho*rho*(4*ay+8*a*a),-6*a*rho*rho*H,2*lam*a*rho*rho],
         [0,3+3*H*H-lam*S,-2*lam*H],[2*a*S,2*H*S,(4-lam)*S]]
    force = [6*rho*uz-24*a*rho*rho*u,0,0]
    U, Ut, Uz, Uzz, Utt = [f,b,d],[ft,bt,dt],[fz,bz,dz],[fzz,bzz,dzz],[ftt,btt,dtt]
    originals = [first_y_equation,second_y_equation,third_y_equation]
    ell = P.var('ell')
    for i in range(3):
        rhs = sum((P(A[i][j])*Ut[j]+P(B[i][j])*Uz[j]+P(D[i][j])*U[j] for j in range(3)),P())+force[i]
        zero('matrix_amplitude_'+str(i), originals[i], Uzz[i]-Utt[i]+rhs)
        scaled = {name: ell**-1*P.var(name+'_scaled') for name in ('ft','bt','dt','fz','bz','dz')}
        for name in ('ftt','btt','dtt','fzz','bzz','dzz'):
            scaled[name] = ell**-2*P.var(name+'_scaled')
        target_time = [P.var(n+'_scaled') for n in ('ft','bt','dt')]
        target_space = [P.var(n+'_scaled') for n in ('fz','bz','dz')]
        scaled_rhs = sum((ell*P(A[i][j])*target_time[j]+
                          ell*P(B[i][j])*target_space[j]+
                          ell*ell*P(D[i][j])*U[j] for j in range(3)),P())+ell*ell*force[i]
        zero('scaled_matrix_amplitude_'+str(i),ell*ell*rhs.sub(scaled),scaled_rhs)

    R, h, w, wz, theta = vars_('R h w wz theta')
    Rz = rho*(Q(9,2)*a-h)*R
    zero('reconstruction_source',6*rho*theta*(Rz*w+R*wz)-24*a*rho*rho*R*theta*w,
         6*rho*R*theta*wz+6*rho*rho*R*theta*(a*Q(1,2)-h)*w,
         (1-theta)*(6*rho*R*wz+6*rho*rho*R*(a*Q(1,2)-h)*w))
    zero('initial_radial_acceleration',force[0].sub({'uz':-rho*uy}),-6*rho*rho*(uy+4*a*u))

    # Covariant lower tangential boundary row and reflected first-order drift.
    v, vy = vars_('v vy')
    zero('boundary_covector_row',(2*a*rho*rho*v+rho*rho*vy)-2*a*rho*rho*v,rho*rho*vy)
    zeta = P.var('zeta')
    zero('reflection_B_parity',(-B[0][0])*(-zeta),B[0][0]*zeta,2*P(B[0][0])*zeta)

    # Exact conditional box consequences. Actual 8ell-collar authentication
    # belongs to the separate coefficient certificate, not this algebra file.
    L,T = Q(1,10**6),Q(3)
    rholo,rhohi,amax,pmax,Hmax=Q(129),Q(131),Q(14,125),Q(7,50000),T*L
    aymax=1/rholo**2+pmax**2/3
    Acap=L*max(3*Hmax+2*amax*rhohi**2,2*amax+3*Hmax,5*Hmax)
    Bcap=L*5*amax*rhohi
    Drows=[L**2*(8+rhohi**2*(4*aymax+8*amax**2)+6*amax*rhohi**2*Hmax+16*amax*rhohi**2),
           L**2*(3+3*Hmax**2+8+16*Hmax),
           L**2*(2*amax+2*Hmax+4)]
    Dcap=max(Drows)
    qcap=T*(Acap+Bcap+T*Dcap)
    claimed=Q(12,1000)
    if injection=='contraction_underbound':
        claimed=Q(1,1000)
    require(qcap<claimed,'IDENTITY_FAILED:contractive_cap')
    require(qcap>Q(1,1000),'ineffective contraction control')
    controls.append('contractive_cap')
    require(Acap<Q(3845,10**6),'A cap exceeded')
    require(Bcap<=Q(7336,10**8),'B cap exceeded')
    require(Dcap<Q(33,10**9),'D cap exceeded')
    return {
      'status':'PASS_EXACT_ALL_POINT_ALGEBRA_AND_CONDITIONAL_RATIONAL_CAPS',
      'scope':'Algebra identities and consequences of stated whole-collar box; no PDE existence or numerical trajectory is asserted by this checker.',
      'identity_count':len(identities),'identities':identities,
      'mutation_witness_count':len(controls),'mutation_witnesses':controls,
      'parameters_exact':{'ell':str(L),'T_scaled':str(T),'lambda':'8','z_collar_required':'1/125000'},
      'conditional_background_box':{'rho':['129','131'],'absolute_a_upper':'14/125','absolute_p_upper':'7/50000','ay_identity':'-rho^-2-p^2/3','H_bound':'abs(tanh(ell*x)) <= 3ell','required_domain':'0<=x<=3, r>=0, x+r<=8'},
      'exact_caps':{'absolute_ay':str(aymax),'A_row_sum':str(Acap),'B_row_sum':str(Bcap),'D_row_sum':str(Dcap),'characteristic_contraction':str(qcap),'sharper_characteristic_contraction':str(T*(Acap+Bcap+T*Dcap/2))},
      'readable_strict_caps':{'A':'3845/1000000','B_nonstrict':'7336/100000000','D':'33/1000000000','q':'3/250'},
      'source_convention':{'R':'alpha/rho^2','h':'U_potential_prime/p','R_z':'rho*(9*a/2-h)*R','theta':'cosh(t)^(-3/2)','u':'R*theta*w','scaled_radial_source':'6*ell*rho*R*theta*w_r+6*ell^2*rho^2*R*theta*(a/2-h)*w'},
      'parity':'J=diag(-1,1,1); A(-r)=J A(r) J, B(-r)=-J B(r) J, D(-r)=J D(r) J, F(-r)=J F(r); jumps allowed',
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject',choices=sorted(MUTATIONS))
    args=parser.parse_args()
    require(not args.output.exists(),'output path already exists')
    result=build(args.inject)
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    data=json.dumps(result,indent=2,sort_keys=True)+'\n'
    with args.output.open('x') as handle:
        handle.write(data)
    print(result['status'])


if __name__=='__main__':
    try:
        main()
    except (ValueError,OSError) as exc:
        print(str(exc),file=sys.stderr)
        sys.exit(2)
