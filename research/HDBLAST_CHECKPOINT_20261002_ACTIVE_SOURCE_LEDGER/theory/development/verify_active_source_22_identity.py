#!/usr/bin/env python3
"""Pure active Ward algebra and artificial polynomial-source moment controls.

This verifier has no file input option and does not read saved physical arrays.
Its numerical examples use explicitly fabricated rational constants only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import mpmath as m
import sympy as s


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def symbolic():
    checks, mutations = [], []
    def zero(name, expr):
        require(s.cancel(s.expand(expr)) == 0, name)
        checks.append(name)
    def reject(name, expr):
        require(s.cancel(s.expand(expr)) != 0, 'Mutation survived: ' + name)
        mutations.append(name)
    L, k, ep = s.symbols('L k epsilon', positive=True)
    X, Y, Z, T, f = s.symbols('X Y Z T f', real=True)
    h = s.symbols('h0:6', real=True)
    g = 4*L**2*h[0]-2*L*h[1]-h[2]
    def D(expr):
        return s.expand(L**2*s.diff(expr,L)+Z*s.diff(expr,X)+T*s.diff(expr,Y)
            +(-2*k*T-ep*f)*s.diff(expr,Z)+2*k*Z*s.diff(expr,T)
            +sum(h[i+1]*s.diff(expr,h[i]) for i in range(5)))
    Rm=((2*k*k+3*L*L)*X-k*T-L*Z)/(2*k*ep)
    Pm=((2*k*k/3-L*L)*X-k*T-L*Z)/(2*k*ep)
    R0=(2*k*k+3*L*L)/(4*k)
    P0=(2*k*k/3-L*L)/(4*k)
    CR=L*h[1]/(2*k)-2*k*h[0]-2*L*L*h[0]/k
    CP=L*h[1]/(2*k)-2*k*h[0]/3
    G=(3*L*L*X-L*Z)/(2*k*ep)
    cR=X-T/(2*k)
    Fm=L*(Rm-3*Pm)
    FC=L*(CR-3*CP)-3*h[1]*(R0+P0)
    R,P=Rm+CR,Pm+CP
    F=L*(R-3*P)-3*h[1]*(R0+P0)
    zero('unrestricted real-forcing invariant',D(cR))
    zero('mode stress constant plus primitive',Rm-k*cR/ep-G)
    zero('mode primitive derivative with forcing work',D(G)-Fm-L*f/(2*k))
    zero('direct metric contact derivative with canonical source work',D(CR)-FC+L*g/(2*k))
    zero('full bare Ward including general source mismatch',D(R)-F-L*(f-g)/(2*k))
    zero('primitive differentiates to independently defined direct F',D(G+CR)-F-L*(f-g)/(2*k))
    zero('on-action full Ward',(D(R)-F).subs(f,g))
    zero('direct background Ward',D(R0)-L*(R0-3*P0))
    canonical_energy=2*k*X-T
    zero('unrestricted free canonical energy variation invariant',D(canonical_energy))
    source_hamiltonian=canonical_energy/(2*ep)+g/(4*k)
    zero('linear canonical Hamiltonian external source work',D(source_hamiltonian)-D(g)/(4*k))
    Md,Ma,J,dc,flow,op,dr,ledger=s.symbols('M_discrete M_analytic J_source delta_contact flow operator deltaR Simpson',real=True)
    I=s.symbols('I_continuous')
    Emom=dc+(Md-Ma)*J
    zero('signed active ledger bookkeeping',(dr-ledger)-((dr-I)+(I-ledger)))
    zero('mixed momentum continuous residual decomposition', (flow+dc+(Md-Ma)*J+op)-(flow+Emom+op))
    reject('omit fixed comoving cutoff baseline work',3*h[1]*(R0+P0))
    reject('force sign inversion',L*g/k)
    reject('drop metric D operator from density',D(L*h[1]/(2*k)))
    reject('drop physical volume contacts',D(-4*h[0]*R0))
    reject('drop source mismatch work',L*(f-g)/(2*k))
    reject('project arbitrary retained real invariant',k*cR/ep)
    # Generic subtraction identities are inherited from the fully proved
    # W0/W2/W4 inventory, rather than assumed to arise from mode dynamics.
    SR, SP, dSR, dSP, SRp, dSRp = s.symbols('SR SP dSR dSP SRp dSRp', real=True)
    CRsub = -dSR + 4*h[0]*SR
    CPsub = -dSP + 4*h[0]*SP
    CRsubprime = -dSRp + 4*h[1]*SR + 4*h[0]*SRp
    # Baseline: SR'=L(SR-3SP). Direction: dSR'=L(dSR-3dSP)+h'(SR-3SP).
    wardrules = {SRp:L*(SR-3*SP), dSRp:L*(dSR-3*dSP)+h[1]*(SR-3*SP)}
    Fsub=L*(CRsub-3*CPsub)+3*h[1]*(SR+SP)
    zero('subtraction contact contribution to fixed-band Ward',CRsubprime.subs(wardrules)-Fsub)
    zero('renormalized contact work cancels forcing work',
        (D(CR)+CRsubprime.subs(wardrules))-(FC+Fsub)+L*g/(2*k))
    # Primitive recurrence verified against an arbitrary polynomial source.
    t, lam = s.symbols('eta lambda', nonzero=True)
    A = s.symbols('a0:6', real=True)
    H = sum(A[j]*t**j for j in range(6))
    forcing = 4*H/t**2+2*s.diff(H,t)/t-s.diff(H,t,2)
    quotient, remainder = s.div(s.expand(t*t*forcing),t*t,t)
    zero('polynomial rational forcing inventory',forcing-quotient-remainder/t**2)
    for n in range(1,6):
        # d(eta^n exp(lambda eta)/lambda - n J_(n-1)/lambda) = eta^n exp.
        jm=s.symbols('Jprev')
        derivative=(n*t**(n-1)+lam*t**n)*s.exp(lam*t)/lam-n*t**(n-1)*s.exp(lam*t)/lam
        zero('polynomial phase moment recurrence degree '+str(n),derivative-t**n*s.exp(lam*t))
    zero('first inverse-power phase primitive',s.diff(s.Ei(lam*t),t)-s.exp(lam*t)/t)
    zero('second inverse-power phase primitive',s.diff(lam*s.Ei(lam*t)-s.exp(lam*t)/t,t)-s.exp(lam*t)/t**2)
    return checks, mutations


def synthetic():
    records=[]
    coeff=s.symbols('a0:6')
    t=s.symbols('eta')
    # Every numeric datum here is fabricated; these are NOT the study profiles.
    h_expr=sum(s.Rational(v,d)*t**j for j,(v,d) in enumerate([(1,7),(-2,11),(3,19),(1,23),(-1,29),(1,41)]))
    g_expr=s.expand(4*h_expr/t**2+2*s.diff(h_expr,t)/t-s.diff(h_expr,t,2))
    q,rem=s.div(s.expand(t*t*g_expr),t*t,t)
    qcoeff=[q.coeff(t,j) for j in range(s.degree(q,t)+1)]
    inv1=rem.coeff(t,1);inv2=rem.coeff(t,0)
    def rational(x):
        return m.mpf(int(x.p))/int(x.q)
    for dps in (80,100):
        with m.workdps(dps):
            a,b=-m.mpf(23)/5,-m.mpf(37)/10
            k=m.mpf(13)/7;om=2*k;ep=m.mpf(1)/10000
            ua=m.mpc(m.mpf(2)/17,-m.mpf(3)/31)
            wa=m.mpc(m.mpf(5)/43,m.mpf(7)/47)
            hc=[rational(h_expr.coeff(t,j)) for j in range(6)]
            qc=[rational(v) for v in qcoeff]
            aa,bb=rational(inv1),rational(inv2)
            def hp(x,n=0):
                return sum(hc[j]*m.factorial(j)/m.factorial(j-n)*x**(j-n) for j in range(n,6))
            def g(x):return 4*hp(x)/x**2+2*hp(x,1)/x-hp(x,2)
            def J(x,lam):
                if lam==0:
                    return sum(qc[j]*x**(j+1)/(j+1) for j in range(len(qc)))+aa*m.log(abs(x))-bb/x
                e=m.exp(lam*x);mom=[e/lam]
                for j in range(1,len(qc)):mom.append(x**j*e/lam-j*mom[-1]/lam)
                return sum(qc[j]*mom[j] for j in range(len(qc)))+aa*m.ei(lam*x)+bb*(lam*m.ei(lam*x)-e/x)
            def moments(x):
                M0=J(x,0)-J(a,0)
                ME=m.exp(1j*om*x)*(J(x,-1j*om)-J(a,-1j*om))
                return M0,ME
            def modes(x):
                M0,ME=moments(x);E=m.exp(1j*om*(x-a))
                return ua+(E-1)*wa/(1j*om)-ep*(ME-M0)/(1j*om),E*wa-ep*ME
            def direct(x):
                u,w=modes(x);L=-1/x;h0,h1=hp(x),hp(x,1)
                rm=((2*k*k+3*L*L)*u.real-k*w.imag-L*w.real)/(2*k*ep)
                pm=((2*k*k/3-L*L)*u.real-k*w.imag-L*w.real)/(2*k*ep)
                cr=L*h1/(2*k)-2*k*h0-2*L*L*h0/k
                cp=L*h1/(2*k)-2*k*h0/3
                r0=(2*k*k+3*L*L)/(4*k);p0=(2*k*k/3-L*L)/(4*k)
                return rm+cr,pm+cp,L*(rm+cr-3*(pm+cp))-3*h1*(r0+p0)
            def primitive(x):
                u,w=modes(x);L=-1/x
                return (3*L*L*u.real-L*w.real)/(2*k*ep)+L*hp(x,1)/(2*k)-2*k*hp(x)-2*L*L*hp(x)/k
            M0,ME=moments(b)
            M0_quad=m.quad(g,[a,(a+b)/2,b])
            ME_quad=m.quad(lambda z:m.exp(1j*om*(b-z))*g(z),[a,(a+b)/2,b])
            integral=m.quad(lambda z:direct(z)[2],[a,(a+b)/2,b])
            canon=primitive(b)-primitive(a)
            endpoint=direct(b)[0]-direct(a)[0]
            u,w=modes(b)
            invariant=(u.real-w.imag/om)-(ua.real-wa.imag/om)
            gaps={'nonoscillatory_moment':abs(M0-M0_quad),'oscillatory_moment':abs(ME-ME_quad),
                  'direct_F_integral_vs_primitive':abs(integral-canon),'direct_density_vs_integrand_primitive':abs(endpoint-canon),
                  'arbitrary_invariant_retained':abs(invariant),
                  'u_ODE_derivative':abs(m.diff(lambda z:modes(z)[0],b)-w),
                  'w_ODE_derivative':abs(m.diff(lambda z:modes(z)[1],b)-(1j*om*w-ep*g(b)))}
            limit=m.mpf('1e-55')
            for key,value in gaps.items():require(value<limit,key)
            records.append({'decimal_digits':dps,'gaps':{key:m.nstr(value,30) for key,value in gaps.items()},'criterion':'absolute 1e-55','fabricated_source':'degree5 rational polynomial, no study-source evaluation'})
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    require(not args.output.exists(),'Fresh evidence required')
    checks,mutations=symbolic()
    fabricated=synthetic()
    receipt={'status':'PASS_PURE_SYMBOLIC_AND_SYNTHETIC','identity_count':len(checks),'identities':checks,
             'mutation_count':len(mutations),'mutations_rejected':mutations,'synthetic_passes':fabricated,
             'physical_evaluations':0,'saved_scientific_values_loaded':0,'python_optimization':sys.flags.optimize,
             'sympy_version':s.__version__,'mpmath_version':m.__version__,
             'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'scope':'Generic action/ODE/primitive identities and fabricated rational polynomial source only. Generic full subtraction Ward identities are inherited pinned inputs, not independently rederived in this verifier.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:receipt[key] for key in ['status','identity_count','mutation_count','physical_evaluations']}))
if __name__=='__main__':main()
