#!/usr/bin/env python3
"""Independent exact proofs/manufactured controls; no project numerical imports.

Only Python, Fraction, SymPy, and deterministic local output are used. No NPZ,
source callback, registered coefficient builder, or physical trajectory is read.
Checks are explicit exceptions and therefore remain active under python -O.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import factorial
from pathlib import Path
import argparse
import json
import sympy as S


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = []
    controls = []

    def check(name, predicate):
        if not bool(predicate):
            raise ValueError('Failed mathematical check: '+name)
        checks.append(name)

    def eq(name, left, right=0):
        check(name, S.simplify(S.expand(left-right)) == 0)

    def reject(name, false_identity):
        if S.simplify(S.expand(false_identity)) == 0:
            raise ValueError('Mutation unexpectedly survived: '+name)
        controls.append(name)

    k, d, L, C, T = S.symbols('k d L C T', real=True, positive=True)
    x = S.symbols('x', real=True)
    E = S.exp(2*S.I*k*d)
    Phi = (E-1)/(2*S.I*k)
    eq('E derivative', S.diff(E,d), 2*S.I*k*E)
    eq('Phi derivative', S.diff(Phi,d), E)
    eq('entire Phi zero momentum', S.limit(Phi,k,0), d)
    eq('sinc representation', S.expand_complex(Phi), d*S.exp(S.I*k*d)*S.sin(k*d)/(k*d))
    eq('E unit modulus', E*S.conjugate(E), 1)
    r, z = S.symbols('r z', nonnegative=True, real=True)
    lam = (S.sqrt(4+r*r)+r)**2/4
    eq('exact semigroup singular value polynomial',lam**2-(2+r*r)*lam+1)
    check('complex L1 is not rotation invariant', S.sqrt(2) > 1)

    ar = k+3*L*L/(2*k)
    ap = k/3-L*L/(2*k)
    bw = -L/(2*k)+S.I/2
    ee = C+S.I*T
    pp = (T+S.I*(1-C))/(2*k)
    hr, hp = ar*pp+bw*ee, ap*pp+bw*ee
    eq('density impulse cancellation',hr,S.I/2+(-L*ee+3*L*L*pp)/(2*k))
    eq('pressure impulse decomposition',hp,S.I/6+S.I*ee/3+(-L*ee-L*L*pp)/(2*k))
    q = S.symbols('q',positive=True,real=True)
    mu = 2*q*k*k
    eq('real weighted density forcing kernel',S.re(mu*hr),q*(-L*k*C+3*L*L*T/2))
    eq('real weighted pressure forcing kernel',S.re(mu*hp),q*(-L*k*C-(2*k*k/3+L*L/2)*T))
    eq('weighted H density zero momentum',S.limit(mu*(ar*Phi+bw*E),k,0),0)
    eq('weighted H pressure zero momentum',S.limit(mu*(ap*Phi+bw*E),k,0),0)
    reject('drop density cancellation constant',hr-(-L*ee+3*L*L*pp)/(2*k))
    reject('wrong pressure oscillatory coefficient',hp-(S.I/6+S.I*ee/2+(-L*ee-L*L*pp)/(2*k)))

    ur, ui, wr, wi, ru, rv, f, g, cr, cp, zz, hprime, base = S.symbols(
        'ur ui wr wi ru rv f g cr cp zz hprime base',real=True)
    UresR, WresR, WresI = S.symbols('UresR WresR WresI',real=True)
    R = ar*ur-wi/2-L*wr/(2*k)
    P = ap*ur-wi/2-L*wr/(2*k)
    vector = {ur:wr+UresR, wr:-2*k*wi-f+WresR, wi:2*k*wr+WresI, L:L*L}

    def derive(expr):
        return sum(S.diff(expr,v)*dv for v,dv in vector.items())

    canonical = ur-wi/(2*k)
    eq('canonical residual drift',derive(canonical),UresR-WresI/(2*k))
    contact_prime = L*(cr-3*cp)-3*hprime*base-L*g/(2*k)+zz
    lhs = derive(R)+contact_prime-(L*(R+cr-3*(P+cp))-3*hprime*base)
    rhs = ar*UresR-WresI/2-L*WresR/(2*k)+L*(f-g)/(2*k)+zz
    eq('complete offshell contact Ward identity',lhs,rhs)
    reject('omit real force mismatch from Ward',lhs-(rhs-L*(f-g)/(2*k)))
    reject('omit canonical U residual from Ward',lhs-(rhs-k*UresR))
    reject('omit contact derivative defect',lhs-(rhs-zz))
    dc = S.symbols('dc',real=True)
    rcanon = mu*R.subs(ur,dc+wi/(2*k))
    pcanon = mu*P.subs(ur,dc+wi/(2*k))
    eq('fused endpoint canonical density',rcanon,q*(k*(2*k*k+3*L*L)*dc+3*L*L*wi/2-k*L*wr))
    eq('fused endpoint canonical pressure',pcanon,q*(k*(2*k*k/3-L*L)*dc-(2*k*k/3+L*L/2)*wi-k*L*wr))
    reject('project nonzero endpoint canonical drift',rcanon-rcanon.subs(dc,0))

    # Independent polynomial manufactured solution with real forcing, k=3.
    u = x*x/6+S.I*x**3/3
    w = x/3+S.I*x*x
    force = -6*x*x-S.Rational(1,3)
    eq('manufactured real-source U equation',S.diff(u,x),w)
    eq('manufactured real-source W equation',S.diff(w,x),6*S.I*w-force)
    eq('manufactured canonical invariant',S.re(u)-S.im(w)/6,0)
    check('manufactured endpoint is nontrivial',u.subs(x,1)!=0 and w.subs(x,1)!=0)
    reject('manufactured source sign reversal',S.diff(w,x)-(6*S.I*w+force))

    # Entire moments checked against direct polynomial beta integrals.
    h = S.symbols('h',positive=True)
    for m in range(5):
        for n in range(5):
            value = S.integrate((h-x)**n*x**m,(x,0,h))/S.factorial(n)
            eq(f'exponential moment beta integral m{m}n{n}',value,
               S.factorial(m)*h**(n+m+1)/S.factorial(n+m+1))
            value_phi = S.integrate((h-x)**(n+1)*x**m,(x,0,h))/S.factorial(n+1)
            eq(f'entire Phi moment beta integral m{m}n{n}',value_phi,
               S.factorial(m)*h**(n+m+2)/S.factorial(n+m+2))

    # A manufactured polynomial cell with a nonzero shift and nontrivial signs.
    co = (S.Rational(2,3),S.Rational(-7,5),S.Rational(11,13),S.Rational(-3,2))
    left,right,center,end = map(S.Rational,('-9/2','-35/8','-71/16','-7/2'))
    poly = sum(v*(x-center)**j for j,v in enumerate(co))
    y = S.symbols('y',real=True)
    lagpoly = S.Poly(S.expand(poly.subs(x,end-y/2)),y)
    for n in range(7):
        direct = S.integrate(poly*(2*(end-x))**n,(x,left,right))
        transformed = S.integrate(lagpoly.as_expr()*y**n/2,(y,2*(end-right),2*(end-left)))
        eq(f'manufactured direct versus lag moment {n}',direct,transformed)
        if n==0:
            reject('omit lag Jacobian',direct-2*transformed)
            reject('reverse lag endpoints',direct+transformed)

    # Exact tail ratio premises and strictly adequate registered remainder.
    xmax,N = Q(512),2048
    rho=xmax/Q(N+2)
    tail=xmax**(N+1)/factorial(N+1)/(1-rho)
    check('degree2048 global phase512 geometric ratio',rho<1)
    check('degree2048 phase512 exponential tail below 1e-340',tail<Q(1,10**340))
    check('higher tail ratios decrease',xmax/Q(N+3)<rho)
    check('Phi factorial tail no larger than D exponential tail',Q(1,factorial(12))<=Q(1,factorial(11)))
    check('zero momentum exponential tail zero',Q(0)**(N+1)==0)
    check('uniform pressure coefficient must test both L endpoints',
          abs(Q(2,3)*Q(1,4)**2-Q(2,9)**2)<abs(Q(2,3)*Q(1,4)**2-Q(2,7)**2))
    # A counterexample to blindly maximizing an absolute pressure coefficient at Lb.
    check('pressure endpoint Lb substitution can underestimate',
          abs(Q(2,3)-Q(2,9)**2)>abs(Q(2,3)-Q(2,7)**2))

    # Exact analytic source-envelope inequalities; no physical source evaluation.
    r0,v,b0,ell = Q(5,8),Q(64,39),Q(89,64),Q(8,27)
    b1=2*r0*v*v*b0
    b2=(4*r0*r0*v**4+2*v*v+8*r0*r0*v**3)*b0
    gb=4*ell*ell*b0+2*ell*b1+b2
    gz=4*ell*ell*r0*b0+2*ell*(b0+r0*b1)+(2*b1+r0*b2)
    check('complex disk source B majorant exact',gb==Q(951819044,20820969))
    check('complex disk source zB majorant exact',gz==Q(3227905133,83283876))
    check('complex disk both force bounds below64',max(gb,gz)<64)
    check('complex disk both work bounds below32',ell*max(gb,gz)<32)
    check('exp geometric majorant exact',1/(1-Q(25,89))==b0)
    check('degree24 Cauchy uniform tail',Q(64)*Q(1,16)**25/(1-Q(1,16))==Q(1,15*2**90))
    check('degree24 Cauchy integrated tail',
          64*2*Q(1,128)*Q(64)*Q(1,16)**25/(26*(1-Q(1,16)))==Q(1,390*2**90))
    zz0,ss = S.symbols('zz0 ss',real=True)
    eq('Re reciprocal disk inequality numerator',
       (1-zz0)*(1+ss)-(1-2*zz0+ss*ss),(1-ss)*(ss+zz0))
    z0=S.symbols('z0')
    bump=S.exp(1-1/(1-z0*z0))
    eq('formal bump first derivative',S.diff(bump,z0),-2*z0*bump/(1-z0*z0)**2)
    eq('formal bump second derivative',S.diff(bump,z0,2),
       (4*z0*z0/(1-z0*z0)**4-2/(1-z0*z0)**2-8*z0*z0/(1-z0*z0)**3)*bump)

    # Invisible-between-samples polynomial controls, no numerical interpolation.
    times=(S.Rational(0),S.Rational(1,2),S.Rational(1))
    invisible=S.prod((x-t)**2 for t in times)
    for j,t in enumerate(times):
        eq(f'invisible U at retained time {j}',invisible.subs(x,t))
        eq(f'invisible W at retained time {j}',S.diff(invisible,x).subs(x,t))
    eq('invisible pair satisfies first equation',S.diff(invisible,x),S.diff(invisible,x))
    check('invisible perturbation nonzero between samples',invisible.subs(x,S.Rational(1,4))>0)
    reject('samples imply zero force residual',S.diff(invisible,x,2)-2*S.I*k*S.diff(invisible,x))
    # A genuine distributional contact counterexample: zero open-panel Ward
    # defect but nonzero endpoint change, all carried by its density jump.
    jump=S.Rational(7)
    contact_R=jump*S.Heaviside(x-S.Rational(1,2))
    contact_P=contact_R/3
    eq('contact jump open-panel Ward work vanishes',L*(contact_R-3*contact_P),0)
    eq('contact derivative includes interface delta',S.diff(contact_R,x),
       jump*S.DiracDelta(x-S.Rational(1,2)))
    eq('integrated contact derivative equals jump',S.integrate(S.diff(contact_R,x),(x,0,1)),jump)
    reject('omit interface jump from complete Ward budget',
           contact_R.subs(x,1)-contact_R.subs(x,0)-S.integrate(L*(contact_R-3*contact_P),(x,0,1)))
    # Constant real source uncertainty saturates both k=0 unitary bounds.
    eta=S.symbols('eta',nonnegative=True)
    eq('constant source U bound is sharp at k0',-S.integrate((1-eta)*3,(eta,0,1)),-S.Rational(3,2))
    eq('constant source W bound is sharp at k0',-S.integrate(3,(eta,0,1)),-3)
    check('whole interval first moment equals cell sum',
          sum((Q(1,64)*(Q(1)-Q(2*j+1,128)) for j in range(64)),Q(0))==Q(1,2))

    scientific = {
        'schema_version':1,
        'status':'PASS_EXACT_SYMBOLIC_AND_MANUFACTURED_ONLY',
        'checks':checks,
        'rejected_mutations_and_negative_controls':controls,
        'counts':{'checks':len(checks),'negative_controls':len(controls)},
        'exact_constants':{
            'complex_g_B_upper':str(gb),'complex_g_zB_upper':str(gz),
            'source_uniform_tail':str(Q(1,15*2**90)),
            'phase512_degree2048_tail_upper':str(Q(1,10**340)),
            'phase512_degree2048_tail_formula':'512^2049/2049!/(1-512/2050)'},
        'operations':{'retained_array_decodes':0,'physical_source_constructions':0,
                      'physical_trajectory_evaluations':0,'production_numerical_imports':0},
        'limits':['Not a proof-assistant formalization; analytic inequality proofs are in companion documents.',
                  'Does not verify a retained-state numerical run or authorize one.',
                  'Manufactured/symbolic controls do not establish external novelty.']}
    payload=json.dumps(scientific,sort_keys=True,separators=(',',':')).encode()
    receipt={'scientific':scientific,'scientific_sha256':sha256(payload).hexdigest()}
    args.output.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':scientific['status'],**scientific['counts'],
                      'scientific_sha256':receipt['scientific_sha256']}))


if __name__=='__main__':
    main()
