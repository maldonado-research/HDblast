#!/usr/bin/env python3
"""Pure algebra: independently derive full-metric bridge and exact primitives.

No physical source evaluation, momentum quadrature or modes occur here.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sympy as S


def derive():
    w, L, v = S.symbols('w L v', positive=True)
    hs = S.symbols('h0:7', real=True)
    h, h1, h2, h3, h4 = hs[:5]
    D = lambda f: S.diff(f,L)*L**2+S.diff(f,w)*2*L**3/w+sum(S.diff(f,hs[j])*hs[j+1] for j in range(6))
    red = lambda f: S.expand(S.cancel(f))
    w1, w2 = 2*L**3/w, 6*L**4/w-4*L**6/w**3
    dw = 2*L**2*h/w
    dw1, dw2 = red(D(dw)), red(D(D(dw)))
    b = -(w*w+4*L*L)/3
    db = -4*L*L*h
    c, dc = 1-b/w**2, -db/w**2+2*b*dw/w**3
    U = -L**2/w-w2/(4*w**2)+3*w1**2/(8*w**3)
    dU = (-h2-2*L*h1)/(2*w)+L**2*dw/w**2-dw2/(4*w**2)+w2*dw/(2*w**3)+3*w1*dw1/(4*w**3)-9*w1**2*dw/(8*w**4)
    U, dU = red(U), red(dU)
    Up, Upp, dUp, dUpp = map(red,(D(U),D(D(U)),D(dU),D(D(dU))))
    U4 = red(-U*U/(2*w)-Upp/(4*w**2)+w2*U/(4*w**3)+3*w1*Up/(4*w**3)-3*w1*w1*U/(4*w**4))
    dU4 = red(-U*dU/w+U*U*dw/(2*w**2)-dUpp/(4*w**2)+Upp*dw/(2*w**3)+(dw2*U+w2*dU)/(4*w**3)-3*w2*U*dw/(4*w**4)+3*(dw1*Up+w1*dUp)/(4*w**3)-9*w1*Up*dw/(4*w**4)-3*(2*w1*dw1*U+w1*w1*dU)/(4*w**4)+3*w1*w1*U*dw/w**5)
    J2 = red(L*w1/w**2+w1*w1/(4*w**3))
    dJ2 = red((h1*w1+L*dw1)/w**2-2*L*w1*dw/w**3+w1*dw1/(2*w**3)-3*w1*w1*dw/(4*w**4))
    J4 = red(L*Up/w**2-2*L*w1*U/w**3+w1*Up/(2*w**3)-3*w1*w1*U/(4*w**4))
    dJ4 = red((h1*Up+L*dUp)/w**2-2*L*Up*dw/w**3-2*(h1*w1*U+L*dw1*U+L*w1*dU)/w**3+6*L*w1*U*dw/w**4+(dw1*Up+w1*dUp)/(2*w**3)-3*w1*Up*dw/(2*w**4)-3*(2*w1*dw1*U+w1*w1*dU)/(4*w**4)+3*w1*w1*U*dw/w**5)
    Qs = red(1/(2*w)-U/(2*w**2))
    dQs = red(-dw/(2*w**2)-dU/(2*w**2)+U*dw/w**3)
    Rs = red(w/2+(L*L/w+J2)/4+(U*U/w-L*L*U/w**2+J4)/4)
    Ps = red((w*w-2*L*L)/(6*w)+(c*U+L*L/w+J2)/4+(c*U4+b*U*U/w**3-L*L*U/w**2+J4)/4)
    dRs = red(dw/2+(2*L*h1/w-L*L*dw/w**2+dJ2)/4+(2*U*dU/w-U*U*dw/w**2-2*L*h1*U/w**2-L*L*dU/w**2+2*L*L*U*dw/w**3+dJ4)/4)
    dPs = red(-(w*w-2*L*L)*dw/(6*w**2)+(dc*U+c*dU+2*L*h1/w-L*L*dw/w**2+dJ2)/4+(dc*U4+c*dU4+db*U*U/w**3+2*b*U*dU/w**3-3*b*U*U*dw/w**4-2*L*h1*U/w**2-L*L*dU/w**2+2*L*L*U*dw/w**3+dJ4)/4)
    g = 4*L*L*h-2*L*h1-h2
    g1, g2 = red(D(g)), red(D(D(g)))
    mU, mUp = g/(2*w), g1/(2*w)-g*w1/(2*w**2)
    mUpp = g2/(2*w)-g1*w1/w**2-g*w2/(2*w**2)+g*w1*w1/w**3
    mU4 = -U*mU/w-mUpp/(4*w**2)+w2*mU/(4*w**3)+3*w1*mUp/(4*w**3)-3*w1*w1*mU/(4*w**4)
    mJ4 = L*mUp/w**2-2*L*w1*mU/w**3+w1*mUp/(2*w**3)-3*w1*w1*mU/(4*w**4)
    mR = (g/w-L*L*mU/w**2+mJ4)/4
    mP = (c*mU-g/w+c*mU4+2*b*U*mU/w**3+g*U/w**2-L*L*mU/w**2+mJ4)/4
    Ar, Ap = (h2+4*L*h1)/2, -h2/2
    A0=red(Rs-3*Ps-2*L*L*Qs+(D(D(Qs))-2*L*D(Qs)-2*L*L*Qs)/2)
    dA=red(dRs-3*dPs-4*L*L*h*Qs-2*L*L*dQs+(D(D(dQs))-2*h1*D(Qs)-2*L*D(dQs)-2*h2*Qs-2*L*L*dQs)/2-4*h*A0)
    contacts = {'q':red(-g/(4*w**3)-dQs), 'rho':red(mR-dRs+Ar*Qs), 'p':red(mP-dPs+Ap*Qs), 'anomaly':dA}
    def inventory(expr):
        terms = {}
        for term in S.Add.make_args(S.expand(expr)):
            exponent = int(term.as_powers_dict().get(w,0))
            terms[-exponent] = terms.get(-exponent,0)+term/w**exponent
        return {n:S.factor(a) for n,a in sorted(terms.items()) if a != 0}
    inventories = {name:inventory(expr) for name,expr in contacts.items()}
    for items in inventories.values():
        if not all(n>=5 and n%2==1 for n in items):
            raise RuntimeError('Nonconvergent bridge contact remains')
    def I(n):
        jmax=(n-5)//2
        return (2*L*L)**S.Rational(3-n,2)*sum((-1)**j*S.binomial(jmax,j)*v**(2*j+3)/S.Integer(2*j+3) for j in range(jmax+1))
    local = {name:S.factor(sum(a*I(n) for n,a in items.items())/(2*S.pi**2)) for name,items in inventories.items()}
    K=S.sqrt(2)*L*v/S.sqrt(1-v*v)
    Wend=S.sqrt(2)*L/S.sqrt(1-v*v)
    A=S.symbols('A',real=True)
    def primitive(n):
        if n==-1: return (K*Wend*(2*K*K+2*L*L)-4*L**4*A)/8
        if n==1: return (K*Wend-2*L*L*A)/2
        if n==3: return A-v
        return I(n)
    Rbase = S.factor((K**4/8+3*L*L*K*K/8-sum(a*primitive(n) for n,a in inventory(Rs).items()))/(2*S.pi**2))
    Pbase = S.factor((K**4/24-L*L*K*K/8-sum(a*primitive(n) for n,a in inventory(Ps).items()))/(2*S.pi**2))
    Rbase, Pbase = map(S.factor,(S.cancel(Rbase),S.cancel(Pbase)))
    if Rbase.has(A) or Pbase.has(A): raise RuntimeError('Residual divergent baseline logarithm')
    Q0=(v*v/(2*(1+v))-v**3/48-v**5/16)/(2*S.pi**2)
    Dv = lambda f: S.diff(f,L)*L**2-S.diff(f,v)*L*v*(1-v*v)+sum(S.diff(f,hs[j])*hs[j+1] for j in range(6))
    qlocal=S.factor(-2*h*L*L*Q0+local['q'])
    Rlocal=S.factor(-4*h*Rbase+Ar*L*L*Q0+local['rho'])
    Plocal=S.factor(-4*h*Pbase+Ap*L*L*Q0+local['p'])
    expressions={'Q0':Q0,'Q0_prime':S.factor(Dv(Q0)),'Q0_second':S.factor(Dv(Dv(Q0))),'rho0':S.factor(Rbase/L**4),'p0':S.factor(Pbase/L**4),'q_local':qlocal,'q_local_prime':S.factor(Dv(qlocal)),'q_local_second':S.factor(Dv(Dv(qlocal))),'rho_local':Rlocal,'p_local':Plocal,'rho_local_prime':S.factor(Dv(Rlocal)),'anomaly':local['anomaly']}
    targets={'q':(12*L*L*h-2*L*h1-h2)/(48*S.pi**2),'rho':-L*(61*L*L*h1+3*L*h2-3*h3)/(240*S.pi**2),'p':-(7*L**3*h1-44*L*L*h2-3*L*h3+3*h4)/(720*S.pi**2),'anomaly':(-116*L**4*h+94*L**3*h1+47*L*L*h2-3*h4)/(240*S.pi**2)}
    checks=[]
    for name,target in targets.items():
        if S.simplify(local[name].subs(v,1)-target)!=0: raise RuntimeError('Continuum contact mismatch '+name)
        checks.append({'name':'continuum_'+name,'passed':True})
    for n in range(5,16,2):
        kk=S.symbols('K',positive=True)
        iv=I(n).subs(v,kk/S.sqrt(kk*kk+2*L*L))
        if S.simplify(S.diff(iv,kk)-kk*kk/(kk*kk+2*L*L)**S.Rational(n,2))!=0: raise RuntimeError('Primitive mismatch')
        checks.append({'name':'primitive_'+str(n),'passed':True})
    baseline_targets={'Q0':(1/(2*K)-Qs.subs(w,Wend))/L**2,
                      'rho0':((2*K*K+3*L*L)/(4*K)-Rs.subs(w,Wend))/L**4,
                      'p0':((2*K*K/3-L*L)/(4*K)-Ps.subs(w,Wend))/L**4}
    for name,target in baseline_targets.items():
        primitive_error=S.factor(S.diff(expressions[name],v)-K*K*target*S.diff(K,v)/(2*S.pi**2))
        if primitive_error != 0: raise RuntimeError('Full baseline primitive mismatch '+name)
        if S.simplify(expressions[name].subs(v,0)) != 0: raise RuntimeError('Baseline lower endpoint mismatch '+name)
        checks.append({'name':'complete_baseline_primitive_'+name,'passed':True})
    ward_local=S.factor(L*L*Q0*(g1-2*L*g)/2+Dv(Rlocal)-L*Rlocal+3*L*Plocal+3*h1*(Rbase+Pbase))
    if ward_local != 0: raise RuntimeError('Finite-K full metric Ward contact bridge mismatch')
    checks.append({'name':'independently_defined_metric_stresses_Ward_bridge','passed':True})
    T=S.symbols('T',real=True)
    shift={hs[j]:S.factorial(j)*L**(j+1)*T for j in range(7)}
    for name,expr in {'g':g,'q':qlocal.subs(v,1),'rho':Rlocal.subs(v,1),'p':Plocal.subs(v,1)}.items():
        if S.simplify(expr.subs(shift)) != 0: raise RuntimeError('Coordinate translation mismatch '+name)
        checks.append({'name':'pure_coordinate_translation_'+name,'passed':True})
    checks.extend([{'name':'bridge_inventory_converges','passed':True},{'name':'baseline_logs_cancel','passed':True}])
    return inventories,expressions,checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise RuntimeError('Refusing to overwrite algebra evidence')
    inv,expr,checks=derive()
    args.output.mkdir(parents=True)
    document={'passed':True,'scope':'Exact symbols only; no physical sources, modes or quadratures','sympy':S.__version__,'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'coefficient_inventory':{name:{str(n):str(x) for n,x in items.items()} for name,items in inv.items()},'closed_expressions':{name:str(x) for name,x in expr.items()},'checks':checks}
    (args.output/'METRIC_CONTACT_ALGEBRA.json').write_text(json.dumps(document,indent=2,sort_keys=True)+'\n')
    lines=['"""Generated exact rational contacts; importing performs no physical evaluation."""','import math','', 'def local_coefficients(L, v, jet):','    h0, h1, h2, h3, h4, h5 = jet','    return {']
    # The anomaly is complementary exact theory, outside the numerical round.
    for name,x in expr.items():
        if name not in ('anomaly','Q0_prime','Q0_second'):
            lines.append('        '+repr(name)+': '+S.pycode(x)+',')
    lines+=['    }','']
    (args.output/'metric_contact_coefficients.py').write_text('\n'.join(lines))
    print('Primary exact metric contact derivation passed:',len(checks),'checks')


if __name__=='__main__': main()
