#!/usr/bin/env python3
"""Exact rational derivative caps; no sampled source or physical evaluation."""
from fractions import Fraction as F
import argparse
import hashlib
import json
from pathlib import Path
import sys
import sympy as s

ROOT_BITS=60
EXP_ORDER=128

def require(c,m):
    if not c: raise RuntimeError(m)

def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):
    q=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]);return min(q),max(q)
def polynomial(cs,x):
    ans=(F(0),F(0))
    for c in cs:ans=add(mul(ans,x),(F(c),F(c)))
    return ans
def root_sqrt(x):
    l,h=F(0),F(1)
    for _ in range(ROOT_BITS):
        m=(l+h)/2
        if m*m<=x:l=m
        else:h=m
    return l,h
def exp_bounds(t):
    require(0<=t<EXP_ORDER+2,'Positive exponential remainder invalid')
    term=total=F(1)
    for n in range(1,EXP_ORDER+1):term*=t/n;total+=term
    rem=term*t/(EXP_ORDER+1)/(1-t/(EXP_ORDER+2))
    return total,total+rem
def ceil(x):return (x.numerator+x.denominator-1)//x.denominator
def ff(x):return F(int(x.p),int(x.q))
def abs_interval(x):return max(abs(x[0]),abs(x[1]))
def serialize(q):return str(q.numerator) if q.denominator==1 else str(q.numerator)+'/'+str(q.denominator)

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    require(not args.output.exists(),'Refusing overwrite')
    u,z=s.symbols('u z');d=1-u*u
    results={};polys={};identities=0;root_count=0
    for source,P0 in [('positive_B',s.Integer(1)),('signed_uB',u)]:
        P=[P0]
        for n in range(8):
            P.append(s.expand(d*d*s.diff(P[-1],u)+(4*n*u*d-2*u)*P[-1]))
        polys[source]={str(n):[int(c) for c in s.Poly(P[n],u).all_coeffs()] for n in range(9)}
        rows=[]
        for n in range(8):
            # Derivative recurrence is independently checked against direct
            # differentiation of P_n exp(1-1/d)/d^(2n), after factoring exp.
            check=s.cancel(s.diff(P[n]/d**(2*n),u)*d**(2*n+2)-2*u*P[n]-P[n+1])
            require(check==0,'Recurrence failed');identities+=1
            powers=[power[0] for power,coef in s.Poly(P[n+1],u).terms()]
            lowest=min(powers);require(lowest in (0,1),'Unexpected zero multiplicity')
            parity=(n+1)%2 if source=='positive_B' else (n+2)%2
            require(all(power%2==parity for power in powers),'Parity failed');identities+=1
            pp=s.Poly(P[n+1]/u**lowest,u)
            critical=s.Poly(sum(coef*z**(power[0]//2) for power,coef in pp.terms()),z)
            require(s.degree(s.gcd(critical,critical.diff()))==0,'Repeated critical root')
            intervals=s.polys.polytools.intervals(critical,eps=s.Rational(1,2**ROOT_BITS),inf=0,sup=1)
            require(len(intervals)==critical.count_roots(0,1),'Exhaustive root count failed')
            vals=[];cert_roots=[]
            if lowest==1:vals.append((F(P[n].subs(u,0)),F(P[n].subs(u,0))))
            for (lohi,multiplicity) in intervals:
                require(multiplicity==1,'Multiple isolated root')
                zl,zr=map(ff,lohi)
                require(0<zl<=zr<1,'Critical root outside pulse')
                ul=root_sqrt(zl)[0];ur=root_sqrt(zr)[1]
                np=polynomial(polys[source][str(n)],(ul,ur))
                tl,tr=zl/(1-zl),zr/(1-zr)
                E=(1/exp_bounds(tr)[1],1/exp_bounds(tl)[0])
                den=((1-zl)**(-2*n),(1-zr)**(-2*n))
                vals.append(mul(mul(np,E),den))
                cert_roots.append([serialize(zl),serialize(zr)])
                root_count+=1
            sup=max([F(0)]+[abs_interval(q) for q in vals])
            # A full total variation bound is useful when bounding L1 norms
            # of h derivatives by Leibniz, though the simple support-length
            # caps also remain available.
            signed=[(-hi,-lo) if ((n%2==1) if source=='positive_B' else (n%2==0)) else (lo,hi) for lo,hi in reversed(vals[-len(cert_roots):]) ] if cert_roots else []
            if lowest==1:signed.append(vals[0])
            signed.extend(vals[-len(cert_roots):] if cert_roots else [])
            variation=F(0);prev=(F(0),F(0))
            for q in signed+[(F(0),F(0))]:
                variation+=abs_interval((q[0]-prev[1],q[1]-prev[0]));prev=q
            rows.append({'derivative_order':n,'sup_upper_integer':ceil(sup),'total_variation_upper_integer':ceil(variation),'positive_critical_root_count':len(cert_roots),'zero_is_critical':lowest==1,'critical_z_intervals':cert_roots,'upper_rational_sha256':hashlib.sha256(serialize(sup).encode()).hexdigest()})
        sup=[F(row['sup_upper_integer']) for row in rows]
        tv=[F(row['total_variation_upper_integer']) for row in rows]
        # L<=1/3 on the source support; L^(m)=m! L^(m+1).
        # g=4L² h-2L h'-h''; exact Leibniz coefficients below.
        gsup=[];gl1=[]
        hl1=[2*sup[0]]+tv[:-1]
        for n in range(6):
            coeffs=[]
            for j in range(n+1):
                coeffs.append((j,F(4*s.binomial(n,j)*s.factorial(n-j+1),3**(n-j+2))))
                coeffs.append((j+1,F(2*s.binomial(n,j)*s.factorial(n-j),3**(n-j+1))))
            coeffs.append((n+2,F(1)))
            gsup.append(sum(c*sup[j] for j,c in coeffs))
            gl1.append(sum(c*hl1[j] for j,c in coeffs))
        N=[gsup[j+2]+gl1[j+3] for j in range(3)]
        results[source]={'h_derivative_caps':rows,'g_sup_upper_rational':list(map(serialize,gsup)),'g_L1_upper_rational':list(map(serialize,gl1)),'g_N_upper_rational':list(map(serialize,N)),'simple_support_length_N_upper_rational':[serialize(gsup[j+2]+2*gsup[j+3]) for j in range(3)],'norm_proof':'N_j(g)<=sup|g^(j+2)|+L1|g^(j+3)|. Leibniz with L<=1/3 bounds terms, using L1|h^m|=TV(h^(m-1)) for m>=1.'}
    report={'status':'PASS','scope':'Pure exact algebra, exhaustive rational polynomial root isolation, rational exponential inequalities. No pulse sampling, response quadrature, modes, physical sources or roots.','optimization_level':sys.flags.optimize,'polynomial_recurrence_order':8,'h_sup_caps_through_order':7,'identities':identities,'positive_critical_root_intervals':root_count,'root_bits':ROOT_BITS,'positive_exponential_Taylor_order':EXP_ORDER,'polynomials':polys,'envelopes':results,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'identities':identities,'positive_critical_root_intervals':root_count,'g_N_upper_rational':{src:results[src]['g_N_upper_rational'] for src in results}}))

if __name__=='__main__':main()
