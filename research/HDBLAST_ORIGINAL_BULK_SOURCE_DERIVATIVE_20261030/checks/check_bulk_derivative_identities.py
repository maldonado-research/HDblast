#!/usr/bin/env python3
"""Independent exact physical-jet and complete source chain-rule audit.

No scientific producer imports. Polynomial identities support the analytic
derivation; arbitrary actual-field/gauge premises remain explicit.
"""
import argparse, hashlib, json, math
from fractions import Fraction as Q
from pathlib import Path

PIN='34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
CHECKS=[]
CONTROLS=[]
def need(ok,why):
    if not ok:raise RuntimeError(why)
    CHECKS.append(why)
def reject(name,wrong,correct):
    need(wrong!=correct,'reject '+name);CONTROLS.append(name)
def add(*parts):
    z={}
    for p in parts:
        for m,c in p.items():z[m]=z.get(m,Q(0))+c
    return {m:c for m,c in z.items() if c}
def scale(p,c):return {m:v*c for m,v in p.items() if v*c}
def mul(a,b):
    z={}
    for ma,ca in a.items():
        for mb,cb in b.items():
            d=dict(ma)
            for n,p in mb:d[n]=d.get(n,0)+p
            key=tuple(sorted((n,p) for n,p in d.items() if p))
            z[key]=z.get(key,Q(0))+ca*cb
    return {m:c for m,c in z.items() if c}
def mono(coef=1,**p):return {tuple(sorted((n,k) for n,k in p.items() if k)):Q(coef)} if coef else {}
def diff(p,name):
    z={}
    for m,c in p.items():
        d=dict(m);n=d.get(name,0)
        if n:
            d[name]=n-1;k=tuple(sorted((v,j) for v,j in d.items() if j))
            z[k]=z.get(k,Q(0))+n*c
    return {m:c for m,c in z.items() if c}
def power(p,n):
    need(n>=0,'polynomial nonnegative exponent')
    z=mono()
    for _ in range(n):z=mul(z,p)
    return z
def decode(terms):
    z={}
    for t in terms:
        need(set(t)=={'coefficient','powers'},'encoded Laurent term shape')
        k=tuple(sorted((s,int(n)) for s,n in t['powers'].items() if n))
        z[k]=z.get(k,Q(0))+Q(t['coefficient'])
    return {m:c for m,c in z.items() if c}
def canonical(p):
    return [[str(c),[[n,e] for n,e in m]] for m,c in sorted(p.items())]
def digest(p):return hashlib.sha256(json.dumps(canonical(p),sort_keys=True,separators=(',',':')).encode()).hexdigest()
def dt_symbol(name):
    if name[0] in 'ABCDELMs' and name[1:] in ('','y','t','yy','yt','tt'):
        return mono(**{name[0]+{'':'t','y':'yt','t':'tt','yy':'yyt','yt':'ytt','tt':'ttt'}[name[1:]]:1})
    if name=='c':return mono(**{'c':1,'theta':1})
    if name=='theta':return add(mono(),mono(-1,theta=2))
    need(name in ('Vp','Vpp','Vppp','a','ay','ell','p','py','rho'),'static background symbol '+name)
    return {}
def dt(p,omit=None):
    out={}
    for m,c in p.items():
        for name,n in m:
            if name==omit:continue
            d=dict(m);d[name]=n-1
            k=tuple(sorted((v,j) for v,j in d.items() if j))
            out=add(out,mul({k:c*n},dt_symbol(name)))
    return out
def evaluate(p,values):
    z=Q(0)
    for m,c in p.items():
        q=c
        for name,n in m:q*=values[name]**n
        z+=q
    return z
class Dual:
    def __init__(self,v,d=0):self.v,self.d=Q(v),Q(d)
    def __mul__(self,b):return Dual(self.v*b.v,self.d*b.v+self.v*b.d)
    def __pow__(self,n):return Dual(self.v**n, n*self.v**(n-1)*self.d if n else 0)
def evaluate_dual(p,values):
    out=Q(0)
    for m,c in p.items():
        z=Dual(c)
        for name,n in m:z=z*values[name]**n
        out+=z.d
    return out

def physical_conversion_audit():
    ell=Q(3,7)
    L=add(mono(2),mono(3,r=1),mono(5,r=2));a=diff(L,'r')
    Dy=lambda p:scale(mul(L,diff(p,'r')),-1)
    Dt=lambda p:scale(diff(p,'x'),1/ell)
    count=0
    for i,j in ((1,2),(2,1),(3,0),(0,3),(4,0),(3,1),(2,2),(1,3),(2,3),(3,2)):
        amp=add(mono(Q(2,3),x=i,r=j),mono(Q(-5,7),x=i+1,r=j+1))
        xrr=diff(diff(diff(amp,'x'),'r'),'r')
        xr=diff(diff(amp,'x'),'r')
        xxr=diff(diff(diff(amp,'x'),'x'),'r')
        xxx=diff(diff(diff(amp,'x'),'x'),'x')
        yyt=scale(add(mul(power(L,2),xrr),mul(mul(a,L),xr)),1/ell)
        ytt=scale(mul(L,xxr),-1/ell**2)
        ttt=scale(xxx,1/ell**3)
        need(yyt==Dt(Dy(Dy(amp))),'ordinary yyt conversion fixture '+str((i,j)))
        need(ytt==Dt(Dt(Dy(amp))),'ordinary ytt conversion fixture '+str((i,j)))
        need(ttt==Dt(Dt(Dt(amp))),'ordinary ttt conversion fixture '+str((i,j)))
        if (i,j)==(1,2):
            reject('omit radial rho product in yyt',scale(mul(power(L,2),xrr),1/ell),yyt)
        if (i,j)==(2,1):
            reject('wrong ytt sign',scale(mul(L,xxr),1/ell**2),ytt)
            reject('missing ell in ytt',scale(mul(L,xxr),-1/ell),ytt)
        count+=3
    for i in range(5):
        for j in range(5-i):
            normalized=Q(2*i+3*j+1,7)
            p=mono(normalized,x=i,r=j)
            for _ in range(i):p=diff(p,'x')
            for _ in range(j):p=diff(p,'r')
            need(p==mono(math.factorial(i)*math.factorial(j)*normalized),'normalized Taylor factorial '+str((i,j)))
    reject('raw xrr factorial omitted',Q(1),Q(math.factorial(1)*math.factorial(2)))
    return count

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-map',type=Path,required=True)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args();raw=args.source_map.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==PIN,'complete map authenticated')
    maps=json.loads(raw);need(maps['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','complete map status')
    need(maps['actual_mode_domain']=='QEred,Qphi require C1=0 as an identity, including its derivative; general alternatives retained','actual harmonic neighbourhood premise')
    physical_count=physical_conversion_audit()
    names=[];polys=[]
    for field,rows in maps['maps']['metric_forcing'].items():
        for i,row in enumerate(rows):names.append('metric_forcing/'+field+'/'+str(i));polys.append(decode(row))
    for i,row in enumerate(maps['maps']['scalar_forcing']):names.append('scalar_forcing/'+str(i));polys.append(decode(row))
    need(len(polys)==23 and len(set(names))==23,'all23 complete forcing rows')
    ordinary=[dt(p) for p in polys]
    derivs=[mul(mono(ell=1),p) for p in ordinary]
    normalizations=[]
    for name,p,dp in zip(names,polys,ordinary):
        parts=name.split('/')
        if parts[0]=='metric_forcing':
            natural=decode(maps['maps']['QEred'][parts[1]][int(parts[2])]);factor=mono(-2,ell=2,rho=2)
        else:
            natural=decode(maps['maps']['Qphi'][int(parts[1])]);factor=mono(1,ell=2,rho=2)
        need(p==mul(factor,natural),'complete forcing normalization '+name)
        need(dp==mul(factor,dt(natural)),'complete differentiated forcing normalization '+name)
        normalizations.append(name)
    allnew=set(n for p in derivs for m in p for n,e in m)-set(maps['input_symbols'])
    expected={ch+suf for ch in 'ABCDELMs' for suf in ('yyt','ytt','ttt')}
    need(allnew<=expected,'all derivative-map new inputs are ordinary third time-containing jets')
    need(allnew==expected,'all24 new metric/scalar jets occur')
    checks=0
    for fixture in range(4):
        symbols=sorted(set(maps['input_symbols'])|expected)
        values={n:Q(3*(i+1)+fixture+1,11+fixture) for i,n in enumerate(symbols)}
        # Keep explicit scales positive and nonsingular; fixture choices need
        # not represent a physical solution to check Laurent chain identities.
        values.update({'ell':Q(fixture+2,17),'rho':Q(fixture+4,3),'c':Q(fixture+3,2),'theta':Q(fixture+1,9)})
        duals={}
        for name in maps['input_symbols']:
            if name[0] in 'ABCDELMs' and name[1:] in ('','y','t','yy','yt','tt'):
                dvalue=values[name+'t']
            elif name=='c':dvalue=values['c']*values['theta']
            elif name=='theta':dvalue=1-values['theta']**2
            else:dvalue=Q(0)
            duals[name]=Dual(values[name],dvalue)
        for name,p,d in zip(names,polys,derivs):
            need(evaluate(d,values)==values['ell']*evaluate_dual(p,duals),'independent dual vs Laurent derivative '+str(fixture)+' '+name)
            checks+=1
        if fixture==0:
            cwrong=[mul(mono(ell=1),dt(p,omit='c')) for p in polys]
            theta_wrong=[mul(mono(ell=1),dt(p,omit='theta')) for p in polys]
            need(any(evaluate(d,values)!=evaluate(w,values) for d,w in zip(derivs,cwrong)),'nonzero generic cosh derivative witness')
            reject('omit cosh derivative map',[evaluate(w,values) for w in cwrong],[evaluate(d,values) for d in derivs])
            reject('omit tanh derivative map',[evaluate(w,values) for w in theta_wrong],[evaluate(d,values) for d in derivs])
            reject('omit scaled-time outer ell',[evaluate(dt(p),values) for p in polys],[evaluate(d,values) for d in derivs])
    report={'status':'PASS_INDEPENDENT_COMPLETE_BULK_TIME_DERIVATIVE_IDENTITIES',
      'source_map_sha256':PIN,'physical_conversion_polynomial_identities':physical_count,
      'complete_rows':names,'source_term_count':sum(len(p) for p in polys),
      'differentiated_term_count':sum(len(p) for p in derivs),'new_third_jet_symbols':sorted(allnew),
      'ordinary_physical_time_derivative_map_sha256':{n:digest(p) for n,p in zip(names,ordinary)},
      'complete_scaled_derivative_map_sha256':{n:digest(p) for n,p in zip(names,derivs)},
      'complete_forcing_normalization_identities':normalizations,
      'dual_chain_rule_fixture_rows':checks,
      'conversion_formulas':{'yyt':'(L^2 j_xrr+a L j_xr)/ell','ytt':'-L j_xxr/ell^2','ttt':'j_xxx/ell^3','L':'1/(ell rho), L_r=a'},
      'source_derivative':'F_x=ell D_t complete stored metric/scalar forcing; no extra forcing normalization',
      'actual_premises':'same actual first field; C1=0 in a neighbourhood; quantitative gauge/master C4 still supplied separately',
      'wrong_formula_controls_rejected':CONTROLS,'distinct_check_count':len(set(CHECKS)),
      'scope':'symbolic full-map and physical-coordinate identity audit, not a standalone numeric F_x bound, target C2 error, nonlinear realization or cosmological evidence'}
    out=json.dumps(report,indent=2,sort_keys=True)+'\n'
    if args.output:
        need(not args.output.exists(),'fresh output');args.output.write_text(out)
    else:print(out,end='')
if __name__=='__main__':main()
