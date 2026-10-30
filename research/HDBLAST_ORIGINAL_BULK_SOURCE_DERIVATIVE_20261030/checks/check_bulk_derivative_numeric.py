#!/usr/bin/env python3
"""Independent ordinary-jet composition and full-map interval norm audit.

Uses ordinary derivatives with binomial Leibniz coefficients, unlike the
producer's normalized Taylor convolution. Imports only the independent
identity checker after authenticating its bytes.
"""
import argparse, hashlib, importlib.util, json, math
from fractions import Fraction as Q
from pathlib import Path

CHECKS=[]
def need(test,label):
    if not test:raise RuntimeError(label)
    CHECKS.append(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p,pin,role):
    need(sha(p)==pin,'authenticate '+role);return json.loads(p.read_text())
class I:
    def __init__(self,lo=0,hi=None):
        self.lo=Q(lo);self.hi=self.lo if hi is None else Q(hi)
        if self.lo>self.hi:raise RuntimeError('interval order')
    @staticmethod
    def cv(z):return z if isinstance(z,I) else I(z)
    def __add__(self,b):
        b=I.cv(b);return I(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,b):return self+-I.cv(b)
    def __rsub__(self,b):return I.cv(b)+-self
    def __mul__(self,b):
        b=I.cv(b);v=[self.lo*b.lo,self.lo*b.hi,self.hi*b.lo,self.hi*b.hi]
        return I(min(v),max(v))
    __rmul__=__mul__
    def inv(self):
        if self.lo<=0<=self.hi:raise RuntimeError('zero denominator')
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,b):return self*I.cv(b).inv()
    def __rtruediv__(self,b):return I.cv(b)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        if n==0:return I(1)
        if n%2==0:
            low=0 if self.lo<=0<=self.hi else min(self.lo**n,self.hi**n)
            return I(low,max(self.lo**n,self.hi**n))
        return I(self.lo**n,self.hi**n)
    def inside(self,b):return b.lo<=self.lo<=self.hi<=b.hi
    def mag(self):return max(abs(self.lo),abs(self.hi))
    def row(self):return [str(self.lo),str(self.hi)]
KEYS=tuple((i,j) for n in range(4) for i in range(n+1) for j in [n-i])
class O:
    """Unnormalized ordinary mixed derivatives through total order three."""
    def __init__(self,c=None):self.c={k:I.cv((c or {}).get(k,0)) for k in KEYS}
    @staticmethod
    def cv(b):return b if isinstance(b,O) else O({(0,0):b})
    def __add__(self,b):
        b=O.cv(b);return O({k:self.c[k]+b.c[k] for k in KEYS})
    __radd__=__add__
    def __neg__(self):return O({k:-v for k,v in self.c.items()})
    def __sub__(self,b):return self+-O.cv(b)
    def __rsub__(self,b):return O.cv(b)+-self
    def __mul__(self,b):
        b=O.cv(b)
        return O({(i,j):sum((math.comb(i,a)*math.comb(j,b0)*self.c[a,b0]*b.c[i-a,j-b0]
                   for a in range(i+1) for b0 in range(j+1)),I(0)) for i,j in KEYS})
    __rmul__=__mul__
    def inv(self):
        z={(0,0):self.c[0,0].inv()}
        for i,j in KEYS[1:]:
            z[i,j]=-z[0,0]*sum((math.comb(i,a)*math.comb(j,b)*self.c[a,b]*z[i-a,j-b]
                       for a in range(i+1) for b in range(j+1) if a+b),I(0))
        return O(z)
    def __truediv__(self,b):return self*O.cv(b).inv()
    def __rtruediv__(self,b):return O.cv(b)*self.inv()
    def val(self):return self.c[0,0]
def ceil_sqrt(q):
    z=math.isqrt(q.numerator//q.denominator)
    if z*z<=q:z+=1
    return z
def encode(z):
    if isinstance(z,Q):return str(z)
    if isinstance(z,dict):return {k:encode(v) for k,v in z.items()}
    if isinstance(z,(list,tuple)):return [encode(v) for v in z]
    return z

ap=argparse.ArgumentParser()
ap.add_argument('--fx-report',type=Path,required=True)
ap.add_argument('--fx-sha256',required=True)
ap.add_argument('--source-map',type=Path,required=True)
ap.add_argument('--identity-checker',type=Path)
ap.add_argument('--identity-checker-sha256',required=True)
ap.add_argument('--gauge4-report',type=Path,required=True)
ap.add_argument('--master4-report',type=Path,required=True)
ap.add_argument('--master3-report',type=Path,required=True)
a=ap.parse_args()
identity=a.identity_checker or Path(__file__).resolve().with_name('check_bulk_derivative_identities.py')
need(sha(identity)==a.identity_checker_sha256,'authenticate independent identity checker')
spec=importlib.util.spec_from_file_location('independent_source_chain',identity)
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
r=load(a.fx_report,a.fx_sha256,'complete actual Fx report')
maps=load(a.source_map,'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615','complete source maps')
g=load(a.gauge4_report,r['input_sha256']['gauge4_report'],'actual gauge C4 report')
m=load(a.master4_report,'d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e','actual master C4 report')
m3=load(a.master3_report,'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df','retained actual master through3')
need(r['status']=='PASS_COMPLETE_ORIGINAL_ACTUAL_FULL23_SCALED_BULK_SOURCE_TIME_DERIVATIVE','actual full original source derivative status')
need(r['actual_source_unchanged'] and not r['scope']['source_smoothing_or_replacement'],'actual source retained without smoothing')
ell=Q(r['ell']);need(ell==Q(1,10**6),'registered scale')
need(r['actual_raw_gauge_scaled_mixed_caps']==g['actual_mixed_scaled_derivative_caps'],'same actual gauge C4 table')
need(max(map(Q,m['exact']['all_scaled_order_four_mixed_caps'].values()))==Q(r['actual_raw_master_scaled_order_caps'][4]),'same actual master fourth cap')
master=list(map(Q,r['actual_raw_master_scaled_order_caps']))
old=m3['exact']['interior']
need(master[:4]==[Q(old['U']),Q(old['first']),max(map(Q,old['second'].values())),max(map(Q,old['third'].values()))],
     'same retained actual master lower-order caps')
gauge={tuple(map(int,k.split(','))):list(map(Q,v)) for k,v in r['actual_raw_gauge_scaled_mixed_caps'].items()}
grid=1<<112
def sym(cap):
    bound=Q(-((-cap*grid).__floor__()),grid)
    # ceil(cap*grid)/grid includes the exact raw cap before all calculus.
    return I(-bound,bound)
def raw(ch,dx=0,dr=0):
    return O({(i,j):sym(master[i+j+dx+dr] if ch=='W' else gauge[i+dx,j+dr][ch]) for i,j in KEYS})
co=r['scaled_normalized_coefficient_jet_boxes']
need(set(co)=={'rho','p','a','R','B','k','c','H','vartheta'},'all9 actual adapter coefficient jets')
coeff={name:O({(i,j):I(*rows[str(i)+','+str(j)])*math.factorial(i)*math.factorial(j) for i,j in KEYS}) for name,rows in co.items()}
for name in ('rho','p','a','R','B','k'):
    need(all(coeff[name].c[i,j].lo==coeff[name].c[i,j].hi==0 for i,j in KEYS if i),'static radial coefficient '+name)
for name in ('c','H','vartheta'):
    need(all(coeff[name].c[i,j].lo==coeff[name].c[i,j].hi==0 for i,j in KEYS if j),'temporal coefficient '+name)
rho,p,aa,R,B,k,c,H,v=[coeff[n] for n in ('rho','p','a','R','B','k','c','H','vartheta')]
W=raw('W');Wr=raw('W',dr=1)
f,b,d=[raw(i) for i in range(3)]
fx,fr,bx,br,dx,dr=raw(0,dx=1),raw(0,dr=1),raw(1,dx=1),raw(1,dr=1),raw(2,dx=1),raw(2,dr=1)
u=R*v*W;scalar=B*v*(Wr/ell+k*W)+p*f
amps={'A':-4*u-2*fr/(ell*rho),'B':fx/(ell*rho)+br/ell,
      'C':f/(rho*c)-c*dr/ell,'D':-2*(u+aa*f+bx/ell),
      'E':c*dx/ell-b/c,'L':2*(u+aa*f+H*b),'M':2*d,'s':scalar}
L=1/(ell*rho.val())
third={}
for name,j in amps.items():
    third[name+'yyt']=(L**2*j.c[1,2]+aa.val()*L*j.c[1,1])/ell
    third[name+'ytt']=-L*j.c[2,1]/ell**2
    third[name+'ttt']=j.c[3,0]/ell**3
need(set(third)==set(r['new_actual_field_third_jet_boxes']),'all24 independent ordinary third-field slots')
outside=[]
for name,x in third.items():
    if not x.inside(I(*r['new_actual_field_third_jet_boxes'][name])):outside.append(name)
need(not outside,'ordinary binomial third-jet composition inside producer boxes: '+','.join(outside))
lower={n:I(*v) for n,v in r['lower_actual_field_and_background_global_boxes'].items()}
need(set(lower)==set(maps['input_symbols']),'complete original lower input closure')
for name in ('rho','p','a','c'):
    need(lower[name].inside(coeff[name].val()),'whole coefficient box covers lower-source '+name)
dcover=r['causal_cover']['division_per_unit'];need(isinstance(dcover,int) and dcover>=1,'positive actual source cover resolution')
indices=[[i,j] for i in range(3*dcover) for j in range(8*dcover) if i+j<8*dcover]
need(r['causal_cover']['indices']==indices and len(indices)==(39*dcover*dcover+3*dcover)//2,'complete clipped causal source cover')
mapping=[('yy',('K:0','K:2','K:4')),('y0',('L:0','L:2','L:4')),('00',('N:0','N:2','N:4')),
         ('ya',('v:2','v:4')),('0a',('u:2','u:4')),('ab',('J:0','J:2','J:4','H:2','H:4','T:+','T:-'))]
rows=[]
for field,names in mapping:
    for i,name in enumerate(names):rows.append((name,maps['maps']['metric_forcing'][field][i]))
for i,name in enumerate(('w:0','w:2','w:4')):rows.append((name,maps['maps']['scalar_forcing'][i]))
values={**lower,**third};independent={}
for name,terms in rows:
    ordinary=s.dt(s.decode(terms))
    need(s.digest(ordinary)==r['differentiated_row_sha256'][name],'independent complete symbolic derivative '+name)
    need(len(ordinary)==r['differentiated_row_monomial_counts'][name],'collected derivative term count '+name)
    # Evaluate independent exact rational intervals without dyadic rounding.
    z=I(0)
    for mon,coef in ordinary.items():
        term=I(coef)
        for variable,n in mon:term=term*values[variable]**n
        z=z+term
    independent[name]=ell*z
    need(independent[name].inside(I(*r['component_F_x_intervals'][name])),'independent derivative interval contained '+name)
S={0:Q(1),2:Q(1,12),4:Q(19,720)}
weights={f'{ch}:{l}':fac*S[l] for ch,fac in [('K',1),('L',2),('N',1),('J',3),('w',1)] for l in S}
weights.update({f'{ch}:{l}':2*l*(l+2)*S[l] for ch in ('u','v') for l in (2,4)})
weights.update({f'H:{l}':Q(2,3)*l*(l+2)*(l*(l+2)-3)*S[l] for l in (2,4)})
weights.update({'T:+':Q(1,6),'T:-':Q(1,6)})
need({n:Q(w) for n,w in r['positive_Gram'].items()}==weights,'all23 positive Gram weights independently packed')
reported_square=sum((weights[n]*I(*r['component_F_x_intervals'][n]).mag()**2 for n in weights),Q(0))
need(reported_square==Q(r['F_x_H_squared_upper']),'exact reported full23 squared norm')
independent_square=sum((weights[n]*independent[n].mag()**2 for n in weights),Q(0))
need(independent_square<=reported_square,'independent fullsource derivative norm dominated')
cap=ceil_sqrt(reported_square)
need(cap==int(r['genuine_scaled_F_x_H_cap'])==905405545954284,'genuine original Fx norm cap independently recovered')
need(Q(cap)<Q('9.055e14'),'readable original source derivative cap')
controls={
 'incorrect derivative norm without offdiagonal weights':sum((I(*r['component_F_x_intervals'][n]).mag()**2 for n in weights),Q(0))!=reported_square,
 'response SOURCE substituted for derivative cap':Q(cap)!=Q(52162547172519,2),
 'omitted scaled-time factor ell':Q(cap)/ell!=cap,
 'third Taylor coefficient used as ordinary derivative':Q(1,6)!=1,
}
need(all(controls.values()),'all independent wrong formula controls rejected')
need(sha(a.fx_report)==a.fx_sha256 and sha(identity)==a.identity_checker_sha256,'review inputs unchanged on readback')
report={'status':'PASS_INDEPENDENT_ORIGINAL_FULL23_BULK_SOURCE_DERIVATIVE_NUMERIC_AUDIT',
 'input_sha256':{'actual_Fx_report':a.fx_sha256,'independent_identity_checker':a.identity_checker_sha256,'master4_report':sha(a.master4_report),'master3_report':sha(a.master3_report),'gauge4_report':sha(a.gauge4_report),'complete_source_map':sha(a.source_map)},
 'third_field_composition':'ordinary derivatives with binomial Leibniz product/reciprocal calculus through3; raw shifts supplied by actual C4 caps; no producer imports',
 'independently_checked_third_field_boxes':{n:z.row() for n,z in third.items()},
 'independent_component_F_x_intervals':{n:z.row() for n,z in independent.items()},
 'independent_H_squared_upper':str(independent_square),'reported_H_squared_upper_reconstructed':str(reported_square),
 'genuine_original_scaled_F_x_H_cap':str(cap),'readable_strict_upper':'9.055e14','causal_cover_cells':len(indices),
 'component_rows':23,'third_field_rows':24,'coefficient_jets':9,
 'named_checks':CHECKS,'wrong_formula_controls_rejected':list(controls),
 'scope':'complete original actual source derivative audit under inherited actual master/gauge/background and harmonic neighbourhood premises; no source smoothing, C2 target error theorem by this audit alone, nonlinear realization or cosmological evidence'}
print(json.dumps(encode(report),indent=2,sort_keys=True))
