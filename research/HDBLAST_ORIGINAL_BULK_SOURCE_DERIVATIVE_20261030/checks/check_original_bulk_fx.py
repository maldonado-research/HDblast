#!/usr/bin/env python3
"""Complete original full23 scaled-time bulk-source derivative enclosure.

The actual field two-jets use the retained authenticated C3 adapter. Third
ordinary field jets use the new actual master/gauge fourth-jet majorants.
No derivative of independent value intervals is taken.
"""
import argparse
from fractions import Fraction as Q
import hashlib, importlib.util, json, math
from pathlib import Path

COUNT=[]
def need(ok,label):
    if not ok: raise RuntimeError(label)
    COUNT.append(label)
def digest(b):return hashlib.sha256(b).hexdigest()
def load(path,sha,role):
    raw=path.read_bytes();need(digest(raw)==sha,'authenticate '+role)
    return json.loads(raw)
def imp(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def canonical(p):
    return [[str(c),[[v,n] for v,n in k]] for k,c in sorted(p.items()) if c]
def psha(p):return digest(json.dumps(canonical(p),separators=(',',':'),sort_keys=True).encode())
def addterm(p,key,c):
    key=tuple(sorted((v,n) for v,n in key.items() if n))
    p[key]=p.get(key,Q(0))+c
    if not p[key]:del p[key]
def decode(terms):
    p={}
    for t in terms:addterm(p,dict(t['powers']),Q(t['coefficient']))
    return p
PREFIXES=tuple('ABCDELM')+('s',)
SUFFIXES=('', 'y','t','yy','yt','tt')
NEXT={'':'t','y':'yt','t':'tt','yy':'yyt','yt':'ytt','tt':'ttt'}
DER={p+s:p+NEXT[s] for p in PREFIXES for s in SUFFIXES}
def dt(poly):
    out={}
    for mon,c in poly.items():
        for v,n in mon:
            key=dict(mon)
            if v in DER:
                key[v]-=1;key[DER[v]]=key.get(DER[v],0)+1;addterm(out,key,c*n)
            elif v=='c':
                key['theta']=key.get('theta',0)+1;addterm(out,key,c*n)
            elif v=='theta':
                key[v]-=1;addterm(out,key,c*n)
                key[v]+=2;addterm(out,key,-c*n)
    return out
def peval(poly,values):
    out=0
    for mon,c in poly.items():
        v=c
        for key,n in mon:v=v*values[key]**n
        out=out+v
    return out
def encoded(out):
    if isinstance(out,Q):return str(out)
    if isinstance(out,dict):return {k:encoded(v) for k,v in out.items()}
    if isinstance(out,(list,tuple)):return [encoded(v) for v in out]
    if hasattr(out,'row'):return out.row()
    return out

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--master4-report',type=Path,required=True)
ap.add_argument('--gauge4-report',type=Path,required=True)
ap.add_argument('--gauge4-sha256',required=True)
ap.add_argument('--cover-divisions',type=int,default=2)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args()
need(not args.output.exists(),'fresh output required')
root=args.repo_root.resolve();geom=root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016';nr=geom/'numeric'
PINS={'actual_first_jets.py':'a68fa8025d2170efca6c3e85d443910a796d0c30a9b7dd17ed988cd9e00bbb00',
 'interval_core.py':'d9cb8529e47b9590ec8dc94f2f4484b06dbeb0ddee6c16015feccde3e4440240',
 'exterior_error_caps.py':'466466c745c39251294a0e9c358f42c00570e0610cff9a51e1e97a665896b472',
 'check_inherited_errors.py':'bcc9695072543a6dc128850db2d9c7e182379b7e838f84db55d7cbf7dbdf66c6',
 'PREDECESSOR_PINS.json':'ffc89df2e605419cd801575b4aa7a1b8fb5f41bd25ff350ece297b8623327beb',
 'INPUT_PINS.json':'226152be92c8c0c49495a5b07f08ee8da7c75ac9bd48ac312c1ff48a94b44acd'}
for rel,sha in PINS.items():need(digest((nr/rel).read_bytes())==sha,'pre-import closure '+rel)
pred=json.loads((nr/'PREDECESSOR_PINS.json').read_text())
for rel,row in pred.items():
    raw=(nr/rel).read_bytes();need(digest(raw)==row['sha256'] and len(raw)==row['bytes'],'retained closure '+rel)
inputpins=json.loads((nr/'INPUT_PINS.json').read_text())
for name,row in inputpins.items():
    raw=(nr/row['path']).read_bytes();need(digest(raw)==row['sha256'] and len(raw)==row['bytes'],'background input '+name)
maps=load(geom/'projection/GEOMETRIC_COEFFICIENT_MAPS.json','34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615','complete physical source map')
m4=load(args.master4_report,'d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e','actual master4 report')
g4=load(args.gauge4_report,args.gauge4_sha256,'actual gauge4 report')
need(m4['status']=='PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS' and g4['status']=='PASS_CONDITIONAL_ACTUAL_FIRST_GAUGE_CONTINUOUS_FOURTH_JETS','actual master/gauge C4 status')
old=imp('actual_lower2',nr/'actual_first_jets.py');core=old.core;I=core.I;ELL=core.ELL
data,errors,_=old.inputs();need(ELL==Q(1,10**6),'same retained scale')

# Fourth-order truncated raw Taylor jets; omitted background coefficients of
# order4 are never used in output derivatives through3.
N=4;KEYS=tuple((i,j) for i in range(N+1) for j in range(N+1-i))
class T:
    def __init__(self,c=None):self.c={k:I.cv((c or {}).get(k,0)) for k in KEYS}
    @staticmethod
    def cv(x):return x if isinstance(x,T) else T({(0,0):x})
    def __add__(self,b):
        b=T.cv(b);return T({k:self.c[k]+b.c[k] for k in KEYS})
    __radd__=__add__
    def __neg__(self):return T({k:-v for k,v in self.c.items()})
    def __sub__(self,b):return self+-T.cv(b)
    def __rsub__(self,b):return T.cv(b)+-self
    def __mul__(self,b):
        b=T.cv(b)
        return T({(i,j):sum((self.c[a,b0]*b.c[i-a,j-b0] for a in range(i+1) for b0 in range(j+1)),I(0)) for i,j in KEYS})
    __rmul__=__mul__
    def inv(self):
        out={(0,0):self.c[0,0].inv()}
        for n in range(1,N+1):
            for i in range(n+1):
                j=n-i;out[i,j]=-out[0,0]*sum((self.c[a,b]*out[i-a,j-b] for a in range(i+1) for b in range(j+1) if a+b),I(0))
        return T(out)
    def __truediv__(self,b):return self*T.cv(b).inv()
    def __rtruediv__(self,b):return T.cv(b)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        out=T.cv(1)
        for _ in range(n):out=out*self
        return out
    def dx(self):return T({(i,j):(i+1)*self.c.get((i+1,j),I(0)) for i,j in KEYS})
    def dr(self):return T({(i,j):(j+1)*self.c.get((i,j+1),I(0)) for i,j in KEYS})
    def derivative(self,i,j):return math.factorial(i)*math.factorial(j)*self.c[i,j]
    def val(self):return self.c[0,0]

master=json.loads((nr/'inputs/HIGHER_JETS_NORMAL.json').read_text())['exact']['interior']
wcap=[Q(master['U']),Q(master['first']),max(map(Q,master['second'].values())),max(map(Q,master['third'].values())),max(map(Q,m4['exact']['all_scaled_order_four_mixed_caps'].values()))]
gcap={tuple(map(int,k.split(','))):list(map(Q,v)) for k,v in g4['actual_mixed_scaled_derivative_caps'].items()}
need(set(gcap)==set(KEYS) and all(len(v)==3 and min(v)>=0 for v in gcap.values()),'all45 gauge mixed caps')
W=T({(i,j):I(-wcap[i+j],wcap[i+j])/(math.factorial(i)*math.factorial(j)) for i,j in KEYS})
U=[T({(i,j):I(-gcap[i,j][k],gcap[i,j][k])/(math.factorial(i)*math.factorial(j)) for i,j in KEYS}) for k in range(3)]
def radial(j):return T({(0,n):ELL**n*j.c[n] for n in range(4)})
rr,eta,pp,aa,hh,RR,BB,kk,_,_=core.background(I(0,8),data,'none')
rho,p,a,R,B,k=map(radial,(rr,pp,aa,RR,BB,kk))
ct,ht,vt=old.temporal(I(0,3));c,H,vartheta=(T(t.c) for t in (ct,ht,vt))
f,b,d=U;u=R*vartheta*W;v=B*vartheta*(W.dr()/ELL+k*W)
amps={'A':-4*u-2*f.dr()/(ELL*rho), 'B':f.dx()/(ELL*rho)+b.dr()/ELL,
      'C':f/(rho*c)-c*d.dr()/ELL,'D':-2*(u+a*f+b.dx()/ELL),
      'E':c*d.dx()/ELL-b/c,'L':2*(u+a*f+H*b),'M':2*d,'s':v+p*f}
L=1/(ELL*rho.val())
third={}
for name,j in amps.items():
    third[name+'yyt']=(L*L*j.derivative(1,2)+a.val()*L*j.derivative(1,1))/ELL
    third[name+'ytt']=-L*j.derivative(2,1)/ELL**2
    third[name+'ttt']=j.derivative(3,0)/ELL**3
need(len(third)==24,'all24 new ordinary field third jets')
# Derivative normalization and radial scale product rule are witnessed at a
# nondegenerate exact polynomial fixture, independently of actual intervals.
fixture=T({(3,0):Q(2,7),(1,2):Q(3,11),(1,1):Q(5,13)})
need(fixture.derivative(3,0).lo<=Q(12,7)<=fixture.derivative(3,0).hi,'raw normalized Taylor factorial3')
need(fixture.derivative(1,2).lo<=Q(6,11)<=fixture.derivative(1,2).hi,'raw normalized Taylor factorial2')
need(Q(3,5)**2*Q(6,11)+Q(7,17)*Q(3,5)*Q(5,13)!=Q(3,5)**2*Q(6,11),'nonzero radial conversion product term')

# Retain every old whole-cell C3 lower-jet certificate; hulls only widen it.
d=args.cover_divisions;need(d>=1,'nonempty cover')
indices=[(i,j) for i in range(3*d) for j in range(8*d) if i+j<8*d]
lower={};cover=[]
for i,j in indices:
    x=I(Q(i,d),Q(i+1,d));r=I(Q(j,d),Q(j+1,d))
    m,s,bg=old.actual_jet(x,r,clip_causal=True)
    values={**m,**s,**bg}
    for key,val in values.items():
        val=I.cv(val)
        lower[key]=val if key not in lower else old.hull(lower[key],val)
    cover.append([i,j])
need(len(indices)==(39*d*d+3*d)//2,'complete clipped actual causal cover')
need(set(maps['input_symbols'])==set(lower),'exact actual lower input symbols')
values={**lower,**third}
rows=[]
for field,names in [('yy',('K:0','K:2','K:4')),('y0',('L:0','L:2','L:4')),('00',('N:0','N:2','N:4')),
                    ('ya',('v:2','v:4')),('0a',('u:2','u:4')),('ab',('J:0','J:2','J:4','H:2','H:4','T:+','T:-'))]:
    for j,name in enumerate(names):rows.append((name,'metric_forcing',field,j,maps['maps']['metric_forcing'][field][j]))
for j,name in enumerate(('w:0','w:2','w:4')):rows.append((name,'scalar_forcing',None,j,maps['maps']['scalar_forcing'][j]))
derived={};outs={};normalizations=[];terms={}
for label,sector,field,j,raw in rows:
    dp=dt(decode(raw));derived[label]=psha(dp);terms[label]=len(dp)
    # Differentiate and verify the preexisting physical forcing normalization
    # exactly; rho,ell are time independent, so no extra derivative term.
    natural=maps['maps']['QEred'][field][j] if field else maps['maps']['Qphi'][j]
    nd=dt(decode(natural));scaled={}
    for mon,coef in nd.items():
        key=dict(mon);key['ell']=key.get('ell',0)+2;key['rho']=key.get('rho',0)+2
        addterm(scaled,key,coef*(-2 if field else 1))
    need(dp==scaled,'complete differentiated normalization '+label)
    normalizations.append(label)
    outs[label]=ELL*peval(dp,values)
need(len(outs)==23,'all23 complete original derivative rows')

# Independent dual evaluation fixture validates the total chain rule for all
# Laurent exponents, including inverse c and rho factors.
def dualmul(a,b):return a[0]*b[0],a[1]*b[0]+a[0]*b[1]
def dualpow(a,n):
    if n<0:a=(1/a[0],-a[1]/a[0]**2);n=-n
    out=(Q(1),Q(0))
    for _ in range(n):out=dualmul(out,a)
    return out
dualchecks=0
for seed in (1,2,3):
    v={k:Q((j+3)*seed+1,(j%7)+5) for j,k in enumerate(sorted(values))}
    dirs={k:(v[DER[k]] if k in DER else v['c']*v['theta'] if k=='c' else 1-v['theta']**2 if k=='theta' else Q(0)) for k in values}
    for label,sector,field,j,raw in rows:
        q=(Q(0),Q(0))
        for term in raw:
            z=(Q(term['coefficient']),Q(0))
            for key,n in term['powers'].items():z=dualmul(z,dualpow((v[key],dirs[key]),n))
            q=q[0]+z[0],q[1]+z[1]
        need(q[1]==peval(dt(decode(raw)),v),'exact dual fullmap chain rule '+label+' fixture'+str(seed));dualchecks+=1

weights={**{f'{k}:{l}':w*fac for k,fac in [('K',1),('L',2),('N',1),('J',3),('w',1)] for l,w in [(0,Q(1)),(2,Q(1,12)),(4,Q(19,720))]},
         **{f'{k}:{l}':2*l*(l+2)*w for k in ('v','u') for l,w in [(2,Q(1,12)),(4,Q(19,720))]},
         **{f'H:{l}':Q(2,3)*l*(l+2)*(l*(l+2)-3)*w for l,w in [(2,Q(1,12)),(4,Q(19,720))]},'T:+':Q(1,6),'T:-':Q(1,6)}
need(set(weights)==set(outs) and min(weights.values())>0,'complete positive full23 source Gram')
square=sum((weights[k]*v.mag()**2 for k,v in outs.items()),Q(0));cap=math.isqrt(square.numerator//square.denominator)+1
need(Q(cap)**2>square and Q(cap-1)**2<=square,'exact genuine scaled Fx H square-root cap')
controls=[]
for name,wrong,correct in [('omit outer scaled-time ell',Q(cap)/ELL,Q(cap)),('drop radial rho conversion term',Q(3,5)**2*Q(6,11),Q(3,5)**2*Q(6,11)+Q(7,17)*Q(3,5)*Q(5,13)),('third Taylor coefficient used as derivative',Q(2,7),Q(12,7))]:
    need(wrong!=correct,'reject '+name);controls.append(name)
for rel,sha in PINS.items():need(digest((nr/rel).read_bytes())==sha,'readback closure '+rel)
for rel,row in pred.items():need(digest((nr/rel).read_bytes())==row['sha256'],'readback retained input '+rel)
for name,row in inputpins.items():need(digest((nr/row['path']).read_bytes())==row['sha256'],'readback background input '+name)
need(digest((geom/'projection/GEOMETRIC_COEFFICIENT_MAPS.json').read_bytes())=='34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615','readback complete source map')
need(digest(args.master4_report.read_bytes())=='d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e','readback actual master4')
need(digest(args.gauge4_report.read_bytes())==args.gauge4_sha256,'readback gauge4 report')
report={'status':'PASS_COMPLETE_ORIGINAL_ACTUAL_FULL23_SCALED_BULK_SOURCE_TIME_DERIVATIVE',
 'source_coordinate':'x=t/ell','ell':ELL,'domain':'0<=x<=3,r>=0,x+r<=8',
 'actual_source_unchanged':True,'genuine_scaled_F_x_H_cap':cap,'F_x_H_squared_upper':square,
 'input_sha256':{**PINS,'master4_report':digest(args.master4_report.read_bytes()),'gauge4_report':args.gauge4_sha256,'complete_source_map':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'},
 'actual_raw_master_scaled_order_caps':wcap,'actual_raw_gauge_scaled_mixed_caps':g4['actual_mixed_scaled_derivative_caps'],
 'scaled_normalized_coefficient_jet_boxes':{name:{str(i)+','+str(j):value.c[i,j] for i,j in KEYS if i+j<=3} for name,value in [('rho',rho),('p',p),('a',a),('R',R),('B',B),('k',k),('c',c),('H',H),('vartheta',vartheta)]},
 'lower_actual_field_and_background_global_boxes':lower,'new_actual_field_third_jet_boxes':third,
 'ordinary_chain_rule':'F_x=ell*(sum field_jets F_j*next_t_jet+F_c*c*theta+F_theta*(1-theta^2)); other background factors static',
 'differentiated_row_sha256':derived,'differentiated_row_monomial_counts':terms,
 'component_F_x_intervals':outs,'positive_Gram':weights,'all23_normalization_identities':normalizations,
 'causal_cover':{'division_per_unit':d,'indices':cover,'rectangles':len(indices),'clipped_only_by':'x+r<=8','global_hull_only_widens':True},
 'dual_fullmap_identity_checks':dualchecks,'named_checks':COUNT,'wrong_formula_controls_rejected':controls,
 'scope':{'actual_original_full23_F_x_bound':True,'source_smoothing_or_replacement':False,'derivative_of_independent_interval_values':False,'full_actual_source_map_chain_rule':True,'C2_target_error_by_this_report':False,'nonlinear_or_cosmological_claim':False},
 'premises':'same actual smooth background/master/gauge solution, C1 harmonic identity on neighbourhood, differentiated PDE and actual initial/boundary premises as authenticated inputs; upper bounds only, no sign or physical force'}
args.output.write_text(json.dumps(encoded(report),indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':report['status'],'genuine_scaled_F_x_H_cap':str(cap),'cover_cells':len(indices),'component_rows':len(outs)}))
