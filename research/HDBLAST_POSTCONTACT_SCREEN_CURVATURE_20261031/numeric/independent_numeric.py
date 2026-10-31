#!/usr/bin/env python3
"""Independent aggregation and evaluation of postcontact TT source intervals.

The authenticated old actual_first_jets module is a shared model input. No new
postcontact producer, source evaluator, shell evaluator or screen evaluator is
imported. This file owns the integer interval arithmetic and all aggregation.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
GRID = 1 << 112
T, R, ELL = F(5,4), F(1,8), F(1,10**6)
U, L = T+R, T-R
LIFT_SIGN = -1
REFLECTED_LOWER_SIGN = 1
SCREEN_GAUGE_FACTOR = 1
FEEDBACK_KEEP_BD = 1

def need(ok, message):
    if not ok: raise ValueError(message)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def floor(q): return q.numerator // q.denominator
def ceil(q): return -floor(-q)

class Box:
    """Outward fixed-grid intervals, implemented using signed integers."""
    def __init__(self, lo=0, hi=None):
        hi = lo if hi is None else hi
        if hasattr(lo,'lo') and hi is lo:
            lo, hi = lo.lo, lo.hi
        lo, hi = F(lo), F(hi)
        need(lo <= hi, 'BOX_ORDER')
        self.a, self.b = floor(lo*GRID), ceil(hi*GRID)
    @classmethod
    def ints(cls,a,b):
        out=object.__new__(cls);out.a,out.b=a,b;return out
    @property
    def lo(self): return F(self.a,GRID)
    @property
    def hi(self): return F(self.b,GRID)
    @staticmethod
    def cv(v): return v if isinstance(v,Box) else Box(v)
    def __add__(self,v):
        v=self.cv(v);return Box.ints(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self): return Box.ints(-self.b,-self.a)
    def __sub__(self,v): return self+-self.cv(v)
    def __rsub__(self,v): return self.cv(v)+-self
    def __mul__(self,v):
        v=self.cv(v);products=[x*y for x in (self.a,self.b) for y in (v.a,v.b)]
        a,b=min(products),max(products)
        return Box.ints(a//GRID,-((-b)//GRID))
    __rmul__=__mul__
    def inv(self):
        need(self.a>0 or self.b<0,'BOX_ZERO_DENOMINATOR')
        return Box.ints((GRID*GRID)//self.b,-((-GRID*GRID)//self.a))
    def __truediv__(self,v): return self*self.cv(v).inv()
    def __rtruediv__(self,v): return self.cv(v)*self.inv()
    def __pow__(self,n):
        need(isinstance(n,int),'INTEGER_POWER')
        if n<0:return self.inv()**(-n)
        if n==0:return Box(1)
        if n%2==0:
            half=self**(n//2)
            low=0 if half.a<=0<=half.b else min(half.a**2,half.b**2)
            high=max(half.a**2,half.b**2)
            return Box.ints(low//GRID,-((-high)//GRID))
        return self*(self**(n-1))
    def mag(self): return max(abs(self.lo),abs(self.hi))
    def row(self): return [str(self.lo),str(self.hi)]

def encoded(v):
    if isinstance(v,Box):return v.row()
    if isinstance(v,F):return str(v)
    if isinstance(v,dict):return {k:encoded(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [encoded(x) for x in v]
    return v

def lift_dx(g): return LIFT_SIGN*g
def endpoint_pair(plus,minus): return (plus+REFLECTED_LOWER_SIGN*minus)/2
def feedback_sigma(a,b,d): return a+FEEDBACK_KEEP_BD*(b+T*d/2)
def quadratic(m,b):
    rho,c,theta=b['rho'],b['c'],b['theta'];rad=rho*c
    return ((m['D']+4*m['M']-m['L'])*m['Mt']-SCREEN_GAUGE_FACTOR*2*m['M']*m['Lt'])/(4*rho)-m['E']*(m['D']+m['L'])/(2*rad)+(m['C']/rad-m['My']/2)*(rad*(m['Ct']-theta*m['C'])/(2*rho)-rad*m['Ey']/2+m['B']/2)

def elementary_guards():
    # These fixtures test live computational functions and yield exact failed
    # witnesses under source mutations. They are independent of scientific data.
    need(lift_dx(F(7,3))==F(-7,3),'LIFT_SIGN_WITNESS:got='+str(lift_dx(F(7,3)))+';expected=-7/3')
    v=T*endpoint_pair(F(1),F(1))
    need(v==T,'REFLECTED_KERNEL_WITNESS:got='+str(v)+';expected=5/4')
    m={k:F(0) for k in 'B C D E L M Ct Ey Lt Mt My'.split()};m.update(M=F(1),Lt=F(1))
    b={'rho':F(1),'c':F(1),'theta':F(0)};v=quadratic(m,b)
    need(v==F(-1,2),'GAUGE_QUADRATIC_WITNESS:got='+str(v)+';expected=-1/2')
    v=feedback_sigma(F(1,7),F(2,9),F(3,11))
    expected=F(1,7)+F(2,9)+T*F(3,11)/2
    need(v==expected,'SCALAR_FEEDBACK_WITNESS:got='+str(v)+';expected='+str(expected))
    # In particular interval inverse and even powers are checked independently.
    for a,b in ((-7,-2),(2,7),(-3,5)):
        x=Box(F(a,3),F(b,3));y=x*x
        need(y.lo<=min(F(a*a,9),F(b*b,9),F(0) if a<=0<=b else F(a*a,9)) and y.hi>=max(F(a*a,9),F(b*b,9)), 'INTERVAL_MULTIPLICATION')

def authenticate(repo,pins):
    inventories={}
    need(repo.is_dir() and not repo.is_symlink(),'ROOT_TYPE')
    for label,pin in pins.items():
        if label.startswith('research/'):
            p=repo/label;need(p.is_file() and not p.is_symlink() and digest(p)==pin,'SELECTED_SOURCE:'+label);continue
        root=repo/'research'/label
        need(root.is_dir() and not root.is_symlink(),'CAPSULE_ROOT:'+label)
        manifest=root/'MANIFEST.json';need(not manifest.is_symlink() and digest(manifest)==pin,'CAPSULE_MANIFEST:'+label)
        rows=json.loads(manifest.read_text())['files'];seen={'MANIFEST.json'}
        for row in rows:
            rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SAFE_CLOSED_MEMBER:'+label)
            seen.add(str(rel));p=root/rel
            need(p.is_file() and not p.is_symlink() and p.stat().st_size==row['bytes'] and digest(p)==row['sha256'],'CLOSED_MEMBER:'+label+'/'+str(rel))
            need(all(not q.is_symlink() for q in p.parents),'MEMBER_ANCESTOR_SYMLINK')
        actual=set()
        for p in root.rglob('*'):
            need(not p.is_symlink(),'TREE_SYMLINK:'+label)
            if p.is_file():actual.add(str(p.relative_to(root)))
        need(actual==seen,'CLOSED_INVENTORY:'+label)
        inventories[label]={'manifest_sha256':pin,'members_including_manifest':len(seen)}
    return inventories

def load_jets(repo):
    p=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py'
    spec=importlib.util.spec_from_file_location('numeric_independent_common_actual_jets',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    need(m.ELL==ELL and m.core.GRID==GRID,'COMMON_INPUT_UNITS')
    return m

class Evaluation:
    def __init__(self,repo,jets):
        self.jets=jets;self.calls=0
        mp=repo/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
        maps=json.loads(mp.read_text())['maps'];rows=maps['metric_forcing']['ab']
        need(rows[5]==rows[6] and len(rows[5])==78,'COMPLETE_IDENTICAL_TT_SOURCE')
        self.terms=rows[5]
        self.variables=sorted({k for row in self.terms for k in row['powers']})
    def fields(self,x,r):
        self.calls+=1;need(self.calls<=140,'FROZEN_MAXIMUM_140_COMMON_JET_CALLS')
        m,s,b=self.jets.actual_jet(self.jets.I(*x),self.jets.I(*r),mode='actual')
        return tuple({k:Box(v.lo,v.hi) if hasattr(v,'lo') else Box(v) for k,v in row.items()} for row in (m,s,b))
    def forcing(self,x,r):
        m,s,b=self.fields(x,r);env={**m,**s,**b};out=Box(0)
        # Literal Laurent monomial evaluation, independent of old encoded
        # evaluator and main producer's dyadic implementation.
        for row in self.terms:
            term=Box(F(row['coefficient']))
            for name,power in sorted(row['powers'].items()):term=term*(env[name]**power)
            out=out+term
        return out
    def shell(self,x):
        m,s,b=self.fields(x,(F(0),F(0)))
        # ab/5 and ab/6 have trace weight0 and Hessian weight2/3.
        return -F(4,3)*ELL*b['rho']*(b['a']*m['A']-F(2,3)*b['p']*s['s'])*m['M']

def operator(jets):
    data,_,_=jets.inputs()
    rho,eta,p,a,*_=jets.core.background(jets.I(0,U),data,'none')
    rho,eta,p,a=[Box(v.val().lo,v.val().hi) for v in (rho,eta,p,a)]
    cc,hh,_=jets.temporal(jets.I(0,T));cc,hh=[Box(v.val().lo,v.val().hi) for v in (cc,hh)]
    ay=-1/rho**2-p**2/3;b=-(ay+a**2);d=1/rho**2-a**2
    V=Box(0)
    for coef in reversed(jets.core.U):V=V*eta+coef
    A=-3*ELL*hh;B=-3*ELL*rho*a
    D=ELL**2*(-6/cc**2+rho**2*(-2*a**2+2*(hh/rho)**2-2*b-8*d+F(4,3)*V))
    simple=-8*ELL**2/cc**2
    caps=[z.mag() for z in (A,B,D)];sigma=feedback_sigma(*caps);q=T*sigma
    need(0<q<1,'FULL_FEEDBACK_CONTRACTION')
    return {'A':A,'B':B,'D_original':D,'D_field_Hamilton_simplification':simple,'caps':dict(zip(('A','B','D'),caps)),'sigma':sigma,'q':q,'simplified_D_within_original_D':D.lo<=simple.lo<=simple.hi<=D.hi}

def source_cover(ev):
    rows=[];cap=F(0);skipped=[]
    # Freeze before computation: original8x8grid, projected intersections with
    # 0<=s<=T, max(0,s-L)<=v<=U-s. No adaptive subdivision.
    for i in range(8):
        a,b=T*i/8,T*(i+1)/8
        for j in range(8):
            c,d=U*j/8,U*(j+1)/8
            xb=min(b,U-c,L+d)
            if xb<a:skipped.append([i,j]);continue
            rlo=max(c,F(0),a-L);rhi=min(d,U-a)
            if rhi<rlo:skipped.append([i,j]);continue
            value=ev.forcing((a,xb),(rlo,rhi));cap=max(cap,value.mag())
            rows.append({'grid':[i,j],'x':[a,xb],'r':[rlo,rhi],'F_T':value})
    return rows,cap,skipped

def campaign(ev):
    op=operator(ev.jets);cover,f,skipped=source_cover(ev)
    print(json.dumps({'stage':'independent_clipped_cover','cells':len(cover)}),flush=True)
    shell=[];g=F(0)
    for i in range(8):
        a,b=L*i/8,L*(i+1)/8;value=ev.shell((a,b));g=max(g,value.mag());shell.append({'x':[a,b],'G':value})
    endpoint=ev.shell((L,L));curve=[];integral=Box(0)
    knots=sorted({T*i/32 for i in range(33)}|{L})
    for a,b in zip(knots,knots[1:]):
        positive=(U-b,U-a);reflected=tuple(sorted((abs(a-L),abs(b-L))))
        vp=ev.forcing((a,b),positive);vm=ev.forcing((a,b),reflected)
        term=(b-a)*endpoint_pair(vp,vm);integral=integral+term
        curve.append({'x':[a,b],'positive_r':positive,'reflected_r':reflected,'F_positive':vp,'F_reflected':vm,'contribution':term})
    seed=integral+lift_dx(endpoint)
    full_norm=(2*g+T*f)/(1-op['q']);tail=op['q']*full_norm
    ux=seed+Box(-tail,tail)
    m,s,bg=ev.fields((T,T),(R,R));fq=quadratic(m,bg)
    scale=F(3,2)/(ELL*bg['rho']**2*bg['c']);quad=4*fq/(bg['rho']*bg['c'])
    linear=scale*ux;h=linear+quad;sign='negative' if h.hi<0 else ('positive' if h.lo>0 else 'inconclusive')
    need(len(curve)==33 and sum((row['x'][1]-row['x'][0] for row in curve),F(0))==T,'COMPLETE_CHARACTERISTIC_COVER')
    need(scale.lo>0,'POSITIVE_SCREEN_SCALE')
    return {'status':'PASS_INDEPENDENT_NUMERICAL_ENCLOSURE' if sign!='inconclusive' else 'INDEPENDENT_VALID_ENCLOSURE_SIGN_INCONCLUSIVE','sign':sign,'point':[T,R],'ell':ELL,'operator':op,'source_cover':cover,'source_cover_omitted_empty_grid_cells':skipped,'source_cap':f,'shell_history':shell,'shell_cap':g,'shell_endpoint':endpoint,'characteristics':curve,'characteristic_integral':integral,'seed_x':seed,'full_C1_norm_upper':full_norm,'complete_feedback_tail_Ux':tail,'actual_U_Tx':ux,'apex_metric':m,'apex_scalar':s,'apex_background':bg,'screen_Fquad':fq,'screen_quadratic':quad,'screen_scale':scale,'screen_linear':linear,'screen_h':h,'screen_opposite_orientation':-h,'actual_jet_calls':ev.calls,'fixed_budget':{'grid':[8,8],'clipped_cover_cells':len(cover),'curve_segments':33,'curve_source_cells':66,'shell_history_cells':8,'shell_endpoint_cells':1,'apex_cells':1,'maximum_calls':140,'adaptive_refinement':False},'independence':'Common authenticated old actual_first_jets only; own integer interval, literal78term map evaluation, projected cone cover, shell formula, endpoint quadrature, feedback and compact fixed quadratic. No new scientific producer import.','limits':['Conditional inherited actual first-solution C3 and background/response regularity premises retained.','D simplification conditional on field and Hamilton constraints; primary tail uses unsimplified D.','Specified slicing observer and screen; no gauge-invariance or finite-amplitude nonlinear-contact claim.','Independent enclosure uses clipped fixed grid, so its numerical interval is not expected to equal the producer overcover interval.']}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--guards-only',action='store_true');args=ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    elementary_guards()
    if args.guards_only:
        args.output.write_text(json.dumps({'status':'PASS_LIVE_EXACT_IMPLEMENTATION_GUARDS'})+'\n');return
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text());repo=args.source_root.resolve();before=authenticate(repo,pins);ownpin=digest(Path(__file__))
    note=(repo/'research/HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015/background/inputs/PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md').read_text()
    need('a_y=-\\rho^{-2}-p^2/3' in note and 'a^2=\\rho^{-2}+p^2/12-U/6' in note,'BACKGROUND_FIELD_HAMILTON_PREMISES')
    need(repo not in args.output.resolve().parents,'OUTPUT_OUTSIDE_SOURCE')
    need(all(not p.is_symlink() for p in args.output.parents),'OUTPUT_ANCESTOR_SYMLINK')
    result=campaign(Evaluation(repo,load_jets(repo)))
    need(authenticate(repo,pins)==before and digest(Path(__file__))==ownpin,'POST_AUTHENTICATION')
    result.update(schema='HDBLAST_INDEPENDENT_POSTCONTACT_TT_NUMERIC_V1',source_authentication=before,source_pins=pins,checker_sha256=ownpin,all_sources_authenticated_before_import_and_after=True)
    args.output.write_text(json.dumps(encoded(result),indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'h_display_only':[float(result['screen_h'].lo),float(result['screen_h'].hi)],'calls':result['actual_jet_calls'],'sha256':digest(args.output)}),flush=True)

if __name__=='__main__':main()
