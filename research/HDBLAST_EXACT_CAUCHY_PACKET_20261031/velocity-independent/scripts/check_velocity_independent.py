#!/usr/bin/env python3
"""Independent exact Cauchy packet velocity algebra and gates, no producer imports."""
import argparse, hashlib, json, math, sys
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
BITS=128
SOURCE_MAP_PIN='ce31d4d0a720cc0779261184494a83aa59349fa0449f8a86cc1245bda1595a4b'
def need(ok,why):
    if not ok: raise RuntimeError(why)
class I:
    def __init__(self,lo=0,hi=None):
        lo=Q(lo);hi=lo if hi is None else Q(hi);need(lo<=hi,'INTERVAL_ORDER')
        self.lo=Q((lo*(1<<BITS)).__floor__(),1<<BITS)
        self.hi=Q((hi*(1<<BITS)).__ceil__(),1<<BITS)
    @staticmethod
    def cv(x):return x if isinstance(x,I) else I(x)
    def __add__(a,b):b=I.cv(b);return I(a.lo+b.lo,a.hi+b.hi)
    __radd__=__add__
    def __neg__(a):return I(-a.hi,-a.lo)
    def __sub__(a,b):return a+-I.cv(b)
    def __rsub__(a,b):return I.cv(b)+-a
    def __mul__(a,b):
        b=I.cv(b);v=[x*y for x in (a.lo,a.hi) for y in (b.lo,b.hi)];return I(min(v),max(v))
    __rmul__=__mul__
    def inv(a):need(a.lo>0 or a.hi<0,'DENOMINATOR_ZERO');return I(1/a.hi,1/a.lo)
    def __truediv__(a,b):return a*I.cv(b).inv()
    def __rtruediv__(a,b):return I.cv(b)*a.inv()
    def __pow__(a,n):
        if n<0:return a.inv()**(-n)
        if n==0:return I(1)
        if n%2==0:
            b=a**(n//2);return I(0 if b.lo<=0<=b.hi else min(b.lo*b.lo,b.hi*b.hi),max(b.lo*b.lo,b.hi*b.hi))
        return a*a**(n-1)
    def mag(a):return max(abs(a.lo),abs(a.hi))
    def minimum(a):return 0 if a.lo<=0<=a.hi else min(abs(a.lo),abs(a.hi))
    def row(a):return [a.lo,a.hi]
class J:
    def __init__(self,c):self.c=list(map(I.cv,c))
    @property
    def n(self):return len(self.c)-1
    def cv(a,b):return b if isinstance(b,J) else J([b]+[0]*a.n)
    def __add__(a,b):b=a.cv(b);return J([x+y for x,y in zip(a.c,b.c)])
    __radd__=__add__
    def __neg__(a):return J([-x for x in a.c])
    def __sub__(a,b):return a+-a.cv(b)
    def __rsub__(a,b):return a.cv(b)+-a
    def __mul__(a,b):b=a.cv(b);return J([sum((a.c[k]*b.c[n-k] for k in range(n+1)),I(0)) for n in range(min(a.n,b.n)+1)])
    __rmul__=__mul__
    def inv(a):
        out=[a.c[0].inv()]
        for n in range(1,a.n+1):out.append(-out[0]*sum((a.c[k]*out[n-k] for k in range(1,n+1)),I(0)))
        return J(out)
    def __truediv__(a,b):return a*a.cv(b).inv()
    def __rtruediv__(a,b):return a.cv(b)*a.inv()
    def __pow__(a,n):
        if n<0:return a.inv()**(-n)
        out=a.cv(1)
        for _ in range(n):out=out*a
        return out
    def dz(a):return J([(n+1)*a.c[n+1] for n in range(a.n)])
class P:
    def __init__(self,a,b=None,d=None):
        self.a=a;self.b=a if b is None else b
        self.d=[C(0,a.n)]*3 if d is None else d
    @property
    def n(self):return self.a.n
    def cv(a,b):return b if isinstance(b,P) else P(C(b,a.n))
    def __add__(a,b):b=a.cv(b);return P(a.a+b.a,a.b+b.b,[x+y for x,y in zip(a.d,b.d)])
    __radd__=__add__
    def __neg__(a):return P(-a.a,-a.b,[-x for x in a.d])
    def __sub__(a,b):return a+-a.cv(b)
    def __rsub__(a,b):return a.cv(b)+-a
    def __mul__(a,b):b=a.cv(b);return P(a.a*b.a,a.b*b.b,[x*b.b+a.a*y for x,y in zip(a.d,b.d)])
    __rmul__=__mul__
    def inv(a):
        aa=a.a.inv();bb=a.b.inv();return P(aa,bb,[-aa*x*bb for x in a.d])
    def __truediv__(a,b):return a*a.cv(b).inv()
    def __rtruediv__(a,b):return a.cv(b)*a.inv()
    def __pow__(a,n):
        if n<0:return a.inv()**(-n)
        out=a.cv(1)
        for _ in range(n):out=out*a
        return out
    def dz(a):return P(a.a.dz(),a.b.dz(),[x.dz() for x in a.d])
ELL=Q(1,10**6)
LABELS=[f'{f}:{l}' for f in ('K','L','N','J','w') for l in (0,2,4)]+[f'{f}:{l}' for f in ('v','u') for l in (2,4)]+['H:2','H:4','T:+','T:-']
GRAM=list(map(Q,['1','1/12','19/720','2','1/6','19/360','1','1/12','19/720','3','1/4','19/240','1','1/12','19/720','4/3','19/15','4/3','19/15','20/9','133/15','1/6','1/6']))
UP0=list(map(Q,['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27']))
UP1=[n*x for n,x in enumerate(UP0) if n];UP2=[n*x for n,x in enumerate(UP1) if n]
def C(x,n=5):return J([x]+[0]*n)
def poly(cs,z):
    ans=z*0
    for c in reversed(cs):ans=ans*z+c
    return ans
def expminus(x):
    need(x>=0,'NONNEGATIVE_EXP_ARGUMENT')
    t=x;k=0
    while t>1:t/=2;k+=1
    n=2*((BITS+31)//8)
    while Q(1,math.factorial(n+1))>Q(1,1<<(BITS+20)):n+=2
    upper=sum(((-t)**j/math.factorial(j) for j in range(n+1)),Q(0));lower=upper+(-t)**(n+1)/math.factorial(n+1)
    need(0<lower<=upper,'ALTERNATING_EXP_REMAINDER')
    ans=I(lower,upper)
    for _ in range(k):ans=ans*ans
    return ans
def expjetminus(t):
    v=I(expminus(t.c[0].hi).lo,expminus(t.c[0].lo).hi);ds=J([0]+t.c[1:]);pw=C(1,t.n);out=C(0,t.n)
    for n in range(t.n+1):out=out+pw*(((-1)**n)*v/math.factorial(n));pw=pw*ds
    return out
def registered(bg,q,e):
    rho,eta,p,a,alpha=[bg[k] for k in ('rho','eta','p','a','alpha')]
    ir=1/rho;ip=1/p;h=poly(UP1,eta)*ip;k=rho*(Q(5,2)*a-h);dy=lambda f:-ir*f.dz()
    u=alpha*ir**2*q;v=3*alpha*ir**3*ip*(q.dz()+k*q);vy=dy(v);ut=e*u;d=3*ir*ip;vt=e*(v-d*u)
    aa=-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192*ir**2)*u*u-(2*poly(UP2,eta)+Q(8,3)*p*p)*v*v-ir**2*(12*ut*ut+2*vt*vt)
    dd=ir**2*(36*u*u-2*v*v)
    fs=[rho**2*(aa/6+Q(4,3)*dd)/6,-rho**2*(aa/2+2*dd)/10,-rho**2*(aa-4*dd)/42]
    fy=[dy(f) for f in fs];zero=q*0;Y2=[Q(1,6),Q(1,2),Q(1)]
    # Direct all23 table; no imported packing or initial-algebra method.
    displacement=[zero]*6+[-2*u*u*c for c in Y2]+fs+[-Q(3,2)*ip*f for f in fy]+[zero]*8
    kk=e*ir/3*(18*u*u-v*v+d*u*v)
    pi=e*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v+(30*a-6*ir)*u*u+(2*ir-a)*v*v)
    velocity=[ELL*(2*rho*kk-8*u*ut)*c for c in Y2]
    velocity += [ELL*rho*((-26*u*dy(u)-2*a*u*u)*c-Q(3,2)*fy[j]-3*a*fs[j]) for j,c in enumerate(Y2)]
    velocity += [ELL*(-8*rho*kk+44*u*ut)*c for c in Y2]
    velocity += [ELL*(2*rho*kk+4*u*ut)*c for c in Y2]
    velocity += [ELL*(rho*pi+2*u*vt)*c for c in Y2]
    velocity += [zero,zero,ELL*(Q(11,2)*u*u-fs[1]/2),ELL*(11*u*u-fs[2]/2)]+[zero]*4
    # Explicit five mean overrides preserve the registered mean formulas.
    means_u=[zero,zero,-u*u/3,fs[0],-Q(3,2)*ip*fy[0]]
    means_v=[ELL*(2*rho*kk-8*u*ut)/6,ELL*rho*((-26*u*dy(u)-2*a*u*u)/6-Q(3,2)*fy[0]-3*a*fs[0]),ELL*(-8*rho*kk+44*u*ut)/6,ELL*(2*rho*kk+4*u*ut)/6,ELL*(rho*pi+2*u*vt)/6]
    for idx,uu,vv in zip((0,3,6,9,12),means_u,means_v):displacement[idx]=uu;velocity[idx]=vv
    need(len(displacement)==len(velocity)==23,'COMPLETE_23')
    return displacement,velocity
def sqrt_bounds(x):
    lo=sum((g*z.minimum()**2 for g,z in zip(GRAM,x)),Q(0));hi=sum((g*z.mag()**2 for g,z in zip(GRAM,x)),Q(0));grid=1<<BITS
    return [Q(math.isqrt((lo*grid*grid).__floor__()),grid),Q(math.isqrt((hi*grid*grid).__ceil__())+1,grid)]
def encode(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,I):return encode(x.row())
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x
def closure(root,pins):
    out={};need(root.is_dir() and not root.is_symlink(),'SOURCE_ROOT_TYPE')
    for rel,pin in pins.items():
        name=Path(rel);need(not name.is_absolute() and '..' not in name.parts,'SOURCE_MEMBER_PATH')
        p=root/name;need(p.is_file(),'SOURCE_MEMBER_MISSING:'+rel)
        for parent in (p,*p.parents):
            if parent==root:break
            need(not parent.is_symlink(),'SOURCE_SYMLINK:'+rel)
        actual=hashlib.sha256(p.read_bytes()).hexdigest();need(actual==pin,'SOURCE_MEMBER_HASH:'+rel);out[rel]=actual
    return out
def tail_proof():
    # q^(n)=exp(4-s) P_n(s) v^(n mod2), s=1/d, v=3-2r.
    # s'=-v*s^2, v'=-2, v^2=1-4/s give these exact recurrences.
    ps=[[Q(1)]]
    def plus(a,b):return [(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))]
    def times(a,b):
        out=[Q(0)]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]+=x*y
        return out
    def minus(a,b):return plus(a,[-x for x in b])
    for n in range(5):
        p=ps[-1];dp=[j*p[j] for j in range(1,len(p))]
        q=times([0,0,1] if n%2==0 else [0,-4,1],minus(p,dp))
        if n%2:q=plus(q,[-2*x for x in p])
        ps.append(q)
    s0=Q(4096,63);checks=[]
    for n,p in enumerate(ps):
        mon=[Q(0)]*(2*n)+[Q(1)]
        for sign in (-1,1):
            pol=plus(mon,[sign*x for x in p])
            shifted=[sum((math.comb(k,j)*c*s0**(k-j) for k,c in enumerate(pol) if k>=j),Q(0)) for j in range(len(pol))]
            need(all(x>=0 for x in shifted),'TAIL_POLYNOMIAL_DOMINANCE')
            checks.append({'order':n,'sign':sign,'shifted_coefficients':shifted})
        need(s0>=2*n,'TAIL_EXP_POLYNOMIAL_MONOTONICITY')
    return {'s_min':s0,'polynomials_ascending':ps,'shifted_nonnegative_coefficients':checks,'identity':'q^(n)=exp(4-s)*P_n(s)*(3-2r)^(n mod2); exact recurrence uses s_prime=-(3-2r)*s^2 and (3-2r)^2=1-4/s','inequality':'|q^(n)|<=s^(2n)*exp(4-s)<=s_min^(2n)*exp(4-s_min), n<=5, r in[1,65/64]','flat_endpoint':'Every derivative tends to zero as s tends to infinity.'}
class Evaluator:
    def __init__(self,root,pinpath):
        self.root=root;self.pinpath=pinpath;self.pinraw=pinpath.read_bytes();need(hashlib.sha256(self.pinraw).hexdigest()==SOURCE_MAP_PIN,'PINNED_SOURCE_MAP')
        self.pins=json.loads(self.pinraw);self.before=closure(root,self.pins)
        self.freeze=json.loads((root/'research/HDBLAST_INPUT_PRECISION_20261023/initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
        self.cs={k:list(map(Q,v)) for k,v in self.freeze['representation_background_coefficients_ascending'].items()}
    def background(self,r):
        z=ELL*J([r,1/ELL]+[0]*4);out={}
        for name in ('rho','eta','p','a','alpha'):
            nom=poly(self.cs[name],z);co=[I(*v) for v in self.freeze['shell_Taylor_coefficients'][name]];fifth=I(*self.freeze['uniform_fifth_normalized_ODE_derivative'][name])
            dc=[];rem=[]
            for n in range(6):
                dc.append(sum((math.comb(j,n)*(co[j]-self.cs[name][j])*z.c[0]**(j-n) for j in range(n,5)),I(0)))
                rem.append(math.comb(5,n)*fifth*z.c[0]**(5-n))
            ds=[J(dc),C(0),J(rem)];out[name]=P(nom+ds[0]+ds[2],nom,ds)
        return out
    def evaluate(self,left,right=None,tail=False):
        r=I(left,right);bg=self.background(r);rr=J([r,1/ELL]+[0]*4)
        if tail:
            smin=1/((r.hi-1)*(2-r.hi));ex=expminus(smin-4);q=[]
            for n in range(6):
                need(smin>=2*n,'TAIL_MONOTONICITY');cap=smin**(2*n)*ex.hi;q.append(I(0,cap) if n==0 else I(-cap,cap))
            q=J([q[n]/(ELL**n*math.factorial(n)) for n in range(6)])
        else:
            den=(rr-1)*(2-rr);den.c[0]=Q(1,4)-(r-Q(3,2))**2;q=expjetminus(1/den-4)
        e=expjetminus(ELL*rr)
        return registered(bg,P(q),P(e)),bg
def local(ev):
    (_,vv),_=ev.evaluate(Q(65,64))
    actual=[ELL*z.a.c[1] for z in vv]
    nominal=[ELL*z.b.c[1] for z in vv]
    channels=[[ELL*z.d[k].c[1] for z in vv] for k in range(3)]
    errors=[sum((ch[i] for ch in channels),I(0)) for i in range(23)]
    den=sqrt_bounds(actual);err=sqrt_bounds(errors)
    need(den[0]>0,'POSITIVE_LOCAL_VR_TARGET_DENOMINATOR')
    path=ev.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_PROFILES.json'
    old=json.loads(path.read_bytes());nodes=list(map(Q,old['nodes']));table=[[Q(x) for x in row] for row in old['profiles']['v']]
    need(old['labels']==LABELS and list(map(Q,old['gram']))==GRAM,'OLD_FULL23_LABELS_GRAM')
    need(nodes[1]==Q(65,64) and len(table)==len(nodes) and all(len(row)==23 for row in table),'COMPLETE_OLD_TABLE')
    delta=Q(1,2048);rc=nodes[1]
    need(nodes[0]<rc-delta<rc<rc+delta<nodes[2],'OLD_KERNEL_WINDOW_ISOLATES_NODE')
    kp=[Q(35,32)*x for x in [1,0,-3,0,3,0,-1]]
    mass=sum((x/Q(k+1) for k,x in enumerate(kp)),Q(0))
    need(mass==Q(1,2),'OLD_KERNEL_EXACT_HALF_MASS')
    left=[(table[1][i]-table[0][i])/(nodes[1]-nodes[0]) for i in range(23)]
    right=[(table[2][i]-table[1][i])/(nodes[2]-nodes[1]) for i in range(23)]
    oldvr=[mass*(a+b) for a,b in zip(left,right)]
    olderr=[a-b for a,b in zip(actual,oldvr)];oldnorm=sqrt_bounds(olderr)
    need(oldnorm[0]>0,'POSITIVE_OLD_LOCAL_VR_ERROR_LOWER')
    need(err[1]<=den[0]/10,'LOCAL_VR_RELATIVE_GATE')
    need(err[1]<oldnorm[0]/2,'LOCAL_VR_TRUE_ERROR_FACTOR_TWO_GATE')
    need(all(z.lo==z.hi==0 for z in channels[1]),'SHARED_Q_E_ZERO_CHANNEL')
    rows=[{'label':name,'Gram':g,'actual_v_r':a,'nominal_v_r':b,'paired_error_v_r':d,'paired_channels':[ch[i] for ch in channels],'old_left_secant':left[i],'old_right_secant':right[i],'old_smoothed_v_r':oldvr[i],'actual_minus_old_smoothed_v_r':olderr[i]} for i,(name,g,a,b,d) in enumerate(zip(LABELS,GRAM,actual,nominal,errors))]
    return {'point':['0','65/64'],'actual_v_r_H_interval':den,'new_error_v_r_H_interval':err,'new_over_target_upper':err[1]/den[0],'old_error_v_r_H_interval':oldnorm,'new_upper_over_old_lower':err[1]/oldnorm[0],'gates_pass':True,'kernel_half_mass':mass,'kernel_polynomial':kp,'kernel_width':delta,'kernel_window':[rc-delta,rc+delta],'adjacent_old_nodes':nodes[:3],'old_all_nodes_and_full_velocity_table_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'rows':rows}

def full_coverage(ev):
    old=json.loads((ev.root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json').read_text())
    thresholds=[Q(old['initial_error_H_caps'][k]) for k in ('v0','v0_r')]
    def leaf(left,right,tail,depth):
        (_,vv),bg=ev.evaluate(left,right,tail)
        ch=[[[ (math.factorial(j)*ELL**j*z.d[k].c[j]).mag() for z in vv] for j in range(2)] for k in range(3)]
        aa=[[(math.factorial(j)*ELL**j*z.b.c[j]).mag() for z in vv] for j in range(2)]
        need(all(v==0 for row in ch[1] for v in row),'GLOBAL_SHARED_Q_E_DIFFERENCE_EXACT_ZERO')
        return {'r':[left,right],'tail':tail,'depth':depth,'channel_component_caps':ch,'nominal_absolute_component_caps':aa,'actual_background_denominators':{k:bg[k].a.c[0] for k in ('rho','p')}}
    pieces=[(Q(1),Q(65,64),True)]+[(Q(k,64),Q(k+1,64),False) for k in range(65,112)]
    leaves=[]
    for idx,(left,right,tail) in enumerate(pieces):
        leaves.append(leaf(left,right,tail,0))
        if idx%8==0:print(json.dumps({'velocity_coverage_leaf':idx,'bits':BITS}),flush=True)
    def totals():
        ch=[[[max(l['channel_component_caps'][k][j][i] for l in leaves) for i in range(23)] for j in range(2)] for k in range(3)]
        aa=[[max(l['nominal_absolute_component_caps'][j][i] for l in leaves) for i in range(23)] for j in range(2)]
        hh=[[sqrt_bounds([I(x) for x in row])[1] for row in chan] for chan in ch]
        caps=[sum((hh[k][j] for k in range(3)),Q(0)) for j in range(2)]
        return ch,aa,hh,caps
    ch,aa,hh,caps=totals()
    while any(a>b for a,b in zip(caps,thresholds)) and len(leaves)<4096:
        # Split the leaf with largest exact normalized sum-of-channel H charge.
        viable=[(max(sum((sqrt_bounds([I(x) for x in l['channel_component_caps'][k][j]])[1] for k in range(3)),Q(0))/thresholds[j] for j in range(2)),idx) for idx,l in enumerate(leaves) if l['depth']<16 and not l['tail']]
        if not viable:break
        _,idx=max(viable);l=leaves.pop(idx);a,b=l['r'];m=(a+b)/2
        leaves[idx:idx]=[leaf(a,m,False,l['depth']+1),leaf(m,b,False,l['depth']+1)]
        ch,aa,hh,caps=totals()
        if len(leaves)%16==0:print(json.dumps({'adaptive_leaves':len(leaves),'error_H':[float(x) for x in caps]}),flush=True)
    need(len(leaves)<=4096 and all(l['depth']<=16 for l in leaves),'FROZEN_BUDGET')
    need(leaves[0]['r'][0]==1 and leaves[-1]['r'][1]==Q(7,4) and all(a['r'][1]==b['r'][0] for a,b in zip(leaves,leaves[1:])),'EXACT_WHOLE_FOOTPRINT_COVERAGE')
    return {'leaves':leaves,'leaves_count':len(leaves),'max_depth':max(l['depth'] for l in leaves),'exact_zero_region':[0,1],'channel_component_suprema':ch,'nominal_absolute_component_suprema':aa,'channel_error_H_caps':hh,'error_H_caps':caps,'nominal_absolute_H_caps':[sqrt_bounds([I(x) for x in row])[1] for row in aa],'retained_error_H_thresholds':thresholds,'initial_gates_pass':[a<=b for a,b in zip(caps,thresholds)],'scope':'Fresh whole [0,7/4] analytic nominal velocity v/v_r absolute and actual-minus-nominal paired bounds; no old P1/convolution bound used for the new candidate.'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);ap.add_argument('--bits',choices=[128,192],type=int,default=128);ap.add_argument('--coverage',action='store_true');ap.add_argument('--source-pins',type=Path,default=Path(__file__).resolve().parent.parent/'source-pins.json');args=ap.parse_args()
    global BITS;BITS=args.bits
    need(not args.out.exists(),'FRESH_OUTPUT_REQUIRED');self_before=Path(__file__).read_bytes();ev=Evaluator(args.source,args.source_pins)
    tail=tail_proof();loc=local(ev);coverage=full_coverage(ev) if args.coverage else None
    need(closure(args.source,ev.pins)==ev.before and args.source_pins.read_bytes()==ev.pinraw,'SOURCE_CLOSURE_AND_MAP_UNCHANGED_AFTER')
    need(Path(__file__).read_bytes()==self_before,'SELF_UNCHANGED_AFTER')
    report={'status':'PASS_INDEPENDENT_EXACT_VELOCITY_LOCAL_GATES','candidate':'W1_exact_Cauchy_packet','bits':BITS,'local':loc,'tail':tail,'source_sha256_before_and_after':ev.before,'source_pins_sha256':SOURCE_MAP_PIN,'executed_helper_sha256':hashlib.sha256(self_before).hexdigest(),'own_direct_all23_formula':True,'producer_scientific_imports':False,'global_coverage_performed':False}
    if coverage:
        report['global_coverage_performed']=True;report['coverage']=coverage
        report['status']='PASS_INDEPENDENT_EXACT_VELOCITY_INITIAL_GATES' if all(coverage['initial_gates_pass']) else 'INCONCLUSIVE_INDEPENDENT_VELOCITY_GLOBAL_GATES'
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(encode(report),indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':report['status'],'target_H':[float(x) for x in loc['actual_v_r_H_interval']],'new_H':[float(x) for x in loc['new_error_v_r_H_interval']],'old_H':[float(x) for x in loc['old_error_v_r_H_interval']],'ratio':float(loc['new_over_target_upper'])}),flush=True)
    if coverage:print(json.dumps({'error_H':[float(x) for x in coverage['error_H_caps']],'absolute_H':[float(x) for x in coverage['nominal_absolute_H_caps']],'leaves':coverage['leaves_count']}),flush=True)
if __name__=='__main__':main()
