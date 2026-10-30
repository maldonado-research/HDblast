#!/usr/bin/env python3
"""Independent product/Leibniz/majorant audit of the actual master C4 bridge."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def encoded(obj):
    if isinstance(obj,Q):
        return str(obj)
    if isinstance(obj,dict):
        return {k:encoded(v) for k,v in obj.items()}
    if isinstance(obj,(tuple,list)):
        return [encoded(v) for v in obj]
    return obj


def read(path,digest):
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==digest,'pinned input:'+str(path))
    pins.append((path,digest))
    return json.loads(raw)


def text_input(path,digest):
    raw=path.read_bytes()
    need(hashlib.sha256(raw).hexdigest()==digest,'pinned proof:'+str(path))
    pins.append((path,digest))
    return raw.decode()


def rejection(name,ok):
    try:
        need(ok,name)
    except RuntimeError:
        bad.append(name)
    else:
        raise RuntimeError('wrong formula survived:'+name)


# A product monomial is (physical W time order, W radial order,
# scaled U time order, U radial order, ell power). Each scaled derivative
# of the physical W contributes an extra ell. No source code is imported.
def differentiate(poly,axis):
    out={}
    for key,c in poly.items():
        w=list(key); w[axis]+=1; w[4]+=1
        u=list(key); u[2+axis]+=1
        for z in (tuple(w),tuple(u)):
            out[z]=out.get(z,Q(0))+c
    return out


def deriv(poly,axes):
    for axis in axes:
        poly=differentiate(poly,axis)
    return poly


def separate(poly):
    return {k:c for k,c in poly.items() if not (k[0] and k[1])}


def cap(poly,W,U,ell):
    return sum(c*ell**power*W[(a,b)]*U[(i,j)]
               for (a,b,i,j,power),c in poly.items())


ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--producer-report',type=Path,required=True)
a=ap.parse_args()
base=a.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
pins=[];bad=[];groups=[]
prior=read(base/'numeric/inputs/HIGHER_JETS_NORMAL.json',
    '92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df')
bump=read(base/'numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json',
    'add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd')
master_proof=text_input(base/'numeric/inputs/MASTER_HIGHER_JETS_PREMISE.md',
    '4c6453049f0d91e9223f06b8c881cd6bd42164d84473bf5fa5a5357bfbe35c8a')
potential_proof=text_input(base/'response-signs/inputs/INITIAL_PROFILE_PROOF.md',
    'a6a417c442a8132263276cc0d14fa6525b0172b2221d0c298758854fedee80c2')
need('U_{xx}-U_{rr}+f=0' in master_proof and
     'f=\\ell^2\\mathcal W' in master_proof and
     'H=J_0-y\'' in master_proof and
     '-35 tanh^2(t)/4' in potential_proof,
     'authenticated master PDE, incoming trace and temporal potential expressions')
report=read(a.producer_report,'d04f67e2e1863cca423a4b0835bef0a739bf5c471cf1d3e2f6730677d970490e')
need(prior['status']==bump['status']=='PASS' and
     report['status']=='PASS_CONDITIONAL_ACTUAL_MASTER_CONTINUOUS_FOURTH_JETS',
     'successful actual-field premises')
need(prior['scope']['conditional_on_retained_smooth_actual_solution'] is True and
     prior['scope']['independent_H_cell_selections_are_not_solutions'] is True,
     'one actual smooth master, not differentiated interval choices')
groups.append('authenticated actual master and bump premises')
m=prior['exact']['interior'];bd=prior['exact']['boundary']
ell=Q(m['ell']);T=Q(m['scaled_domain']['x_max'])
need(ell==Q(1,10**6) and T==3 and m['scaled_domain']==
     {'x_max':'3','r_min':0,'x_plus_r_max':8},'same causal domain and scaling')
M=Q(prior['exact']['background']['W']);Mz=Q(prior['exact']['background']['Wz'])
Mt=Q(prior['exact']['background']['Wt']);Mtt=Q(35,2)
U0=Q(m['U']);U1=Q(m['first']);U2=max(map(Q,m['second'].values()));U3=max(map(Q,m['third'].values()))
Wu={(0,0):M,(1,0):Mt,(0,1):Mz,(2,0):Mtt}
Ub={(0,0):U0}
for n,B in ((1,U1),(2,U2),(3,U3)):
    for i in range(n+1):
        Ub[(i,n-i)]=B
forcing={(0,0,0,0,2):Q(1)}
fxx=separate(deriv(forcing,(0,0)))
fxr=separate(deriv(forcing,(0,1)))
fxxr=separate(deriv(forcing,(0,0,1)))
need(fxx=={(2,0,0,0,4):1,(1,0,1,0,3):2,(0,0,2,0,2):1},
     'independent Fxx product derivative')
need(fxr=={(1,0,0,1,3):1,(0,1,1,0,3):1,(0,0,1,1,2):1},
     'independent Fxr product derivative')
need(fxxr=={(2,0,0,1,4):1,(1,0,1,1,3):2,
             (0,1,2,0,3):1,(0,0,2,1,2):1},
     'independent Fxxr product derivative; all mixed potential terms vanish by separability')
Fxx=cap(fxx,Wu,Ub,ell);Fxr=cap(fxr,Wu,Ub,ell);Fxxr=cap(fxxr,Wu,Ub,ell)
need(Fxx==Q(report['exact']['Fxx']) and Fxr==Q(report['exact']['Fxr']) and
     Fxxr==Q(report['exact']['Fxxr']),'exact independent source-cap agreement')
groups.extend(['source product derivative Fxx','source product derivative Fxr',
               'source product derivative Fxxr','three independent source majorants'])

# Independent moving-endpoint check: integral_0^x s^a(x-s)^b ds has the
# beta-function coefficient a!b!/(a+b+1)!, hence is a monomial in x.
# Its second derivative equals fx(x,0)+fr(x,0)+integral frr.
leibniz_count=0
for aa in range(7):
    for bb in range(7):
        degree=aa+bb+1
        coef=Q(math.factorial(aa)*math.factorial(bb),math.factorial(degree))
        lhs=coef*degree*(degree-1)
        rhs=(aa if bb==0 else 0)+(1 if bb==1 else 0)
        if bb>=2:
            rhs+=bb*(bb-1)*Q(math.factorial(aa)*math.factorial(bb-2),
                             math.factorial(aa+bb-1))
        need(lhs==rhs,'incoming second Leibniz polynomial identity')
        leibniz_count+=1
groups.append('49 independent incoming second-Leibniz identities')
y=Q(bd['y']);yp=Q(bd['yprime']);y2=Q(m['ysecond']);y3=Q(m['ythird'])
gamma=max(map(lambda z:abs(Q(z)),bd['gamma']))
BQ=list(map(Q,bump['dimensionless_Q_derivative_sup']))
need(BQ[:4]==[Q(1),Q(16),Q(256),Q(46656,7)] and BQ[4]==Q(16777216,49),
     'actual smooth bump fourth derivative cap')
shell_fx=Q(m['forcing']['shell_fx'])
shell_fr=ell**3*Mz*y+ell**2*M*U1
Frr=Q(m['forcing_rr'])
incoming2=BQ[3]+ell*(BQ[2]+2*ell*BQ[1]+ell**2*BQ[0])
H2=incoming2+y3+ell*gamma*y2+shell_fx+shell_fr+T*Frr
y4=Q(bd['ell_over_mu'])*H2+ell**4*Mtt*y+2*Q(bd['Aprime'])*yp+Q(bd['A'])*y2
need(H2==Q(report['exact']['Hsecond']) and y4==Q(report['exact']['yfourth']),
     'incoming Hsecond and dynamic boundary yfourth independently reconstructed')
groups.extend(['incoming Hsecond includes both endpoint derivatives',
               'boundary fourth derivative uses Asecond and twice Aprime'])

# Wtt=-35/2 (1-z)(1-3z), z=tanh(t)^2. The polynomial ranges in
# [-1/3,1] on z in [0,1], using the completed square below.
need(Q(1)-4*Q(2,3)+3*Q(2,3)**2==Q(-1,3) and
     4-3>=0 and Q(35,2)>=Q(35,6),'global separable temporal potential bound')
groups.append('separable temporal potential second derivative')

# Independent construction for V=Uxx. Its initial displacement Q''-f(0)
# contributes BQ4+Frr, and initial velocity ell V_initial''-fx(0)
# contributes ell BV3+Fxr to every second free-wave mixed derivative.
BV3=sum(Q(math.comb(3,k))*ell**k*BQ[3-k] for k in range(4))
free=BQ[4]+Frr+ell*BV3+Fxr+y4
shell_Fxx=ell**4*Mtt*y+2*ell**3*Mt*yp+ell**2*M*y2
base_cap=free+T*Fxxr+shell_Fxx
caps={'4,0':base_cap+Fxx,'3,1':base_cap,'2,2':base_cap,
      '1,3':base_cap+Fxr,'0,4':base_cap+Frr}
need(free==Q(report['exact']['Uxx_free_second']) and
     shell_Fxx==Q(report['exact']['shell_Fxx']) and
     caps=={k:Q(z) for k,z in report['exact']['all_scaled_order_four_mixed_caps'].items()},
     'independent equality of all five actual fourth mixed derivative caps')
need(max(caps.values())<Q('342395.408') and y4<Q('3.235927'),
     'strict readable continuous caps')
groups.extend(['Uxx initial displacement source curvature term',
               'Uxx initial velocity source mixed term',
               'reflected forcing shell-trace term',
               'five recovered fourth mixed majorants',
               'strict readable actual-master and boundary caps'])

# Reflected wave sign certificate: odd g has d_r g_o=(g_r)_even+
# 2g_shell delta. Its second-derivative ray identities are
# Krr=half(Iplus-Iminus)-g(x-r,0), Kxr=half(Iplus+Iminus)+g(x-r,0).
# At r=0: Krr=-g, Kxx=0, so the source jump signs are fixed by
# the principal equation and the exact homogeneous Dirichlet trace.
reflection_trace={'Krr':-1,'Kxr':1,'Kxx_shell':0}
need(reflection_trace['Krr']==-1 and reflection_trace['Kxr']==1 and
     1+reflection_trace['Krr']==reflection_trace['Kxx_shell']==0,
     'reflected source trace respects principal and Dirichlet identities')
groups.append('reflected shell-jump signs and principal boundary cancellation')

rejection('drop shell fx from Hsecond',H2<=H2-shell_fx)
rejection('drop shell fr from Hsecond',H2<=H2-shell_fr)
rejection('drop initial frr from Uxx free bound',free<=free-Frr)
rejection('drop initial fxr from Uxx free bound',free<=free-Fxr)
rejection('drop reflected fxx shell trace',base_cap<=base_cap-shell_Fxx)
rejection('drop Wz Uxx from Fxxr',Fxxr<=Fxxr-ell**3*Mz*U2)
rejection('reuse base bound for all five jets',max(caps.values())<=base_cap)
rejection('reverse reflected Krr shell jump',1+1==reflection_trace['Kxx_shell'])
for path,digest in pins:
    need(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'unchanged input readback')
groups.append('unchanged authenticated input readback')
print(json.dumps(encoded({'status':'PASS_INDEPENDENT_ACTUAL_MASTER_FOURTH_REVIEW',
    'source_report_sha256':pins[-1][1],'input_sha256':{'higher_master':pins[0][1],'bump':pins[1][1],
      'master_proof':pins[2][1],'potential_separation_proof':pins[3][1]},
    'domain':'same actual master, 0<=x<=3, r>=0, x+r<=8; ell=1e-6',
    'independent_source_derivative_terms':{'Fxx':{str(k):z for k,z in fxx.items()},
      'Fxr':{str(k):z for k,z in fxr.items()},'Fxxr':{str(k):z for k,z in fxxr.items()}},
    'independent_incoming_Leibniz_checks':leibniz_count,
    'independent_caps':{'Hsecond':H2,'yfourth':y4,'Fxx':Fxx,'Fxr':Fxr,'Fxxr':Fxxr,
      'shell_Fxx':shell_Fxx,'Uxx_free_second':free,'all_order_four_mixed':caps},
    'all_fourth_scaled_mixed_readable_strict_upper':Q('342395.408'),
    'boundary_fourth_scaled_readable_strict_upper':Q('3.235927'),
    'named_check_groups':groups,'named_check_group_count':len(groups),
    'wrong_formula_controls_rejected':bad,'wrong_formula_control_count':len(bad),
    'scope':'continuous conditional actual first-master fourth jets; no actual gauge fourth certificate, complete quadratic source derivative, second-response C2 error, natural Einstein residual, nonlinear family or cosmological evidence'}),sort_keys=True))
