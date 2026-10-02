#!/usr/bin/env python3
"""Pure symbolic audit: no physical source, spectrum, modes, or quadrature."""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import io
import json
from math import factorial
from pathlib import Path
import platform
import sympy as s
import metric_wkb as mw

p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
checks=[];mutations=[]
def require(condition,message):
    if not condition: raise RuntimeError(message)
def eq(label,left,right=0):
    difference=s.expand(left-right)
    require(difference==0,label+': '+str(s.factor(difference)))
    checks.append(label)
def reject(label,residual):
    require(s.expand(residual)!=0,'Mutation was not detected: '+label)
    mutations.append(label)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
# Definitions are the author's separately expanded directional derivation.
namespace={'__file__':str(Path(__file__).with_name('explore_contacts.py'))}
with contextlib.redirect_stdout(io.StringIO()):
    exec(compile(Path(__file__).with_name('explore_contacts.py').read_text(),namespace['__file__'],'exec'),namespace)
L,M,a,w,h,h1,h2,h3,h4,h5=[namespace[n] for n in ('L','M','a','w','h','h1','h2','h3','h4','h5')]
D=namespace['D'];R=namespace['R'];P=namespace['P'];S=namespace['S']
dR=namespace['dR'];dP=namespace['dP'];dS=namespace['dS'];U=namespace['U'];V=namespace['V'];u=namespace['u'];v=namespace['v'];sigma=namespace['sigma']
special=lambda e:s.expand(e.subs(M,2*L**2))
# Full generic-jet recurrence audited against an independently expanded chain rule.
mw.ONE=s.Integer(1);mw.HALF=s.Rational(1,2)
Lp,wp=s.symbols('Lp wp',positive=True)
hj=s.symbols('hj0:5',real=True)
ad=tuple(s.factorial(n)*Lp**(n+1) for n in range(5))
dad=tuple(sum(s.binomial(n,j)*ad[n-j]*hj[j] for j in range(n+1)) for n in range(5))
result=mw.generic_subtractions(wp**2-2*Lp**2,ad,dad,r=2)
translate={Lp:L,wp:w,**dict(zip(hj,[h,h1,h2,h3,h4]))}
reference={'R':R,'P':P,'S':S,'W2':U,'W4':V}
dreference={'R':dR,'P':dP,'S':dS,'W2':u,'W4':v}
for name in reference:
    eq('generic evaluator baseline '+name,result[name]['value'].subs(translate),special(reference[name]))
    eq('generic evaluator directional '+name,result[name]['delta'].subs(translate),special(dreference[name]))
for n,label in [(1,'first'),(2,'second')]:
    q=S;dq=dS
    for _ in range(n):q=D(q);dq=D(dq)
    eq('generic evaluator S '+label,result['S'][label].subs(translate),special(q))
    eq('generic evaluator delta S '+label,result['S']['delta_'+label].subs(translate),special(dq))
for n,label in [(1,'first'),(2,'second')]:
    baseline=U;direction=u
    for _ in range(n):baseline=D(baseline);direction=D(direction)
    eq('generic evaluator W2 '+label,result['W2'][label].subs(translate),special(baseline))
    eq('generic evaluator delta W2 '+label,result['W2']['delta_'+label].subs(translate),special(direction))
# Fixed-r graded conservation, before specializing M=2 L^2.
for label,rho,pressure in [('0',namespace['w']/2,namespace['k2']/(6*w)),('2',(L**2/w+namespace['J2'])/4,(namespace['c']*U+L**2/w+namespace['J2'])/4),('4',(U**2/w-L**2*U/w**2+namespace['J4'])/4,(namespace['c']*V+namespace['b']*U**2/w**3-L**2*U/w**2+namespace['J4'])/4)]:
    eq('exact subtraction Ward grade '+label,D(rho)-L*rho+3*L*pressure)
eq('exact full metric directional subtraction Ward',D(dR)-L*dR+3*L*dP-h1*R+3*h1*P)
# Full minimal bare-mode metric Ward, with fixed-state canonical response.
k,A,B=s.symbols('k A B',nonzero=True,real=True)
def Db(e):
    return D(e)+s.diff(e,A)*B+s.diff(e,B)*(-4*k*k*A-special(sigma)/k)
Rb=(2*k*k+3*L*L)/(4*k);Pb=(2*k*k/3-L*L)/(4*k)
dRb=(3*L*L*A-L*B+(L*h1+2*L*L*h)/k)/2-4*h*Rb
dPb=(-(4*k*k/3+L*L)*A-L*B+(L*h1-2*L*L*h)/k)/2-4*h*Pb
eq('bare baseline Ward',D(Rb)-L*Rb+3*L*Pb)
eq('bare fixed-state metric Ward',Db(dRb)-L*dRb+3*L*dPb+3*h1*(Rb+Pb))
# Exact convergent bridge contacts and continuum primitive integration.
contact_expect={'Q':(12*L**2*h-2*L*h1-h2)/24,'rho':-L*(61*L**2*h1+3*L*h2-3*h3)/120,'p':-(7*L**3*h1-44*L**2*h2-3*L*h3+3*h4)/360}
contact_inventory={}
for label,expr in namespace['contacts'].items():
    expr=special(expr);integ=0;terms=[]
    for piece in s.Add.make_args(expr):
        power=piece.as_powers_dict().get(w,0)
        require(power<=-5 and power%2==1,'Nonconvergent contact term: '+str(piece))
        n=-power;coef=piece/w**power
        primitive=s.sqrt(s.pi)*s.gamma((n-3)/2)/(4*s.gamma(n/2))*(2*L**2)**((3-n)/2)
        integ+=coef*primitive
    eq('continuum contact primitive '+label,s.simplify(integ),contact_expect[label])
    contact_inventory[label]={str(n):str(s.factor(expr.coeff(w,-n))) for n in range(5,16,2) if expr.coeff(w,-n)!=0}
# Finite-K rational primitives: differential and endpoint checks, no evaluation.
vvar,K,mass=s.symbols('v K mass',positive=True)
for n in range(5,16,2):
    polynomial=sum((-1)**ell*s.binomial((n-5)//2,ell)*vvar**(2*ell+3)/s.Integer(2*ell+3) for ell in range((n-5)//2+1))
    eq('finite-K primitive transformed derivative n'+str(n),s.diff(polynomial,vvar),vvar**2*(1-vvar**2)**((n-5)//2))
    eq('finite-K primitive zero endpoint n'+str(n),polynomial.subs(vvar,0))
    eq('finite-K primitive continuum endpoint n'+str(n),polynomial.subs(vvar,1),s.sqrt(s.pi)*s.gamma(s.Rational(n-3,2))/(4*s.gamma(s.Rational(n,2))))
# Pure coordinate shift: full subtraction varies as T times eta derivative.
T=s.symbols('T',real=True)
gauge={h:L*T,h1:L**2*T,h2:2*L**3*T,h3:6*L**4*T,h4:24*L**5*T,h5:120*L**6*T}
eq('coordinate-shift canonical forcing',special(sigma).subs(gauge))
eq('metric forcing is curvature response',special(sigma),-L**2*(6*(h2+2*L*h1-4*L**2*h)/L**2)/6)
for label,delta,baseline in [('R',dR,R),('P',dP,P),('S',dS,S)]:
    eq('fixed-comoving-band coordinate-shift '+label,special(delta).subs(gauge),T*special(D(baseline)))
pi=s.pi;Q0=1/(12*pi*pi);rho0=11/(960*pi*pi);p0=-rho0
cq=contact_expect['Q']/(2*pi*pi*L**2)
cr=contact_expect['rho']/(2*pi*pi*L**4)
cp=contact_expect['p']/(2*pi*pi*L**4)
eq('continuum physical Q gauge null',(-2*h*Q0+cq).subs(gauge))
eq('continuum physical rho gauge null',(-4*h*rho0+(h2+4*L*h1)*Q0/(2*L**2)+cr).subs(gauge))
eq('continuum physical p gauge null',(-4*h*p0-h2*Q0/(2*L**2)+cp).subs(gauge))
# Independently varied local anomaly and complete finite-band subtraction trace.
C=(D(D(S))-2*L*D(S)-2*L**2*S)/2-M*S+R-3*P
dC=(D(D(dS))-2*L*D(dS)-2*h1*D(S)-2*L**2*dS-2*h2*S)/2-2*M*h*S-M*dS+dR-3*dP
anomaly_direction=special(dC-4*h*C)
integ=0
for piece in s.Add.make_args(anomaly_direction):
    power=piece.as_powers_dict().get(w,0);require(power<=-5 and power%2==1,'Anomaly is not convergent')
    n=-power;integ+=piece/w**power*s.sqrt(s.pi)*s.gamma((n-3)/2)/(4*s.gamma(n/2))*(2*L**2)**((3-n)/2)
expected_anomaly=(-116*L**4*h+94*L**3*h1+47*L**2*h2-3*h4)/120
eq('full metric local anomaly continuum integral',s.simplify(integ),expected_anomaly)
# Separate curvature expression confirms the same metric local trace term.
deltaH=h1/L-h
deltad=D(deltaH)/L;deltae=D(deltad)/L;deltaf=D(deltae)/L
eq('cosmic-curvature anomaly directional variation',(232*deltaH-14*deltad-42*deltae-6*deltaf)*L**4/240,expected_anomaly)
# Mutation witnesses are symbolic nonidentities, not physical alternate runs.
reject('drop metric operator density prefactor',4*h*Rb)
reject('drop metric operator pressure prefactor',4*h*Pb)
reject('drop deltaL stress contact',L*h1/k)
reject('drop explicit metric mass contact',2*L**2*h/k)
reject('use scalar Q without metric contacts',-2*h*Q0+cq)
reject('drop W4 pressure subtraction term',special(metric_direction:=dP-namespace['metric'](namespace['PP']-namespace['CC']*namespace['VV']/4)))
reject('replace finite-band primitive by continuum',vvar**3/3-s.Rational(1,3))
reject('omit finite-band baseline rho+p Ward contact',3*h1*(Rb+Pb-special(R+P)))
reject('move positive reference with geometry',M*h/w)
reject('wrong directional trace fourth derivative sign',2*h4)
# Record symbolic inventory and pin every inherited source actually used.
root=Path(__file__).resolve().parent
inputs=[root/'reference_inputs'/name for name in ('forced_stress.py','COMMON_ACTION_SMOOTH_FRW.md','MATCHED_STRESS_RESPONSE_PROTOCOL.md','DESITTER_COMMON_ACTION_DERIVATION.md','derive_local_trace.py')]
inputs.append(root/'explore_contacts.py')
receipt={'status':'PASS','scope':'Pure exact analytic proof. No new physical source, mode, spectrum, stress, quadrature, or numerical response evaluation.','identity_count':len(checks),'identities':checks,'mutation_count':len(mutations),'mutations':mutations,'python':platform.python_version(),'sympy':s.__version__,'input_sha256':{str(path.relative_to(root)):sha(path) for path in inputs},'module_sha256':sha(Path(__file__).with_name('metric_wkb.py')),'source_sha256':sha(Path(__file__)),'contact_inventory':contact_inventory,'finite_trace_metric_inventory':{str(n):str(s.factor(anomaly_direction.coeff(w,-n))) for n in range(5,16,2) if anomaly_direction.coeff(w,-n)!=0},'physical_results_computed':False}
Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['status','identity_count','mutation_count','python','sympy','module_sha256','physical_results_computed']}))
