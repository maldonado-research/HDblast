#!/usr/bin/env python3
"""Bounded, whole-footprint initial/free-wave approximation gate.

No source integration or bulk iteration. Rational coefficients specify a single
C1 initial function; cell boxes are only certificates for that function.
"""
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audit',HERE/'authenticate_inputs.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);out=importlib.util.module_from_spec(spec);spec.loader.exec_module(out);return out
R=Path('/workspace/hdblast-research-work')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--cells',type=int,default=128);ap.add_argument('--dependency-map',type=Path,help='JSON object with P,N,G,C capsule-root paths; embedded expected hashes authenticate them');a=ap.parse_args()
 dependency_map=json.loads(a.dependency_map.read_text()) if a.dependency_map else None
 roots={label:Path(dependency_map[label]).resolve() if dependency_map is not None else audit.ROOT/spec[0] for label,spec in audit.SPECS.items()}
 before=audit.run(dependency_map) # complete closures authenticated before numerical imports
 import math
 from fractions import Fraction as F
 init=module('initial_profiles',roots['N']/'numeric/evaluate_initial_profiles.py')
 adap=module('registered_adapter',roots['G']/'projection/inputs/source_adapters.py')
 op=module('registered_linear_operator',roots['G']/'projection/inputs/linear_operator.py')
 I,J,ell=init.I,init.J,init.ELL;data,pins=init.load();npow=40;power=1<<npow
 def midpoint(x):return (x.lo+x.hi)/2
 # Rational polynomial coefficients, selected once, shared in every profile.
 state={k:J([midpoint(I(*data['shell']['boxes'][k+'_b']['exact']))]+[0]*4) for k in ('rho','eta','p','a')};state['alpha']=J([1]+[0]*4)
 for order in range(4):
  rho,e,p,aa,alpha=[state[k] for k in ('rho','eta','p','a','alpha')];h=init.poly(init.UP[1],e)/p
  rhs={'rho':-rho*rho*aa,'eta':-rho*p,'p':-rho*(init.poly(init.UP[1],e)-4*aa*p),'a':1/rho+rho*p*p/3,'alpha':rho*(F(5,2)*aa-h)*alpha}
  for k in state:state[k].c[order+1]=I(midpoint(rhs[k].c[order]/(order+1)))
 coefficients={k:[midpoint(v) for v in value.c] for k,value in state.items()}
 def constant(v):return J([v]+[0]*4)
 def jetpoly(cs,z):
  out=constant(0)
  for c in reversed(cs):out=out*z+c
  return out
 def smoothstep(t):return 126*t**5-420*t**6+540*t**7-315*t**8+70*t**9
 def square_power(j):
  for _ in range(npow):j=j*j
  return j
 def nominal(r):
  rr=J([r,1/ell,0,0,0]);z=ell*rr
  bg={k:jetpoly(cs,z) for k,cs in coefficients.items()};rho,e,p,aa,alpha=[bg[k] for k in ('rho','eta','p','a','alpha')]
  h=init.poly(init.UP[1],e)/p;k=rho*(F(5,2)*aa-h);C=alpha/rho**2;B=3*alpha/(rho**3*p)
  if r.hi<=F(65,64):q=constant(0)
  else:
   d=(rr-1)*(2-rr);base=1-(1-4*d)/(1+(power-4)*d);q=square_power(base)
   if r.lo<F(66,64):q=q*smoothstep(64*rr-65)
  E=square_power(1/(1+z/power))
  return rho,e,p,aa,h,C,B,k,init.poly(init.UP[2],e),E,q
 def registered(bg):
  rho,e,p,aa,h,C,B,k,U2,E,q=bg;ir=1/rho;ip=1/p;dy=lambda j:-ir*j.dz()
  u=C*q;v=B*(q.dz()+k*q);vy=dy(v);ut=E*u;d=3*ir*ip;vt=E*(v-d*u)
  A=-2*vy**2-8*p*u*vy+8*p*(h+2*aa)*u*v+(16*p*p-192*ir**2)*u*u-(2*U2+F(8,3)*p*p)*v*v-ir**2*(12*ut*ut+2*vt*vt)
  D=ir**2*(36*u*u-2*v*v);sectors=[A/6+F(4,3)*D,A/2+2*D,A-4*D]
  fs=[-rho**2*sectors[j]/den for j,den in enumerate((-6,10,42))];fy=[dy(x) for x in fs]
  K=E*ir/3*(18*u*u-v*v+d*u*v)
  Pi=E*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v+(30*aa-6*ir)*u*u+(2*ir-aa)*v*v)
  # Proved exact Fmom=3K and Gmom-Fmom_y=p Pi. Only the latter
  # difference is used by the adapter; supplying Fmom_y=0,Gmom=pPi
  # evaluates that exact difference without numerical cancellation.
  source=adap.initial_frame_data(rho=rho,a=aa,p=p,u=u,uy=dy(u),ut=ut,v=v,vt=vt,Fmom=3*K,Gmom=p*Pi,Fmom_y=0,f=fs,fy=fy)
  us,vs=op.pack_initial(source,ell);ix=[op.IDX[k+':0'] for k in ('K','L','N','J','w')]
  us=[us[i] if isinstance(us[i],J) else constant(us[i]) for i in ix];vs=[vs[i] if isinstance(vs[i],J) else constant(vs[i]) for i in ix]
  return [x.val() for x in us],[(ell*x.dz()).val() for x in us],[x.val() for x in vs]
 def actual(r,tail=False):return registered((*init.background(r,data,'none'),init.bump(r,data,tail)))
 def sqrtup(q):
  grid=1<<112;n=math.isqrt((q*grid*grid).__ceil__())+1;return F(n,grid)
 def norm(caps):return sqrtup(sum(w*c*c for w,c in zip((1,2,1,3,1),caps)))
 caps=[[F(0)]*5 for _ in range(3)];rows=[]
 if a.cells%64 or a.cells<64:raise ValueError('cells multiple64')
 # Breakpoints of the cutoff are exact partition endpoints.
 for k in range(a.cells,7*a.cells//4):
  r=I(F(k,a.cells),F(k+1,a.cells));tail=r.hi<=F(65,64)
  av=actual(r,tail);nv=registered(nominal(r));errs=[[(x-y).mag() for x,y in zip(ar,nr)] for ar,nr in zip(av,nv)]
  for j in range(3):caps[j]=[max(x,y) for x,y in zip(caps[j],errs[j])]
  rows.append({'r':r.row(),'analytic_actual_flat_tail':tail,'component_error_caps':[[str(x) for x in rr] for rr in errs]})
 delta=[norm(c) for c in caps];efree=max((delta[0]+F(3,2)*delta[2])/3,delta[1]+delta[2])
 point=[]
 for r0 in (F(5,4),F(3,2),F(7,4)):
  av=actual(I(r0));nv=registered(nominal(I(r0)))
  point.append({'r':str(r0),'actual_initial_intervals':[[x.row() for x in row] for row in av],'nominal_initial_intervals':[[x.row() for x in row] for row in nv],'scalar_value_oracle_radius':str((av[0][-1].hi-av[0][-1].lo)/2),'scalar_radial_oracle_radius':str((av[1][-1].hi-av[1][-1].lo)/2)})
 after=audit.run(dependency_map)
 if before!=after:raise ValueError('capsules changed during evaluation')
 report={'status':'FAIL_PROVISIONAL_INITIAL_FREE_WAVE_ACCURACY_GATE_BOUND_INSUFFICIENT','target_initial_free_wave_error':'10000000000000000000000000000','actual_error_lower_bound_proved':False,'mathematical_impossibility_claimed':False,
 'representation':{'mean_order':['K0','L0','N0','J0','w0'],'gram_diagonal':['1','2','1','3','1'],'r_footprint':['0','7/4'],'zero_on':['0','65/64'],'ell':str(ell),'background_polynomial_variable':'z=ell*r','background_coefficients_ascending':{k:[str(x) for x in c] for k,c in coefficients.items()},'packet':'d=(r-1)(2-r); Qhat=[1-(1-4d)/(1+(2^40-4)d)]^(2^40) times cutoff; zero for r<=65/64','cutoff':'S(64r-65) on [65/64,66/64], S=126t^5-420t^6+540t^7-315t^8+70t^9; 1 above66/64','temporal_initial_factor':'Ehat=(1+z/(2^40))^(-2^40)','construction':'all shared u,v, Hamiltonian f, scalar zeta, momentum K/Pi, then registered initial_frame_data and pack_initial; no independent band midpoints','regularity':'Qhat C4 across cutoff, rational on each patch; Uhat0 C1, Vhat0 continuous; exact zero extension/reflection; L0 value identically zero','free_wave_candidate':'half of reflected Uhat0(r+x)+Uhat0(r-x), plus half integral of reflected Vhat0 from r-x to r+x','whole_cone_C1':'free wave from this single initial representation, D+ target cone; L0 uses odd extension, other four means even extension; zero for x+r<=1'},
 'error_budget':{'labels':['delta0_H','deltar_H','deltav_H'],'component_error_caps':[[str(x) for x in row] for row in caps],'H_error_caps':[str(x) for x in delta],'H_error_caps_display':[float(x) for x in delta],'free_wave_X_error_cap':str(efree),'free_wave_X_error_cap_display':float(efree),'formula':'max((delta0+(3/2)deltav)/3,deltar+deltav)','passes':efree<=10**28,'included_once':['retained shell/background interval uncertainty','actual background ODE integral/jet enclosures','rational packet and cutoff approximation','nominal background polynomial approximation','ordinary t velocity and scaled radial conversions','whole-cell variation and outward rational arithmetic'],'additional_error_terms':'none for this exact rational candidate; no bulk/forcing/iteration approximation claimed'},
 'point_preflight':point,'partition':{'step':str(F(1,a.cells)),'whole_cells':len(rows),'zero_region':'[0,1] actual and candidate exactly zero','rows':rows},
 'input_numeric_pins':pins,'manifest_closures':{k:{q:v[q] for q in ('root','manifest_sha256','external_pin_source','external_pin_source_sha256','files_verified','declared_bytes_verified')} for k,v in before['closures'].items()},'reauthenticated_after_evaluation':True,
 'limitations':['Failed upper error gate is not a lower bound on actual approximation error.','Point oracle radii retain input enclosure width at zero spatial width; finer spatial cells cannot remove that fixed part of this oracle.','No evaluated total response U, B1, source reintegration, invariant observable sign, or Einstein residual certificate is produced.']}
 a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'status':report['status'],'cells':len(rows),'caps':[float(x) for x in delta],'free_wave_error':float(efree)}))
if __name__=='__main__':main()
