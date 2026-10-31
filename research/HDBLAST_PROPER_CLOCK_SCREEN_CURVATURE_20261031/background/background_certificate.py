#!/usr/bin/env python3
"""Exact-rational background shoot/clock/Jacobian certificate; no float solve.

Public API: certify(repo_root, injection='none') returns canonical bound rows.
Every repository path is relative. authenticate() completes before source import.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
NUM=Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric')
POT=Path('research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential')

def need(ok,msg):
 if not ok: raise ValueError(msg)

def digest(path):
 raw=path.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def authenticate(repo_root,pins_path=None):
 root=Path(repo_root).resolve();path=HERE/'SOURCE_PINS.json' if pins_path is None else Path(pins_path)
 pins=json.loads(path.read_text())
 need(pins['paths_are_relative_to_repo_root'] is True,'relative repository pin contract')
 for rel,pin in pins['files'].items():
  p=Path(rel);need(not p.is_absolute() and '..' not in p.parts,'unsafe repository pin path')
  need(digest(root/p)==pin,'repository authentication before import: '+rel)
 # The source closure is also tied to the original read-only capsule manifests.
 for capsule in (NUM.parent,POT.parent):
  manifest=json.loads((root/capsule/'MANIFEST.json').read_text())
  indexed={r['path']:r for r in manifest['files']}
  for rel,pin in pins['files'].items():
   try: sub=str(Path(rel).relative_to(capsule))
   except ValueError: continue
   if sub=='MANIFEST.json':continue
   need(sub in indexed,'pinned source missing in retained manifest: '+rel)
   need(indexed[sub]['bytes']==pin['bytes'] and indexed[sub]['sha256']==pin['sha256'],'retained manifest pin mismatch: '+rel)
 return root,pins

def load_background(repo_root,pins_path=None):
 root,pins=authenticate(repo_root,pins_path)
 spec=importlib.util.spec_from_file_location('proper_clock_authenticated_core',root/NUM/'interval_core.py')
 core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
 data,_=core.load()
 return core,data,pins

def row(lo,hi=None):return [str(F(lo)),str(F(lo if hi is None else hi))]
def record(iv):return row(iv.lo,iv.hi)

def certify(repo_root,injection='none'):
 need(injection in ('none','wrong_acceleration_sign','wrong_clock_unit','corrupt_domain_margin'),'injection mode')
 freeze=json.loads((HERE/'STAGE_FREEZE.json').read_text())
 core,data,pins=load_background(repo_root);I=core.I;ell=core.ELL
 need(ell==F(freeze['ell']),'frozen length unit')
 T=F(freeze['target']['t']);Z=F(freeze['target']['z']);H=F(freeze['budgets']['coordinate_time_horizon'])
 need(T==5*ell/4 and Z==ell/8,'frozen target coordinates')
 # Includes all potentially successful nonnegative launches, in addition to
 # the strict compact Gaussian label tube, so uniqueness is not only local.
 bg=core.background(I(0,F(1,4)),data,'none')
 rho,eta,p,a=bg[:4];rho,p,a=rho.val(),p.val(),a.val()
 b=rho*a; bz=1-b*b+rho*rho*p*p/3
 pot=json.loads((Path(repo_root)/POT/'ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
 need(pot['status']=='PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ','actual potential/background premise status')
 for k,iv in (('rho',rho),('p',p),('a',a)):
  outer=list(map(F,pot['background_intervals'][k]['exact']))
  need(outer[0]<=iv.lo<=iv.hi<=outer[1],'local actual background remains in authenticated wider collar: '+k)
 m=b.lo;M=b.hi;K=F(freeze['budgets']['abs_b_z'])
 need(F(freeze['budgets']['b_lower'])<m<M<F(freeze['budgets']['b_upper']),'frozen strict acceleration bounds')
 need(bz.hi<0 and bz.mag()<K,'strict signed Jacobi coefficient bound')
 acceleration_sign=1 if injection=='wrong_acceleration_sign' else -1
 chi=I(acceleration_sign*M*T,acceleration_sign*m*T) if acceleration_sign<0 else I(m*T,M*T)
 need(chi.hi<0,'inward launch must have negative rapidity; wrong-sign control')
 # tanh u <= u, tanh u >= u-u^3/3. These follow by differentiating
 # u-tanh u, then using tanh^2 u <= u^2; u>=0.
 shift_lo=m*T*T/2-M**3*T**4/12;shift_hi=M*T*T/2
 need(0<shift_lo<shift_hi<Z,'strict positive launch and radial displacement')
 launch=I(Z-shift_hi,Z-shift_lo)
 need(launch.lo>0 and ell/16<launch.lo<launch.hi<3*ell/16,'root inside compact label tube')
 # Since b bounds are outward, the endpoint inequalities are weak; enlarge
 # the bracket by one rational output-grid unit to make the IVT signs strict.
 q=F(1,core.GRID);shoot_bracket=I(launch.lo-q,launch.hi+q)
 endpoint_residuals={'lower_upper':str(-q),'upper_lower':str(q)}
 jlo=1-K*H*H
 need(K*H*H<1 and jlo>0,'Jacobi first-exit closure |j|<2 and positive j')
 # sech u >= 1-u^2/2, established by cosh u*(1-u^2/2)<=1
 # (or its derivative). Integrated against positive rho bounds.
 clock_lo=rho.lo*(T-M*M*T**3/6);clock_hi=rho.hi*T
 if injection=='wrong_clock_unit':clock_lo=(T-M*M*T**3/6)/rho.hi;clock_hi=T/rho.lo
 need(128*T<clock_lo<=clock_hi<132*T,'physical proper clock unit dτ/dt=ρ sechχ')
 tau=I(clock_lo,clock_hi)
 tau_max=F(freeze['compact_gaussian_tube']['proper_tau'][1]);label_lo=ell/16;label_hi=3*ell/16
 clock_H_lower=rho.lo*(H-M*M*H**3/6)
 need(clock_H_lower>tau_max>tau.hi,'compact proper-time tube reaches beyond target with strict clock headroom')
 radial_upper=label_hi+M*H*H/2
 rbox=list(map(F,freeze['enclosing_chart_tube']['scaled_z']));xbox=list(map(F,freeze['enclosing_chart_tube']['scaled_t']))
 if injection=='corrupt_domain_margin':rbox[1]=F(1,8)
 need(xbox==[0,H/ell] and 0<rbox[0]<label_lo/ell and radial_upper/ell<rbox[1],'declared chart tube must contain all connecting trajectories')
 need(xbox[1]<3 and xbox[1]+rbox[1]<8,'actual first-field tube strict domain margin')
 # rho sech chi is positive throughout, hence coordinate t is a valid
 # parameter. Determinant at fixed proper time is t_tau*j, not just j.
 detlo=jlo/rho.hi;radial_gram_lo=rho.lo*rho.lo*jlo*jlo
 angular_gram_lo=rho.lo*rho.lo
 tube_chi=I(-M*H,0);sechlo=1-M*M*H*H/2
 need(sechlo>0 and detlo>0 and radial_gram_lo>0 and angular_gram_lo>0,'strict timelike/inverse-map/screen margins')
 # Same background event Pdagger: no event shift is inserted into the TT API.
 target_bg=core.background(I(F(1,8)),data,'none');target_rho=target_bg[0].val()
 en=core.expbox(I(T));ep=1/en;cosh=(ep+en)/2
 epchi=core.expbox(-chi);ce=1/epchi
 coshchi=(epchi+ce)/2;coshchi=I(max(F(1),coshchi.lo),coshchi.hi);sinhchi=(epchi-ce)/2
 A_target=target_rho*I(1-K*T*T,1)*coshchi
 A_tube=I(rho.lo*jlo,rho.hi/(1-M*M*H*H/2))
 # Clock label derivative at fixed t follows from Gauss lemma:
 # tau_zi = rho*j*sinh(chi), hence t_zi|tau=-j*sinh(chi)*cosh(chi).
 clock_label_target=target_rho*I(1-K*T*T,1)*sinhchi
 # expbox(v) means exp(-v), per authenticated interval_core implementation.
 need(sinhchi.hi<0 and coshchi.lo>0,'hyperbolic sign and orientation')
 return {'schema':'proper-clock-background-certificate-v1','status':'PASS_ACTUAL_BACKGROUND_PROPER_CLOCK_SHOOT_AND_GAUSSIAN_TUBE',
  'source_pins_sha256':digest(HERE/'SOURCE_PINS.json')['sha256'],'stage_freeze_sha256':digest(HERE/'STAGE_FREEZE.json')['sha256'],
  'ell':str(ell),'target':{'t':str(T),'z':str(Z),'scaled_x':'5/4','scaled_r':'1/8','same_background_event_as_retained_Pdagger':True},
  'arithmetic':{'decisions':'exact rational only','outward_grid_denominator':str(core.GRID),'ODE_method':'analytic integral inequalities, no numerical integration','root_freeze':'exact implicit unique actual-background root; interval endpoints and floats are not launch data','clock_freeze':'exact integral along that implicit root; interval endpoints and floats are not clock data'},
  'actual_background_local':{'scaled_z':row(0,F(1,4)),'rho':record(rho),'p':record(p),'a':record(a),'rho_a':record(b),'rho_a_z':record(bz)},
  'shoot':{'t_dagger':str(T),'z_dagger':str(Z),'z_i':record(shoot_bracket),'scaled_z_i':record(shoot_bracket/ell),'radial_shift':row(shift_lo,shift_hi),'endpoint_residual_bounds':endpoint_residuals,'chi_dagger':record(chi),'tau_dagger':record(tau),'shooting_derivative_at_fixed_t':row(1-K*T*T,1),'uniqueness':'all positive launches on the retained inward regular branch that reach the target time; launches above z_dagger cannot hit because a_z>0 and a>0 imply z is increasing'},
  'boost':{'chi_uniform_all_shooting_roots':True,'chi_dagger':record(chi),'cosh_chi_dagger':record(coshchi),'sinh_chi_dagger':record(sinhchi),'exp_chi_dagger':record(epchi),'exp_minus_chi_dagger':record(ce),'rho_at_Pdagger':record(target_rho),'cosh_t_dagger':record(cosh)},
  'compact_tube':{'proper_tau':row(0,tau_max),'scaled_initial_z':row(F(1,16),F(3,16)),'whole_sphere':True,'enclosing_scaled_x':list(map(str,xbox)),'enclosing_scaled_r':list(map(str,rbox)),'trajectory_scaled_r':row(F(1,16),radial_upper/ell),'rapidity':record(tube_chi),'clock_at_horizon_lower':str(clock_H_lower),'proper_clock_headroom':str(tau_max-tau.hi),'clock_horizon_margin':str(clock_H_lower-tau_max),'initial_radial_margin':str(label_lo),'chart_radial_margin':str(ell*rbox[1]-radial_upper),'x_lt_3_margin':str(3-xbox[1]),'x_plus_r_lt_8_margin':str(8-xbox[1]-rbox[1])},
  'gaussian':{'shooting_jacobian_j':row(jlo,1),'chi_label_derivative_k':row(0,2*K*H),'coordinate_determinant_dt_dz_over_dtau_dzi_lower':str(detlo),'radial_gram_lower':str(radial_gram_lo),'angular_gram_relative_to_kappa_lower':str(angular_gram_lo),'clock_derivative_lower':str(rho.lo*sechlo),'radial_jacobi_norm_A_at_target':record(A_target),'radial_jacobi_norm_A_on_tube':record(A_tube),'clock_label_derivative_at_fixed_t_target':record(clock_label_target),'gaussian_cross_terms':'zero by Gauss lemma','identities':['det d(t,z)/d(tau,zi)=j cosh(chi)/rho','partial_zi z at fixed tau=j cosh(chi)^2','gamma_zi_zi=rho^2 j^2 cosh(chi)^2','gamma_AB=(rho cosh(t))^2 kappa_AB'],'injectivity':'t_tau>0 and j>0 imply global injectivity on the stated radial-proper-time rectangle; angular labels are unchanged'},
  'family_persistence':{'statement':'For the supplied exact C2 parameter family with piecewise spacetime C4 regularity, continuous ODE and label-Jacobi dependence on this compact tube preserves positive clock, radial and angular Gram, coordinate determinant, and strict domain margins on a common sufficiently small parameter neighborhood. A compact finite sphere atlas and inverse-map/C1 closeness preserve injectivity. The neighborhood may be ineffective.','effective_parameter_radius':False,'premises':['same initial physical hypersurface and labels','future unit-normal launches','exact-family existence and regularity premise','whole family uniformly defined on an open neighborhood of the compact tube']},
  'scope':{'background_uncertainty_retained':True,'coherent_actual_background_root':True,'ordinary_float_integration_used_for_certification':False,'parameter_variation_ODE_certified_here':False,'second_response_or_final_sign_certified_here':False},
  'input_sha256':{k:v['sha256'] for k,v in pins['files'].items()}}

def canonical_bytes(report):return (json.dumps(report,sort_keys=True,indent=2)+'\n').encode()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path);p.add_argument('--inject-error',choices=['none','wrong_acceleration_sign','wrong_clock_unit','corrupt_domain_margin'],default='none');args=p.parse_args()
 report=certify(args.repo_root,args.inject_error);raw=canonical_bytes(report)
 if args.output is not None:
  need(not args.output.exists(),'refuse overwrite');args.output.write_bytes(raw)
 print(json.dumps({'status':report['status'],'report_sha256':hashlib.sha256(raw).hexdigest(),'z_i':report['shoot']['z_i'],'tau_dagger':report['shoot']['tau_dagger'],'chi_dagger':report['shoot']['chi_dagger']},sort_keys=True))

if __name__=='__main__':main()
