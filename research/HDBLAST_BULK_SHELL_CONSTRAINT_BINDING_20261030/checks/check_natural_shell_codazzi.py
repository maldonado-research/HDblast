#!/usr/bin/env python3
"""Exact natural QE bulk/shell Codazzi source-map binding; standard library.

The implementation side evaluates authenticated exported natural geometric
maps. The independent expected side differentiates the physical shell QR
profiles and contracts the intrinsic Lorentzian shell divergence. No first
bulk field equation or harmonic normal row is imposed.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode=True
ENGINE_PIN='178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2'
PINS={
 'map':('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json','34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'),
 'evaluator':('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/source_coefficients.py','bcdfb1de30355699241fd68b0499831d9a805a2405148d1b13a572ef6a6ecf3d'),
 'shell_adapter':('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py','0233ff25b9873dfc4060ba0d4963f50d0b198dbe71368625017840d225b19f8a'),
 'geometric_derivation':('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/geometry/GEOMETRIC_SOURCE_DERIVATION.md','599c59092ac15c0a375b758e6462fe73aa7969e11965fcef4de6c1b74c60bfe0'),
 'physical_shell_proof':('research/HDBLAST_COMMON_GAUGE_FORCING_20261016/boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md','581bcadecd44c408484d20632d7d3e524fb0bd098370d1fd24da2744bd6c1cf1'),
}
CONTROLS={
 'reduced_instead_of_natural':'OFFSHELL_CODAZZI_y0:0',
 'forcing_instead_of_natural':'OFFSHELL_CODAZZI_y0:0',
 'half_second_variation':'OFFSHELL_CODAZZI_y0:0',
 'omit_offshell_scalar_term':'OFFSHELL_CODAZZI_y0:0',
 'positive_time_trace':'OFFSHELL_CODAZZI_y0:0',
 'omit_temporal_FRW_connection':'OFFSHELL_CODAZZI_y0:0',
 'omit_angular_FRW_connection':'OFFSHELL_CODAZZI_ya:0',
 'coordinate_derivative_as_frame':'OFFSHELL_CODAZZI_y0:0',
 'wrong_scalar_shell_sign':'PHYSICAL_CODAZZI_y0:0',
 'wrong_background_tension_sign':'OFFSHELL_CODAZZI_y0:0',
 'missing_rho_time_scale':'OFFSHELL_CODAZZI_y0:0',
}

def need(ok,why):
 if not ok:raise ValueError(why)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def auth(root):
 out={}
 for key,(rel,pin) in PINS.items():
  path=root/rel;need(path.is_file() and not path.is_symlink(),'SOURCE_TYPE:'+key)
  need(sha(path)==pin,'SOURCE_PIN:'+key)
  out[key]={'path':rel,'bytes':path.stat().st_size,'sha256':pin}
 engine=Path(__file__).resolve().with_name('codazzi_laurent.py')
 need(engine.is_file() and not engine.is_symlink() and sha(engine)==ENGINE_PIN,'ENGINE_PIN')
 return out

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def calculate(root,mutation=None):
 before=auth(root)
 alg=load('codazzi_exact_ring',Path(__file__).resolve().with_name('codazzi_laurent.py'))
 P,dt=alg.P,alg.dt
 adapter=load('codazzi_actual_shell_adapter',root/PINS['shell_adapter'][0])
 evaluator=load('codazzi_registered_map_evaluator',root/PINS['evaluator'][0])
 report=json.loads((root/PINS['map'][0]).read_text())
 need(report['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','EXPORTED_MAP_STATUS')
 checks=[]
 def eq(name,left,right=0):
  left,right=P.cv(left),P.cv(right);res=left-right
  need(not res,'NONZERO:'+name+':'+json.dumps(res.encode(),separators=(',',':')))
  checks.append({'name':name,'left_collected_terms':len(left.t),'right_collected_terms':len(right.t),'residual_collected_terms':0})
 vals,g,beta=alg.shell_values(report,False)
 h,X,V,L,d,v,rho,c,theta,a,p,sigma2=[g[n] for n in g]
 sig1=-2*p
 if mutation=='wrong_background_tension_sign':
  sig1=2*p;beta=a*h+sig1*v/3
  vals.update(Dy=-beta,Ly=beta)
  vals['Dyt']=dt(vals['Dy']);vals['Lyt']=dt(vals['Ly'])
 if mutation=='coordinate_derivative_as_frame':
  vals['Dy']+=2*a*vals['D'];vals['Ly']+=2*a*vals['L']
  vals['Dyt']=dt(vals['Dy']);vals['Lyt']=dt(vals['Ly'])
 # Independent physical QR formula from the affine unit-normal Hessian after
 # FIRST physical shell rows. Neither the geometric bulk maps nor the shell
 # adapter defines this expected table.
 A0=a*h*h/4-sig1*h*v/6-sigma2*v*v/6
 q00=-A0+2*beta*X;q0=-beta*V;qa=A0-2*beta*L;qb=-2*beta*d
 T=3*qa-8*qb
 ys=(F(1,6),F(1,2),F(1));gs=(F(4,3),F(2),F(-4))
 vs=(F(1,4),F(1,2));hs=(F(1),F(-2))
 qrad=adapter.shell_sources(a=a,p=p,sigmap=sig1,sigmapp=sigma2,sigmappp=P.var('sigma3'),h=h,X=X,V=V,L=L,d=d,v=v)
 YHESS=(F(-4,9),F(-4,3),F(-8,3),F(-1,4),F(1,6),F(2,3),F(2,3))
 for j,x in enumerate(ys):eq('SAME_ACTUAL_QR_00:'+str(j),qrad['QK00'][j],q00*x)
 for j,x in enumerate(vs):eq('SAME_ACTUAL_QR_0a:'+str(j),qrad['QK0a'][j],q0*x)
 for j in range(7):eq('SAME_ACTUAL_QR_ab:'+str(j),qrad['QKab'][j],qa*(ys[j] if j<3 else 0)+qb*YHESS[j])
 for j,x in enumerate(ys):eq('LORENTZIAN_QR_TRACE:'+str(j),-qrad['QK00'][j]+3*qrad['QKab'][j],(-q00+T)*x)
 # Contracted intrinsic divergence in the physical frame. The coordinate
 # connection derivation is independently checked in the companion review.
 timefactor=1/rho
 if mutation=='missing_rho_time_scale':timefactor=P(1)
 temporal=(-timefactor*dt(T)-theta/rho*(3*q00+T)-8*q0/(rho*c),q0/(rho*c))
 angular=((2*q00-4*qa+10*qb)/(rho*c)-timefactor*dt(q0)-4*theta*q0/rho,qb/(rho*c))
 if mutation=='positive_time_trace':
  # Change only tr(QR)=-q00+T into +q00+T in div(QR)-d tr(QR).
  # This adds -2 d(q00 Y²), distinct from reversing T_t by hand.
  temporal=(temporal[0]-2*dt(q00)/rho,temporal[1])
  angular=(angular[0]-4*q00/(rho*c),angular[1])
 if mutation=='omit_temporal_FRW_connection':temporal=(temporal[0]+theta/rho*(3*q00+T),temporal[1])
 if mutation=='omit_angular_FRW_connection':angular=(angular[0]+4*theta*q0/rho,angular[1])
 expected={'y0':tuple(temporal[0]*x+temporal[1]*y for x,y in zip(ys,gs)),
           'ya':tuple(angular[0]*x+angular[1]*y for x,y in zip(vs,hs))}
 quantities='QE'
 if mutation=='reduced_instead_of_natural':quantities='QEred'
 if mutation=='forcing_instead_of_natural':quantities='metric_forcing'
 def maps(values):
  out={field:tuple(alg.encoded(terms,values) for terms in report['maps'][quantities][field]) for field in ('y0','ya')}
  if mutation=='half_second_variation':out={k:tuple(x/2 for x in v) for k,v in out.items()}
  return out
 # All ordinary normal amplitude derivatives not fixed by physical boundary
 # rows remain free. In particular Ay,By,Cy and yy derivatives are untouched.
 off=maps(vals)
 P1=vals['sy']-p*h/2+sigma2*v/2
 offmatter={'y0':tuple(-2*P1*P.var('vt')/rho*x for x in ys),
            'ya':tuple(-2*P1*v/(rho*c)*x for x in vs)}
 if mutation=='omit_offshell_scalar_term':offmatter={k:tuple(P() for x in v) for k,v in offmatter.items()}
 for field in ('y0','ya'):
  for j,(actual,want,matter) in enumerate(zip(off[field],expected[field],offmatter[field])):
   eq('OFFSHELL_CODAZZI_'+field+':'+str(j),actual,want+matter)
 # Same shell jet class, now impose actual scalar physical row P1=0.
 physical=dict(vals);physical['sy']=p*h/2-sigma2*v/2
 if mutation=='wrong_scalar_shell_sign':physical['sy']=p*h/2+sigma2*v/2
 closed=maps(physical)
 for field in ('y0','ya'):
  for j,(actual,want) in enumerate(zip(closed[field],expected[field])):
   eq('PHYSICAL_CODAZZI_'+field+':'+str(j),actual,want)
 # Independent monomial decoder is checked against the actual implementation
 # evaluator on the SAME exported map and coherent physical jet table.
 mkeys={field+suffix for field in ('A','B','C','D','E','L','M') for suffix in ('','t','y','tt','yt','yy')}
 skeys={'s','st','sy','stt','syt','syy'}
 def registered(values):
  m={k:values[k] for k in mkeys};ss={k:values[k] for k in skeys};b={k:v for k,v in values.items() if k not in mkeys|skeys}
  return evaluator.evaluate_sources(m,ss,b,('QE',),path=root/PINS['map'][0])['QE']
 for tag,values,independent in [('OFFSHELL',vals,off),('PHYSICAL',physical,closed)]:
  got=registered(values)
  for field in ('y0','ya'):
   for j,(x,y) in enumerate(zip(got[field],independent[field])):eq('REGISTERED_MAP_EVALUATOR_'+tag+'_'+field+':'+str(j),x,y)
 free_normal=sorted(n for n in vals if ('y' in n and n not in ('Dy','Ly','Ey','My','Dyt','Lyt','Eyt','Myt','sy','syt')))
 for field in ('y0','ya'):
  for j,expr in enumerate(closed[field]):
   need(not(set(expr.names())&set(free_normal)),'UNEXPECTED_NORMAL_JET:'+field+':'+str(j))
 # Check the external inputs and arithmetic source again after all use.
 need(auth(root)==before,'SOURCE_READBACK_CHANGED')
 return {'status':'PASS_EXACT_NATURAL_BULK_SHELL_CODAZZI_BINDING',
  'source_sha256':before,'arithmetic_helper_sha256':ENGINE_PIN,
  'source_quantity':'natural QE second parameter derivatives; not QEred or principal forcing',
  'independent_expected_QR_primitives':{k:x.encode() for k,x in {'q00':q00,'q0':q0,'qa':qa,'qb':qb}.items()},
  'independent_intrinsic_divergence_primitives':{'temporal':[x.encode() for x in temporal],'angular':[x.encode() for x in angular]},
  'actual_natural_QE_physical_modal_coefficients':{k:[x.encode() for x in v] for k,v in closed.items()},
  'actual_natural_QE_offshell_modal_coefficients':{k:[x.encode() for x in v] for k,v in off.items()},
  'first_scalar_physical_residual':P1.encode(),
  'offscalar_scope':'First physical metric rows R1=0 and background R0=P0=0 are retained; only the first scalar residual P1 is left free. This is not an arbitrary R1-free off-shell source identity.',
  'modal_rows':{'natural_Codazzi_with_P1_free':5,'natural_Codazzi_with_P1_zero':5},
  'named_identity_rows':len(checks),'checks':checks,
  'normal_jets_retained_formal':free_normal,
  'unneeded_background_quantities':['ay','py','Vp','Vpp','Vppp'],
  'premises':['fixed physical shell with B=C=0 and all tangential derivatives of these traces zero',
   'background R0=0 with sigma0=6a and P0=0 with sigma1=-2p',
   'first physical metric rows Dy=-beta,Ly=beta,Ey=My=0 and their tangential derivatives only',
   'ordinary physical y,t derivatives, radially parallel frame, c_t=c theta',
   'nonzero rho,c and static shell background parameters',
   'for final5rows only: sy=p h/2-sigma2 v/2'],
  'not_imposed':['harmonic Ay row','any first bulk PDE or reduced acceleration','normal derivative of a shell-only identity','background bulk ODE','cosh identity 1/c^2+theta^2=1'],
  'scope':'Exact source-map binding to existing Codazzi and physical-shell derivations. No new evolution/uniqueness theorem, physical accuracy, finite-amplitude nonlinear family or cosmological-origin result.'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--mutation',choices=CONTROLS);ap.add_argument('--run-controls',action='store_true');a=ap.parse_args()
 need(not a.output.exists(),'REFUSE_OVERWRITE')
 result=calculate(a.repo_root.resolve(),a.mutation)
 if a.run_controls:
  runs=[]
  with tempfile.TemporaryDirectory() as temp:
   for mutation,reason in CONTROLS.items():
    path=Path(temp)/(mutation+'.json')
    cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(Path(__file__).resolve()),'--repo-root',str(a.repo_root.resolve()),'--output',str(path),'--mutation',mutation]
    run=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
    need(run.returncode!=0 and 'NONZERO:'+reason+':' in run.stderr and not path.exists(),'CONTROL_NOT_REJECTED:'+mutation+':'+run.stderr[-1500:])
    runs.append({'mutation':mutation,'expected_nonzero_identity':reason,'returncode':run.returncode,'successful_report_absent':True})
  result['rejected_corrupted_executions']=runs
 result['checker_sha256']=sha(Path(__file__))
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'status':result['status'],'identities':result['named_identity_rows'],'controls':len(result.get('rejected_corrupted_executions',[]))}))
if __name__=='__main__':main()
