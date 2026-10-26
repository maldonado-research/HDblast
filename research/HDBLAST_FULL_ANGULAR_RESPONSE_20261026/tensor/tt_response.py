#!/usr/bin/env python3
"""New bounded actual total TT response: full source plus complete scalar tail."""
import argparse, hashlib, importlib.util, json, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
LABELS={'current':'HDBLAST_COUPLED_MEAN_PILOT_20261025','I':'HDBLAST_INPUT_PRECISION_20261023','P':'HDBLAST_COMMON_GAUGE_FORCING_20261016','N':'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016','G':'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016','C':'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023','T':'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023'}
PINS={'current':'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de','I':'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61','P':'8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6','N':'928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44','G':'3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2','C':'72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76','T':'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'}
HELPER_PIN='c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619'
HELPER_REL=Path('research/HDBLAST_SHELL_FEEDBACK_PRECISION_20261025/balanced/exact_boundary.py')
X0,R0=F(1,4),F(5,4)
DECLARATION={'first_n':4,'first_cells':20,'max_n':8,'max_cells':72,'refinement_gate':'N8 only if complete N4 total interval contains zero','no_further_refinement':True,'point_probes':3,'whole_cell_probes':2,'acceptance':'exact rational total TT coefficient interval excludes zero'}

def need(ok,message):
 if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def encode(value):
 if hasattr(value,'lo') and hasattr(value,'hi'):return value.row()
 if isinstance(value,F):return str(value)
 if isinstance(value,dict):return {key:encode(item) for key,item in value.items()}
 if isinstance(value,(list,tuple)):return [encode(item) for item in value]
 return value
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def authenticate(roots):
 # Copied complete-tree validation precedes loading even the frozen helper.
 result={}
 for key,root in roots.items():
  need(not root.is_symlink(),'ROOT_SYMLINK:'+key);manifest=root/'MANIFEST.json'
  need(sha(manifest)==PINS[key],'MANIFEST:'+key);seen={'MANIFEST.json'}
  for entry in json.loads(manifest.read_text())['files']:
   rel=Path(entry['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'MEMBER_PATH:'+key);seen.add(str(rel));path=root/rel
   need(path.is_file() and not path.is_symlink(),'MEMBER_TYPE:'+str(path))
   need(path.stat().st_size==entry['bytes'] and sha(path)==entry['sha256'],'MEMBER_BYTES:'+str(path))
  for path in root.rglob('*'):
   need(not path.is_symlink(),'TREE_SYMLINK:'+str(path))
   if path.is_file():need(str(path.relative_to(root)) in seen,'UNLISTED_FILE:'+str(path))
  result[key]={'manifest_sha256':PINS[key],'files_including_manifest':len(seen)}
 return result

def guard_output(repo_root,output):
 repo=Path(repo_root).resolve();target=Path(output).resolve()
 protected=[(repo/'research'/label).resolve() for label in LABELS.values()]+[(repo/HELPER_REL.parents[1]).resolve()]
 need(not output.exists() and not output.is_symlink(),'FRESH_OUTPUT_REQUIRED')
 need(all(target!=root and root not in target.parents for root in protected),'OUTPUT_INSIDE_FROZEN_CAPSULE')

class Context:
 def __init__(self,repo_root):
  self.repo=Path(repo_root).resolve()
  self.roots={key:self.repo/'research'/label for key,label in LABELS.items()};self.before=authenticate(self.roots)
  self.helper_path=self.repo/HELPER_REL
  need(sha(self.helper_path)==HELPER_PIN,'LATEST_HELPER_PIN');self.helper=load('tt_frozen_boundary_authenticator',self.helper_path)
  need(self.helper.PINS==PINS and self.helper.LABELS==LABELS,'HELPER_TABLE');need(self.helper.authenticate(self.roots)==self.before,'HELPER_COMPLETE_AUTHENTICATION')
  sys.path.insert(0,str(self.roots['G']/'numeric'));sys.path.insert(0,str(self.roots['G']/'projection'))
  import actual_first_jets as jets
  import fast_dyadic
  import source_coefficients
  self.jets,self.dyadic,self.coeff=jets,fast_dyadic,source_coefficients;self.data,_,_=jets.inputs()
  self.maps_path=self.roots['G']/'projection/GEOMETRIC_COEFFICIENT_MAPS.json';self.maps=source_coefficients.load_maps(self.maps_path)
  need(self.maps['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','FULL_SOURCE_MAP_PREMISE')
  self.terms=self.maps['maps']['metric_forcing']['ab'];need(len(self.terms)==7 and self.terms[5]==self.terms[6],'FULL_TT_SOURCES_IDENTICAL')
  self.I,self.ell=jets.I,jets.ELL;need(self.ell==F(1,10**6) and jets.core.GRID==1<<112,'ARITHMETIC_CONVENTION');self.calls=0
 def forcing(self,x,r,audit=False):
  m,s,b=self.jets.actual_jet(x,r,mode='actual');raw={**m,**s,**b}
  converted={key:self.dyadic.D(value) if hasattr(value,'lo') and hasattr(value,'hi') else value for key,value in raw.items()}
  # All terms of the full authenticated metric_forcing TT row are evaluated.
  v=self.coeff.evaluate_encoded(self.terms[5],converted);self.calls+=1
  if audit:
   reference=self.coeff.evaluate_encoded(self.terms[5],raw);need(v.row()==reference.row(),'DYADIC_REFERENCE_ENDPOINT_IDENTITY')
   full=self.coeff.evaluate_sources(m,s,b,names=('metric_forcing',),path=self.maps_path)
   need(v.row()==full['metric_forcing']['ab'][5].row()==full['metric_forcing']['ab'][6].row(),'SELECTED_ROW_EQUALS_FULL_SOURCE')
  return v
 def operator(self):
  x,r=self.I(0,X0),self.I(R0-X0,R0+X0)
  rho,eta,p,a,_,_,_,_,_,_=self.jets.core.background(r,self.data,'none');rho,eta,p,a=rho.val(),eta.val(),p.val(),a.val()
  c,Ht,_=self.jets.temporal(x);c,Ht=c.val(),Ht.val();ay=-1/rho**2-p**2/3;b=-(ay+a**2);d=1/rho**2-a**2
  V=self.jets.core.poly(self.jets.core.U,eta);h=Ht/rho
  A=-3*self.ell*Ht;B=-3*self.ell*rho*a;D=self.ell**2*(-6/c**2+rho**2*(-2*a**2+2*h**2-2*b-8*d+F(4,3)*V))
  caps={'A':A.mag(),'B':B.mag(),'D':D.mag()};sigma=caps['A']+caps['B']+(X0/2)*caps['D'];q=X0*sigma;need(0<=q<1,'SCALAR_C1_RESOLVENT_GATE')
  return {'box':{'x':x,'r':r},'background':{'rho':rho,'eta':eta,'p':p,'a':a,'ay':ay,'b':b,'d':d,'V':V,'c':c,'Ht':Ht,'h':h},'coefficient_intervals':{'A':A,'B':B,'D':D},'absolute_caps':caps,'sigma':sigma,'q':q,'norm':'max(2 sup|u|/T,sup|u_x|,sup|u_r|), T=1/4','proof':'||Kf||X<=T||f||inf. |A ux+B ur+D u|<=sigma||u||X with sigma=a+b+(T/2)d. Complete scalar feedback summed by 1/(1-q).'}
 def finish(self):
  need(authenticate(self.roots)==self.before,'FROZEN_CLOSURES_CHANGED');need(self.helper.authenticate(self.roots)==self.before,'HELPER_CLOSURES_CHANGED');need(sha(self.helper_path)==HELPER_PIN,'HELPER_CHANGED')

def probes(ctx):
 rows=[]
 for label,x,r,w in [('initial_packet_point',F(0),F(5,4),None),('middle_cone_point',F(1,8),F(5,4),None),('apex_point',X0,R0,None),('initial_left_N4_cell',F(0),F(1),F(1,16)),('middle_N4_cell',F(1,8),F(19,16),F(1,16))]:
  xx=ctx.I(x) if w is None else ctx.I(x,x+w);rr=ctx.I(r) if w is None else ctx.I(r,r+w);v=ctx.forcing(xx,rr,audit=(not rows))
  row={'label':label,'x':xx,'r':rr,'forcing_Tplus':v,'Tminus_identical':True,'strict_positive':v.lo>0,'strict_negative':v.hi<0,'display_only':[float(v.lo),float(v.hi)]};rows.append(row)
  print(json.dumps({'probe':label,'display_only':row['display_only']}),flush=True)
 return rows

def partition(ctx,n,op):
 need(n in (4,8),'DECLARED_BUDGET_ONLY');rows=[];area_sum=lo=hi=Fcap=F(0);h=X0/n
 for i in range(n):
  for k in range(i,2*n-i):
   xl,xh=i*h,(i+1)*h;rl,rh=R0-X0+k*h,R0-X0+(k+1)*h;area=h*h/2 if k in (i,2*n-i-1) else h*h
   need(xl>=0 and xh<=X0 and rl-xh>=F(3,4),'EXTERIOR_BOUNDING_CELL');v=ctx.forcing(ctx.I(xl,xh),ctx.I(rl,rh));lo+=area*v.lo/2;hi+=area*v.hi/2;Fcap=max(Fcap,v.mag());area_sum+=area
   rows.append({'i':i,'k':k,'x':[xl,xh],'r':[rl,rh],'intersection':'|r-5/4|<=1/4-x','exact_area':area,'forcing_Tplus':v,'Tminus_identical':True})
  print(json.dumps({'n':n,'row':i+1,'rows':n,'cells':len(rows)}),flush=True)
 need(len(rows)==n*(n+1) and area_sum==X0**2,'EXACT_CONE_COVER')
 sigma,q=op['sigma'],op['q'];norm=X0*Fcap/(1-q);error=X0**2*sigma*norm/2;total=[lo-error,hi+error]
 result={'partition':{'n':n,'h':h,'cells':len(rows),'area':area_sum,'proof':'For i=0..N-1 take k=i..2N-i-1; outer triangles area h²/2 and inner cells area h². Exact cone intersections have disjoint interiors and cover the whole cone.'},'actual_source_cell_enclosures':rows,'TT_forcing_absolute_cap':Fcap,'zero_data_Duhamel_interval':[lo,hi],'full_scalar_C1_norm_upper':norm,'all_feedback_point_value_error_upper':error,'total_Tplus_response_interval':total,'Tminus_total_identical':True,'strict_positive':total[0]>0,'strict_negative':total[1]<0,'passes_sign_gate':total[0]>0 or total[1]<0,'display_only':{'duhamel_interval':[float(lo),float(hi)],'feedback_error_upper':float(error),'total_response_interval':list(map(float,total)),'forcing_absolute_cap':float(Fcap)}}
 print(json.dumps({'n':n,**result['display_only'],'passes_sign_gate':result['passes_sign_gate']}),flush=True);return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--probes-only',action='store_true');p.add_argument('--output',type=Path,required=True);args=p.parse_args();guard_output(args.repo_root,args.output);own_pin=sha(Path(__file__));started=time.monotonic();ctx=Context(args.repo_root);sample_rows=probes(ctx);op=ctx.operator();runs=[]
 if not args.probes_only:
  runs.append(partition(ctx,4,op))
  if not runs[-1]['passes_sign_gate']:runs.append(partition(ctx,8,op))
 ctx.finish();need(sha(Path(__file__))==own_pin,'EVALUATOR_CHANGED')
 result={'schema':'HDBLAST_BOUNDED_ACTUAL_TOTAL_TT_RESPONSE_V1','status':'PROBED_ONLY' if not runs else ('CERTIFIED_NONZERO_TOTAL_TT_COEFFICIENT' if runs[-1]['passes_sign_gate'] else 'TOTAL_TT_SIGN_GATE_FAILED_AT_DECLARED_MAXIMUM'),'evaluator_sha256':own_pin,'authentication':ctx.before,'latest_authentication_helper':{'path':ctx.helper_path.relative_to(ctx.repo).as_posix(),'sha256':HELPER_PIN},'all_seven_complete_closures_authenticated_before_project_imports_and_after_run':True,'budget_declaration':DECLARATION,'source':{'path':ctx.maps_path.relative_to(ctx.repo).as_posix(),'sha256':sha(ctx.maps_path),'row':'metric_forcing/ab/5 = metric_forcing/ab/6','full_geometric_terms_per_TT_coefficient':len(ctx.terms[5]),'actual_oracle':'actual_first_jets.actual_jet(mode=actual), physical jets through2 with authenticated C3 derivative errors retained','point_reference_audit':'selected dyadic endpoints equal Fraction and full metric_forcing endpoints'},'structural_premises':{'TT_A':'-3 ell tanh(ell x)','TT_B':'-3 ell rho a','TT_D':'ell²[-6/c²+rho²(-2a²+2h²-2b-8d+4V/3)]','TT_M':'0','TT_initial_value':'exactly0','TT_initial_velocity':'exactly0','TT_rows_and_columns_decoupled':True,'both_TT_equations_and_sources_identical':True,'scope':'Exact full23 operator, initial and source identities separately certified; no TT initial uncertainty.'},'domain':{'apex':[X0,R0],'cone':'0<=x<=1/4, 1+x<=r<=3/2-x','shell_reached':False,'reflection_reached':False,'shell_lift':'exactly0 by backward domain of dependence','free_TT':'exactly0'},'probes':sample_rows,'independent_scalar_operator_bounds':op,'partition_runs':runs,'actual_TT_oracle_calls':ctx.calls,'limits':['Conditional registered-gauge angular metric coefficient of full total second response.','All actual geometric forcing terms and complete variable-coefficient scalar feedback retained.','No source derivatives or comparator substitution.','No physical gravitational-wave, gauge-invariant force, nonlinear family, cosmology, origin or novelty claim.','No old full-campaign rerun, frozen edit or GitHub mutation; authenticated source/error helper functions were used.']}
 with args.output.open('x') as handle:json.dump(encode(result),handle,indent=2,sort_keys=True);handle.write('\n')
 print(json.dumps({'status':result['status'],'calls':ctx.calls,'wall_seconds':round(time.monotonic()-started,3)}),flush=True)
if __name__=='__main__':main()
