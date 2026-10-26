#!/usr/bin/env python3
"""Independent TT report review: authentication, geometry, sums, oracle spots.

Does not import or run the TT evaluator.  All decisions use exact Fractions.
"""
import argparse,copy,hashlib,importlib.util,json,sys
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
LABELS={'current':'HDBLAST_COUPLED_MEAN_PILOT_20261025','I':'HDBLAST_INPUT_PRECISION_20261023','P':'HDBLAST_COMMON_GAUGE_FORCING_20261016','N':'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016','G':'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016','C':'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023','T':'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023'}
HELPER_CAPSULE='HDBLAST_SHELL_FEEDBACK_PRECISION_20261025'
PINS={'current':'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de','I':'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61','P':'8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6','N':'928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44','G':'3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2','C':'72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76','T':'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'}
HELPER_PIN='c01a8f5f95737dcacebdc26379ae0ee5863af9c1977890b3901ad6113103e619'
T=Q(1,4);R=Q(5,4);GRID=1<<112

def need(ok,label):
 if not ok:raise ValueError(label)
def guard_output(output,repo_root):
 target=output.resolve()
 protected=[(repo_root/'research'/name).resolve() for name in (*LABELS.values(),HELPER_CAPSULE)]
 need(not output.exists() and not output.is_symlink(),'FRESH_OUTPUT')
 need(all(target!=root and root not in target.parents for root in protected),'OUTPUT_INSIDE_PROTECTED_CAPSULE')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def complete_auth(roots):
 out={}
 for key,root in roots.items():
  need(not root.is_symlink(),'ROOT_SYMLINK')
  manifest=root/'MANIFEST.json';need(sha(manifest)==PINS[key],'MANIFEST:'+key)
  seen={'MANIFEST.json'}
  for entry in json.loads(manifest.read_text())['files']:
   rel=Path(entry['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'MEMBER_PATH');seen.add(str(rel));p=root/rel
   need(p.is_file() and not p.is_symlink() and p.stat().st_size==entry['bytes'] and sha(p)==entry['sha256'],'MEMBER:'+str(p))
  for p in root.rglob('*'):
   need(not p.is_symlink(),'TREE_SYMLINK')
   if p.is_file():need(str(p.relative_to(root)) in seen,'UNLISTED_FILE')
  out[key]={'manifest_sha256':PINS[key],'files_including_manifest':len(seen)}
 return out

# A separate outward rational interval class for operator reconstruction.
class Box:
 def __init__(self,lo,hi=None):
  lo=Q(lo);hi=lo if hi is None else Q(hi);need(lo<=hi,'INTERVAL_ORDER')
  self.lo=Q((lo*GRID).__floor__(),GRID);self.hi=Q((hi*GRID).__ceil__(),GRID)
 @classmethod
 def cv(cls,v):return v if isinstance(v,cls) else cls(v)
 def __add__(self,v):v=self.cv(v);return Box(self.lo+v.lo,self.hi+v.hi)
 __radd__=__add__
 def __neg__(self):return Box(-self.hi,-self.lo)
 def __sub__(self,v):return self+-self.cv(v)
 def __rsub__(self,v):return self.cv(v)+-self
 def __mul__(self,v):
  v=self.cv(v);p=[x*y for x in (self.lo,self.hi) for y in (v.lo,v.hi)];return Box(min(p),max(p))
 __rmul__=__mul__
 def inv(self):need(self.lo>0 or self.hi<0,'ZERO_DENOMINATOR');return Box(1/self.hi,1/self.lo)
 def __truediv__(self,v):return self*self.cv(v).inv()
 def __rtruediv__(self,v):return self.cv(v)*self.inv()
 def __pow__(self,n):
  need(n in (-2,-1,0,1,2),'REVIEW_POWER_RANGE')
  if n<0:return self.inv()**(-n)
  if n==0:return Box(1)
  if n==1:return self
  return Box(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
 def row(self):return [str(self.lo),str(self.hi)]
 def mag(self):return max(abs(self.lo),abs(self.hi))

def clip(poly,axis,edge,greater):
 if not poly:return []
 out=[];prev=poly[-1];pin=prev[axis]>=edge if greater else prev[axis]<=edge
 for point in poly:
  qin=point[axis]>=edge if greater else point[axis]<=edge
  if pin!=qin:
   ratio=(edge-prev[axis])/(point[axis]-prev[axis]);out.append(tuple(prev[j]+ratio*(point[j]-prev[j]) for j in (0,1)))
  if qin:out.append(point)
  prev,pin=point,qin
 return out

def area(box):
 xl,xh,rl,rh=box;poly=[(Q(0),R-T),(Q(0),R+T),(T,R)]
 for axis,edge,greater in ((0,xl,True),(0,xh,False),(1,rl,True),(1,rh,False)):poly=clip(poly,axis,edge,greater)
 return abs(sum((p[0]*q[1]-q[0]*p[1] for p,q in zip(poly,poly[1:]+poly[:1])),Q(0)))/2 if poly else Q(0)

def verify_run(run,op):
 p=run['partition'];n=p['n'];need(n in (4,8),'N_BUDGET');h=T/n
 need(Q(p['h'])==h and p['cells']==n*(n+1),'PARTITION_METADATA')
 rows=run['actual_source_cell_enclosures'];expected=[(i,k) for i in range(n) for k in range(i,2*n-i)]
 need([(v['i'],v['k']) for v in rows]==expected,'EXACT_CELL_ROSTER')
 lo=hi=total=f=Q(0)
 for row in rows:
  i,k=row['i'],row['k'];box=tuple(map(Q,row['x']+row['r']));need(box==(i*h,(i+1)*h,R-T+k*h,R-T+(k+1)*h),'BOX_GEOMETRY')
  ar=area(box);need(ar>0 and ar==Q(row['exact_area']),'EXACT_POLYGON_AREA');total+=ar
  l,u=map(Q,row['forcing_Tplus']);need(l<=u and row['Tminus_identical'] is True,'SOURCE_INTERVAL')
  lo+=ar*l/2;hi+=ar*u/2;f=max(f,abs(l),abs(u))
 need(total==T*T==Q(p['area']),'CONE_TOTAL_AREA')
 need([str(lo),str(hi)]==run['zero_data_Duhamel_interval'],'DUHAMEL_SUM')
 need(f==Q(run['TT_forcing_absolute_cap']),'FORCING_CAP')
 a,b,d=(Q(op['absolute_caps'][key]) for key in ('A','B','D'));sigma=a+b+T*d/2;q=T*sigma
 need(0<=q<1 and sigma==Q(op['sigma']) and q==Q(op['q']),'SCALAR_CONTRACTION')
 norm=T*f/(1-q);e=T*T*f*q/(2*(1-q));interval=[lo-e,hi+e]
 need(norm==Q(run['full_scalar_C1_norm_upper']),'RESOLVENT_NORM');need(e>0 and e==Q(run['all_feedback_point_value_error_upper']),'COMPLETE_FEEDBACK')
 need(list(map(str,interval))==run['total_Tplus_response_interval'],'TOTAL_INTERVAL')
 need(run['strict_positive']==(interval[0]>0) and run['strict_negative']==(interval[1]<0),'SIGN_FLAGS')
 need(run['passes_sign_gate']==(interval[0]>0 or interval[1]<0) and run['Tminus_total_identical'] is True,'SIGN_GATE')
 return {'n':n,'cells':len(rows),'area':str(total),'Duhamel_interval':list(map(str,(lo,hi))),'TT_forcing_absolute_cap':str(f),'sigma':str(sigma),'q':str(q),'complete_feedback_error':str(e),'total_interval':list(map(str,interval)),'display_only':{'total_interval':list(map(float,interval)),'feedback_error':float(e)}}

def rejected(thunk,label):
 try:thunk()
 except ValueError:return label
 raise ValueError('MUTATION_NOT_REJECTED:'+label)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--evaluator',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();repo_root=args.repo_root.resolve();guard_output(args.output,repo_root)
 own_pin=sha(Path(__file__));report_pin=sha(args.report);evaluator_pin=sha(args.evaluator);d=json.loads(args.report.read_text());need(d['evaluator_sha256']==evaluator_pin,'EVALUATOR_BINDS_REPORT')
 need(d['source']['path']=='research/'+LABELS['G']+'/projection/GEOMETRIC_COEFFICIENT_MAPS.json','LOGICAL_SOURCE_PATH')
 need(d['latest_authentication_helper']=={'path':'research/'+HELPER_CAPSULE+'/balanced/exact_boundary.py','sha256':HELPER_PIN},'LOGICAL_HELPER_PATH_AND_PIN')
 need('wall_seconds' not in d,'CANONICAL_REPORT_EXCLUDES_TIMING')
 roots={key:repo_root/'research'/value for key,value in LABELS.items()};before=complete_auth(roots)
 helper_path=repo_root/'research'/HELPER_CAPSULE/'balanced/exact_boundary.py';need(sha(helper_path)==HELPER_PIN,'LATEST_HELPER_PIN')
 helper=load('tt_review_frozen_boundary_auth',helper_path);need(helper.PINS==PINS and helper.LABELS==LABELS,'LATEST_HELPER_TABLES');need(helper.authenticate(roots)==before,'LATEST_HELPER_AUTH');need(d['authentication']==before and d['all_seven_complete_closures_authenticated_before_project_imports_and_after_run'],'REPORTED_AUTH')
 sys.path.insert(0,str(roots['G']/'numeric'));sys.path.insert(0,str(roots['G']/'projection'))
 import actual_first_jets as jets
 import source_coefficients as coefficient_api
 maps_path=roots['G']/'projection/GEOMETRIC_COEFFICIENT_MAPS.json';maps=coefficient_api.load_maps(maps_path);need(sha(maps_path)==d['source']['sha256'],'SOURCE_MAP_REPORT_BINDING')
 terms=maps['maps']['metric_forcing']['ab'];need(len(terms)==7 and terms[5]==terms[6] and len(terms[5])==d['source']['full_geometric_terms_per_TT_coefficient']==78,'FULL_TT_MAP_ROWS')
 data,_,_=jets.inputs();need(jets.ELL==Q(1,10**6) and jets.core.GRID==GRID,'ORACLE_ARITHMETIC')
 # Reconstruct through authenticated G background oracle, then separate Box arithmetic.
 bg=jets.core.background(jets.I(1,Q(3,2)),data,'none');cv=lambda value:Box(value.lo,value.hi)
 rho,eta,p,a=(cv(bg[k].val()) for k in (0,1,2,3));tc,th,_=jets.temporal(jets.I(0,T));c,Ht=cv(tc.val()),cv(th.val());ell=jets.ELL
 ay=-1/rho**2-p**2/3;b=-(ay+a**2);dd=1/rho**2-a**2;V=Box(0)
 for item in reversed(jets.core.U):V=V*eta+item
 hh=Ht/rho
 A=-3*ell*Ht;B=-3*ell*rho*a;D=ell**2*(-6/c**2+rho**2*(-2*a**2+2*hh**2-2*b-8*dd+Q(4,3)*V))
 op=d['independent_scalar_operator_bounds'];background={'rho':rho,'eta':eta,'p':p,'a':a,'ay':ay,'b':b,'d':dd,'V':V,'c':c,'Ht':Ht,'h':hh}
 need(op['box']=={'x':['0','1/4'],'r':['1','3/2']},'OPERATOR_CONE_BOUNDING_BOX')
 for key,val in background.items():need(val.row()==op['background'][key],'BACKGROUND_RECONSTRUCTION:'+key)
 for key,val in {'A':A,'B':B,'D':D}.items():need(val.row()==op['coefficient_intervals'][key] and val.mag()==Q(op['absolute_caps'][key]),'OPERATOR_CAP:'+key)
 run_reviews=[verify_run(run,op) for run in d['partition_runs']];need(len(run_reviews)>=1 and run_reviews[0]['n']==4,'FIRST_N4')
 need(len(run_reviews)==1 if d['partition_runs'][0]['passes_sign_gate'] else len(run_reviews)==2 and run_reviews[1]['n']==8,'DECLARED_REFINEMENT_GATE')
 need(d['status']=='CERTIFIED_NONZERO_TOTAL_TT_COEFFICIENT' and d['partition_runs'][-1]['passes_sign_gate'],'FINAL_STATUS')
 # Full geometric API spot replay, actual solution errors retained.
 point=d['probes'][0];rows=d['partition_runs'][0]['actual_source_cell_enclosures'];clipped=next(row for row in rows if row['i']==0 and row['k']==0);full=next(row for row in rows if row['i']==0 and row['k']==3)
 replay=[];point_raw=None
 for label,row in [('initial_point',point),('clipped_cell',clipped),('full_cell',full)]:
  xx=jets.I(*map(Q,row['x']));rr=jets.I(*map(Q,row['r']));m,s,bg0=jets.actual_jet(xx,rr,mode='actual');cm,_,_=jets.actual_jet(xx,rr,mode='comparator')
  need(any(m[key].hi-m[key].lo>cm[key].hi-cm[key].lo for key in m),'ACTUAL_C3_UNCERTAINTY_RETAINED')
  values=coefficient_api.evaluate_sources(m,s,bg0,names=('metric_forcing',),path=maps_path);plus,minus=values['metric_forcing']['ab'][5:7]
  need(plus.row()==minus.row()==row['forcing_Tplus'],'FULL_SOURCE_API_SPOT:'+label)
  replay.append({'label':label,'x':row['x'],'r':row['r'],'full_metric_forcing_Tplus':plus.row(),'Tminus_identical':True,'actual_C3_uncertainty_retained':True})
  if label=='initial_point':point_raw={**m,**s,**bg0}
 # Negative controls must be rejected by the independent review checks.
 original=d['partition_runs'][0];mutations=[]
 bad=copy.deepcopy(original);bad['actual_source_cell_enclosures'].pop();mutations.append(rejected(lambda:verify_run(bad,op),'drop_cell'))
 bad=copy.deepcopy(original);bad['actual_source_cell_enclosures'][0]['exact_area']=str(2*Q(bad['actual_source_cell_enclosures'][0]['exact_area']));mutations.append(rejected(lambda:verify_run(bad,op),'double_clipped_cell_area'))
 bad=copy.deepcopy(original);bad['all_feedback_point_value_error_upper']='0';bad['total_Tplus_response_interval']=bad['zero_data_Duhamel_interval'];mutations.append(rejected(lambda:verify_run(bad,op),'omit_complete_feedback'))
 dropped=coefficient_api.evaluate_encoded(terms[5][1:],point_raw);need(dropped.row()!=point['forcing_Tplus'],'SOURCE_TERM_DROP_CONTROL');mutations.append('drop_first_full_geometric_TT_source_term')
 need(complete_auth(roots)==before and helper.authenticate(roots)==before,'SEVEN_CLOSURES_UNCHANGED');need(sha(helper_path)==HELPER_PIN and sha(Path(__file__))==own_pin and sha(args.report)==report_pin and sha(args.evaluator)==evaluator_pin,'REVIEW_INPUTS_UNCHANGED')
 out={'schema':'HDBLAST_INDEPENDENT_TOTAL_TT_REPORT_REVIEW_V1','status':'PASS_INDEPENDENT_TOTAL_TT_NUMERIC_REVIEW','report_sha256':report_pin,'evaluator_sha256':evaluator_pin,'checker_sha256':own_pin,'authentication':before,'latest_helper_sha256':HELPER_PIN,'all_seven_closures_authenticated_before_project_imports_and_after_run':True,'operator_background_and_coefficients_reconstructed_with_independent_grid_interval_engine':True,'all_recorded_cell_geometry_source_sums_and_scalar_resolvent_arithmetic_checked':True,'run_reviews':run_reviews,'full_source_api_spot_replays':replay,'negative_controls_rejected':mutations,'limits':['Structural decoupling, exact zero TT initial data, and formal-source/gauge existence remain separately authenticated premises.','Every recorded source interval participates in the sum; only three source intervals were independently replayed.','No physical gravitational-wave, gauge-invariant force, cosmological or nonlinear-family claim.']}
 with args.output.open('x') as handle:json.dump(out,handle,indent=2,sort_keys=True);handle.write('\n')
 print(json.dumps({'status':out['status'],'source_replays':len(replay),'negative_controls':len(mutations),'total_interval_display_only':run_reviews[-1]['display_only']['total_interval']},indent=2))
if __name__=='__main__':main()
