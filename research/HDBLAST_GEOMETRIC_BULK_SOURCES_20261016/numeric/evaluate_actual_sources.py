#!/usr/bin/env python3
"""Continuous actual evolved geometric quadratic source cell evaluations.

The coefficient maps are exact algebraically collected projections.  Actual
first solution jet enclosures, not the free comparisons, are substituted into
those maps.  Every reported interval covers its entire stated cell, including
the whole S^3 represented by the exact angular basis.
"""
import argparse,hashlib,importlib.util,json,math,sys
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT.parent/'projection'))
import actual_first_jets as j
import shell_source_values as shell
import source_coefficients as coefficients
import angular_generic as angular
import fast_dyadic
I=j.I;need=j.need
NAMES=('QE','QEred','Qphi','QC','metric_forcing','scalar_forcing')

def encode(x):
 if hasattr(x,'lo') and hasattr(x,'hi'):return x.row()
 if isinstance(x,dict):return {key:encode(value) for key,value in x.items()}
 if isinstance(x,(tuple,list)):return [encode(value) for value in x]
 return str(x)

def leaves(tree,prefix=()):
 if isinstance(tree,dict):
  for key,value in tree.items():yield from leaves(value,prefix+(str(key),))
 elif isinstance(tree,(tuple,list)):
  for key,value in enumerate(tree):yield from leaves(value,prefix+(str(key),))
 else:yield '/'.join(prefix),tree

def cell(x,r,maps,mutation='none',clipped=False,audit_actual=False):
 actual_m,s,b=j.actual_jet(x,r,mode='comparator' if mutation=='comparator_only' else 'actual',clip_causal=clipped)
 if audit_actual or mutation=='comparator_only':
  cm,_,_=j.actual_jet(x,r,mode='comparator',clip_causal=clipped)
  need(any(actual_m[k].hi-actual_m[k].lo>cm[k].hi-cm[k].lo for k in actual_m),'actual solution derivative errors retained rather than comparator substitution')
 convert=lambda row:{key:fast_dyadic.D(value) if hasattr(value,'lo') and hasattr(value,'hi') else value for key,value in row.items()}
 values=coefficients.evaluate_sources(convert(actual_m),convert(s),convert(b),names=NAMES,path=maps)
 if mutation=='halve_hessian':
  values={name:({field:tuple(v/2 for v in row) for field,row in fields.items()} if isinstance(fields,dict) else tuple(v/2 for v in fields)) for name,fields in values.items()}
 if mutation=='reverse_metric_forcing':values['metric_forcing']={field:tuple(-v for v in row) for field,row in values['metric_forcing'].items()}
 return values

def summarize(rows):
 caps={};strict={}
 for row in rows:
  for path,value in leaves(row['raw_sources']):
   caps[path]=max(caps.get(path,F(0)),value.mag())
   if value.lo>0:strict[path]=strict.get(path,0)+1
   if value.hi<0:strict[path]=strict.get(path,0)+1
 return {'absolute_coefficient_caps':{path:str(value) for path,value in caps.items()},
         'strict_nonzero_cell_counts':strict}

def forcing_H_bound(summary):
 """Probability sphere L2, positive metric/scalar Frobenius bundle norm."""
 caps={key:F(value) for key,value in summary['absolute_coefficient_caps'].items()}
 square=F(0);blocks={}
 for field,basis,multiplicity in (('yy',angular.SCALARS,1),('y0',angular.SCALARS,2),('00',angular.SCALARS,1),('ya',angular.GRADS,2),('0a',angular.GRADS,2),('ab',angular.TENSORS,1),('scalar',angular.SCALARS,1)):
  name='scalar_forcing' if field=='scalar' else 'metric_forcing'
  prefix=name if field=='scalar' else name+'/'+field
  cc=[caps[prefix+'/'+str(i)] for i in range(len(basis))]
  gram=[[angular.inner(a,b) for b in basis] for a in basis]
  bound=multiplicity*sum((abs(gram[i][k])*cc[i]*cc[k] for i in range(len(basis)) for k in range(len(basis))),F(0))
  blocks[field]={'multiplicity':multiplicity,'gram':[[str(v) for v in row] for row in gram],'squared_norm_upper':str(bound)};square+=bound
 upper=math.isqrt(square.numerator//square.denominator)+1
 need(F(upper)**2>=square,'exact positive bundle norm square-root roundup')
 parent=F(171610068644000)
 return {'norm':'Probability-normalized S3 L2 of the positive Frobenius metric-plus-scalar bundle; mixed tensor blocks counted twice.',
         'squared_upper':str(square),'F_H_upper':str(upper),'block_Gram_certificates':blocks,
         'retained_parent_F_H':str(parent),'strictly_improves_parent_bulk_forcing_bound':F(upper)<parent,
         'comparison_scope':'Bulk C0 forcing only. No operator, shell Green-map, response quadrature, or full continuous second-response residual is updated by this number.'}

def run(maps,profile_cells=64,cover_divisions=2,mutation='none'):
 need(profile_cells>=8 and cover_divisions>=1,'nonempty source profile and full causal cover')
 j.inputs();map_report=coefficients.load_maps(maps)
 dyadic_checks=fast_dyadic.verify_equivalence(I)
 need(map_report['status']=='PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION','exact all-point source map premise')
 raw_maps=Path(maps).read_bytes()
 witness_x=I(F(1,4),F(1,4)+F(1,4096));witness_r=I(F(5,4),F(5,4)+F(1,4096))
 witness=cell(witness_x,witness_r,maps,mutation,audit_actual=True)
 need(witness['QE']['yy'][0].hi < -5*10**17,'actual full-cell Einstein Hessian negative witness and second-derivative normalization')
 need(witness['QEred']['yy'][0].hi < -5*10**17,'actual reduced Einstein Hessian negative witness')
 need(witness['QE']['y0'][0].lo > 5*10**17,'actual mixed Einstein Hessian positive witness')
 need(witness['metric_forcing']['yy'][0].lo>16*10**9,'registered reduced metric forcing sign and scale')
 im,is_,ib=j.actual_jet(witness_x,witness_r)
 reference=coefficients.evaluate_sources(im,is_,ib,names=NAMES,path=maps)
 need(encode(reference)==encode(witness),'integer dyadic and rational source-map endpoint identity')
 shell.algebra_check('omit_normal_normalization' if mutation=='omit_shell_normalization' else None)
 reflected_x=I(F(3,2),F(3,2)+F(1,4096));reflected_r=I(F(1,4),F(1,4)+F(1,4096))
 reflected=cell(reflected_x,reflected_r,maps)
 shell_x=I(F(5,4),F(5,4)+F(1,4096));shell_values=shell.actual_shell_sources(shell_x)
 need(shell_values['one_copy_outward_hessian']['scalar'][0].lo>4*10**8,'actual full-cell positive scalar shell Hessian')
 # A finite rectangular partition of [0,3]x[0,8], clipped only by x+r<=8.
 # Every D point belongs to one retained rectangle.  The diagonal parts
 # are explicitly included, rather than silently discarded as small cells.
 d=cover_divisions
 expected=[(i,k) for i in range(3*d) for k in range(8*d) if i+k<8*d]
 indices=expected[:-1] if mutation=='drop_corner_cell' else expected
 need(indices==expected,'complete clipped causal triangle box cover including diagonal cells')
 global_rows=[]
 for index,(i,k) in enumerate(indices):
  x=I(F(i,d),F(i+1,d));r=I(F(k,d),F(k+1,d));values=cell(x,r,maps,clipped=True)
  global_rows.append({'index':[i,k],'x':x.row(),'r':r.row(),'intersection':'x+r<=8','raw_sources':values})
 profile=[]
 for k in range(profile_cells):
  r=I(1+F(k,profile_cells),1+F(k+1,profile_cells));values=cell(witness_x,r,maps)
  profile.append({'index':k,'x':witness_x.row(),'r':r.row(),'raw_sources':values})
 # Shell whole-time interval cover, with all anisotropic tangential rows.
 shell_rows=[]
 for k in range(3*d):
  x=I(F(k,d),F(k+1,d));shell_rows.append({'x':x.row(),'sources':shell.actual_shell_sources(x)})
 summary=summarize(global_rows);summary['positive_bundle_bulk_forcing_H_bound']=forcing_H_bound(summary);profile_summary=summarize(profile)
 for rows in (global_rows,profile):
  for row in rows:row['sources']=encode(row.pop('raw_sources'))
 return {'status':'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES',
         'domain':{'ell':str(j.ELL),'causal_triangle':'0<=x<=3,r>=0,x+r<=8','angular':'entire unit S3; exact unnormalized registered scalar/vector/tensor coefficient basis'},
         'source_map':{'bytes':len(raw_maps),'sha256':hashlib.sha256(raw_maps).hexdigest(),'status':map_report['status']},
         'normalization':'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.',
         'actual_first_jet_method':'Explicit common correlated Wbar,(P,0,0), plus authenticated actual C3 derivative-error caps and actual background ODE jets. Exterior r>=x uses the separately replayed causal-lift restriction. This is evaluated enclosure of the actual evolved source, not substitution of the comparisons.',
         'arithmetic':{'outward_grid_denominator':str(j.core.GRID),'decisions':'exact rational only','source_interval_engine':'integer grid-tick operations with exactly the reference Fraction rounding order','reference_arithmetic_differential_checks':dyadic_checks,'complete_witness_source_endpoints_identical_to_reference_Fraction_arithmetic':True,'source_cancellations':'exact independent-jet Laurent maps collected before interval evaluation'},
         'bulk_witness':{'x':witness_x.row(),'r':witness_r.row(),'sources':encode(witness),'strict_signs':{'QE/yy/0':'negative','QEred/yy/0':'negative','QE/y0/0':'positive','metric_forcing/yy/0':'positive'}},
         'reflected_witness':{'x':reflected_x.row(),'r':reflected_r.row(),'sources':encode(reflected),'scope':'Actual enclosures after reflection; no source sign is asserted unless the retained interval excludes zero.'},
         'global_cover':{'division_per_unit':d,'cells':len(global_rows),'proof':'Rectangles [i/d,(i+1)/d] x [k/d,(k+1)/d] with 0<=i<3d,0<=k<8d,i+k<8d, each intersected with x+r<=8. Their union is exactly D. Interval comparisons on the rectangular superset only widen actual-D enclosures.','summary':summary,'cell_enclosures':global_rows},
         'evolved_radial_profile':{'x':witness_x.row(),'r':['1','2'],'cells':profile_cells,'gap_free_equal_partition':True,'summary':profile_summary,'cell_enclosures':profile},
         'shell_witness':{'x':shell_x.row(),'r':['0','0'],'sources':encode(shell_values),'one_copy_Qp_mean_strict_positive':True},
         'shell_whole_time_cover':{'x':['0','3'],'cells':3*d,'cell_enclosures':encode(shell_rows)},
         'scope':{'actual_evolved_common_gauge_bulk_sources_evaluated':True,'all_registered_bulk_and_physical_shell_rows_retained':True,
                  'numerical_first_Picard_refinement':False,'numerical_second_response_evaluated':False,'continuous_second_response_residual_certified':False,
                  'actual_smooth_background_master_gauge_existence':'authenticated retained conditional premise','invariant_force_or_finite_amplitude_claim':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--maps',type=Path,default=ROOT.parent/'projection/GEOMETRIC_COEFFICIENT_MAPS.json');p.add_argument('--profile-cells',type=int,default=64);p.add_argument('--cover-divisions',type=int,default=2)
 p.add_argument('--inject-error',default='none',choices=['none','comparator_only','halve_hessian','reverse_metric_forcing','drop_corner_cell','omit_shell_normalization'])
 args=p.parse_args();need(not args.output.exists(),'refuse overwrite');r=run(args.maps,args.profile_cells,args.cover_divisions,args.inject_error)
 args.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'global_cells':r['global_cover']['cells'],'profile_cells':r['evolved_radial_profile']['cells'],'shell_cells':r['shell_whole_time_cover']['cells']},indent=2))
if __name__=='__main__':main()
