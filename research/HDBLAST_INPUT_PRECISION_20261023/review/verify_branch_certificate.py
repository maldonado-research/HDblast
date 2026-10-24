#!/usr/bin/env python3
"""Independent branch certificate verifier. Does not import the producer."""
import argparse,hashlib,json,math
from fractions import Fraction as Q
from pathlib import Path
checks=[]
def need(ok,why):
 if not ok:raise ValueError(why)
 checks.append(why)
def plus(p,q):
 r=dict(p)
 for k,v in q.items():r[k]=r.get(k,Q(0))+v
 return {k:v for k,v in r.items() if v}
def times(p,q):
 r={}
 for (i,j),v in p.items():
  for (k,l),w in q.items():r[(i+k,j+l)]=r.get((i+k,j+l),Q(0))+v*w
 return {k:v for k,v in r.items() if v}
def sc(p,a):return {k:v*a for k,v in p.items() if v*a}
def der(p,axis):
 r={}
 for k,v in p.items():
  n=k[axis]
  if n:
   kk=list(k);kk[axis]-=1;r[tuple(kk)]=v*n
 return r
def point(p,x,h):return sum((v*x**i*h**j for (i,j),v in p.items()),Q(0))
def mag(p,x,h):return sum((abs(v)*x**i*h**j for (i,j),v in p.items()),Q(0))
def bounds(p,x,h):
 c=p.get((0,0),Q(0));r=mag(plus(p,{(0,0):-c}),x,h)
 return c-r,c+r
# Intervals use plain pairs and independent arithmetic rather than producer I.
def ia(a,b):return a[0]+b[0],a[1]+b[1]
def im(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def ip(a,n):
 if not n:return Q(1),Q(1)
 if n%2:return a[0]**n,a[1]**n
 return (Q(0) if a[0]<=0<=a[1] else min(a[0]**n,a[1]**n),max(a[0]**n,a[1]**n))
def box(p,x,h):
 r=(Q(0),Q(0))
 for (i,j),v in p.items():r=ia(r,im((v,v),im(ip(x,i),ip(h,j))))
 return r
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--branch-dir',type=Path,default=Path('/tmp/hdblast-input-precision-20261023/branch'));ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 root=args.branch_dir;raw=(root/'GRAPH_POLYNOMIALS.json').read_bytes();g=json.loads(raw);report=json.loads((root/'BRANCH_COMPARISON_REPORT.json').read_text())
 need(hashlib.sha256(raw).hexdigest()==report['graph_sha256'],'graph_sha256')
 P={s:{(int(t['i']),int(t['j'])):Q(int(t['num']),int(t['den'])) for t in g[s]} for s in ('U','M','Z','A','B','R','D')}
 U,M,Z,A,B,R,D=[P[k] for k in ('U','M','Z','A','B','R','D')]
 one={(0,0):Q(1)};x={(1,0):Q(1)};h={(0,1):Q(1)}
 W={(0,0):Q(1,3),(2,0):Q(1),(3,0):Q(1,3)};Wp=der(W,0)
 need(U==plus(sc(times(Wp,Wp),Q(1,2)),sc(times(W,W),-Q(2,3))),'potential_identity')
 need(times(M,x)==der(U,0),'divided_potential_identity')
 need(B==plus(sc(times(h,der(Z,1)),2),sc(Z,-4)),'graph_transport_B_identity')
 Rcomputed=plus(plus(plus(M,sc(times(Z,Z),-1)),sc(times(times(x,der(Z,0)),Z),-1)),times(A,B))
 Dcomputed=plus(plus(plus(h,sc(times(times(x,x),times(Z,Z)),Q(1,12))),sc(U,-Q(1,6))),sc(times(A,A),-1))
 need(R==Rcomputed,'graph_residual_exact_identity');need(D==Dcomputed,'metric_square_residual_exact_identity')
 need(min(i+j for i,j in R)>=5 and min(i+j for i,j in D)>=5,'residual_vanishing_order')
 e=Q(9,100000);h0=Q(1,100);amin=Q(1,9)-Q(1,1000000);amax=Q(3,20);kap=Q(37,20)
 zlo,zhi=bounds(Z,e,h0);alo,ahi=bounds(A,e,h0);blo,bhi=bounds(B,e,h0);xz=mag(times(x,der(Z,0)),e,h0)
 need(0<zlo<zhi<2,'entry_graph_between_zero_and_two');need(alo>0 and bhi<0,'square_root_and_secant_signs')
 damping=4*amin+zlo-xz
 need(damping>=kap and kap>12*amax,'uniform_damping_and_memory_margin')
 need(Q(report['exact_checks']['secant_damping_lower']['exact'])==damping,'reported_damping_exact')
 C=[Q(0) for _ in range(6)];fac=max(abs(blo),abs(bhi))/(amin+alo)
 for poly,mult in ((R,Q(1)),(D,fac)):
  for (i,j),v in poly.items():
   degree=min(j,5);C[degree]+=abs(v)*mult*e**i*h0**(j-degree)
 Ej=[C[j]/(kap-2*j*amax) for j in range(6)]
 need(C==list(map(Q,report['residual_majorant_Cj'])),'reported_residual_majorant_exact')
 need(Ej==list(map(Q,report['barrier_Ej'])),'reported_barrier_coefficients_exact')
 c=Q(96714321425692,161839258930731);delta=Q(1,1000)
 f={(0,0):1+c,(1,0):c}
 A_shell=plus(sc(times(W,f),Q(1,9)),sc(Wp,-c/Q(12)))
 B_shell=plus(sc(times(f,f),Q(1,36)),sc(one,-c*c/Q(48)))
 H=plus(sc(A_shell,delta),sc(B_shell,delta*delta));Hp=der(H,0)
 center=-Q(9,64)*c*delta;oldrad=Q(799,250)*delta**2;old=(center-oldrad,center+oldrad)
 need(old==tuple(map(Q,report['old_beta'])),'retained_beta_interval_exact')
 Hb=box(H,old,(Q(0),Q(0)));need(Hb==tuple(map(Q,report['old_H2'])),'shell_curvature_box_exact')
 need(0<Hb[0]<Hb[1]<h0 and max(abs(v) for v in old)<e,'comparison_domain_reached_before_shell')
 E=2*(Hb[1]/h0)**6+sum((Ej[j]*Hb[1]**j for j in range(6)),Q(0));scalar=e*E
 need(E==Q(report['endpoint_error_z']['exact']),'endpoint_slope_error_exact')
 need(scalar==Q(report['endpoint_error_scalar_residual']['exact']),'endpoint_scalar_error_exact')
 gp=ia(ia(box(Z,old,Hb),im(old,ia(box(der(Z,0),old,Hb),im(box(der(Z,1),old,Hb),box(Hp,old,(Q(0),Q(0))))))),box(der(Wp,0),old,(Q(0),Q(0))))
 need(gp==tuple(map(Q,report['scalar_polynomial_derivative'])) and gp[0]>0,'scalar_derivative_interval_exact_and_positive')
 seed=Q(report['seed_steps'][-1]['beta_seed']);hs=point(H,seed,Q(0));gres=seed*point(Z,seed,hs)+point(Wp,seed,Q(0))+delta*c/2
 rad=(abs(gres)+scalar)/gp[0]
 need(rad==Q(report['pre_rounding_radius']['exact']),'mean_value_radius_exact')
 beta=tuple(map(Q,report['beta']))
 need(old[0]<=beta[0]<=seed-rad<=seed+rad<=beta[1]<=old[1],'outward_beta_enclosure_inside_retained_interval')
 need((beta[1]-beta[0])/2==Q(report['beta_radius']['exact']) and (beta[1]-beta[0])/2<Q(1,10**11),'target_radius_exact')
 # Independent source-tree verification against the manifest hash already
 # frozen externally in the current source selection. No source execution.
 source=Path('/workspace/hdblast-research-work/publish-HDblast-20261014/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011');manifest_raw=(source/'MANIFEST.json').read_bytes();manifest=json.loads(manifest_raw)
 need(hashlib.sha256(manifest_raw).hexdigest()=='917480172ae98b2bf39f29d37fe0ad3cf279ad0d8508c98efd469efdba26c0ec','selected_complete_source_manifest_pin')
 expected={'MANIFEST.json'}
 for r in manifest['files']:
  q=Path(r['path']);need(not q.is_absolute() and '..' not in q.parts and str(q) not in expected,'source_path_'+str(q));expected.add(str(q));p=source/q;b=p.read_bytes();need(not p.is_symlink() and len(b)==r['bytes'] and hashlib.sha256(b).hexdigest()==r['sha256'],'source_payload_'+str(q))
 need({str(p.relative_to(source)) for p in source.rglob('*') if p.is_file()}==expected,'complete_source_inventory')
 out={'status':'PASS_INDEPENDENT_BRANCH_ARITHMETIC_AND_PROOF_REVIEW','checks_passed':len(checks),'checks':checks,'producer_sha256':hashlib.sha256((root/'validate_branch.py').read_bytes()).hexdigest(),'graph_sha256':hashlib.sha256(raw).hexdigest(),'report_sha256':hashlib.sha256((root/'BRANCH_COMPARISON_REPORT.json').read_bytes()).hexdigest(),'verified_beta':[str(v) for v in beta],'verified_radius':str((beta[1]-beta[0])/2),'source_payload_count':len(manifest['files']),'review_scope':'Independent exact graph identities, finite-domain signs, coefficientwise defect/barrier, actual shell map, mean-value contraction, source tree, and accompanying scalar-comparison argument. No producer import or scientific source import.','premises':['Retained theorem existence/uniqueness in its specified Dirichlet tube, convexity, regular-cone signs and positive a.','Elementary scalar ODE comparison; regular cone gives z->0, M<4 and a>0 imply 0<=z<=2.','Proof reviewed internally; no external peer review or formal proof assistant.'],'limitations':['Necessary enclosure of the retained existent root; not an independently solved endpoint IVP.','No local initial-function allocation, total U, B1 or physical observable sign follows from this review.']}
 args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['checks_passed'])
if __name__=='__main__':main()
