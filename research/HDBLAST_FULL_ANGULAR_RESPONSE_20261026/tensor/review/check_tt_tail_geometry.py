#!/usr/bin/env python3
"""Independent exact cone geometry and scalar resolvent arithmetic.

No project modules or source oracle are imported.  The triangle is clipped
against each rectangular cell and its area is computed by the shoelace rule.
"""
from fractions import Fraction as Q
import argparse,json
from pathlib import Path

T=Q(1,4); R=Q(5,4)
PROTECTED_CAPSULES=('HDBLAST_COUPLED_MEAN_PILOT_20261025','HDBLAST_INPUT_PRECISION_20261023','HDBLAST_COMMON_GAUGE_FORCING_20261016','HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016','HDBLAST_GEOMETRIC_BULK_SOURCES_20261016','HDBLAST_FORCED_RESPONSE_BOUNDS_20261023','HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023','HDBLAST_SHELL_FEEDBACK_PRECISION_20261025')

def need(ok):
    if not ok:
        raise ValueError("EXACT_REVIEW_CHECK_FAILED")

def guard_output(output,repo_root):
    target=output.resolve()
    protected=[(repo_root/'research'/name).resolve() for name in PROTECTED_CAPSULES]
    if output.exists() or output.is_symlink():
        raise ValueError('FRESH_OUTPUT_REQUIRED')
    if any(target==root or root in target.parents for root in protected):
        raise ValueError('OUTPUT_INSIDE_PROTECTED_CAPSULE')

def clip(poly, axis, edge, keep_greater):
    if not poly:
        return []
    out=[]
    prev=poly[-1]
    pin=prev[axis]>=edge if keep_greater else prev[axis]<=edge
    for point in poly:
        qin=point[axis]>=edge if keep_greater else point[axis]<=edge
        if pin!=qin:
            ratio=(edge-prev[axis])/(point[axis]-prev[axis])
            out.append(tuple(prev[j]+ratio*(point[j]-prev[j]) for j in (0,1)))
        if qin:
            out.append(point)
        prev,pin=point,qin
    return out

def triangle_intersection_area(xlo,xhi,rlo,rhi):
    poly=[(Q(0),R-T),(Q(0),R+T),(T,R)]
    for axis,edge,greater in ((0,xlo,True),(0,xhi,False),(1,rlo,True),(1,rhi,False)):
        poly=clip(poly,axis,edge,greater)
    if not poly:
        return Q(0)
    return abs(sum((p[0]*q[1]-q[0]*p[1] for p,q in zip(poly,poly[1:]+poly[:1])),Q(0)))/2

def partition(n):
    h=T/n; rows=[]
    for i in range(n):
        for k in range(i,2*n-i):
            box=(i*h,(i+1)*h,R-T+k*h,R-T+(k+1)*h)
            actual=triangle_intersection_area(*box)
            declared=h*h/(2 if k in (i,2*n-i-1) else 1)
            need(actual==declared and actual>0)
            rows.append({'i':i,'k':k,'exact_area':str(actual)})
    need(len(rows)==n*(n+1))
    need(sum((Q(row['exact_area']) for row in rows),Q(0))==T*T)
    # Every bounding rectangle lies in the exterior actual-oracle domain.
    need(all(R-T+k*h-(i+1)*h>=Q(3,4) for i in range(n) for k in range(i,2*n-i)))
    # Geometric lower boundary of the true cone: r-x>=1.
    need(R-T==1 and R+T<Q(8))
    # Mutation controls: a dropped cell and doubled edge triangle fail area.
    need(sum((Q(row['exact_area']) for row in rows[:-1]),Q(0))!=T*T)
    need(sum((Q(row['exact_area']) for row in rows),Q(0))+h*h/2!=T*T)
    return {'n':n,'cells':len(rows),'area':str(T*T),'method':'triangle clipped by four halfplanes; exact rational shoelace area','rows':rows,'mutation_controls_passed':2}

def scalar_tail(a,b,d,f):
    """Output exact bound for nonnegative scalar coefficient/source caps."""
    need(min(a,b,d,f)>=0)
    sigma=a+b+T*d/2
    q=T*sigma
    need(q<1)
    xbound=T*f/(1-q)
    feedback=T*T*sigma*xbound/2
    need(feedback==T*T*f*q/(2*(1-q)))
    return {'A_upper':str(a),'B_upper':str(b),'D_upper':str(d),'F_upper':str(f),'sigma':str(sigma),'q':str(q),'X_upper':str(xbound),'apex_feedback_upper':str(feedback)}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    guard_output(args.output,args.repo_root.resolve())
    result={'status':'PASS_INDEPENDENT_EXACT_TT_TAIL_AND_CONE_GEOMETRY','partitions':[partition(4),partition(8)],'tail_formula':'q=T*(a+b+T*d/2); X<=T*f/(1-q); E<=T^2*f*q/(2*(1-q))','arithmetic_identity_control':scalar_tail(Q(1,10**12),Q(1,20000),Q(1,10**8),Q(1000)),'limitations':['No source oracle executed.','Numerical forcing and coefficient caps must be authenticated separately.','The norm proof requires zero-data decoupled scalar TT dynamics and exterior domain of dependence.']}
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True)
        handle.write('\n')
    print(json.dumps({key:result[key] for key in ('status','tail_formula')},indent=2))
