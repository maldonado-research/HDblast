#!/usr/bin/env python3
"""Independent symbolic assembly, norm, ledger and six-sample audit."""
import argparse
import ast
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True
BASE = Path(__file__).resolve().parents[2]

def need(ok, marker):
    if not ok: raise ValueError(marker)
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    out = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(out)
    return out
def plus(a,b): return (a[0]+b[0],a[1]+b[1])
def minus(a,b): return (a[0]-b[1],a[1]-b[0])
def times(a,c): return (a[0]*c,a[1]*c) if c>=0 else (a[1]*c,a[0]*c)
def row(a): return tuple(map(F,a))
def contains(a,b): return a[0]<=b[0] and b[1]<=a[1]

def main():
    global BASE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--capsule-root',type=Path,default=BASE)
    parser.add_argument('--normal-results',type=Path)
    parser.add_argument('--optimized-results',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    BASE=args.capsule_root.resolve()
    normal_results=(args.normal_results or BASE/'compute/response-normal').resolve()
    need(not args.output.exists(),'FRESH_OUTPUT_REQUIRED')
    need(BASE not in args.output.resolve().parents,'OUTPUT_INSIDE_CAPSULE')
    paths=[BASE/'compute/full23_response.py',BASE/'compute/angular_ridge.py',
           BASE/'initial/full23_initial.py',BASE/'initial/results/INITIAL_REPORT.json',
           BASE/'initial/results/INITIAL_PROFILES.json']
    paths += [normal_results/name for name in
              ('RESPONSE_REPORT.json','RESPONSE_PROFILES.json','EVALUATED_CANDIDATE.json')]
    input_names=[str(path.relative_to(BASE)) for path in paths[:5]]+['normal-results/'+path.name for path in paths[-3:]]
    before={name:sha(path) for name,path in zip(input_names,paths)}
    report,profiles,evaluated=[json.loads(path.read_text()) for path in paths[-3:]]
    if args.optimized_results is not None:
        for path in paths[-3:]:
            need(path.read_bytes()==(args.optimized_results/path.name).read_bytes(),'NORMAL_OPTIMIZED_OUTPUT_EQUALITY')
    initial=json.loads(paths[3].read_text())
    structural=load('response_review_structure',BASE/'structure/exact_angular_structure.py')
    closure=structural.authenticate(args.repo_root)
    op=load('response_review_frozen_operator',args.repo_root/'research'/structural.ROOTS['G'][0]/'projection/inputs/linear_operator.py')
    # Exact polynomial identity of the NEW assembly body against the retained
    # full operator, with background equations substituted symbolically.
    L=structural.L
    def power(self,n):
        need(isinstance(n,int) and n>=0,'NONNEGATIVE_POWER')
        out=L(1)
        for unused in range(n): out=out*self
        return out
    L.__pow__=power
    z=structural.symbols('ell rho a p eta c Ht V V1 V2')
    rho,a,p,c,Ht,eta=[z[key] for key in ('rho','a','p','c','Ht','eta')]
    h=Ht/rho; d=1/(rho*rho)-a*a; b=d+p*p/3
    prim=op.matrices()
    def matrix_sum(*terms):
        return [[sum((coefficient*matrix[i][j] for coefficient,matrix in terms),L(0))
                 for j in range(23)] for i in range(23)]
    tree=ast.parse(paths[0].read_text())
    cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='Response')
    assemble=next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='assemble')
    nodes=[node for node in assemble.body if isinstance(node,ast.Assign)
           and isinstance(node.targets[0],ast.Name) and node.targets[0].id in ('A','B','C','D','M')]
    need(len(nodes)==5,'ASSEMBLY_AST_TARGETS')
    ns={**z,'F':F,'ht':Ht,'invc':1/c,'h':h,'d':d,'b':b,'m':prim,
        'self':SimpleNamespace(ell=z['ell'],matrix_sum=matrix_sum)}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),'<new-assembly-exact-audit>','exec'),ns)
    expected=op.assemble(ell=z['ell'],rho=rho,a=a,h=h,c=c,Ht=Ht,p=p,
                         py=z['V1']-4*a*p,b=b,d=d,V=z['V'],V1=z['V1'],V2=z['V2'],
                         sigma1=-2*p,sigma2=4*(1+eta))
    need(all(ns[key]==expected[key] for key in ('A','B','D','M')),'FULL_SYMBOLIC_COEFFICIENT_IDENTITY')
    need(max(-prim['Delta'][i][i] for i in range(23))==24,'ANGULAR_ROUGH_SPECTRAL_CAP')

    # Reconstruct the fresh coefficient boxes only; no old producer is run.
    module=load('new_response_reviewed',paths[0])
    response=module.Response(args.repo_root)
    need(response.local_before==report['local_dependency_hashes'],'LOCAL_DEPENDENCIES_RECORDED')
    for key in ('A','B','D','M'):
        need([[str(v) for v in r] for r in response.matrices[key]]==report['constant_matrices'][key], 'CONSTANT_MATRIX_REPLAY_'+key)
        need([[v.row() for v in r] for r in response.actual_boxes[key]]==report['actual_coefficient_boxes'][key], 'ACTUAL_BOX_REPLAY_'+key)
    weights=list(map(F,report['gram']))
    need(weights==op.gram(),'GRAM_CONVENTION')
    blocks=list(report['orthogonal_blocks'].values())
    ids=[[op.IDX[label] for label in block] for block in blocks]
    norm_checks=0
    for table_name,block_name,global_name in (
      ('constant_matrices','constant_matrix_block_H_operator_upper','constant_matrix_H_operator_upper'),
      ('coefficient_entry_error_caps','coefficient_error_block_H_operator_upper','coefficient_error_H_operator_upper')):
        for key in ('A','B','D','M'):
            matrix=[list(map(F,r)) for r in report[table_name][key]]
            membership={i:k for k,block in enumerate(ids) for i in block}
            need(all(not matrix[i][j] or membership[i]==membership[j] for i in range(23) for j in range(23)), 'ORTHOGONAL_BLOCK_SUPPORT')
            caps=[]
            for name,block in zip(report['orthogonal_blocks'],ids):
                cap=F(report[block_name][key][name])
                squared=sum((weights[i]*matrix[i][j]**2/weights[j] for i in block for j in block),F(0))
                need(cap>=0 and cap*cap>=squared,'WEIGHTED_FROBENIUS_UPPER')
                caps.append(cap);norm_checks+=1
            need(F(report[global_name][key])==max(caps),'ORTHOGONAL_MAX_NORM')
    for key in ('A','B','D','M'):
        for i in range(23):
            for j in range(23):
                lo,hi=map(F,report['actual_coefficient_boxes'][key][i][j]); mid=F(report['constant_matrices'][key][i][j])
                need(F(report['coefficient_entry_error_caps'][key][i][j])>=max(abs(lo-mid),abs(hi-mid)), 'COEFFICIENT_ENTRY_ERROR_COVERAGE')

    ledger={k:F(v) for k,v in report['ledger'].items()}
    nu,ndu,nv=map(F,initial['nominal_H_caps'])
    eu,edu,ev=map(F,initial['total_H_caps'])
    pu,pdu,pv=map(F,initial['P1_H_error_caps'])
    x0=max((nu+F(3,2)*(nv+pv))/3,ndu+nv+pv)
    efree=max((eu+F(3,2)*ev)/3,edu+ev)
    q=F(3409,8000000); source=F(52162547172519,2)
    e0=efree+pv+source
    coeff={k:F(v) for k,v in report['coefficient_error_H_operator_upper'].items()}
    caps={k:F(v) for k,v in report['constant_matrix_H_operator_upper'].items()}
    dfc=(coeff['A']+coeff['B']+F(3,4)*coeff['D'])*x0
    dfp=(caps['A']+caps['B'])*pdu+caps['D']*pu
    dgc=F(1,4)*coeff['M']*x0; dgp=caps['M']*pu
    e1=F(3,2)*(dfc+dfp)+dgc+dgp
    ucap=(x0+e0)/(1-q); error=q*q*ucap+(1+q)*e0+e1
    comparisons={'B0_X_upper':x0,'free_input_error_X_upper':efree,'e0':e0,
                 'delta_f_coefficient':dfc,'delta_f_profile':dfp,
                 'delta_g_coefficient':dgc,'delta_g_profile':dgp,'e1':e1,
                 'actual_response_X_upper':ucap,'total_X_error_upper':error}
    need(all(ledger[k]==v for k,v in comparisons.items()),'INDEPENDENT_FULL_LEDGER')

    # Direct integration of each velocity segment, independent of hinge sums.
    nodes=list(map(F,profiles['nodes'])); values=[list(map(F,r)) for r in profiles['profiles']['v']]
    parity=response.kernel.PARITY
    def direct_v(t,primitive=False):
        negative=t<0; t=abs(t); out=[F(0)]*23
        for left,right,u,v in zip(nodes,nodes[1:],values,values[1:]):
            if t<left: break
            z=min(t,right)-left
            slope=[(v[i]-u[i])/(right-left) for i in range(23)]
            if primitive:
                out=[out[i]+z*u[i]+z*z*slope[i]/2 for i in range(23)]
            elif left<=t<=right:
                out=[u[i]+z*slope[i] for i in range(23)];break
            if t<=right: break
        if negative: out=[v*parity[i]*(-1 if primitive else 1) for i,v in enumerate(out)]
        return out
    candidate=response.kernel.from_data(profiles)
    sample_checks=0; ordinary_time_checks=0
    for sample in evaluated['samples']:
        x,r=F(sample['x']),F(sample['r']); aa,bb=r+x,r-x
        def initial_at(t):
            data,_=response.context.profiles(response.I(abs(t)))
            return [[times(row(box.row()),parity[i] if t<0 else 1) for i,box in enumerate(data[0])],
                    [times(row(box.row()),-parity[i] if t<0 else 1) for i,box in enumerate(data[1])]]
        ua,da=initial_at(aa); ub,db=initial_at(bb)
        va,vb=direct_v(aa),direct_v(bb); pa,pb=direct_v(aa,True),direct_v(bb,True)
        free={'value':[times(plus(plus(ua[i],ub[i]),(pa[i]-pb[i],)*2),F(1,2)) for i in range(23)],
              'dx':[times(plus(minus(da[i],db[i]),(va[i]+vb[i],)*2),F(1,2)) for i in range(23)],
              'dr':[times(plus(plus(da[i],db[i]),(va[i]-vb[i],)*2),F(1,2)) for i in range(23)]}
        feedback=candidate.feedback(x,r,response.matrices['M'])
        for kind in ('value','dx','dr'):
            need([str(v) for v in feedback[kind]]==sample['feedback_exact'][kind],'SAVED_FEEDBACK_REPLAY')
            for i in range(23):
                saved_free=row(sample['free_wave_interval'][kind][i])
                need(contains(saved_free,free[kind][i]),'DIRECT_SEGMENT_FREE_WAVE_ENCLOSURE')
                complete=row(sample['complete_candidate_interval_evaluation'][kind][i])
                need(contains(complete,plus(free[kind][i],(feedback[kind][i],)*2)),'COMPLETE_SAMPLE_ENCLOSURE')
                lower,upper=row(sample['conditional_exact_response_component_enclosure'][kind][i])
                radius=min(complete[0]-lower,upper-complete[1]);factor=3 if kind=='value' else 1
                need(radius>=0 and radius*radius*weights[i]>=factor*factor*error*error,'COMPONENT_ERROR_RADIUS')
                sample_checks+=1
        for i in range(23):
            direct_dx=plus(free['dx'][i],(feedback['dx'][i],)*2)
            saved=row(sample['complete_candidate_interval_evaluation']['ordinary_t_velocity'][i])
            need(contains(saved,times(direct_dx,1/response.ell)),'ORDINARY_TIME_VELOCITY_CONVERSION')
            ordinary_time_checks+=1
    need(sample_checks==6*3*23,'ALL_SIX_FULL23_SAMPLES')
    need(structural.authenticate(args.repo_root)==closure,'FROZEN_CLOSURES_CHANGED')
    need({name:sha(path) for name,path in zip(input_names,paths)}==before,'REVIEWED_FILES_CHANGED')
    result={'status':'PASS_INDEPENDENT_FULL23_COEFFICIENT_LEDGER_SAMPLE_REVIEW',
            'reviewed_files':before,'symbolic_entries_checked':4*23*23,
            'symbolic_background_substitutions':['b=d+p^2/3','py=V1-4ap','sigma1=-2p','sigma2=4(1+eta)','h=Ht/rho'],
            'weighted_block_norm_bounds_checked':norm_checks,'all_entry_error_caps_cover_coefficient_boxes':True,
            'all_ledger_components_recomputed_from_initial_report':True,
            'direct_velocity_segment_free_wave_component_checks':sample_checks,
            'ordinary_time_velocity_checks':ordinary_time_checks,
            'normal_optimized_all_output_files_identical':True if args.optimized_results is not None else None,
            'feedback_sample_replay_with_independently_reviewed_kernel':True,
            'angular_rough_spectral_cap':24,'total_X_error_upper':str(error),
            'limitations':['Inherited actual-background interval premises retained','Source omission bound and contraction theorem retained conditionally','Kernel formulas rely on separate universal polynomial review','No frozen producer rerun','No pointwise C2 Einstein residual or nonlinear realization']}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'sample_checks':sample_checks,'error_display':float(error)}))

if __name__=='__main__':main()
