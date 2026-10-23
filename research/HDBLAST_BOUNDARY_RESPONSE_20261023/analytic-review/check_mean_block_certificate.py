#!/usr/bin/env python3
"""Exact five-mean-block certificate, conditional on pinned actual source inputs.

No source/PDE evaluations. Run python -I -S -B, optionally -O.
Source interval integration is separately audited and joined by an exact pin.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, importlib.util, json, sys

ROOT=Path('/workspace/hdblast-research-work/geometric-bulk-source-20261016/capsule-source')
WORK=Path('/workspace/hdblast-research-work/boundary-response-20261023')
COMPUTE=WORK/'compute'
N16_PIN='6caf16c304fb6db940899d58d09b648066eb2356385189ccb79eb2b0a6879675'
PINS={
    'projection/inputs/linear_operator.py':'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189',
    'projection/angular_generic.py':'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9',
    'numeric/ACTUAL_SOURCES_NORMAL.json':'45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
    'geometry/inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md':'5de45b21236cfe8c6d99fb403bacc4624069c8d96d6cf0ad2ec89b8e2d3617e1',
}

def need(ok,s):
    if not ok: raise ValueError(s)

def pin(p):
    b=p.read_bytes()
    if p==Path(__file__): label='analytic-review/check_mean_block_certificate.py'
    elif p.is_relative_to(ROOT): label='old/'+p.relative_to(ROOT).as_posix()
    elif p.is_relative_to(COMPUTE): label='compute/'+p.relative_to(COMPUTE).as_posix()
    else: raise ValueError('unregistered pin path')
    return {'path':label,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()}

def add(*args):
    return [[sum((m[i][j] for m in args),Q(0)) for j in range(5)] for i in range(5)]

def mul(c,m): return [[c*v for v in row] for row in m]
def absm(m): return [[abs(v) for v in row] for row in m]
def mv(m,v): return [sum((x*y for x,y in zip(row,v)),Q(0)) for row in m]

def solve(m,v):
    w=[row[:]+[x] for row,x in zip(m,v)]
    for j in range(5):
        k=next(i for i in range(j,5) if w[i][j])
        w[j],w[k]=w[k],w[j]
        den=w[j][j];w[j]=[x/den for x in w[j]]
        for i in range(5):
            if i!=j:
                c=w[i][j];w[i]=[x-c*y for x,y in zip(w[i],w[j])]
    return [row[-1] for row in w]

def decimal_bracket(value,places=6):
    scale=10**places
    lo=(value*scale).__floor__();hi=(value*scale).__ceil__()
    def fmt(n):
        sign='-' if n<0 else '';whole,frac=divmod(abs(n),scale)
        return sign+str(whole)+'.'+str(frac).zfill(places)
    return [fmt(lo),fmt(hi)]

def calculate():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,'run python -I -S -B')
    paths=[ROOT/name for name in PINS]
    pins=[pin(p) for p in paths]
    for row in pins:
        name=row['path'].removeprefix('old/')
        need(row['sha256']==PINS[name],'authenticated input before import '+name)
    own=pin(Path(__file__))
    spec=importlib.util.spec_from_file_location('mean_block_readonly_operator',paths[0])
    op=importlib.util.module_from_spec(spec);spec.loader.exec_module(op)
    labels=['K:0','L:0','N:0','J:0','w:0'];ii=[op.IDX[k] for k in labels]
    full=op.matrices()
    need(set(full)=={'I','O','Sa','Sh','Sah','Rb','Rd','O2','Delta','Ca','Ch','Aphi','Bphi','Ephi','Pmet','Pg','Pww','W_K','W_L','By_a','By_p','By_sigma1','By_sigma2','PD','PN'},'every retained primitive included')
    for name,m in full.items():
        need(all(m[i][j]==m[j][i]==0 for i in ii for j in range(23) if j not in ii),'exact block closure '+name)
    mats={name:[[m[i][j] for j in ii] for i in ii] for name,m in full.items()}
    need([op.gram()[i] for i in ii]==[1,2,1,3,1],'positive Gram block')
    need([mats['PD'][i][i] for i in range(5)]==[0,1,0,0,0],'mean boundary projector')
    angular_spec=importlib.util.spec_from_file_location('mean_block_readonly_angular',ROOT/'projection/angular_generic.py')
    angular=importlib.util.module_from_spec(angular_spec);angular_spec.loader.exec_module(angular)
    need([angular.inner(v,angular.SCALARS[0]) for v in angular.SCALARS]==[1,0,0],'probability scalar means')
    need([angular.inner(v,angular.TENSORS[0]) for v in angular.TENSORS]==[3,0,0,0,0,0,0],'sphere metric mean projection and tensor orthogonality')
    for name in ['Delta','Ca','Ch','Ephi']:
        need(all(v==0 for row in mats[name] for v in row),'vanishing angular primitive '+name)

    ell,R,alpha,P=Q(1,10**6),Q(131),Q(14,125),Q(7,50000)
    Ht=Q(3,10**6);h=Ht/129
    bc=alpha**2+Q(1,10000);dc=Q(1,129**2)+alpha**2
    def term(c,name): return mul(c,absm(mats[name]))
    AA=add(term(3*ell*Ht,'I'),term(2*ell*R*alpha,'O'),term(ell*R*P,'Aphi'))
    BB=add(term(3*ell*R*alpha,'I'),term(ell*R*P,'Bphi'))
    CC=add(mul(alpha**2,absm(add(mats['Sa'],mul(-1,mats['O2'])))),
        mul(alpha*h,absm(add(mats['Sah'],mul(-3,mats['O'])))),term(h*h,'Sh'),
        term(2*bc,'Rb'),term(2*dc,'Rd'),term(Q(4,3),'Pmet'),term(Q(4,3000),'Pg'),
        term(4,'Pww'),term(Q(1,500)+4*alpha*P,'W_K'),term(3*P*h,'W_L'))
    DD=mul(ell**2*R**2,CC)
    MM=mul(ell*R,add(term(alpha,'By_a'),term(P,'By_p'),term(2*P,'By_sigma1'),term(4,'By_sigma2')))
    sigma=add(AA,BB,mul(3,DD))
    contraction=add(mul(Q(3,2),sigma),mul(3,MM))
    q=max(sum(row) for row in contraction)
    need(q<1,'positive matrix resolvent convergence')
    report=json.loads((ROOT/'numeric/ACTUAL_SOURCES_NORMAL.json').read_text())
    need(report['status']=='PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES','actual source input status')
    need(report['normalization']=='Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.','source conventions')
    caps=report['global_cover']['summary']['absolute_coefficient_caps']
    fc=[Q(caps[k]) for k in ['metric_forcing/yy/0','metric_forcing/y0/0','metric_forcing/00/0','metric_forcing/ab/0','scalar_forcing/0']]
    shell=report['shell_whole_time_cover']
    need(shell['x']==['0','3'],'whole boundary-time cover')
    shell_intervals=[tuple(map(Q,cell['x'])) for cell in shell['cell_enclosures']]
    need(shell_intervals[0][0]==0 and shell_intervals[-1][1]==3 and all(lo<hi for lo,hi in shell_intervals) and all(left[1]==right[0] for left,right in zip(shell_intervals,shell_intervals[1:])),'continuous gap-free inherited shell source cover')
    def cap(key): return max(abs(Q(v)) for cell in shell['cell_enclosures'] for v in cell['sources']['one_copy_outward_hessian'][key][0])
    q00,qab,qp,qcy=cap('00'),cap('ab'),cap('scalar'),Q(caps['QC/y/0'])
    gc=[2*ell*R*(q00+3*qab+qcy),Q(0),2*ell*R*q00,2*ell*R*qab,ell*R*qp]
    sc=[Q(3,2)*f+g for f,g in zip(fc,gc)]
    identity=[[Q(i==j) for j in range(5)] for i in range(5)]
    bound=solve(add(identity,mul(-1,contraction)),sc)
    need(all(x>=0 for x in bound),'nonnegative component bounds')
    need(bound==[s+x for s,x in zip(sc,mv(contraction,bound))],'exact vector supersolution')
    ebulk=Q(5,8)*mv(sigma,bound)[0]
    erobin=Q(1,32)*mv(MM,bound)[0]
    total=ebulk+erobin
    # Positive weighted maximum norm with source-sensitive weights bound.
    weighted_q=max(mv(contraction,bound)[i]/bound[i] for i in range(5))
    need(weighted_q<1,'weighted maximum norm contraction')
    old=Q(14866325944167915,3992782)
    n16path=COMPUTE/'BOUNDARY_N16_NORMAL.json'
    n16opt=COMPUTE/'BOUNDARY_N16_OPTIMIZED.json'
    need(pin(n16path)['sha256']==N16_PIN and n16path.read_bytes()==n16opt.read_bytes(),'pinned identical normal and optimized source integrals')
    n16=json.loads(n16path.read_text())
    need(n16['status']=='CONDITIONAL_ACTUAL_BOUNDARY_AUXILIARY_MEAN_ENCLOSURE','source interval status')
    need(n16['partition']['n']==16,'source interval partition')
    need(n16['domain']['a']=='3/2' and n16['domain']['b']=='1/4' and n16['domain']['d']=='5/4','same apex and history')
    slo,shi=map(Q,n16['signed_integral_ledger']['S_equals_KF_plus_LG_mean'])
    need(slo<=shi,'source interval order')
    wlo,whi=slo-total,shi+total
    need(wlo>0,'strict positive mean margin after complete block feedback')
    need(pins==[pin(p) for p in paths] and own==pin(Path(__file__)),'read-only mathematical inputs and checker unchanged')
    result={'status':'PASS_CONDITIONAL_MEAN_BLOCK_POSITIVE_POINT_CERTIFICATE',
        'source_evaluations':0,'labels':labels,'gram':['1','2','1','3','1'],
        'exact_two_way_primitive_block_closure':True,'input_pins':pins,
        'primitive_names':sorted(full),'source_integral_pins':[pin(n16path),pin(n16opt)],'checker_sha256':own['sha256'],
        'background_majorants':{k:str(v) for k,v in {'ell':ell,'rho_abs':R,'a_abs':alpha,'p_abs':P,'Ht_abs':Ht,'h_abs':h,'curvature_b_abs':bc,'curvature_d_abs':dc,'V_abs':Q(1),'V1_abs':Q(1,1000),'V2_abs':Q(4),'py_abs':Q(1,500),'sigma1_abs':2*P,'sigma2_abs':Q(4)}.items()},
        'component_source_caps':{'F':[str(v) for v in fc],'G':[str(v) for v in gc]},
        'component_operator_majorants':{k:[[str(v) for v in row] for row in m] for k,m in {'A':AA,'B':BB,'D0':DD,'M':MM,'X_contraction':contraction}.items()},
        'component_X_bounds':[str(v) for v in bound],
        'unweighted_max_row_contraction':str(q),'source_weighted_max_row_contraction':str(weighted_q),
        'bulk_feedback_value_bound':str(ebulk),'Robin_feedback_value_bound':str(erobin),'total_feedback_value_bound':str(total),
        'old_feedback_value_bound':str(old),'improvement_factor':str(old/total),
        'evaluated_source_integral_interval':[str(slo),str(shi)],
        'actual_Wyy_probability_mean_interval':[str(wlo),str(whi)],
        'strict_positive':wlo>0,'strict_positive_margin':str(wlo),
        'outward_decimal_brackets':{k:decimal_bracket(v) for k,v in {'source_lower':slo,'source_upper':shi,'bulk_feedback':ebulk,'Robin_feedback':erobin,'total_feedback':total,'W_lower':wlo,'W_upper':whi}.items()},
        'limits':['Conditional on retained smooth actual first solution, full linear operator, source caps, and separate source-integral enclosure review.',
            'Quadratic sources include all first-field harmonics before projection.',
            'Exact block closure retains every possible contribution to the mean; the complement remains in the full response but cannot feed these five means linearly.',
            'Componentwise comparison uses a positive weighted max norm; it does not drop positive Gram factors from an H norm.',
            'This is one auxiliary forced-increment mean, not total physical U, a full waveform, or a standalone sourced physical-constraint solution.',
            'No source recalculation, optimization search, or numerical feedback field.']}
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source-root',type=Path,default=ROOT)
    parser.add_argument('--compute-root',type=Path,default=COMPUTE)
    args=parser.parse_args();ROOT=args.source_root.resolve();COMPUTE=args.compute_root.resolve()
    need(not args.output.exists(),'refuse overwrite output')
    result=calculate();args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status']);print(json.dumps(result['outward_decimal_brackets'],indent=2))
