#!/usr/bin/env python3
"""Independent exact one-point geometry/arithmetic audit; no source evaluation.

Synthetic controls are exact operator witnesses, not actual response samples.
All retained inputs are opened read-only. Explicit checks survive python -O.
Input roots are configurable; pin paths are stable role-relative names.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True

BASE = Path('/workspace/hdblast-research-work')
FROZEN = BASE/'forced-response-20261023/package-work/attempt-2/source'
FOLLOW = BASE/'forced-response-20261023/followon-boundary-formulation'
OLD = BASE/'geometric-bulk-source-20261016/capsule-source'

def need(ok, message):
    if not ok:
        raise ValueError(message)

def pin(path, label=None):
    raw = path.read_bytes()
    return {'path':str(path) if label is None else label,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}

def area(poly):
    return abs(sum((x*v-y*u for (x,y),(u,v) in zip(poly,poly[1:]+poly[:1])),Q(0)))/2

def run(frozen_root=FROZEN, follow_root=FOLLOW, old_root=OLD):
    FROZEN, FOLLOW, OLD = map(Path, (frozen_root, follow_root, old_root))
    def input_pin(path):
        for role, root in (('frozen',FROZEN),('follow',FOLLOW),('old',OLD)):
            try:
                relative = path.relative_to(root)
            except ValueError:
                continue
            return pin(path,role+'/'+relative.as_posix())
        raise ValueError('Input outside declared roots: '+str(path))
    inputs = [FOLLOW/name for name in ('NEXT_BOUNDARY_CALCULATION.md',
        'ONE_POINT_FEEDBACK_ENCLOSURE_DRAFT.md','ONE_POINT_FEEDBACK_AUDIT.md')]
    inputs += [FROZEN/'inputs'/name for name in ('FULL_SECOND_RESPONSE_CONSTRUCTION.md',
        'FORCED_INCREMENT_BOUND.md','ACTUAL_SOURCES_NORMAL.json','REFINED_RESPONSE_NORMAL.json')]
    inputs += [OLD/name for name in ('projection/inputs/linear_operator.py',
        'projection/angular_generic.py','projection/source_coefficients.py',
        'projection/GEOMETRIC_COEFFICIENT_MAPS.json','numeric/actual_first_jets.py',
        'numeric/shell_source_values.py','numeric/fast_dyadic.py',
        'numeric/interval_core.py','GAUGE_COVECTOR_BOOKKEEPING.md')]
    pins = [input_pin(p) for p in inputs]
    known = {
        'FULL_SECOND_RESPONSE_CONSTRUCTION.md':'5de45b21236cfe8c6d99fb403bacc4624069c8d96d6cf0ad2ec89b8e2d3617e1',
        'FORCED_INCREMENT_BOUND.md':'317262b4581037f00e6f65cba6771ed49a0bc8c89bf5d6b5bf7759b758eb6d7f',
        'GEOMETRIC_COEFFICIENT_MAPS.json':'34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
    }
    for row in pins:
        name=Path(row['path']).name
        if name in known:
            need(row['sha256']==known[name], 'retained known input pin: '+name)

    a,b=Q(3,2),Q(1,4)
    d,c=a-b,a+b
    # Cartesian vertices directly from |u-b|<=a-s, independently of any mesh.
    plus=[(Q(0),Q(0)),(Q(0),c),(a,b),(d,Q(0))]
    minus=[(Q(0),Q(0)),(Q(0),d),(d,Q(0))]
    exclusive=[(Q(0),d),(Q(0),c),(a,b),(d,Q(0))]
    zero=[(Q(0),Q(0)),(Q(0),Q(1)),(Q(1),Q(0))]
    Aplus,Aminus,Aexclusive,Azero=map(area,(plus,minus,exclusive,zero))
    need((Aplus,Aminus,Aexclusive,Azero)==(Q(47,32),Q(25,32),Q(11,16),Q(1,2)), 'independent polygon areas')
    need(Aplus-Aminus==Aexclusive and Aplus+Aminus==a*a,'fold partition')
    active_weight=(Aplus-Azero+Aminus-Azero)/2
    need(active_weight==Q(5,8),'both zero triangles removed')

    actual=json.loads((FROZEN/'inputs/ACTUAL_SOURCES_NORMAL.json').read_text())
    refined=json.loads((FROZEN/'inputs/REFINED_RESPONSE_NORMAL.json').read_text())
    FH=Q(actual['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    GH=Q(refined['retained_join']['G_H'])
    need((FH,GH)==(17387513156573,3851400),'retained full-state source caps')
    sigma=Q(1,10000)+Q(1,10000)+3*Q(1,1000000)
    m=Q(1,2000)
    KX=max(a,a*a/6)
    LX=max(Q(1),d/3)
    q=KX*sigma+LX*3*m
    S=KX*FH+LX*GH
    BD=S/(1-q)
    Ecoarse=(a*a*sigma/2+3*d*m)*BD
    Ezero=(active_weight*sigma+3*(d-1)*m)*BD
    Etrace=(active_weight*sigma+m*(d-1)**2/2)*BD
    need(q==Q(3609,2000000) and q<1,'restricted contraction')
    need(S==Q(52162547172519,2),'restricted source bound')
    need(BD==Q(52162547172519000000,1996391),'restricted resolvent')
    need(Ecoarse==Q(877739181271977213,15971128),'unrestricted point allowance')
    need(Ezero==Q(209432626897663785,15971128),'zero-wedge allowance')
    need(Etrace==Q(14866325944167915,3992782),'trace-improved allowance')
    need(Etrace<Ezero<Ecoarse,'strict improvement chain')

    spec=importlib.util.spec_from_file_location('readonly_operator',OLD/'projection/inputs/linear_operator.py')
    operator=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(operator)
    mats=operator.matrices()
    pd=[operator.LABELS[i] for i in range(23) if mats['PD'][i][i]]
    pn=[operator.LABELS[i] for i in range(23) if mats['PN'][i][i]]
    need(pd==['L:0','L:2','L:4','v:2','v:4'] and len(pn)==18,'five selected Dirichlet and eighteen Robin amplitudes')
    need('K:0' in pn and operator.gram()[operator.IDX['K:0']]==1,'mean yy H-dual normalization')

    # These unequal exact witnesses identify each false formula specifically.
    controls=[]
    def reject(name, correct, wrong):
        need(correct!=wrong,'negative control indistinguishable: '+name)
        controls.append({'name':name,'expected':str(correct),'false_formula':str(wrong),'rejected':True})
    robin=(Aplus+Aminus)/2
    dirichlet=(Aplus-Aminus)/2
    reject('omit_reflected_Robin_cone',robin,Aplus/2)
    reject('extend_Dirichlet_components_evenly',dirichlet,robin)
    reject('drop_wave_kernel_half',robin,Aplus+Aminus)
    reject('remove_zero_wedge_from_only_one_cone',active_weight,(Aplus+Aminus-Azero)/2)
    # g(s)=s has compatible g(0)=0; Lg=-d^2/2, radial derivative +d.
    reject('wrong_boundary_lift_sign',-d*d/2,d*d/2)
    reject('wrong_boundary_radial_derivative_sign',d,-d)
    q00,qab,qcy=Q(2),Q(3),Q(5)
    correct_trace=-q00+3*qab+qcy
    reject('Euclidean_time_trace',correct_trace,q00+3*qab+qcy)
    reject('four_sphere_dimensions',correct_trace,-q00+4*qab+qcy)
    reject('omit_QCy',correct_trace,-q00+3*qab)
    reject('two_copy_QK_reused_as_one_copy',correct_trace,-2*q00+6*qab+qcy)
    ell,rho=Q(1,1000000),Q(130)
    gn=-2*(-q00+3*qab)-2*qcy
    got=operator.boundary_forcing(ell,rho,{'GN':(gn,0,0),
        'QK00':(q00,0,0),'QK0a':(0,0),'QKab':(qab,0,0,0,0,0,0),'Qp':(0,0,0)})
    need(got[operator.IDX['K:0']]==2*ell*rho*correct_trace,'retained boundary adapter agrees with independent normal projection')
    reject('multiply_value_allowance_by_three',Etrace,3*Etrace)
    need(pins==[input_pin(p) for p in inputs],'retained inputs unchanged during audit')
    return {'status':'PASS_INDEPENDENT_CONDITIONAL_ONE_POINT_ANALYTIC_AUDIT',
        'scope':'Exact geometry, rational constants and synthetic convention witnesses; no actual source integral evaluated; accepted smooth first/second construction remains a premise.',
        'input_pins':pins,
        'geometry':{'a':str(a),'b':str(b),'d':str(d),'direct_area':str(Aplus),'reflected_area':str(Aminus),'exclusive_area':str(Aexclusive),'zero_triangle_area_each':str(Azero),'active_direct_area':str(Aplus-Azero),'active_reflected_area':str(Aminus-Azero),'active_Robin_kernel_weight':str(active_weight)},
        'full_state':{'selected_Dirichlet_amplitudes':pd,'selected_Robin_count':len(pn),'all_bulk_and_Robin_couplings_enclosed':True,'mean_H_dual_norm':'1'},
        'bounds':{key:str(value) for key,value in {'FH':FH,'GH':GH,'sigma':sigma,'m':m,'K_X':KX,'L_X':LX,'q_D':q,'S_D_cap':S,'B_D':BD,'E_coarse':Ecoarse,'E_zero':Ezero,'E_trace':Etrace}.items()},
        'synthetic_false_formula_controls':controls,
        'not_claimed':['new source evaluation','sign separation','standalone physical constraint solution for W','full numerical waveform','independent replay of inherited PDE/source certificates'],
        'checker_sha256':pin(Path(__file__))['sha256']}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--frozen-root',type=Path,default=FROZEN)
    parser.add_argument('--follow-root',type=Path,default=FOLLOW)
    parser.add_argument('--old-root',type=Path,default=OLD)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    need(not args.output.exists(),'refuse output overwrite')
    result=run(args.frozen_root,args.follow_root,args.old_root)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
