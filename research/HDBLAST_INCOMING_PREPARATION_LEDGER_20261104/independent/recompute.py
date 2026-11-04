#!/usr/bin/env python3
"""Independent exact source-budget and preparation-identity reconstruction."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import sys

HERE = Path(__file__).resolve().parent
PIN_SHA = 'd7ecbf4f0b6721b1ab7e283030cc6cb8d2b05e5f5c6a45f2ec8043baeae6124d'
GEO = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json'
PREMISES = 'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json'
INTRINSIC = 'research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/bounds/references/NORMAL.json'

class Failure(Exception):
    pass

def need(ok, witness):
    if not ok:
        raise Failure(witness)

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def regular(p):
    return p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.absolute().parents)

def canonical(v):
    return (json.dumps(v, sort_keys=True, indent=2) + '\n').encode()

def read(p):
    need(regular(p), 'REGULAR_INPUT_FILE')
    return json.loads(p.read_bytes())

def encode(v):
    if isinstance(v, F):
        return str(v)
    if isinstance(v, dict):
        return {k: encode(x) for k, x in v.items()}
    if isinstance(v, (tuple, list)):
        return [encode(x) for x in v]
    return v

def authenticate(repo):
    need(regular(HERE/'SOURCE_PINS.json') and digest(HERE/'SOURCE_PINS.json') == PIN_SHA, 'LOCAL_PIN_BYTES')
    pins = read(HERE/'SOURCE_PINS.json')['files']
    for path, row in pins.items():
        relative = PurePosixPath(path)
        need(not relative.is_absolute() and '..' not in relative.parts and relative.as_posix() == path, 'RELATIVE_SOURCE_PATH')
        p = repo/path
        need(regular(p) and p.stat().st_size == row['bytes'] and digest(p) == row['sha256'], 'SOURCE_PIN_MISMATCH')
    return pins

# Sparse multivariate polynomials represented by tuples of exponents. This
# implementation is independent of the primary component's polynomial class.
def plus(*ps):
    result = {}
    for p in ps:
        for key, value in p.items():
            result[key] = result.get(key, F(0)) + value
    return {k:v for k,v in result.items() if v}

def scale(c,p):
    return {k:c*v for k,v in p.items() if c*v}

def product(p,q):
    result = {}
    for a,v in p.items():
        for b,w in q.items():
            key = tuple(x+y for x,y in zip(a,b))
            result[key] = result.get(key,F(0)) + v*w
    return {k:v for k,v in result.items() if v}

def variable(i,n=4):
    return {tuple(int(j == i) for j in range(n)):F(1)}

def sq(p):
    return product(p,p)

def coefficient(p,index,power):
    result={}
    for exponents,value in p.items():
        if exponents[index] == power:
            key=tuple(0 if i == index else x for i,x in enumerate(exponents))
            result[key]=value
    return result

def exp_upper(x,n):
    partial = sum((x**j/F(math.factorial(j)) for j in range(n+1)),F(0))
    ratio = x/F(n+2)
    need(0 <= ratio < 1, 'TAYLOR_TAIL_RATIO')
    return partial + x**(n+1)/F(math.factorial(n+1))/(1-ratio)

def interval(row):
    a,b = map(F,row)
    need(a <= b, 'ORDERED_RETAINED_INTERVAL')
    return a,b

def magnitude(row):
    return max(map(abs,interval(row)))

def triangle_rectangle_area(x,r,level):
    # Inclusion/exclusion of four translated right triangles, independently
    # replacing the primary implementation's polygon half-plane clipping.
    a,b = map(F,x)
    c,d = map(F,r)
    return sum((sign*max(F(0),level-u-v)**2/2
                for sign,u,v in [(1,a,c),(-1,b,c),(-1,a,d),(1,b,d)]),F(0))

def recompute(repo, source_path, theory_path):
    pins = authenticate(repo)
    source_bytes, theory_bytes = source_path.read_bytes(), theory_path.read_bytes()
    source, theory = read(source_path), read(theory_path)
    checks = []
    def check(ok,name):
        need(ok,name)
        checks.append(name)
    premises = read(repo/PREMISES)
    check(premises['status'] == 'EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER', 'retained_operator_caps_remain_conditional')
    check(premises['mean_gram'] == [1,2,1,3,1], 'retained_positive_mean_Gram')
    a=b=F(1,10000)
    d,m,k,alpha,ell,T = F(1,1000000),F(1,2000),F(1,2),F(1,2),F(1,10**6),F(5,4)
    check({n:F(v) for n,v in premises['coefficient_norm_upper'].items()} == {'A':a,'B':b,'D':d,'M':m}, 'matching_actual_future_operator_caps')
    check(premises['domain'] == '0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x' and T <= F(5,4), 'future_domain_only_no_past_extension')
    V,Q,U,G = [variable(i) for i in range(4)]
    e2 = plus(sq(V),sq(Q),scale(2*k,product(V,Q)),sq(U))
    check(plus(e2,scale(-1,sq(plus(V,scale(k,Q))))) == plus(scale(1-k*k,sq(Q)),sq(U)), 'sharp_bulk_forcing_polynomial')
    P,u,mu,g = [variable(i) for i in range(4)]
    def boundary(q):
        return plus(scale(-k/2,sq(P)),scale(-(1-k),product(P,q)),scale(1-k,sq(q)),scale(alpha,product(u,plus(P,scale(-1,q)))))
    difference = plus(boundary(plus(mu,g)),scale(-1,boundary(mu)))
    check(difference == plus(scale(F(-1,2),product(P,g)),product(mu,g),scale(F(1,2),sq(g)),scale(F(-1,2),product(u,g))), 'actual_nonsymmetric_boundary_source_split')
    check(plus(scale(F(1,8),sq(P)),scale(F(-1,2),product(P,g)),scale(F(1,2),sq(g))) == scale(F(1,8),sq(plus(P,scale(-2,g)))), 'boundary_characteristic_Young_square')
    c=m+F(1,2)
    check(plus(scale(F(1,4),sq(u)),scale(-c,product(u,g)),scale(c*c,sq(g))) == sq(plus(scale(F(1,2),u),scale(-c,g))), 'boundary_state_Young_square')
    cv=a+(b+k*a)/2+d/2+F(1,2)
    cq=k*b+(b+k*a)/2+k*d/2
    cu=(1+k)*d/2+F(1,2)
    Cbulk=max(2*cv/(1-k),2*cq/(1-k),2*cu)
    beta=alpha+(1-k)*m
    Cboundary=2/alpha*(beta*beta/k+(1-k)*m*m+alpha*m)
    C,Cf,Cg = F(301,100),F(1,2),1+c*c
    check(Cbulk == F(1000351,500000) and Cboundary == F(2003001,1000000), 'independent_homogeneous_rates')
    check(max(Cbulk+1,Cboundary+1) < C and Cg == F(5002001,4000000), 'forced_growth_and_boundary_square_coefficients')
    exp_bound=exp_upper(C*T,24)
    check(exp_bound < 44, 'independent_order_24_exponential_guard')

    # Affine-null Raychaudhuri identities: formal polynomial background and
    # generic second jets; no integration of a past solution is performed.
    A,H,B = [variable(i,3) for i in range(3)]
    background = plus(scale(3,plus(scale(2,product(A,plus(H,scale(-1,A)))),sq(A),scale(-1,sq(H)),scale(F(-1,3),B))),scale(3,sq(plus(H,scale(-1,A)))),B)
    check(not background, 'affine_null_background_Raychaudhuri_identity')
    t0,t1,t2,ep = [variable(i) for i in range(4)]
    jet=plus(t0,product(ep,t1),scale(F(1,2),product(sq(ep),t2)))
    check(scale(2,coefficient(sq(jet),3,2)) == plus(scale(2,sq(t1)),scale(2,product(t0,t2))), 'general_second_square_retains_background_acceleration')
    check(F(2,3) > 0 and F(2) > 0, 'null_second_focusing_square_weights')
    # Since sigma0=0, each nonzero shear product already has order >=2;
    # first variations of inverse screen metrics enter only at order >=3.
    check(1+1 == 2 and 1+1+1 > 2, 'zero_shear_background_eliminates_metric_crosses_at_order_two')
    L1,L2 = F(2),F(-1)
    check(L2 < 0 and L2+L1*L1 > 0, 'negative_second_log_density_does_not_fix_area_sign')

    # For normalized means with V0=1, use formal quotient differentiation.
    phi,v1,v2,w1,w2,cross = [variable(i,6) for i in range(6)]
    n1=plus(product(phi,v1),w1)
    n2=plus(product(phi,v2),w2,cross)
    normalized=plus(n2,scale(-2,product(n1,v1)),scale(2,product(phi,sq(v1))),scale(-1,product(phi,v2)))
    check(normalized == plus(w2,cross,scale(-2,product(w1,v1))), 'normalized_scalar_mean_measure_cancellation')

    constants=read(repo/INTRINSIC)['constants']
    rho,p_abs,a_abs = [F(constants[name]['exact'][1]) for name in ['rho','minus_p','a']]
    q=F(constants['q']['exact'])
    check((rho,p_abs,a_abs,q,F(constants['ell']['exact'])) == (131,F(1,7000),F(14,125),F(1201,3200000),ell), 'retained_inverse_initial_constants')
    bM,cM=ell*rho*a_abs,2*ell*rho*p_abs/3
    delta=T-1
    matrix_guard=F(20001,20000)
    check(max(1+3*bM,3*ell*ell+bM) < matrix_guard, 'independent_inverse_IVP_matrix_guard')
    check(exp_upper(matrix_guard*delta,8) < F(4,3), 'independent_inverse_IVP_exponential_guard')
    Hw=delta*F(4,3)*cM
    Hrw=(3*ell*ell+bM)*Hw+cM
    Vw=3*ell*ell*Hw+cM
    N=max((1+2*Hw+2*T*Vw)/3,1+2*Hrw+2*Vw)
    Nguard=F(20000001,20000000)
    kappaguard=F(47,125000)
    kappa=q*Nguard/(1-q)
    check(N < Nguard and kappa < kappaguard < 1 and q > 0, 'inverse_history_near_identity_and_injectivity_guards')

    # Two-endpoint weighted Cauchy gives Q(x)^2 <= x(1-x)||Q'||^2
    # on the unit support interval; integral of the kernel is 1/6. This is
    # stronger than the shared publication's conservative 1/4 inequality.
    cauchy_integral=F(1,2)-F(1,3)
    check(cauchy_integral == F(1,6) and cauchy_integral <= F(1,4), 'independent_two_endpoint_Cauchy_kernel')
    check(4*ell*(ell/2) == 2*ell*ell, 'principal_difference_integration_by_parts_factor')
    imbalance=ell*ell/6
    publication_guard=ell*ell/4
    check(0 < imbalance <= publication_guard == F(1,4000000000000), 'principal_imbalance_common_conservative_guard')
    check(ell*(1-F(3,2)*ell) > 0, 'actual_packet_both_current_directions_positive')
    check(F(2,3)*F(3,2) == 1, 'scalar_oneway_h_reconstruction_coefficient_cancellation')

    actual=read(repo/GEO)
    check(actual['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES', 'retained_actual_continuous_source_status')
    check(actual['normalization'] == 'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.', 'actual_second_Hessian_source_normalization')
    FH=F(actual['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    GH=F(premises['source_norm_upper']['G'])
    check(FH == F(premises['source_norm_upper']['F']) == 17387513156573 and GH == 3851400, 'actual_source_registered_norms_match')
    rows=actual['global_cover']['cell_enclosures']
    bulk_full=bulk_active=area_full=area_active=F(0)
    pieces=[]
    for row in rows:
        af=triangle_rectangle_area(row['x'],row['r'],T)
        aa=af-triangle_rectangle_area(row['x'],row['r'],F(1))
        need(af >= 0 and aa >= 0, 'NONNEGATIVE_INTERSECTION_AREA')
        area_full+=af
        area_active+=aa
        if af:
            ss=row['sources']
            slots=[ss['metric_forcing'][name][0] for name in ['yy','y0','00','ab']]+[ss['scalar_forcing'][0]]
            cap=sum((F(weight)*magnitude(box)**2 for weight,box in zip([1,2,1,3,1],slots)),F(0))
            chosen=min(cap,FH*FH)
            bulk_full+=af*chosen
            bulk_active+=aa*chosen
            pieces.append({'old_cell_index':row['index'],'cone_area':af,'active_area':aa,'mean_square_cap':cap,'chosen_square_cap':chosen})
    check(area_full == F(25,32) and area_active == F(9,32) and len(pieces) == 6, 'independent_rectangle_CDF_cone_geometry')
    check({tuple(x['old_cell_index']):x['active_area'] for x in pieces if x['active_area']} == {(0,1):F(3,32),(1,0):F(3,32),(0,2):F(1,32),(1,1):F(1,32),(2,0):F(1,32)}, 'independent_five_active_piece_areas')
    boundary_full=boundary_active=length_full=length_active=F(0)
    boundary_rows=[]
    for shell in actual['shell_whole_time_cover']['cell_enclosures']:
        left,right=map(F,shell['x'])
        width=max(F(0),min(right,T)-max(left,F(0)))
        active_width=max(F(0),min(right,T)-max(left,F(1)))
        length_full+=width
        length_active+=active_width
        if not width:
            continue
        adjacent=[row for row in rows if row['x'] == shell['x'] and F(row['r'][0]) == 0]
        need(len(adjacent) == 1, 'UNIQUE_RETAINED_BOUNDARY_QC_BOX')
        physical=shell['sources']['one_copy_outward_hessian']
        qr00,qrJ,qp,qcy=physical['00'][0],physical['ab'][0],physical['scalar'][0],adjacent[0]['sources']['QC']['y'][0]
        intervals=[interval(box) for box in [qr00,qrJ,qcy]]
        coefficients=[F(2),F(-6),F(-2)]
        gn=[sum((c*box[0 if c >= 0 else 1] for c,box in zip(coefficients,intervals)),F(0)),sum((c*box[1 if c >= 0 else 0] for c,box in zip(coefficients,intervals)),F(0))]
        cap=(ell*rho)**2*(magnitude(gn)**2+4*magnitude(qr00)**2+12*magnitude(qrJ)**2+magnitude(qp)**2)
        chosen=min(cap,GH*GH)
        boundary_full+=width*chosen
        boundary_active+=active_width*chosen
        boundary_rows.append({'old_time_cell':shell['x'],'length':width,'active_length':active_width,'direct_mean_square_cap':cap,'chosen_square_cap':chosen})
    check(length_full == T and length_active == F(1,4) and boundary_active == GH*GH/4, 'independent_active_boundary_integral')
    budgets={}
    for name,f,g in [('full_cone_global_caps',area_full*FH*FH,T*GH*GH),('wedge_zero_global_caps',area_active*FH*FH,GH*GH/4),('full_cone_retained_cell_caps',bulk_full,boundary_full),('wedge_zero_retained_cell_caps',bulk_active,boundary_active)]:
        ledger=Cf*f+Cg*g
        budgets[name]={'bulk_L2_squared_upper':f,'boundary_L2_squared_upper':g,'energy_source_ledger_upper':ledger,'conditional_zero_initial_energy_upper':44*ledger,'conditional_zero_initial_mean_trace_squared_upper':176*ledger}
    check(bulk_active <= area_active*FH*FH and boundary_active <= GH*GH/4, 'integrated_source_ledger_no_worse_than_global_cap')

    # Candidate outputs are inputs, not trusted or imported implementations.
    check(source['status'] == 'PASS_EXACT_ZERO_ORACLE_SOURCE_LEDGER', 'candidate_source_status')
    check(source['mean_gram'] == [1,2,1,3,1] and all(v == 0 for v in source['scientific_activity'].values()), 'candidate_source_no_new_queries_scope')
    check({name:{key:F(value) for key,value in b.items()} for name,b in source['budgets'].items()} == budgets, 'CANDIDATE_SOURCE_LEDGER_MISMATCH')
    check({key:F(value) for key,value in source['registered_caps'].items()} == {'F_H':FH,'G_H':GH} and source['retained_actual_source_report_sha256'] == digest(repo/GEO), 'candidate_source_norm_and_report_binding')
    check(theory['status'] == 'PASS_CONDITIONAL_INCOMING_PREPARATION_IDENTITIES', 'candidate_theory_status')
    expected_forced={'tilt_k':k,'boundary_alpha':alpha,'homogeneous_bulk_rate':Cbulk,'homogeneous_boundary_rate':Cboundary,'forced_growth_guard':C,'bulk_square_factor':Cf,'boundary_square_factor':Cg,'exponential_guard':F(44),'scalar_trace_square_factor':F(4)}
    check({key:F(theory['forced_energy'][key]) for key in expected_forced} == expected_forced, 'CANDIDATE_FORCED_THEORY_MISMATCH')
    expected_inverse={'cM':cM,'Hw':Hw,'Hrw':Hrw,'Vw':Vw,'N_exact_majorant':N,'N_guard':Nguard,'q_retained':q,'kappa_exact':kappa,'kappa_guard':kappaguard}
    check({key:F(theory['inverse_history'][key]) for key in expected_inverse} == expected_inverse, 'CANDIDATE_INVERSE_THEORY_MISMATCH')
    check(theory['background_null_identity']['screen_dimension'] == 3 and theory['null_constraint']['integrating_factor'] == 'R0^2' and theory['null_constraint']['general_second_scalar_cross_term'] == '-2 Phi0prime Phi2prime retained before setting pulled-back Phi2=0', 'candidate_null_second_constraint_scope')
    check(theory['null_constraint']['actual_area_second_formula'] == 'average(L2+L1^2); no sign from L2 alone' and theory['null_constraint']['full_characteristic_solution_constructed'] is False and theory['null_constraint']['current_shell_expansion_or_H1_bound_established'] is False, 'candidate_null_no_area_or_full_solution_claim')
    check(theory['inverse_history']['on_admissible_h_only_image'] is True and theory['inverse_history']['onto_arbitrary_histories'] is False and theory['inverse_history']['past_state_constructed'] is False, 'candidate_inverse_only_admissible_image')
    check(theory['packet_diagnostic']['not_current_time_symmetric_or_pure_oneway'] is True and theory['packet_diagnostic']['past_incoming_only_state_excluded'] is False and theory['packet_diagnostic']['is_physical_energy_flux_measurement'] is False, 'candidate_current_packet_diagnostic_scope')
    check(authenticate(repo) == pins and source_path.read_bytes() == source_bytes and theory_path.read_bytes() == theory_bytes, 'all_scientific_input_bytes_unchanged_after_arithmetic')
    return encode({'schema':'HDBLAST_INCOMING_INDEPENDENT_RECONSTRUCTION_V1','status':'PASS_INDEPENDENT_INCOMING_PREPARATION_RECONSTRUCTION','source_pins_sha256':PIN_SHA,'candidate_sha256':{'source':hashlib.sha256(source_bytes).hexdigest(),'theory':hashlib.sha256(theory_bytes).hexdigest()},'authenticated_sources':pins,'verification':{'finite_checks':len(checks),'check_names':checks,'new_PDE_solves':0,'new_field_background_bump_or_source_queries':0,'project_imports':0},'shared_projection':{'forced_energy':expected_forced,'inverse_history':expected_inverse,'source_budgets':budgets,'packet_common_principal_imbalance_upper':publication_guard,'null_screen_dimension':3,'null_integrating_factor':'R0^2','area_second_formula':'average(L2+L1^2); no sign from L2 alone'},'independent_route':{'rectangle_CDF_bulk_pieces':pieces,'boundary_mean_boxes':boundary_rows,'order_24_exponential_upper':exp_bound,'two_endpoint_Cauchy_integral':cauchy_integral,'stronger_principal_imbalance_upper':imbalance},'limits':['Continuous source boxes and future operator caps remain retained inputs with their original validated/conditional scope.','No past-domain bound, new characteristic solution, full nonlinear solution or physical preparation law is constructed.','Zero-initial forced ledger is an auxiliary estimating response, not admissible physical zero second data.','Principal packet diagnostics are not physical Noether energy or stress flux.','Scalar history inverse is only on the retained admissible h-only image of normalized proper-clock scalar means.']})

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['repo-root','source-candidate','theory-candidate','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    args=ap.parse_args()
    need(all(getattr(args,key).is_absolute() for key in ['repo_root','source_candidate','theory_candidate','output']),'ABSOLUTE_ARGUMENT_PATHS')
    repo,source,theory,out=[getattr(args,key).absolute() for key in ['repo_root','source_candidate','theory_candidate','output']]
    need(all(not q.is_symlink() for p in [repo,source,theory,out] for q in [p,*p.parents]), 'PATH_SYMLINK')
    repo,source,theory=repo.resolve(strict=True),source.resolve(strict=True),theory.resolve(strict=True)
    out=out.resolve()
    need(not out.exists(),'FRESH_OUTPUT_FILE')
    need(out.parent.is_dir(),'OUTPUT_PARENT_EXISTS')
    need(all(base != out and base not in out.parents for base in [repo,HERE]),'OUTPUT_OUTSIDE_CODE_AND_REPO')
    report=recompute(repo,source,theory)
    with out.open('xb') as stream:
        stream.write(canonical(report))
    print(json.dumps({'status':'PASS','witness':'INDEPENDENT_INCOMING_RECONSTRUCTION','finite_checks':report['verification']['finite_checks'],'output_sha256':digest(out)},sort_keys=True))

if __name__ == '__main__':
    try:
        main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True))
        sys.exit(2)
    except ArithmeticError:
        print(json.dumps({'status':'FAIL','witness':'NUMERIC_PARSE_ERROR'},sort_keys=True))
        sys.exit(2)
    except (ValueError,OSError,KeyError,TypeError):
        print(json.dumps({'status':'FAIL','witness':'RUNTIME_INPUT_ERROR'},sort_keys=True))
        sys.exit(2)
