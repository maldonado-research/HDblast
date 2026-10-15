#!/usr/bin/env python3
"""Independent rational 5D metric-jet checks; imports no producer science code."""
import argparse
import csv
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import random

N = 5
I = range(N)

def arr(shape, value=F(0)):
    return [arr(shape[1:], value) for _ in range(shape[0])] if shape else value

def mm(a, b):
    return [[sum(a[i][k]*b[k][j] for k in I) for j in I] for i in I]

def neg(a):
    return [[-x for x in row] for row in a]

def eq(a, b, label):
    if a != b:
        raise ValueError(label + ': nonzero rational residual ' + str(a-b))

def jets(seed, injection):
    rng = random.Random(seed)
    rnd = lambda: F(rng.randint(-7, 7), rng.randint(1, 6))
    g = [[F(0) for _ in I] for _ in I]
    for i in I:
        g[i][i] = F((-1 if i == 0 else 1)*(i+1))
    inv0 = [[(1/g[i][i] if i == j else F(0)) for j in I] for i in I]
    h, dg, dh, ddg, ddh = arr([N,N]), arr([N,N,N]), arr([N,N,N]), arr([N,N,N,N]), arr([N,N,N,N])
    for a in I:
        for b in range(a, N):
            h[a][b] = h[b][a] = rnd()
            for k in I:
                dg[k][a][b] = dg[k][b][a] = rnd()
                dh[k][a][b] = dh[k][b][a] = rnd()
                for l in range(k, N):
                    for target in (ddg, ddh):
                        v = rnd()
                        target[k][l][a][b] = target[l][k][a][b] = v
                        target[k][l][b][a] = target[l][k][b][a] = v
    inv = [inv0, neg(mm(mm(inv0,h),inv0)), mm(mm(mm(mm(inv0,h),inv0),h),inv0)]
    gd = [dg, dh, arr([N,N,N])]
    gdd = [ddg, ddh, arr([N,N,N,N])]
    dinv = arr([3,N,N,N])
    for k,l,a,b in itertools.product(range(3), I, I, I):
        dinv[k][l][a][b] = -sum(inv[r][a][c]*gd[s][l][c][d]*inv[k-r-s][d][b]
            for r in range(k+1) for s in range(k-r+1) for c in I for d in I)
    gamma, dgamma = arr([3,N,N,N]), arr([3,N,N,N,N])
    for k,a,b,c in itertools.product(range(3),I,I,I):
        gamma[k][a][b][c] = sum(inv[r][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
            for r in range(k+1) for d in I)/2
        for l in I:
            dgamma[k][l][a][b][c] = sum(
                dinv[r][l][a][d]*(gd[k-r][b][c][d]+gd[k-r][c][b][d]-gd[k-r][d][b][c])
                + inv[r][a][d]*(gdd[k-r][l][b][c][d]+gdd[k-r][l][c][b][d]-gdd[k-r][l][d][b][c])
                for r in range(k+1) for d in I)/2
    # Independent covariant A formula, using the background connection.
    nabh = arr([N,N,N])
    for a,b,c in itertools.product(I,I,I):
        nabh[a][b][c] = dh[a][b][c]-sum(gamma[0][d][a][b]*h[d][c]+gamma[0][d][a][c]*h[b][d] for d in I)
    A = [[[sum(inv0[a][d]*(nabh[b][c][d]+nabh[c][b][d]-nabh[d][b][c]) for d in I)/2 for c in I] for b in I] for a in I]
    H = mm(inv0,h)
    B = [[[-2*sum(H[a][d]*A[d][b][c] for d in I) for c in I] for b in I] for a in I]
    dH = [mm(dinv[0][l],h) for l in I]
    for l,a,b in itertools.product(I,I,I):
        dH[l][a][b] += sum(inv0[a][d]*dh[l][d][b] for d in I)
    dB = [[[[-2*sum(dH[l][a][d]*A[d][b][c]+H[a][d]*dgamma[1][l][d][b][c] for d in I) for c in I] for b in I] for a in I] for l in I]
    nabB = arr([N,N,N,N])
    checks = 0
    for l,a,b,c in itertools.product(I,I,I,I):
        nabB[l][a][b][c] = dB[l][a][b][c]+sum(gamma[0][a][l][d]*B[d][b][c]-gamma[0][d][l][b]*B[a][d][c]-gamma[0][d][l][c]*B[a][b][d] for d in I)
    for a,b,c in itertools.product(I,I,I):
        eq(A[a][b][c],gamma[1][a][b][c],'first connection')
        eq(B[a][b][c],2*gamma[2][a][b][c],'second connection')
        checks += 2
    p, chi1 = [rnd() for _ in I], [rnd() for _ in I]
    phi2, chi2 = arr([N,N]), arr([N,N])
    for a in I:
        for b in range(a,N):
            phi2[a][b] = phi2[b][a] = rnd()
            chi2[a][b] = chi2[b][a] = rnd()
    chi, Up, Upp, Uppp = rnd(), rnd(), rnd(), rnd()
    for a in I:
        for b in range(a,N):
            directR = 2*(sum(dgamma[2][c][c][a][b]-dgamma[2][b][c][a][c] for c in I)
                +sum(gamma[r][c][c][d]*gamma[2-r][d][a][b]-gamma[r][c][b][d]*gamma[2-r][d][a][c] for r in range(3) for c in I for d in I))
            sourceR = sum(nabB[c][c][a][b]-nabB[b][c][a][c] for c in I)+2*sum(A[c][c][d]*A[d][a][b]-A[c][b][d]*A[d][a][c] for c in I for d in I)
            if injection == 'half_ricci':
                sourceR /= 2
            matter = -2*chi1[a]*chi1[b]-F(2,3)*(Upp*chi**2*g[a][b]+2*Up*chi*h[a][b])
            candidate = sourceR+matter
            if injection == 'drop_potential_cross':
                candidate += F(4,3)*Up*chi*h[a][b]
            eq(candidate,directR+matter,'Einstein Hessian '+str((a,b)))
            checks += 1
    nabphi = [[phi2[a][b]-sum(gamma[0][c][a][b]*p[c] for c in I) for b in I] for a in I]
    nabchi = [[chi2[a][b]-sum(gamma[0][c][a][b]*chi1[c] for c in I) for b in I] for a in I]
    hUU = neg(inv[1])
    sourceS = sum(2*inv[2][a][b]*nabphi[a][b]-2*hUU[a][b]*nabchi[a][b]
        +sum(2*hUU[a][b]*A[c][a][b]*p[c]-inv0[a][b]*B[c][a][b]*p[c]-2*inv0[a][b]*A[c][a][b]*chi1[c] for c in I) for a in I for b in I)-Uppp*chi**2
    directS = 2*sum(inv[2][a][b]*phi2[a][b]+inv[1][a][b]*chi2[a][b]
        -sum(inv[r][a][b]*gamma[2-r][c][a][b]*p[c] for r in range(3) for c in I)
        -sum(inv[r][a][b]*gamma[1-r][c][a][b]*chi1[c] for r in range(2) for c in I) for a in I for b in I)-Uppp*chi**2
    eq(sourceS,directS,'scalar Hessian')
    checks += 1
    C1 = [sum(inv0[b][c]*A[a][b][c] for b in I for c in I) for a in I]
    QC = [sum(inv0[b][c]*B[a][b][c]-2*hUU[b][c]*A[a][b][c] for b in I for c in I) for a in I]
    for a in I:
        cand = sum(g[a][b]*QC[b]+2*h[a][b]*C1[b] for b in I)
        if injection == 'drop_lower_gauge_cross':
            cand -= 2*sum(h[a][b]*C1[b] for b in I)
        exact = 2*sum(g[a][d]*(inv0[b][c]*gamma[2][d][b][c]+inv[1][b][c]*gamma[1][d][b][c])+h[a][d]*inv0[b][c]*gamma[1][d][b][c] for d in I for b in I for c in I)
        eq(cand,exact,'lower wave-map Hessian')
        checks += 1
    return checks

def run(args):
    named=[]
    mark=lambda name: named.append(name)
    pins=json.loads((args.root/'inputs/INPUT_PINS.json').read_text())
    for item in pins['files']:
        raw=(args.root/item['path']).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=item['sha256'] or len(raw)!=item['bytes']:
            raise ValueError('portable input authentication failed: '+item['path'])
    mark('authenticate_portable_input_inventory')
    report=json.loads((args.report or args.root/'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json').read_text())
    scope=report['scope']
    if (report['status']!='PASS' or scope['full_U2_evaluated'] is not False
        or scope['nonlinear_remainder_established'] is not False
        or scope['coordinate_source_sign_is_gauge_invariant'] is not False
        or scope['EN_scope']!='Retained degree-four reduced Robin forcing only.'):
        raise ValueError('coordinate report scientific scope failed')
    if any(report['exact']['norms'][key]['angular_field']!='H4' for key in ('EN','dEN','EN_material')):
        raise ValueError('degree-four Robin norm projection scope failed')
    mark('coordinate_source_scope_and_degree_four_only_Robin_norms')
    boxpath=args.root/'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'
    box=json.loads(boxpath.read_text())['boxes']
    get=lambda k: tuple(F(v) for v in box[k]['exact'])
    a,s,p=get('a_b'),get('sigma_second_b'),get('p_b')
    pm=max(map(abs,p))
    detlo=24*a[0]*s[0]-16*pm**2
    lambdalo=min(6*a[0]-4*pm,4*s[0]-4*pm)
    if args.inject_error=='wrong_robin_sign': detlo=-detlo
    if not detlo>0 or not lambdalo>0: raise ValueError('actual symmetrized Robin matrix positivity failed')
    mark('actual_shell_Robin_symmetrizer_positive_definite')
    count=sum(jets(seed,args.inject_error) for seed in (271,1927))
    named.extend(['generic_5D_first_connection','generic_5D_second_connection_factorial',
                  'all_fifteen_Einstein_Hessian_components','scalar_Hessian_all_terms',
                  'lower_wave_map_Hessian_cross_term'])
    # Boundary normal-normalization and eliminated first shell jets.
    for aa,pp,ss,u,v in [(F(1,7),F(-2,11),F(13,5),F(3,8),F(-5,9)),(F(2,9),F(-1,19),F(7,3),F(-4,7),F(2,5))]:
        uy=-2*aa*u-pp*v/3; vy=-2*pp*u-ss*v/2
        sk=ss*v*v/6-12*aa*u*u
        sp=-4*pp*u*u+2*ss*u*v-2*v*v
        eq(sk,-4*pp*u*v/3+ss*v*v/6-(4*u*uy+20*aa*u*u),'coordinate Israel source')
        eq(sp,-2*v*v-(4*u*vy+12*pp*u*u),'coordinate scalar source')
        eq(32*u*uy-F(8,3)*v*vy-4*sk,-16*aa*u*u-F(16,3)*pp*u*v+F(2,3)*ss*v*v,'degree-four Robin source')
        eq(-4*v*vy-2*sk,24*aa*u*u+8*pp*u*v+F(5,3)*ss*v*v,'degree-two Robin source')
        # Ambient time derivative equals de Sitter derivative divided by r.
        b,r,base=F(13,2),F(11,3),F(7,5)
        eq((r/b)*base/r,base/b,'proportional ambient velocity')
        count += 5
    named.extend(['coordinate_scalar_and_Israel_normalization',
                  'degree_four_Robin_join','distinct_degree_two_Robin_join',
                  'ambient_clock_proportional_initial_profile'])
    ell=F(1,1000000)
    for rho,c,y,yp,th in [(F(7,3),F(4,5),F(2,7),F(-3,11),F(1,13)),(F(17,2),F(3,5),F(-1,3),F(5,2),F(2,9))]:
        direct=(c*yp-F(3,2)*ell*th*c*y)/(rho**3*ell)
        candidate=c/(rho**3*ell)*(yp-F(3,2)*ell*th*y)
        if args.inject_error=='wrong_metric_clock': candidate*=rho
        eq(candidate,direct,'background proper-time metric derivative')
        count+=1
    mark('background_proper_time_metric_derivative_chain_rule')
    T=F(3,1000000); M=F(102000)
    bgpath=args.root/'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json'
    bg=json.loads(bgpath.read_text())
    state=bg['candidate_state_exact']
    if F(state['ell'])!=ell or F(state['T'])!=T or F(state['Z'])!=8*ell or bg['first_exit']['status']!='PASS':
        raise ValueError('actual background collar state mismatch')
    for component in bg['first_exit']['component_margins'].values():
        if F(component['strict_lower_margin'])<=0 or F(component['strict_upper_margin'])<=0:
            raise ValueError('actual background first-exit margin failed')
    rb=tuple(F(v) for v in bg['first_exit']['component_margins']['rho']['outer']['exact'])
    qb=tuple(F(v) for v in bg['outer_derived']['Q']['exact'])
    if rb[0]<=0 or qb[0]<=0: raise ValueError('positive actual background Q and radius required')
    wpotential=(rb[0]**2*qb[0]-F(35,4)*T*T,rb[1]**2*qb[1])
    if max(map(abs,wpotential))>=M: raise ValueError('actual background potential exceeds contraction input')
    mark('actual_wider_collar_first_exit_and_potential_bound')
    r=M*T*T/2
    if args.inject_error=='drop_bulk_half_factor': r=M*T*T
    eq(r,F(459,1000000000),'interior Volterra contraction')
    if not 0<r<1: raise ValueError('no interior contraction')
    mark('interior_Dirichlet_Volterra_operator_bound')
    priorpath=args.root/'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json'
    prior=json.loads(priorpath.read_text())
    grid=F(prior['cell_data']['time_grid_denominator'])
    den=F(prior['arithmetic']['value_denominator'])
    csvpath=args.root/'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv'
    BH=F(0); last=F(0); ncell=0
    for row in csv.DictReader(csvpath.read_text().splitlines()):
        lo,hi=F(row['left_grid'])/grid,F(row['right_grid'])/grid
        if lo!=last or not lo<hi: raise ValueError('actual-H partition failed')
        hv=max(abs(F(row['H_lo'])/den),abs(F(row['H_hi'])/den))
        BH+=hv*((3-lo)**2-(3-hi)**2)/2
        last=hi; ncell+=1
    if last!=3 or ncell!=4353: raise ValueError('actual-H complete 4353-cell coverage failed')
    q=F(5,2)*ell**2*9/2
    By=ell*BH/(get('mu')[0]*(1-q))
    w0bound=1+T+By; wbound=w0bound/(1-r); first_error=r*wbound
    mark('independent_actual_H_boundary_sup_and_interior_value_error')
    # z+t<=8ell is backward characteristic closed and contains first support.
    Z=8*ell
    if not 2*ell+T<=Z-T: raise ValueError('audited collar misses full first support')
    mark('backward_closed_causal_domain_contains_complete_first_support')
    data={'status':'PASS','scope':'Exact local generic 5D jet tests support analytic full-source identities; no numerical interior evolution or full U2 claim.',
      'check_names':named,'check_count':len(named),'primitive_check_count':count+3,'generic_rational_jet_seeds':[271,1927],
      'shell_box_sha256':hashlib.sha256(boxpath.read_bytes()).hexdigest(),
      'background_derivative_bounds_sha256':hashlib.sha256(bgpath.read_bytes()).hexdigest(),
      'symmetrizer':{'D':[3,8],'determinant_lower_exact':str(detlo),'determinant_lower_display':float(detlo),'euclidean_eigenvalue_lower_exact':str(lambdalo),'euclidean_eigenvalue_lower_display':float(lambdalo),'scope':'D M only; no proof of unconstrained master kinetic matrix'},
      'interior_value_volterra':{'T_exact':str(T),'M_exact':str(M),'Z_exact':str(Z),'domain':'0<=t<=T, 0<=z, z+t<=Z','operator_bound_exact':str(r),'residual_multiplier_exact':str(1/(1-r)),
        'actual_potential_interval_exact':[str(v) for v in wpotential],
        'actual_H_weighted_absolute_integral_exact':str(BH),'boundary_master_sup_exact':str(By),'free_Dirichlet_plus_boundary_sup_exact':str(w0bound),'actual_interior_sup_exact':str(wbound),'one_step_error_sup_exact':str(first_error),
        'display_only':{'boundary_master_sup':float(By),'free_plus_boundary_sup':float(w0bound),'actual_interior_sup':float(wbound),'one_step_error_sup':float(first_error)},'derivative_jets_certified':False},
      'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'PASS','check_count':len(named),'primitive_check_count':count+3}))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument('--report',type=Path)
    ap.add_argument('--inject-error',choices=['half_ricci','drop_potential_cross','drop_lower_gauge_cross','wrong_robin_sign','wrong_metric_clock','drop_bulk_half_factor'])
    args=ap.parse_args()
    try:
        run(args)
    except Exception as exc:
        data={'status':'FAIL','injection':args.inject_error,'error_type':type(exc).__name__,'reason':str(exc)}
        args.output.write_text(json.dumps(data,indent=2)+'\n')
        print(json.dumps(data))
        raise SystemExit(1)

if __name__=='__main__': main()
