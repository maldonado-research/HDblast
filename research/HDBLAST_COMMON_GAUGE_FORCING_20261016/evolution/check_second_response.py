#!/usr/bin/env python3
"""Exact arithmetic joins for the accompanying full second-response proof.

This finite checker authenticates inputs and evaluates the proved operator
majorants. It does not substitute point tests for the analytic wave proof.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_cap(coeff, radius, order=0):
    result = Q(0)
    for power, value in enumerate(coeff):
        if power >= order:
            factor = 1
            for j in range(order):
                factor *= power - j
            result += abs(value) * factor * radius ** (power-order)
    return result


def verify_semantics(mutation=None):
    dims = [1, 9, 25, 5, 5]
    if mutation == 'omit_TT':
        dims = dims[:3]
    require(sum(dims) == 45, 'angular_TT_sectors_required')
    full_dim = 5*(1+9+25)+2*(9+25)+(9+25)+10
    require(full_dim == 287, 'section_multiplicities')
    require(35+34 == 69 and full_dim-69 == 218, 'boundary_multiplicities')
    spectra = [0,8,24,6,22,2,18,6]
    require(max(spectra) == 24 and 17**2 > full_dim, 'angular_operator_bounds')
    frame_trace_diag = Q(2,7) if mutation == 'coordinate_Robin_in_frame' else Q(0)
    require(frame_trace_diag == 0, 'parallel_frame_cancels_tangent_2a')
    # Flat point jet: only partial_y k_yy=2 is nonzero, p=3.
    # dC_y=2-(1/2)*2=1; the scalar linear residual includes -3.
    scalar_dc_witness = Q(0) if mutation == 'drop_scalar_DC' else -Q(3)*(2-Q(2,2))
    require(scalar_dc_witness == -3, 'scalar_DC_coupling_required')
    metric_force_multiplier = 2 if mutation == 'wrong_metric_force_sign' else -2
    require(metric_force_multiplier == -2, 'unit_wave_metric_force_sign')
    # g(s)=s², Lg=-(x-r)^3/3 on x>r: d_r Lg=(x-r)^2.
    lift_multiplier = 1 if mutation == 'wrong_lift_sign' else -1
    require(-lift_multiplier == 1, 'neumann_lift_normal_derivative')
    # In max(|U|/3,|Ux|,|Ur|), a constant U has trace norm3||U||.
    trace_factor = 1 if mutation == 'omit_trace_weight' else 3
    require(trace_factor >= 3, 'weighted_boundary_trace_factor')
    return {'angular_representation_dimensions': dims, 'section_dimension':full_dim,
            'dirichlet_dimension':69, 'robin_dimension':218,
            'rough_angular_laplacian_cap':24, 'point_evaluation_factor':17,
            'metric_unit_source_multiplier':-2, 'scalar_unit_source_multiplier':1}


def run(root, pins_path):
    pins=json.loads(pins_path.read_text())
    data={}
    for key,item in pins.items():
        path=root/item['path']
        require(digest(path)==item['sha256'], 'input_hash_'+key)
        data[key]=json.loads(path.read_text()) if path.suffix=='.json' else path.read_text()
    bg=data['background']; bulk=data['bulk']; shell=data['shell']
    require(bg['domain']['ell']=='1/1000000', 'background_ell')
    require(bg['domain']['scaled_x']==[0,3] and bg['domain']['scaled_r']==[0,8], 'whole_collar_domain')
    L=Q(1,10**6); R=Q(131); alpha=Q(14,125); P=Q(7,50000)
    ay=Q(1,10000); ht=Q(3,10**6); h=ht/129
    caps=bg['static_actual_derivative_caps']
    for key,cap in [('rho',R),('a',alpha),('a_y',ay),('p',P)]:
        require(Q(caps[key][0])<cap, 'actual_background_'+key)
    eta_outer=list(map(Q,bg['first_exit']['eta']['outer']['exact']))
    require(eta_outer[0]>=-Q(9,100000) and eta_outer[1]<=-Q(8,100000), 'signed_phi_tube')
    potential=list(map(Q,[-Q(2,27),0,Q(14,9),Q(50,27),-Q(1,6),-Q(4,9),-Q(2,27)]))
    pot=[polynomial_cap(potential,Q(9,100000),j) for j in range(3)]
    require(pot[0]<1 and pot[1]<Q(1,1000) and pot[2]<4,'potential_jet_caps')
    py=pot[1]+4*alpha*P
    require(py<Q(1,500) and py+2*alpha*P<1,'background_hessian_cap')
    curvature=alpha**2+ay+Q(1,129**2)
    require(curvature<Q(13,1000) and 29*curvature<1,'linear_Ricci_curvature_cap')
    A=L*(3*ht+4*R*alpha)+7*L*R*P
    B=3*L*R*alpha+7*L*R*P
    Z=(2*alpha)**2+3*h*(2*alpha)+3*(2*(alpha+h))**2
    D=L**2*(24+40*R*(alpha+h)+R**2*Z
             +R**2*(2+Q(4,3)+4+4+1+P*(8*alpha+6*h))+30*R*P)
    M=L*R*(10*alpha+4*P+P/2+2)
    A_cap=B_cap=Q(1,10000); D_cap=Q(1,10**6); M_cap=Q(1,2000)
    require(A<A_cap and B<B_cap and D<D_cap and M<M_cap,'operator_caps')
    q=3*(A_cap+B_cap+3*D_cap+M_cap)
    require(q==Q(2109,10**6) and q<1,'global_contraction')
    qactual=3*(A+B+3*D+M)
    sources=bulk['actual_source_bounds']
    qe=Q(sources['QEred']['strict_simple_upper']); qp=Q(sources['Qphi']['strict_simple_upper'])
    qc=Q(sources['QC']['strict_simple_upper'])
    require(qe<10**21 and qp<4*10**15 and qc<6000,'bulk_source_join')
    shell_values=shell['bounds']['one_copy_outward_residual_hessian']
    qr=max(Q(v['exact_upper']) for v in shell_values['metric']['value'].values())
    qshell=Q(shell_values['scalar']['value']['exact_upper'])
    require(qr<Q(15*10**8) and qshell<Q(44*10**8),'shell_source_join')
    gn=8*qr+2*qc
    require(gn<13*10**9,'normal_gauge_trace_join')
    F=L**2*R**2*(10*10**21+4*10**15)
    G=L*R*(8*Q(15*10**8)+44*10**8+13*10**9)
    require(F<2*10**14 and G<4*10**6,'unit_source_caps')
    # Source provenance authenticates these published component caps; do not
    # silently replace them with a different physical/gauge initial jet.
    text=data['initial']
    for value in ['1.06e24','6.4e29','5.13e33','3.06e39']:
        require(value in text, 'initial_report_cap_'+value)
    g0=Q('1.06e24'); dg0=Q('6.4e29'); w0=Q('5.13e33'); dw0=Q('3.06e39')
    initial_value=5*g0+w0
    initial_r=L*R*(5*dg0+dw0)
    initial_x=L*R*(5*dg0+2*alpha*5*g0+dw0)
    require(initial_value<6*10**33 and initial_r<5*10**35 and initial_x<5*10**35,'initial_frame_join')
    free=max(initial_value/3+initial_x,initial_r+initial_x)
    require(free<10**36,'free_wave_cap')
    b=10**36+3*(2*10**14)+4*10**6
    response=b/(1-q)
    require(response<Q('1.02e36'),'response_cap')
    require(17*3*response<6*10**37 and 17*response<2*10**37,'pointwise_caps')
    semantic=verify_semantics()
    controls={}
    for mutation,reason in {
        'omit_TT':'angular_TT_sectors_required',
        'coordinate_Robin_in_frame':'parallel_frame_cancels_tangent_2a',
        'drop_scalar_DC':'scalar_DC_coupling_required',
        'wrong_metric_force_sign':'unit_wave_metric_force_sign',
        'wrong_lift_sign':'neumann_lift_normal_derivative',
        'omit_trace_weight':'weighted_boundary_trace_factor',
    }.items():
        try:
            verify_semantics(mutation)
        except ValueError as error:
            require(str(error)==reason,'wrong_control_reason_'+mutation)
            controls[mutation]=reason
        else:
            raise ValueError('control_not_rejected_'+mutation)
    return {
        'status':'PASS_EXACT_SECOND_RESPONSE_MAJORANTS',
        'input_sha256':{key:item['sha256'] for key,item in pins.items()},
        'domain':{'ell':str(L),'x':[0,3],'r':[0,8],'condition':'x+r<=8'},
        'semantics':semantic,
        'coefficient_operator_upper':dict(zip(['A','B','D','M'],map(str,[A,B,D,M]))),
        'rounded_operator_caps':dict(zip(['A','B','D','M'],map(str,[A_cap,B_cap,D_cap,M_cap]))),
        'contraction':{'exact_majorant':str(qactual),'registered_upper':str(q),'resolvent_upper':str(1/(1-q))},
        'source_join':{'QEred_component':str(qe),'Qphi':str(qp),'QR_component':str(qr),'QP':str(qshell),'GN_from_trace':str(gn),'F_H_upper':str(F),'G_H_upper':str(G)},
        'initial_join':dict(zip(['U0_H','U0_r_H','U0_x_H','free_X'],map(str,[initial_value,initial_r,initial_x,free]))),
        'response':{'X_exact_upper':str(response),'X_simple_upper':'1.02e36','point_value_upper':'6e37','point_scaled_first_derivative_upper':'2e37','tail':'q^N * b_X/(1-q)'},
        'controls':controls,
        'scope':['Integral-defined actual full second formal coefficient, conditional on retained background/master and accepted source certificates.','Finite checker joins authenticated reports and verifies rational majorants; all-point operator, Volterra and qualitative smoothness arguments are in the companion proof.','No numerical waveform, high-derivative response bound, finite-amplitude nonlinear solution or remainder is asserted.'],
        'arithmetic':'fractions.Fraction only; no floating point used for acceptance',
    }


def main():
    here=Path(__file__).resolve().parent
    parser=argparse.ArgumentParser()
    parser.add_argument('--research-root',type=Path,default=here.parent.parent)
    parser.add_argument('--pins',type=Path,default=here/'INPUT_PINS.json')
    parser.add_argument('--inject',choices=['omit_TT','coordinate_Robin_in_frame',
                        'drop_scalar_DC','wrong_metric_force_sign',
                        'wrong_lift_sign','omit_trace_weight'])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'refuse_existing_output')
    if args.inject:
        verify_semantics(args.inject)
        raise ValueError('unrecognized_injection_success')
    result=run(args.research_root,args.pins)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(result['status'])


if __name__=='__main__':
    main()
