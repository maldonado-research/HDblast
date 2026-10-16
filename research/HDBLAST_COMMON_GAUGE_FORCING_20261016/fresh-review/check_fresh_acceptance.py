#!/usr/bin/env python3
"""Separate exact joins for the reviewed conditional second-response theorem.

No scientific producer module is imported. The companion review supplies the
all-point arguments. This bounded program authenticates that scientific corpus
and independently checks selected norm, source, and response inequalities.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def read(root, path):
    return json.loads((root/path).read_text())


def mag(interval):
    lo, hi = map(F, interval)
    need(lo <= hi, 'ordered interval')
    return max(abs(lo), abs(hi))


def audit(root):
    background = read(root, 'bulk-producer/inputs/GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json')
    bulk = read(root, 'bulk-producer/COMMON_GAUGE_BULK_NORMAL.json')
    shell = read(root, 'shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json')
    evolution = read(root, 'evolution/SECOND_RESPONSE_FINAL_NORMAL.json')
    initial = read(root, 'evolution-review/inputs/INITIAL_SPACETIME_FRAME_NORMAL.json')
    ell, R, a, p = F(1, 10**6), F(131), F(14,125), F(7,50000)
    eta = tuple(map(F, background['first_exit']['eta']['outer']['exact']))
    need(-F(9,100000) <= eta[0] <= eta[1] <= -F(8,100000),
         'signed actual whole-collar scalar tube')
    rho_lo,rho_hi = map(F, background['first_exit']['rho']['outer']['exact'])
    need(129 <= rho_lo <= rho_hi <= R, 'whole-collar lower and upper warp bounds')
    for name, upper in [('rho',R), ('a',a), ('p',p), ('a_y',F(1,10000))]:
        need(F(background['static_actual_derivative_caps'][name][0]) < upper,
             'whole-collar coefficient cap: '+name)
    need(background['domain']['ell'] == str(ell), 'same actual length scale')
    need(background['domain']['scaled_x'] == [0,3]
         and background['domain']['scaled_r'] == [0,8], 'whole collar input')

    # Check every contraction in the direct source envelopes independently
    # from the certified field jets. These are component maxima, not energies.
    h0,h1,h2 = map(F,bulk['field_covariant_caps']['metric'])
    e0,e1,e2 = map(F,bulk['field_covariant_caps']['scalar'])
    c = {k:F(v) for k,v in bulk['connection_and_background_caps'].items()}
    P0,P1 = F(3,2)*h1, F(3,2)*h2
    J1 = 10*(h1*P0+h0*P1)
    qe = 10*J1+100*P0**2+2*e1**2+F(2,3)*(c['U2']*e0**2+2*c['U1']*e0*h0)
    qc = 50*h0*P0
    dqc = 50*(h1*P0+h0*P1)
    qphi = 50*(h0*h0*c['Phi2']+h0*e2+h0*P0*c['p'])+c['U3']*e0**2
    direct = {'QE':qe, 'QC':qc, 'nabla_QC':dqc, 'Qphi':qphi, 'QEred':qe+dqc}
    for key, value in direct.items():
        row = bulk['actual_source_bounds'][key]
        need(value == F(row['exact_upper']), 'direct full source contraction: '+key)
        need(value < F(row['strict_simple_upper']), 'strict displayed source cap: '+key)

    # The full section count distinguishes representation dimensions from
    # scalar, vector and tensor section multiplicities.
    scalar_dims = (1,9,25)
    Hdim = 5*sum(scalar_dims)+3*sum(scalar_dims[1:])+5+5
    Ddim = sum(scalar_dims)+sum(scalar_dims[1:])
    need(Hdim == 287 and Ddim == 69 and Hdim-Ddim == 218,
         'all metric/scalar section and boundary multiplicities')
    spectra = [0,8,24, 8-2,24-2, 8-6,24-6, 6,6]
    need(max(spectra) == 24 and 17**2 > Hdim,
         'full-section intrinsic derivative and evaluation constants')

    ht = 3*ell
    beta = ht/129
    A = ell*(3*ht+4*R*a)+7*ell*R*p
    B = 3*ell*R*a+7*ell*R*p
    Z = 4*a*a+6*beta*a+12*(a+beta)**2
    D = ell**2*(24+40*R*(a+beta)+R*R*Z
                +R*R*(2+F(4,3)+4+4+1+p*(8*a+6*beta))+30*R*p)
    for key, value in [('A',A), ('B',B), ('D',D)]:
        need(value == F(evolution['coefficient_operator_upper'][key]),
             'separate exact coefficient majorant: '+key)
    need(A < F(1,10000) and B < F(1,10000) and D < F(1,10**6),
         'full bulk operator caps')

    # Independent Robin estimate: use the squared Frobenius norm of its
    # normalized (trace/2, k_yy, scalar) block. The tangent STF block is zero.
    params = shell['bounds']['parameters']
    sigma1 = F(params['sigma_prime_abs'])
    sigma2 = mag(params['sigma_second_interval'])
    need(sigma1 <= 2*p and sigma2 < 4, 'actual signed shell-tension cap')
    robin_sq = (ell*R)**2*(36*a*a+F(20,9)*sigma1*sigma1+p*p/4+sigma2*sigma2/4)
    need(robin_sq < F(1,2000)**2, 'independent Robin block Frobenius bound')

    q = 3*(F(1,10000)+F(1,10000)+3*F(1,10**6)+F(1,2000))
    need(q == F(2109,10**6) and q == F(evolution['contraction']['registered_upper']),
         'weighted source and boundary trace contraction')
    need(q < 1, 'strict second-response contraction')

    # All ten lower tensor entries, including the six off-diagonal ones,
    # must be supplied. The metric equation uses twice the one-copy Hessian.
    physical = shell['bounds']['one_copy_outward_residual_hessian']
    rows = physical['metric']['value']
    need(set(rows) == {'00','01','02','03','11','12','13','22','23','33'},
         'complete physical tensor source rows')
    qr = max(F(v['exact_upper']) for v in rows.values())
    qp = F(physical['scalar']['value']['exact_upper'])
    forcing = shell['bounds']['coordinate_physical_forcing']
    for key,row in rows.items():
        need(F(forcing['metric']['value'][key]['exact_upper']) == 2*F(row['exact_upper']),
             'metric equation one/two-copy adapter: '+key)
    need(F(forcing['scalar']['value']['exact_upper']) == qp,
         'scalar equation retains one-copy normalization')

    # Sharper independent response bound using actual source envelopes and
    # the explicit initial component caps, rather than the producer's rounded b.
    qred_cap = F(bulk['actual_source_bounds']['QEred']['strict_simple_upper'])
    qphi_cap = F(bulk['actual_source_bounds']['Qphi']['strict_simple_upper'])
    qc_cap = F(bulk['actual_source_bounds']['QC']['strict_simple_upper'])
    gn = 8*qr+2*qc_cap
    Fcap = ell**2*R**2*(10*qred_cap+qphi_cap)
    Gcap = ell*R*(8*qr+qp+gn)
    initial_keys = ('metric2_components',
                    'metric2_full_initial_covariant_derivative_components',
                    'scalar2', 'scalar2_full_initial_covariant_derivative_components')
    initial_caps = [F(initial['convenient_strict_caps'][key]) for key in initial_keys]
    need(initial_caps == list(map(F,('1.06e24','6.4e29','5.13e33','3.06e39'))),
         'actual initial spacetime report caps')
    for key,cap in zip(initial_keys,initial_caps):
        need(F(initial['bounds'][key]['exact_upper']) < cap,
             'strict inherited initial cap: '+key)
    k0,dk0,w0,dw0 = initial_caps
    U0 = 5*k0+w0
    Ur0 = ell*R*(5*dk0+dw0)
    Ux0 = ell*R*(5*dk0+10*a*k0+dw0)
    free = max(U0/3+Ux0,Ur0+Ux0)
    for key,value in [('U0_H',U0),('U0_r_H',Ur0),('U0_x_H',Ux0),('free_X',free)]:
        need(value == F(evolution['initial_join'][key]), 'actual initial/frame join: '+key)
    response = (free+3*Fcap+Gcap)/(1-q)
    need(response < F(evolution['response']['X_simple_upper']), 'independent response norm')
    need(51*response < F(evolution['response']['point_value_upper']), 'independent point values')
    need(17*response < F(evolution['response']['point_scaled_first_derivative_upper']),
         'independent point first derivatives')
    return {'status':'PASS_FRESH_CONDITIONAL_SECOND_RESPONSE_REVIEW',
            'source_contractions':{k:str(v) for k,v in direct.items()},
            'section_dimension':Hdim, 'Dirichlet_dimension':Ddim, 'Robin_dimension':Hdim-Ddim,
            'independent_Robin_operator_squared_upper':str(robin_sq),
            'registered_contraction':str(q),
            'independent_source_join':{'F_H':str(Fcap),'G_H':str(Gcap),'GN':str(gn)},
            'independent_response_X_upper':str(response),
            'scope':'Selected exact norm joins for the analytically reviewed integral-defined formal second coefficient. No producer imports, numerical waveform, nonlinear family/remainder, high-derivative response cap or physical-origin claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args = parser.parse_args()
    need(not args.output.exists(),'refuse to overwrite output')
    pins = json.loads(Path(__file__).with_name('FRESH_REVIEW_PINS.json').read_text())
    for relative,entry in pins['reviewed_files'].items():
        need(not relative.startswith('fresh-review/') and relative not in
             ('MANIFEST.json','REPLAY_ROSTER.json','reproduce.py'),'nonrecursive scientific pins')
        raw = (args.root/relative).read_bytes()
        need(len(raw)==entry['bytes'] and sha256(raw).hexdigest()==entry['sha256'],
             'reviewed scientific file changed: '+relative)
    result = audit(args.root)
    result['reviewed_files'] = pins['reviewed_files']
    with args.output.open('x') as out:
        json.dump(result,out,indent=2,sort_keys=True)
        out.write('\n')


if __name__ == '__main__':
    main()
