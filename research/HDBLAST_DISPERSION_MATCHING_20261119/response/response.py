"""Euclidean finite-bath dispersion with two zero-frequency matching conditions.

The spectral weights and original 64-node bath rule are those of the preceding
finite Yukawa Hamiltonian. No nonzero probe is used to set a matching parameter.
"""
import argparse
import csv
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

HERE = Path(__file__).resolve().parent
PI8 = 8*np.pi**2


@lru_cache(maxsize=None)
def quadrature_rule(order):
    return roots_legendre(order)


def analytic(q, cutoff, z):
    """Continuous spectral integral, with stable twice-subtracted remainder."""
    q, z = np.asarray(q), np.asarray(z)
    c = (cutoff**2-q*q-q*q*np.log(cutoff**2/(q*q)))/PI8
    d = (np.log(cutoff**2/(q*q))+q*q/cutoff**2-1)/PI8
    infinite = ((q*q+z)*np.log1p(z/(q*q))-z)/PI8
    missing_tail = ((q*q+z)*np.log1p(z/cutoff**2)-z*q*q/cutoff**2)/PI8
    matched = infinite-missing_tail
    static = -z*d+matched
    raw = c+static
    return raw,c,-d,static,matched,infinite,missing_tail


def finite_bath(q, cutoff, z, order):
    """Exact response of the quadrature-discretized oscillator Hamiltonian."""
    nodes, quadrature = quadrature_rule(order)
    omega = q+(cutoff-q)*(nodes+1)/2
    dw = (cutoff-q)*quadrature/2
    g2 = omega*(omega*omega-q*q)*dw/(4*np.pi**2)
    om2 = omega*omega
    c = np.sum(g2/om2)
    d = np.sum(g2/om2**2)
    z = np.asarray(z)
    raw = np.sum(g2/(om2+z[...,None]),axis=-1)
    # Subtract under the summation to avoid cancelling O(Lambda^2) terms.
    static = -z*np.sum(g2/(om2*(om2+z[...,None])),axis=-1)
    matched = z*z*np.sum(g2/(om2**2*(om2+z[...,None])),axis=-1)
    return raw,c,-d,static,matched


def stable_integral(q,cutoff,nu):
    """Independent smooth-variable continuous integral of the remainder."""
    ratio2 = (nu/q)**2
    integral,error = quad(lambda t:(1-t)/(1+ratio2*t),q*q/cutoff**2,1,
                          epsabs=1e-14,epsrel=1e-13)
    prefactor = nu**4/(PI8*q*q)
    return prefactor*integral,prefactor*error


def calculate(modes_path,protocol_path,output):
    modes_path,protocol_path,output = map(Path,(modes_path,protocol_path,output))
    modes,protocol = json.loads(modes_path.read_bytes()),json.loads(protocol_path.read_bytes())
    if hashlib.sha256(modes_path.read_bytes()).hexdigest()!=protocol['native_modes_sha256']:
        raise ValueError('Mode catalog differs from the predeclared protocol')
    q,weights = np.array(modes['q']),np.array(modes['weight'])
    cutoffs = np.array(protocol['cutoffs'])
    ratios = np.array(protocol['frequency_ratios_nu_over_q'])
    orders = np.array(protocol['bath_orders'])
    if len(q)!=24 or np.min(q)<=0 or min(cutoffs)<=max(q):
        raise ValueError('Require all24 positive native modes below each cutoff')
    nodes,dw = np.polynomial.legendre.leggauss(24)
    np.testing.assert_allclose(q,4*(nodes+1),rtol=1e-13)
    np.testing.assert_allclose(weights,4*dw*q*q/(2*np.pi**2),rtol=2e-12)
    nu = q[:,None]*ratios
    z = nu**2
    shape = (len(orders),len(cutoffs),len(q),len(ratios))
    fields = {name:np.empty(shape) for name in ('raw','static','matched')}
    sigma0 = np.empty(shape[:-1])
    slope0 = np.empty(shape[:-1])
    exact = {name:np.empty(shape[1:]) for name in ('raw','static','matched','tail','smooth_integral')}
    exact_sigma0 = np.empty(shape[1:-1])
    exact_slope0 = np.empty(shape[1:-1])
    infinite = analytic(q[:,None],cutoffs[-1],z)[5]
    for ci,cutoff in enumerate(cutoffs):
        raw,c,cp,static,matched,_,tail = analytic(q[:,None],cutoff,z)
        exact['raw'][ci],exact['static'][ci],exact['matched'][ci],exact['tail'][ci] = raw,static,matched,tail
        exact_sigma0[ci],exact_slope0[ci] = c[:,0],cp[:,0]
        for qi,k in enumerate(q):
            exact['smooth_integral'][ci,qi] = [stable_integral(k,cutoff,v)[0] for v in nu[qi]]
            for ni,order in enumerate(orders):
                raw,c,cp,static,matched = finite_bath(k,cutoff,z[qi],int(order))
                fields['raw'][ni,ci,qi],fields['static'][ni,ci,qi],fields['matched'][ni,ci,qi] = raw,static,matched
                sigma0[ni,ci,qi],slope0[ni,ci,qi] = c,cp
    quadrature_relative = (fields['matched']-exact['matched'][None])/exact['matched'][None]
    cutoff_relative = exact['tail']/infinite[None]
    tail_to_leading = exact['tail']/(z[None]**2/(PI8*cutoffs[:,None,None]**2))
    tail_exponent = np.log(exact['tail'][:-1]/exact['tail'][1:])/np.log(cutoffs[1:,None,None]/cutoffs[:-1,None,None])
    adjacent_cutoff_change = (exact['matched'][1:]-exact['matched'][:-1])/exact['matched'][1:]
    smooth_relative = abs(exact['smooth_integral']-exact['matched'])/exact['matched']
    np.testing.assert_allclose(exact['smooth_integral'],exact['matched'],rtol=2e-12,atol=1e-15)
    assert np.max(abs(quadrature_relative[-1]))<3e-9
    assert np.all(fields['matched']>0)
    assert np.all(exact['tail']>0)
    assert np.all(np.diff(exact['matched'],axis=0)>0)
    # Ratios are dimensionless; native momentum weights enter only the
    # dimensionless weighted error below, never the one-mode kernel itself.
    weighted_quad = np.sum(weights[None,None,:,None]*abs(fields['matched']-exact['matched'][None]),axis=2)/np.sum(weights[None,:,None]*exact['matched'],axis=1)[None]
    weighted_cutoff = np.sum(weights[None,:,None]*exact['tail'],axis=1)/np.sum(weights[:,None]*infinite,axis=0)[None]
    output.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output/'TABLES.npz',q=q,weight=weights,cutoffs=cutoffs,ratios=ratios,orders=orders,
        nu=nu,z=z,sigma0=sigma0,slope0=slope0,
        finite_raw=fields['raw'],finite_static=fields['static'],finite_matched=fields['matched'],
        analytic_sigma0=exact_sigma0,analytic_slope0=exact_slope0,
        analytic_raw=exact['raw'],analytic_static=exact['static'],analytic_matched=exact['matched'],
        infinite_matched=infinite,missing_tail=exact['tail'],smooth_integral=exact['smooth_integral'],
        quadrature_relative_error=quadrature_relative,cutoff_relative_error=cutoff_relative,
        weighted_quadrature_relative_l1=weighted_quad,weighted_cutoff_relative=weighted_cutoff,
        missing_tail_to_leading_asymptotic=tail_to_leading,missing_tail_doubling_exponent=tail_exponent,
        adjacent_cutoff_relative_change=adjacent_cutoff_change)
    report = dict(schema=1,experiment='fixed-background finite oscillator bath, imaginary-frequency response',
        inputs={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (modes_path,protocol_path)},
        axes=dict(finite=['order','cutoff','mode','frequency_ratio'],analytic=['cutoff','mode','frequency_ratio']),
        units=dict(q='comoving H0',nu='conformal H0',z='conformal H0^2',response='H0^2',sigma0='H0^2',slope0='dimensionless',weight='native comoving H0^3'),
        q_range=[float(min(q)),float(max(q))],native_mode_count=len(q),
        anchors=dict(target_value=0.,target_slope=0.,point='z=0',per_mode=True,nonzero_probes_used=False),
        maximum_smooth_integral_relative_error=float(np.max(smooth_relative)),
        quadrature=[],cutoff=[],matching_coefficients=[])
    for ni,order in enumerate(orders):
        for ci,cutoff in enumerate(cutoffs):
            where = np.unravel_index(np.argmax(abs(quadrature_relative[ni,ci])),quadrature_relative[ni,ci].shape)
            report['quadrature'].append(dict(order=int(order),cutoff=int(cutoff),
                max_relative_error=float(np.max(abs(quadrature_relative[ni,ci]))),
                worst_q=float(q[where[0]]),worst_nu_over_q=float(ratios[where[1]]),
                weighted_relative_l1=weighted_quad[ni,ci].tolist()))
    for ci,cutoff in enumerate(cutoffs):
        report['matching_coefficients'].append(dict(cutoff=int(cutoff),
            analytic_C_range=[float(np.min(exact_sigma0[ci])),float(np.max(exact_sigma0[ci]))],
            analytic_D_range=[float(np.min(-exact_slope0[ci])),float(np.max(-exact_slope0[ci]))],
            C_units='H0^2, conformal',D_units='dimensionless'))
        report['cutoff'].append(dict(cutoff=int(cutoff),
            max_relative_missing_tail_by_ratio=np.max(cutoff_relative[ci],axis=0).tolist(),
            weighted_relative_missing_tail_by_ratio=weighted_cutoff[ci].tolist(),
            overall_max_relative_missing_tail=float(np.max(cutoff_relative[ci])),
            min_relative_missing_tail=float(np.min(cutoff_relative[ci])),
            tail_to_leading_asymptotic_range=[float(np.min(tail_to_leading[ci])),float(np.max(tail_to_leading[ci]))]))
    report['adjacent_cutoffs'] = [dict(lower=int(cutoffs[i]),upper=int(cutoffs[i+1]),
        max_relative_heldout_change=float(np.max(adjacent_cutoff_change[i])),
        missing_tail_doubling_exponent_range=[float(np.min(tail_exponent[i])),float(np.max(tail_exponent[i]))])
        for i in range(len(cutoffs)-1)]
    report['continuum_limits'] = dict(
        exact='R_inf=[(q^2+z)*log1p(z/q^2)-z]/(8*pi^2)',
        leading_missing_tail='nu^4/(8*pi^2*Lambda^2)',
        low_frequency='R_inf=nu^4/(16*pi^2*q^2)+O(nu^6/q^4)',
        proper_time='R_phys(nu_proper)=H0^2/a^2*R_conf(a*nu_proper/H0)',
        dimensional_tail='R_phys missing tail ~ nu_proper^4/(8*pi^2*Lambda_proper^2)',
        frozen_gate='multiply one-mode self-energy by y^2*F0^2')
    with (output/'MODES.csv').open('w',newline='') as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(['mode','q','weight','cutoff','nu_over_q','nu','z','sigma0_analytic','slope0_analytic','R_N64','R_N2048','R_analytic_finite','R_infinite','relative_cutoff_error','relative_N64_error','relative_N2048_error'])
        for ci,cutoff in enumerate(cutoffs):
            for qi,k in enumerate(q):
                for ri,ratio in enumerate(ratios):
                    writer.writerow([qi,k,weights[qi],cutoff,ratio,nu[qi,ri],z[qi,ri],exact_sigma0[ci,qi],exact_slope0[ci,qi],
                        fields['matched'][0,ci,qi,ri],fields['matched'][-1,ci,qi,ri],exact['matched'][ci,qi,ri],infinite[qi,ri],
                        cutoff_relative[ci,qi,ri],quadrature_relative[0,ci,qi,ri],quadrature_relative[-1,ci,qi,ri]])
    report['tables_sha256'] = hashlib.sha256((output/'TABLES.npz').read_bytes()).hexdigest()
    (output/'RESPONSE.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes',type=Path,default=HERE/'NATIVE_MODES.json')
    parser.add_argument('--protocol',type=Path,default=HERE/'PROTOCOL.json')
    parser.add_argument('--output',type=Path,default=HERE)
    args = parser.parse_args()
    report = calculate(args.modes,args.protocol,args.output)
    print(json.dumps({'cutoff':report['cutoff'],
        'N64_and_N2048':[r for r in report['quadrature'] if r['order'] in (64,2048)]},indent=2))
