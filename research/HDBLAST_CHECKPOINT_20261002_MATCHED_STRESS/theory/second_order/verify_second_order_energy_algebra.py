#!/usr/bin/env python3
"""Exact postrun algebra only: no pulse, mode, spectral, or work sampling.

Run from any directory with Python and SymPy. Local reference copies are
hash-checked; optional --repo also checks their original checkout paths.
This is not an implementation or registration of the proposed experiment.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import sympy as S

HERE = Path(__file__).resolve().parent
checks = []
mutations = []


def zero(name, expression):
    reduced = S.simplify(expression)
    if reduced != 0:
        raise RuntimeError(f"Identity failed: {name}: {reduced}")
    checks.append(name)


def reject(name, expression):
    reduced = S.simplify(expression)
    if reduced == 0:
        raise RuntimeError(f"Mutation escaped: {name}")
    mutations.append(name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    pin_data = json.loads((HERE / 'SOURCE_PINS.json').read_text())
    for item in pin_data['sources']:
        local = HERE / item['local_copy']
        if hashlib.sha256(local.read_bytes()).hexdigest() != item['sha256']:
            raise RuntimeError(f"Reference hash mismatch: {local}")
        if args.repo is not None:
            original = args.repo / item['repo_path']
            if hashlib.sha256(original.read_bytes()).hexdigest() != item['sha256']:
                raise RuntimeError(f"Checkout source hash mismatch: {original}")

    k, K, M, h, a, eps = S.symbols('k K M h a epsilon', positive=True)
    f, fp, q, qp, qpp, qppp, E = S.symbols('f fp q qp qpp qppp E', real=True)
    F0, F1, eta = S.symbols('F0 F1 eta', real=True)
    Fr, Fi = S.symbols('Fr Fi', real=True)
    F = Fr + S.I*Fi
    pi = S.pi
    beta = S.I*F/(2*k)
    alpha = -S.I*F0/(2*k)

    zero('canonical beta energy normalization',
         k**3*beta*S.conjugate(beta)/(2*pi**2)
         -k*(Fr**2+Fi**2)/(8*pi**2))
    omega = S.symbols('omega', positive=True)
    zero('omega equals two k Jacobian',
         (k/(8*pi**2)).subs(k, omega/2)/2-omega/(32*pi**2))
    w = -S.exp(2*S.I*k*eta)*F
    zero('forced mode beta extraction', S.exp(-2*S.I*k*eta)*w/(2*S.I*k)-beta)
    zero('bounded complex amplitude sign', -S.exp(-2*S.I*k*eta)*w-F)
    zero('linear Wronskian', alpha+S.conjugate(alpha))
    zero('infrared bounded u1 limit',
         S.limit(alpha + S.I*(F0-2*S.I*k*F1)*S.exp(2*S.I*k*eta)/(2*k),k,0)
         -(F1-eta*F0))
    zero('occupation coherence leading cancellation',
         F0**2/(4*k**2)+(-S.I*F0/(2*k))*(-S.I*F0/(2*k)))

    alpha2_real, beta1_abs2, alpha1_abs2 = S.symbols('alpha2_real beta1_abs2 alpha1_abs2',real=True)
    wronskian2 = 2*alpha2_real+alpha1_abs2-beta1_abs2
    zero('second-order Wronskian coefficient',
         ((1+eps*alpha+eps**2*alpha2_real)
          *(1+eps*S.conjugate(alpha)+eps**2*alpha2_real)
          -eps**2*beta*S.conjugate(beta)).expand().coeff(eps,2)
         -wronskian2.subs({alpha1_abs2:F0**2/(4*k**2), beta1_abs2:(Fr**2+Fi**2)/(4*k**2)}))

    AK = S.asinh(K/M)-K/S.sqrt(K**2+M**2)
    V = K/S.sqrt(K**2+M**2)
    AKprime = S.diff(AK,M)*h*M
    zero('fixed-comoving cutoff moving-reference derivative', AKprime+h*V**3)
    zero('finite-cutoff work drift coefficient', -AKprime*f**2/(32*pi**2)-h*V**3*f**2/(32*pi**2))
    zero('continuum work drift limit', S.limit(-AKprime,K,S.oo)-h)
    g = S.symbols('log_M_over_mu',real=True)
    local = -g*f**2/(32*pi**2)
    local_total_prime = S.diff(local,g)*h+S.diff(local,f)*fp
    zero('local contact exact derivative and drift',
         -g*f*fp/(16*pi**2)-local_total_prime-h*f**2/(32*pi**2))
    zero('physical source work versus canonical source work',
         a**4*((fp-2*h*f)/a**2)*(q/a**2)/2-(fp*q/2-h*f*q))

    R = E-h*qp/2+3*h**2*q/2
    P = E/3+qpp/6-h*qp/2-h**2*q/2
    def derivative(expr):
        return S.diff(expr,h)*h**2+S.diff(expr,q)*qp+S.diff(expr,qp)*qpp+S.diff(expr,qpp)*qppp
    zero('postpulse minimal stress trace', -R+3*P-(qpp/2-h*qp-3*h**2*q))
    zero('postpulse minimal stress Ward identity',derivative(R)-h*(R-3*P))
    zero('canonical/minimal density operator reduction',
         E-h*qp/2+(h**2+2*h**2)*q/2-R)

    L1, Ln, kap = S.symbols('L1 Ln kappa',positive=True)
    zero('infrared spectral energy bound coefficient',
         S.integrate(k*L1**2/(8*pi**2),(k,0,kap))-L1**2*kap**2/(16*pi**2))
    for n in (2,3):
        integral = S.integrate(k*Ln**2/(8*pi**2*(2*k)**(2*n)),(k,K,S.oo))
        general = Ln**2/(8*pi**2*2**(2*n)*(2*n-2)*K**(2*n-2))
        zero(f'UV spectral energy bound n={n}',integral-general)
    zero('third-derivative UV coefficient',
         Ln**2/(8*pi**2*2**6*4*K**4)-Ln**2/(2048*pi**2*K**4))
    zero('finite work kernel removable endpoint',
         S.limit((1-S.cos(2*K*k))/k,k,0))

    reject('wrong beta global phase is caught by complex amplitude', -S.exp(-2*S.I*k*eta)*(-w)-F)
    reject('missing source-work one-half', 2*E-E)
    reject('omitted finite-cutoff drift',h*V**3*f**2/(32*pi**2))
    reject('reversed finite-cutoff drift',2*h*V**3*f**2/(32*pi**2))
    reject('continuum drift substituted at finite K', h*(1-V**3)*f**2/(32*pi**2))
    reject('physical minimal stress replaced by radiation equation of state',P-R/3)
    reject('occupation-only leading variance',F0**2/(4*k**2))
    reject('first-order truncation claimed normalized through quadratic order',wronskian2.subs(alpha2_real,0))

    report = {
        'status':'PASS',
        'scope':'Exact algebra and analytic endpoint coefficients only; no physical response, source, mode, spectral or work numerical evaluation',
        'proposal_status':'UNREGISTERED_AND_NUMERICALLY_UNEXECUTED',
        'python':platform.python_version(), 'sympy':S.__version__,
        'checks_count':len(checks), 'checks':checks,
        'mutations_rejected_count':len(mutations), 'mutations_rejected':mutations,
        'source_pins_checked':len(pin_data['sources']),
        'checkout_sources_checked':args.repo is not None,
        'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    encoded = json.dumps(report,indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded,end='')


if __name__ == '__main__':
    main()
