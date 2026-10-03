"""Prospective independent FRW mode control; no primary implementation imports.

Execution requires an explicitly supplied, hash-checked registration file. Merely
importing this module does not execute field evolution or write a research result.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

R = 4.0
T = 1.0
AMPLITUDES = (0.0, 0.2)
MOMENTA = (0.0, 0.5, 2.0, 8.0, 24.0, 48.0)
TIMES = np.array([0., .125, .25, .375, .5, .625, .75, .875, 1., 1.25])
RTOL, ATOL = 2e-12, 2e-14


def background(eta, amplitude):
    """Analytic logistic derivatives, with exact static exterior branches."""
    s = eta / T
    if s <= 0.0:
        b, bs, bss = 0.0, 0.0, 0.0
    elif s >= 1.0:
        b, bs, bss = 1.0, 0.0, 0.0
    elif min(s, 1.0-s) < 1e-6:
        # Every derivative used here is far below the smallest float64.
        b, bs, bss = (0.0 if s < .5 else 1.0), 0.0, 0.0
    else:
        g = -1.0/s + 1.0/(1.0-s)
        gp = 1.0/s**2 + 1.0/(1.0-s)**2
        gpp = -2.0/s**3 + 2.0/(1.0-s)**3
        q = math.exp(-abs(g))
        b = q/(1.0+q) if g < 0.0 else 1.0/(1.0+q)
        logistic_slope = q/(1.0+q)**2
        # Use q directly: b can round to 1 while its derivatives are nonzero.
        one_minus_twob = (1.0-q)/(1.0+q) if g < 0.0 else (q-1.0)/(q+1.0)
        bs = logistic_slope*gp
        bss = logistic_slope*(gpp + one_minus_twob*gp**2)
    a = math.exp(amplitude*b)
    hh = amplitude*bs/T
    app_over_a = amplitude*bss/T**2 + hh**2
    mass2 = R*(2.0*b-1.0)**2
    mass2_prime = 4.0*R*(2.0*b-1.0)*bs/T
    return a, hh, app_over_a, mass2, mass2_prime


def observables(eta, u, v, momentum, amplitude):
    a, hh, app, mass2, mass2_prime = background(eta, amplitude)
    z = v-hh*u
    uu, zz = abs(u)**2, abs(z)**2
    q = momentum**2+a*a*mass2
    rho = (zz+q*uu)/(2*a**4)
    pressure = (zz-(momentum**2/3+a*a*mass2)*uu)/(2*a**4)
    current = uu/a**2
    # Differentiate the unreduced quadratic expression using the physical ODE.
    vprime = -(q-app)*u
    hprime = app-hh**2
    zprime = vprime-hprime*u-hh*v
    qprime = a*a*(2*hh*mass2+mass2_prime)
    rhoprime = (2*(z.conjugate()*zprime).real + qprime*uu
                + 2*q*(u.conjugate()*v).real)/(2*a**4)-4*hh*rho
    ward = rhoprime+3*hh*(rho+pressure)-mass2_prime*current/2
    # These controls deliberately change one member of the paired observables.
    broken_current = rhoprime+3*hh*(rho+pressure)
    broken_pressure = rhoprime+3*hh*rho-mass2_prime*current/2
    return dict(rho=rho, p=pressure, Q=current, rho_prime=rhoprime,
                exchange=ward, broken_current_exchange=broken_current,
                broken_pressure_exchange=broken_pressure)


def solve_case(momentum, amplitude):
    wi = math.sqrt(momentum**2+R)
    wf = math.sqrt(momentum**2+math.exp(2*amplitude)*R)
    maxstep = min(T/100, .1/wf)
    y0 = [1/math.sqrt(2*wi), 0., 0., -math.sqrt(wi/2)]

    def rhs(eta, y):
        a, _, app, mass2, _ = background(eta, amplitude)
        potential = momentum**2+a*a*mass2-app
        return [y[2], y[3], -potential*y[0], -potential*y[1]]

    def jacobian(eta, y):
        a, _, app, mass2, _ = background(eta, amplitude)
        potential = momentum**2+a*a*mass2-app
        return np.array([[0.,0.,1.,0.],[0.,0.,0.,1.],
                         [-potential,0.,0.,0.],[0.,-potential,0.,0.]])

    # Evolve the actual background on [0,T]; propagate the future analytically.
    sol = solve_ivp(rhs, (0., T), y0, method="Radau", jac=jacobian,
                    rtol=RTOL, atol=ATOL, max_step=maxstep, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    end = sol.y[:, -1]
    ue, ve = complex(end[0],end[1]), complex(end[2],end[3])
    alpha = math.sqrt(wf/2)*ue+1j*ve/math.sqrt(2*wf)
    beta = math.sqrt(wf/2)*ue-1j*ve/math.sqrt(2*wf)
    # Here the future phase origin is eta=T; |beta|^2 is phase invariant.
    samples = []
    for eta in TIMES:
        if eta <= T:
            yy = sol.sol(eta)
            u,v = complex(yy[0],yy[1]), complex(yy[2],yy[3])
        else:
            phase = wf*(eta-T)
            c,s = math.cos(phase), math.sin(phase)
            u,v = ue*c+ve*s/wf, ve*c-ue*wf*s
        wronskian = u*v.conjugate()-v*u.conjugate()
        record = dict(eta=float(eta), u=[u.real,u.imag], v=[v.real,v.imag],
                      wronskian_error=abs(wronskian/1j-1),
                      **observables(eta,u,v,momentum,amplitude))
        if eta == 0. or eta >= T:
            aa = 1. if eta == 0. else math.exp(amplitude)
            ww = wi if eta == 0. else wf
            record["static_subtracted"] = dict(
                rho=record["rho"]-ww/(2*aa**4),
                p=record["p"]-momentum**2/(6*aa**4*ww),
                Q=record["Q"]-1/(2*aa**2*ww))
        if eta >= T:
            phase = complex(math.cos(-2*wf*(eta-T)),math.sin(-2*wf*(eta-T)))
            coherent = (alpha*beta.conjugate()*phase).real
            number = abs(beta)**2
            record["future_from_bogoliubov"] = dict(
                rho=wf*number/aa**4,
                p=(momentum**2*number/3-(2*momentum**2/3+aa*aa*R)*coherent)/(aa**4*wf),
                Q=(number+coherent)/(aa**2*wf))
            record["future_diagonal_only"] = dict(
                rho=wf*number/aa**4,
                p=momentum**2*number/(3*aa**4*wf),Q=number/(aa**2*wf))
        samples.append(record)
    y = sol.y
    uv = y[0]+1j*y[1]
    vv = y[2]+1j*y[3]
    maxwr = float(np.max(np.abs((uv*vv.conjugate()-vv*uv.conjugate())/1j-1)))
    return dict(A=amplitude,k=momentum,rtol=RTOL,atol=ATOL,max_step=maxstep,
                occupation=float(abs(beta)**2),alpha=[alpha.real,alpha.imag],
                beta=[beta.real,beta.imag],bogoliubov_error=abs(abs(alpha)**2-abs(beta)**2-1),
                wronskian_max_internal=maxwr, nfev=sol.nfev,njev=sol.njev,nlu=sol.nlu,
                samples=samples)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    actual = hashlib.sha256(args.protocol.read_bytes()).hexdigest()
    if actual != args.protocol_sha256:
        raise ValueError("Protocol hash mismatch; field evolution was not started")
    cases = [solve_case(k,a) for a in AMPLITUDES for k in MOMENTA]
    zero_controls = all(max(abs(s["broken_current_exchange"]) for s in c["samples"])>1e-5
                        for c in cases)
    pressure_controls = all(max(abs(s["broken_pressure_exchange"]) for s in c["samples"])>1e-5
                            for c in cases if c["A"] != 0.)
    wronskian_pass = all(c["wronskian_max_internal"]<1e-8 for c in cases)
    exchange_pass = all(abs(s["exchange"]) <= 1e-10*(1+abs(s["rho_prime"])
                       +abs(s["rho"])+abs(s["p"])+abs(s["Q"]))
                       for c in cases for s in c["samples"])
    static_pass = all(abs(v)<1e-12*(1+c["k"]) for c in cases
                      for v in c["samples"][0]["static_subtracted"].values())
    future_pass = all(abs(s["static_subtracted"][key]-s["future_from_bogoliubov"][key])
                      <2e-10*(1+abs(s[key]))
                      for c in cases for s in c["samples"] if s["eta"]>=T
                      for key in ("rho","p","Q"))
    result = dict(scope="Independent prescribed-background physical mode check only",
                  solver="SciPy Radau, four real components, analytic Jacobian",
                  protocol_sha256=actual,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  parameters=dict(r=R,T=T,A=AMPLITUDES,k=MOMENTA,times=TIMES.tolist()),
                  gates=dict(wronskian=wronskian_pass,bare_exchange=exchange_pass,
                             missing_current_detected=zero_controls,
                             missing_pressure_detected=pressure_controls,
                             incoming_static_vacuum=static_pass,
                             future_coherent_decomposition=future_pass),cases=cases)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"output":str(args.output),"gates":result["gates"]}))
    return 0 if all(result["gates"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
