"""Energy-neutral uniform canonical shear and its finite-time oscillator control."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def run(checkpoint, normalization):
    with np.load(checkpoint, allow_pickle=False) as snapshot:
        if snapshot["mode"].item() != "full":
            raise ValueError("Expected native full-field variables at fork")
        log_scale, _, phi = snapshot["bulk"][:3, -1]
        X, P, C = np.split(snapshot["matter"][3:], 3)
    event = json.loads(normalization.read_bytes())["matter"]["events"][0]
    a = np.exp(log_scale - event["A_star"])
    mass2 = (100 * (phi - .5))**2
    nodes, weight = np.polynomial.legendre.leggauss(len(X))
    q = 4 * (nodes + 1)
    weight *= 4*q*q/(2*np.pi**2)
    energy2 = mass2 + (q/a)**2
    scale = a**3
    shear = -2 * (weight @ C)/(weight @ X)
    target = (X, P + 2*shear*C + shear**2*X, C + shear*X)
    def sources(state):
        dx, dp = state[0]-X, (state[1]-P)/scale**2
        return dict(rho=float(weight @ (dp + energy2*dx)/2),
                    pressure=float(weight @ (dp-(mass2+(q/a)**2/3)*dx)/2),
                    force=float(1e4*(phi-.5)*(weight @ dx)))
    def occupation(state):
        return .5*(scale*np.sqrt(energy2)*state[0]+state[1]/(scale*np.sqrt(energy2)) - 1)
    def finite(duration, strength):
        shifted = energy2 - strength/(scale*duration)
        if np.min(shifted) <= 0:
            raise ValueError("Finite test uses a positive oscillator frequency squared")
        w = np.sqrt(shifted)
        angle = duration*w
        co, sn = np.cos(angle), np.sin(angle)
        b, d = sn/(scale*w), -scale*w*sn
        dx = -sn**2*X + 2*co*b*C + b*b*P
        return (X+dx, d*d*X+2*co*d*C+co*co*P, co*d*X+(co*co+b*d)*C+b*co*P), dx
    rows=[]
    for duration in (.01,.005,.0025,.00125,.000625):
        def objective(value):
            return weight @ finite(duration,value)[1] / duration
        # The nontrivial zero-work branch sets the weighted field variance change to zero.
        # Its negative sign is specific to the measured covariance and these durations.
        bound=2*abs(shear)
        for _ in range(10):
            if objective(-bound)*objective(0)<=0:
                break
            bound*=2
        else:
            raise ValueError("Could not bracket the nonzero-strength pulse")
        strength=brentq(objective,-bound,0,xtol=1e-14)
        state,dx=finite(duration,strength)
        src=sources(state)
        amplitude=-strength/(scale*duration)
        work=-.5*scale*amplitude*(weight@dx)
        # Also compare to free evolution through the same duration for the field moment.
        # The conserved free energy (and hence the zero-work criterion) is common to
        # the initial state and the freely propagated state on a static background.
        free,_=finite(duration,0.)
        dist=np.sqrt(weight @ ((state[2]-target[2])**2)/(weight @ ((target[2]-C)**2)))
        if abs(work-scale*src['rho'])>3e-12 or abs(src['force'])>1e-10:
            raise ValueError("Finite-pulse work or variance check")
        rows.append(dict(duration=duration,integrated_shear=strength,delta_mass_squared=amplitude,
                         external_work=work,source_changes_vs_initial=src,
                         pressure_change_vs_free=src['pressure']-sources(free)['pressure'],
                         relative_cross_covariance_error_vs_ideal=float(dist)))
    delta_n=occupation(target)-occupation((X,P,C))
    return dict(a=float(a),mass=float(np.sqrt(mass2)),uniform_shear=float(shear),
                ideal_source_changes=sources(target),
                integrated_C=float(weight@C),integrated_C_after=float(weight@target[2]),
                max_extractable_density_uniform=float((weight@C)**2/(2*scale**2*(weight@X))),
                max_extractable_density_modewise=float(weight@(C*C/X)/(2*scale**2)),
                max_determinant_drift=float(np.max(abs(target[0]*target[1]-target[2]**2-(X*P-C*C)))),
                integrated_occupation_change=float(weight@delta_n),max_occupation_change=float(max(abs(delta_n))),
                finite_pulses=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("coarse_checkpoint",type=Path)
    p.add_argument("fine_checkpoint",type=Path)
    p.add_argument("repository",type=Path)
    args=p.parse_args()
    root=args.repository/"research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs"
    report={g:run(getattr(args,g+"_checkpoint"),root/g/"quiet/rolling_summary.json")
            for g in ("coarse","fine")}
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
