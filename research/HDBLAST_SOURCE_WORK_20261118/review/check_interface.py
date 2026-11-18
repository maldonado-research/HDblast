"""Check the band action derivatives, junction conversion and instantaneous loss generator."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


def module(file):
    spec = importlib.util.spec_from_file_location("sector", file)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def five(fn, step=1e-5):
    return (8 * (fn(step) - fn(-step)) - (fn(2*step) - fn(-2*step))) / (12*step)


def review(path):
    sector = module(path)
    rng = np.random.default_rng(19840504)
    defects = dict(action_pressure=0., action_force=0., native_derivatives=0., ward=0., loss=0.,
                   ground_derivative=0., symmetric_point_force=0., physical_number_units=0.,
                   damping_at_mass_zero=0., damping_at_zero_coupling=0., inward_pure_covariance_defect=0.,
                   native_unit_conversion=0., passive_junctions=0., passive_matter_rhs=0.)
    def measure(key, actual, reference):
        scaled = np.max(abs(np.array(actual)-np.array(reference)) / np.maximum(1.,np.abs(reference)))
        defects[key] = float(max(defects[key], scaled))
    for i in range(48):
        rb = 7.84
        shell = sector.QuantumShell(rb, 3.3, 32, True, i%2)
        A = rng.uniform(2.7, 3.7)
        scale = np.exp(A-shell.ref)
        vol = scale**3
        phi = .5 if i%3==0 else rng.uniform(.02, .91)
        AT, phiT = rng.uniform(-.2, 2.), rng.uniform(-1., 1.)
        eB = rb*rng.uniform(.5, 2.)
        lapse = eB/rb
        m = shell.prepare(A, phi)
        x, pp, c = shell.view(m)
        nbar = rng.uniform(0., .8, len(x)) if i%4>=2 else np.zeros(len(x))
        # Apply a real symplectic squeeze and rotation to thermal covariances.
        scale_k = vol*np.hypot(shell.k/scale,100*(phi-.5))
        rr = rng.uniform(-.6, .6, len(x))
        t = rng.uniform(0, 2*np.pi, len(x))
        var1 = (nbar+.5)*np.exp(2*rr)
        var2 = (nbar+.5)*np.exp(-2*rr)
        x[:] = (var1*np.cos(t)**2 + var2*np.sin(t)**2)/scale_k
        pp[:] = (var1*np.sin(t)**2 + var2*np.cos(t)**2)*scale_k
        c[:] = (var1-var2)*np.sin(t)*np.cos(t)
        m[0] = rng.uniform(0.,3.)
        rho, pr, force = shell.quantities(m, A, phi)[0]
        extent = lambda u, v: np.exp(3*(u-shell.ref))*shell.quantities(m,u,v)[0][0]
        measure("action_pressure", five(lambda s:extent(A+s,phi)), -3*vol*pr)
        measure("action_force", five(lambda s:extent(A,phi+s),step=1e-6), vol*force)
        dm = shell.rhs(m, A, phi, AT, phiT, eB)
        measure("native_unit_conversion", np.array(shell.totals(m, A, phi))*rb/6.49e-6,
                (rho+m[0], pr+m[0]/3, force))
        passive = sector.QuantumShell(rb, shell.ref, 32, False, shell.decay)
        measure("passive_junctions", passive.totals(m,A,phi), np.zeros(3))
        measure("passive_junctions", passive.totals_T(m,dm,A,phi,AT,phiT,eB), np.zeros(3))
        measure("passive_matter_rhs", passive.rhs(m,A,phi,AT,phiT,eB), dm)
        deriv = np.array(shell.totals_T(m,dm,A,phi,AT,phiT,eB))
        path_totals=lambda s:np.array(shell.totals(m+s*dm,A+s*AT,phi+s*phiT))
        measure("native_derivatives",five(path_totals,step=4e-7)/shell.mult,deriv/shell.mult)
        measure("ward",deriv[0]/shell.mult,
                phiT*force-3*AT*(rho+pr+4*m[0]/3))
        ex=shell.extras(m,A,phi,AT,phiT,eB)
        E=scale_k/vol
        occ=(scale_k*x+pp/scale_k-1)/2
        measure("physical_number_units",ex['n_hat'],(shell.weight@occ)/vol)
        gamma=(shell.decay*100*(phi-.5))**2/(8*np.pi*E)
        expected_loss=shell.weight@(gamma*E*occ)/vol
        measure("loss",(dm[0]+4*AT*m[0])/lapse, expected_loss)
        if phi==.5:
            measure("symmetric_point_force",force,0.)
            measure("damping_at_mass_zero",expected_loss,0.)
        if not shell.decay:measure("damping_at_zero_coupling",expected_loss,0.)
        if i%4<2:
            # Pure initial covariance. D(det V)/dT is nonnegative under this channel.
            xd, pd, cd = shell.view(dm)
            ddet = xd*pp+pd*x-2*c*cd
            measure("inward_pure_covariance_defect",np.minimum(ddet,0),0)
        prime = shell.prepare(A,phi)
        dp=shell.rhs(prime,A,phi,AT,phiT,eB)
        d0=shell.totals_T(prime,dp,A,phi,AT,phiT,eB)
        measure("ground_derivative",d0[0]/shell.mult,0.)
    if any(x>2e-5 for x in defects.values()):
        raise ValueError(defects)
    return {"deterministic_state_tests":48,"max_scaled_defects":defects}

if __name__=="__main__":
    print(json.dumps(review(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/"coupled/quantum_shell.py"),indent=2))
