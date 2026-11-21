"""Check a necessary partial-wave unitarity condition for the local crossed-vector benchmark."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


MASS = 8300.
G_SQUARED = 12.*np.pi*.30


def partial_magnitude(s):
    """Magnitude of the real diagonal J=0 wave in this helicity convention, m=1."""
    s = np.asarray(s, dtype=float)
    if np.any(s < 4.):
        raise ValueError("Pair energy is below threshold")
    r = (s-4.)/MASS**2
    series = r < 1.e-4
    # Integrals of 1 and u against 1/(1+r u), for 0 <= u <= 1.
    safe = np.where(series, 1., r)
    h = np.where(series, 1.+r*(-.5+r*(1./3.+r*(-.25+r*.2))), np.log1p(safe)/safe)
    j = np.where(series, .5+r*(-1./3.+r*(.25+r*(-.2+r/6.))), (1.-h)/safe)
    beta = np.sqrt((s-4.)/s)
    return beta*G_SQUARED/(16.*np.pi*MASS**2)*((2.*s-4.)*h-4.*j)


def spinor(p, pauli, anti, mass=1.):
    normal = np.sqrt(p[0]+mass)
    px,py,pz = p[1:]
    lower = np.array([[pz,px-1j*py],[px+1j*py,-pz]])@pauli/normal
    large = normal*pauli
    return np.r_[lower,large] if anti else np.r_[large,lower]


def spin_numerators(s, z, mass=1.):
    # Both external pair states transform with zero net helicity about the
    # collision axis. Rotate the outgoing pair together about y.
    momentum = np.sqrt(s/4.-mass**2)
    energy = np.sqrt(s)/2.
    half = np.sqrt(np.array([1.+z,1.-z])/2.)
    start = np.array([0.,0.,momentum])
    end = momentum*np.array([2.*half.prod(),0.,z])
    pauli = np.array([1.,0.])
    u1,v2,u3,v4 = [spinor(np.r_[energy,p],spin,anti,mass) for p,spin,anti in
                   ((start,pauli,False),(-start,pauli,True),(end,half,False),(-end,half,True))]
    identity = np.eye(2)
    gamma0 = np.diag([1.,1.,-1.,-1.])
    paulis = (np.array([[0.,1.],[1.,0.]]), np.array([[0.,-1j],[1j,0.]]), np.diag([1.,-1.]))
    zeros = np.zeros_like(identity)
    gamma = np.array([gamma0]+[np.block([[zeros,k],[-k,zeros]]) for k in paulis])
    metric = np.array([1.,-1.,-1.,-1.])
    def current(a,b):
        return (a.conj()@gamma0)@gamma@b
    u2 = spinor(np.r_[energy,-start],np.array([0.,1.]),False,mass)
    u4 = spinor(np.r_[energy,-end],np.array([-half[1],half[0]]),False,mass)
    return (np.sum(current(u3,u1)*metric*current(v2,v4)),
            np.sum(current(v2,u1)*metric*current(u3,v4)),
            np.sum(current(u3,u1)*metric*current(u4,u2)),
            np.sum(current(u4,u1)*metric*current(u3,u2)))


def relative_minimum_change(a):
    # Distance from a real amplitude to the disk |A-i/2|<=1/2, divided
    # by its magnitude. Rationalized to stay accurate for small amplitudes.
    return abs(a)/(np.hypot(a,.5)+.5)


def controls():
    error = []
    for s in (4.001,18.,MASS**2,50.*MASS**2):
        for mass in (0.,1.):
            for z in (-.999,-.4,0.,.6,.9999):
                unlike_t, unlike_s, like_t, like_u = spin_numerators(s,z,mass)
                analytic_t = 2.*s-6.*mass**2+2.*mass**2*z
                analytic_s = -4.*mass**2*z
                analytic_u_identical = -(2.*s-6.*mass**2-2.*mass**2*z)
                error.extend(abs(exact-prediction)/s for exact,prediction in (
                    (unlike_t,analytic_t),(unlike_s,analytic_s),
                    (like_t,analytic_t),(like_u,analytic_u_identical)))
    if max(error)>3.e-14:
        raise ValueError("Dirac contraction disagrees with the helicity numerators")
    integration_error = []
    for s in (4.001,16.,1.e4,.1*MASS**2,MASS**2,100.*MASS**2,1.e4*MASS**2):
        integral = quad(lambda z:(2.*s-6.+2.*z)/(MASS**2+.5*(s-4.)*(1.-z)),
                        -1.,1.,epsabs=1.e-15,epsrel=5.e-12,limit=200)[0]
        expected = np.sqrt((s-4.)/s)*G_SQUARED*integral/(32.*np.pi)
        integration_error.append(abs(float(partial_magnitude(s))/expected-1.))
    if max(integration_error)>2.e-10:
        raise ValueError("Angular integral disagrees with the closed form")
    s0 = brentq(lambda s: float(partial_magnitude(s))-.5,MASS**2,4.*MASS**2)
    for val in (.01,.5,1.,2.):
        closest = .5j + .5*(val-.5j)/abs(val-.5j)
        if abs(abs(val-closest)/val-relative_minimum_change(val))>1.e-13:
            raise ValueError("Argand disk projection failed")
    return s0, dict(max_normalized_spinor_error=max(error),
                    max_quadrature_relative_error=max(integration_error))


def measure(path, onset):
    with np.load(path,allow_pickle=False) as data:
        s = data["s"]
        incident = data["a"]*data["E_in"].sum(axis=1)
        invalid = s>onset
        amplitude = partial_magnitude(s)
        correction = relative_minimum_change(amplitude)
        changes = []
        for minimum in (.1,.5,.75):
            selected = correction>minimum
            changes.append(dict(minimum=minimum, fraction_events=float(np.mean(selected)),
                           incident_energy_fraction=float(incident[selected].sum()/incident.sum())))
        return dict(label=str(data["label"]),events=len(s),minimum_required_amplitude_changes=changes,
                    fraction_events_above_limit=float(np.mean(invalid)),
                    fraction_incident_conformal_energy_above_limit=float(incident[invalid].sum()/incident.sum()),
                    largest_abs_a0=float(partial_magnitude(s.max())))


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--events",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    args = cli.parse_args()
    onset,checks = controls()
    records = sorted(args.events.glob("*.npz"))
    if not records:
        raise ValueError("Supply the event tapes from the collision-photon study")
    out = dict(s_threshold=onset,root_s_over_M=np.sqrt(onset)/MASS,
               massless_ratio=np.sqrt(np.expm1(4.*np.pi/G_SQUARED)),
               samples=[dict(root_s_over_M=x,abs_a0=float(partial_magnitude((MASS*x)**2)),
                             relative_minimum_change=float(relative_minimum_change(partial_magnitude((MASS*x)**2))))
                        for x in (.5,1.,1.5,2.,3.,7.,10.)],
               checks=checks,runs=[measure(p,onset) for p in records])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
