"""Match a compact scalar pulse to quiet B1 collars using the planar mass constraint.

Run with the archived B1 summary path as argument. Requires NumPy and SciPy.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp
from scipy.optimize import brentq

SUMMARY_SHA256 = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"
CENTER, WIDTH = -0.12, 0.02
PATCH_LEFT = -0.20
PULSE_LEFT, PULSE_RIGHT = CENTER - WIDTH, CENTER + WIDTH


def potential(phi):
    w = 1 - phi + phi**3/3
    wp = phi**2 - 1
    return wp**2/2-2*w**2/3


def potential_prime(phi):
    w = 1 - phi + phi**3/3
    return 2*(phi**2-1)*(phi - 2*w/3)


def potential_difference(phi, h):
    """U(phi+h)-U(phi) in factored form, including its h=0 limit."""
    wp = phi**2 - 1
    dw = h*(wp+phi*h+h*h/3)
    dwp = h*(2*phi+h)
    w = 1 - phi + phi**3/3
    return dwp*(2*wp+dwp)/2 - (2/3)*dw*(2*w+dw)


def bump(z):
    x = (z-CENTER)/WIDTH
    if abs(x) >= 1:
        return 0.,0.
    d = 1-x*x
    f = math.exp(1-1/d)
    return f, -2*x*f/(WIDTH*d*d)


class MatchedPulse:
    """The finite identity-chart patch A=A0, A_T=1, dc=0."""
    def __init__(self, source, max_step=0.0005, method="DOP853", rtol=2e-12):
        payload = Path(source).read_bytes()
        if hashlib.sha256(payload).hexdigest() != SUMMARY_SHA256:
            raise ValueError("Unexpected reference summary")
        cfg = json.loads(payload)
        self.archived_radius = cfg["rho_b"]
        phi = cfg["phi_b_static"]
        delta, c, d = (cfg["params"][k] for k in ("delta","c","d"))
        w = 1 - phi + phi**3/3
        sigma = 2*w+delta*(1+c*phi+d*phi*phi/2)
        sigmap = 2*(phi**2-1)+delta*(c+d*phi)
        self.radius = 1/math.sqrt(sigma**2/36 - sigmap**2/48 + potential(phi)/6)
        initial = [math.log(self.radius), self.radius*sigma/6, phi, -self.radius*sigmap/2]
        self.solve_options = dict(method=method, max_step=max_step, rtol=rtol,
                                  atol=rtol*0.0005, dense_output=True)

        def static_rhs(z, state):
            b, q, phi, r = state
            return (q, 3-3*q*q-(2/3)*math.exp(2*b)*potential(phi), r,
                    math.exp(2*b)*potential_prime(phi)-3*q*r)

        static = solve_ivp(static_rhs,(0,PATCH_LEFT),initial,**self.solve_options)
        if not static.success:
            raise RuntimeError(static.message)
        self.static = static.sol

    def mass_solution(self, amplitude):
        """Return delta-mass and its amplitude derivative, starting on the quiet wall collar."""
        def rhs(z, state):
            delta_mass, sensitivity = state
            b,q,phi,r = self.static(z)
            f,df = bump(z)
            d = amplitude*df
            h = amplitude*f
            radius2 = math.exp(2*b)
            reference_mass = radius2*(1-q*q)
            a0 = q*r*r/(3*(1-q*q))
            da = -2*d*(r+d)/(3*(q-1))
            da_deta = -2*df*(r+2*d)/(3*(q-1))
            volume_factor = (2/3)*radius2**2*q
            forcing = da*reference_mass+volume_factor*potential_difference(phi,h)
            dforcing = da_deta*(reference_mass+delta_mass)+volume_factor*f*potential_prime(phi+h)
            return ((a0+da)*delta_mass+forcing, (a0+da)*sensitivity+dforcing)

        result = solve_ivp(rhs,(PULSE_RIGHT,PULSE_LEFT),[0.,0.],**self.solve_options)
        if not result.success:
            raise RuntimeError(result.message)
        return result.sol

    def mismatch(self, amplitude):
        return float(self.mass_solution(amplitude)(PULSE_LEFT)[0])

    def root(self):
        return brentq(self.mismatch,0.001,0.003,xtol=2e-17,rtol=1e-13)

    def fields(self, z, amplitude, mass_solution):
        """Point values on glued support after shooting; exact quiet data on both collars."""
        if not PATCH_LEFT<=z<=0:
            raise ValueError("Exterior points use the original quiet solution")
        b0,q,phi0,r0 = self.static(z)
        if PULSE_LEFT < z < PULSE_RIGHT:
            f,df=bump(z)
            dm = float(mass_solution(z)[0])
        else:
            f,df,dm = 0.,0.,0.
        h,d = amplitude*f,amplitude*df
        p,r = -d, r0+d
        mu0 = math.exp(2*b0)*(1-q*q)
        db = -0.5*math.log1p(dm/mu0)
        b=b0+db
        du=potential_difference(phi0,h)
        d_e2u = math.exp(2*b0)*(du+math.expm1(2*db)*(potential(phi0)+du))
        db_prime = (2*d_e2u+2*d*(r0+d)-2*p*r/q)/(6*(q-1/q))
        bt = (p*r/3-db_prime)/q
        return dict(A=b0,Q=q,u=1.,B=b,Bz=q+db_prime,BT=bt,phi=phi0+h,
                    p=p,r=r,delta_mass=dm,delta_B=db)

    def checks(self, amplitude, dense):
        results=[]
        for n in (401,801,1601):
            z = np.linspace(PATCH_LEFT,0,n)
            h = z[1]-z[0]
            data=[self.fields(x,amplitude,dense) for x in z]
            b,az,bt,phi,p = (np.array([x[k] for x in data]) for k in ("B","Q","BT","phi","p"))
            def d1(arr):
                return (arr[:-4]-8*arr[1:-3]+8*arr[3:-1]-arr[4:])/(12*h)
            a_zz,b_z,r=d1(az),d1(b),d1(phi)
            mid=slice(2,-2)
            b,az,bt,phi,p=(v[mid] for v in (b,az,bt,phi,p))
            h_terms = (-2*np.exp(2*b)*potential(phi),np.ones_like(p)*6,
                       6*bt, -12*az*az, 6*az*b_z, -6*a_zz, -p*p, -r*r)
            m_terms = (-3*az,3*b_z,3*az*bt,-p*r)
            results.append(dict(nodes=n,
               max_scaled_h=float(max(abs(sum(h_terms))/sum(abs(t) for t in h_terms))),
               max_scaled_m=float(max(abs(sum(m_terms))/sum(abs(t) for t in m_terms)))))
        return results

    def radial_mass_terms(self, amplitude, dense):
        """Quiet-subtracted radial stress integral per normalized planar orbit."""
        def contributions(z):
            b0, q, phi0, r0 = self.static(z)
            field = self.fields(z,amplitude,dense)
            b,p,r = (field[key] for key in ("B","p","r"))
            weight = (2/3)*math.exp(4*b0)
            potential_term = q*potential_difference(phi0, amplitude*bump(z)[0])
            kinetic_term = q*((p*p+r*r)*math.exp(-2*b)-r0*r0*math.exp(-2*b0))/2
            mixed_term = -p*r*math.exp(-2*b)
            return weight*np.array([potential_term,kinetic_term,mixed_term])

        parts,_ = quad_vec(contributions,PULSE_LEFT,PULSE_RIGHT,epsabs=1e-12,epsrel=3e-12)
        total=math.fsum(parts)
        end_difference = -float(dense(PULSE_LEFT)[0])
        return dict(integration="increasing Z from the inner edge to the wall-facing edge",
                    potential=float(parts[0]), kinetic_and_gradient=float(parts[1]),
                    mixed_momentum=float(parts[2]), sum=total,
                    endpoint_delta_mass_difference=end_difference,
                    sum_minus_endpoint=total-end_difference)

    def report(self):
        eta = self.root()
        sol = self.mass_solution(eta)
        z = np.linspace(PULSE_LEFT,PULSE_RIGHT,501)
        values = [self.fields(x,eta,sol) for x in z]
        x_values=[]
        for position in (CENTER-WIDTH/2,CENTER,CENTER+WIDTH/2):
            data=self.fields(position,eta,sol)
            same_phi=brentq(lambda t:self.static(t)[2]-data["phi"],PATCH_LEFT,0,xtol=2e-15)
            b,_,_,r=self.static(same_phi)
            x_quiet=math.exp(-2*b)*r*r
            x_pulse=math.exp(-2*data["B"])*(data["r"]**2-data["p"]**2)
            x_values.append(dict(z=position,phi=float(data["phi"]),quiet_z_at_equal_phi=same_phi,
                                 X=x_pulse,quiet_X=float(x_quiet),difference=float(x_pulse-x_quiet)))
        b,q,_,r=self.static(np.linspace(PATCH_LEFT,0,501))
        return dict(source_sha256=SUMMARY_SHA256, patch=[PATCH_LEFT,0],
            support=[PULSE_LEFT,PULSE_RIGHT], analytic_wall_radius=self.radius,
            archived_wall_radius=self.archived_radius, min_sampled_Q=float(min(q)),
            min_sampled_static_phi_prime=float(min(r)), shooting_amplitude=eta,
            endpoint_mass_defects={str(v):self.mismatch(v) for v in (0.,.001,.003,eta)},
            endpoint_slopes={str(v):float(self.mass_solution(v)(PULSE_LEFT)[1]) for v in (0.,eta)},
            max_sampled_delta_mass=float(max(v["delta_mass"] for v in values)),
            max_sampled_abs_delta_B=float(max(abs(v["delta_B"]) for v in values)),
            max_sampled_abs_BT=float(max(abs(v["BT"]) for v in values)),
            equal_phi_invariants=x_values, radial_mass_terms=self.radial_mass_terms(eta,sol),
            constraint_checks=self.checks(eta,sol))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_summary",type=Path)
    parser.add_argument("--max-step",type=float,default=.0005)
    parser.add_argument("--method",choices=["DOP853","RK45"],default="DOP853")
    opts=parser.parse_args()
    patch=MatchedPulse(opts.source_summary,opts.max_step,opts.method)
    print(json.dumps(patch.report(),indent=2,allow_nan=False))
