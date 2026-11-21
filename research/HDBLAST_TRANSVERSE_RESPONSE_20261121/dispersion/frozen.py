"""Frozen, local, cold-cohort filamentation spectrum on the recorded early wall."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import brentq
from scipy.special import roots_legendre


GAUGE = 2.54
COMOVING_WAVENUMBERS = (32, 128, 512, 2048, 8192, 32768)


def intervals(trace, t, injection=0.):
    """Birth weights and the two axial momenta bounding each linear birth cell."""
    ts = trace["t"]
    end = np.searchsorted(ts, t, side="left")
    a = np.interp(t, ts, trace["a"])
    imp = np.interp(t, ts, trace["I"])
    values = {}
    for key in ("D", "I", "a"):
        values[key] = np.append(trace[key][:end], np.interp(t, ts, trace[key]))
    rest = np.diff(values["D"])/a**3
    momenta = (imp-values["I"] + injection*values["a"])/a
    alive = rest > 0.
    return rest[alive], momenta[:-1][alive], momenta[1:][alive]


def kernel_average(p, q, primitive, value):
    """Average a function over a uniform interval of wall-frame axial momentum."""
    dp = q-p
    short = abs(dp) < 1.e-5 * np.maximum(1., np.maximum(abs(p), abs(q)))
    out = value((p+q)/2.)
    regular = ~short
    out[regular] = (primitive(q[regular])-primitive(p[regular]))/dp[regular]
    return out


def cold_moments(trace, time, charge, injection=0.):
    mass, p, q = intervals(trace, time, injection)
    gamma = lambda z: np.hypot(1., z)
    v = lambda z: z/gamma(z)
    # Rationalize (v(q)-v(p))/(q-p) to keep the ultrarelativistic tail.
    gp, gq = gamma(p), gamma(q)
    inverse3 = np.full_like(p, 1.)
    same = p*q >= 0.
    inverse3[same] = (1. + (1.+p[same]**2+q[same]**2)/(gp[same]*gq[same]+p[same]*q[same]))/(gp[same]*gq[same]*(gp[same]+gq[same]))
    inverse3[~same] = (v(q[~same])-v(p[~same]))/(q[~same]-p[~same])
    vel = (q+p)/(gq+gp)
    stress = kernel_average(p, q, lambda z: .5*(z*gamma(z)-np.arcsinh(z)), lambda z: z*v(z))
    gamma_bound = kernel_average(p, q, lambda z: np.arcsinh(z)-v(z), lambda z: v(z)**2/gamma(z))
    return {"rho": float(np.sum(mass)), "current": float(charge*(mass @ vel)),
            "parallel_pressure": float(mass @ stress), "parallel_plasma_squared":float(charge**2/GAUGE*(mass @ inverse3)),
            "growth_ceiling": float(np.sqrt(max(0.,charge**2/GAUGE*(mass @ gamma_bound))))}


def quadrature(trace, t, order, charge, injection=0.):
    mass, p, q = intervals(trace, t, injection)
    u, w = roots_legendre(order)
    r, z = np.arcsinh(p), np.arcsinh(q)
    # Integrating dp/gamma = dr resolves the youngest particles even when the
    # electric impulse changes by many rest masses across one recorded step.
    rapidity = .5*((z+r)[:,None] + (z-r)[:,None]*u)
    span = q-p
    scale = np.divide(z-r, span, out=1./np.hypot(1.,p), where=span!=0.)
    alpha = (charge**2/(2.*GAUGE) * (mass*scale)[:,None] * w).ravel()
    v = np.tanh(rapidity).ravel()
    inv_gamma = (1./np.cosh(rapidity)).ravel()
    B = np.interp(t, trace["t"], trace["B"])
    return alpha, charge*B*inv_gamma, v, inv_gamma**2


def eigen_growth(alpha, omega, v, inverse2, k):
    """Schur complement of a real symmetric frozen stiffness matrix."""
    a2 = float(alpha @ inverse2)
    omega2 = omega**2
    kv = k*v
    def scalar(s):
        den = 1./(s+omega2)
        yy = s+k*k + s*(alpha @ den)
        yz = alpha @ (omega*kv*den)
        zz = s+k*k+a2-alpha @ (kv**2*den)
        return zz - yz*yz/yy
    upper = float(alpha @ (v*v))
    if upper <= 0.: return 0.
    tiny = max(1.e-22, upper*1.e-17)
    if scalar(tiny)>=0.:return 0.
    high = max(2.*upper,tiny*2.)
    return float(np.sqrt(brentq(scalar,tiny,high,xtol=1.e-15,rtol=1.e-13)))


def snapshot(trace, t, charge, order):
    a = float(np.interp(t,trace["t"],trace["a"]))
    B = float(np.interp(t,trace["t"],trace["B"]))
    moments = cold_moments(trace,t,charge)
    p = moments["parallel_pressure"]
    om = moments["parallel_plasma_squared"]
    current = moments["current"]
    alpha, omega, v, g2 = quadrature(trace,t,order,charge)
    roots = {}
    for ka in COMOVING_WAVENUMBERS:
        k = ka/a
        sign = k*k*(GAUGE*B**2-p) + GAUGE*B**2*om - current**2/GAUGE
        rate = eigen_growth(alpha,omega,v,g2,k)
        if (rate>0.) != (sign<0.):
            raise ValueError("Spectral and static sign tests disagree")
        roots[str(ka)] = rate
    # Independently replace the full longitudinal velocity distribution by a
    # single velocity that gives its exact wall-frame current at fixed density.
    mean_velocity = current/(charge*moments["rho"])
    inverse_gamma = np.sqrt((1.-mean_velocity)*(1.+mean_velocity))
    stream = np.array([mean_velocity])
    inertia = np.array([charge**2*moments["rho"]*inverse_gamma/GAUGE])
    guide = np.array([charge*B*inverse_gamma])
    inverse_sq = np.array([inverse_gamma**2])
    current_matched = {str(ka):eigen_growth(inertia,guide,stream,inverse_sq,ka/a)
                       for ka in COMOVING_WAVENUMBERS}
    before = trace["t"] < t
    times = np.append(trace["t"][before],t)
    radii = np.append(trace["a"][before],a)
    conformal_age = float(trapezoid(1./radii,times))
    err = max(abs(float(alpha@g2)/om-1.),abs(np.sqrt(float(alpha@(v*v)))/moments["growth_ceiling"]-1.))
    if order>=16 and err>1.e-7:raise ValueError("Rapidity quadrature failed the primitive-moment checks")
    return dict(**moments, parallel_pressure_over_guide_tension=p/(GAUGE*B**2),
                charge_current_frequency_squared_over_parallel_plasma_squared=(current/(GAUGE*B))**2/om,
                B=B,a=a,t=t,rates_for_comoving_k=roots,
                current_matched_monoenergetic_gamma=float(1./inverse_gamma),
                current_matched_monoenergetic_rates=current_matched,conformal_age=conformal_age,
                comoving_k_for_one_wavelength_in_light_distance=2.*np.pi/conformal_age,rate_quadrature_moment_check=err)


def growth_clock(trace,t,charge,n,injection=0.):
    times=np.linspace(trace["t"][0],t,n+1)
    rates=np.array([cold_moments(trace,x,charge,injection)["growth_ceiling"] for x in times])
    return float(trapezoid(rates,times))


def controls():
    rng=np.random.default_rng(12891)
    err=0.
    for _ in range(75):
        alpha=np.exp(rng.uniform(-2.,4.,size=7))
        p=rng.normal(0.,8.,size=7)
        g=1./np.hypot(1.,p)
        v=p*g
        om=rng.uniform(.01,4.)*g
        k=float(np.exp(rng.uniform(-1.,4.)))
        mat=np.diag(np.r_[k*k+sum(alpha),k*k+alpha@(g*g),om**2])
        mat[0,2:]=mat[2:,0]=np.sqrt(alpha)*om
        mat[1,2:]=mat[2:,1]=np.sqrt(alpha)*k*v
        ref=float(np.sqrt(max(0.,-np.linalg.eigvalsh(mat)[0])))
        test=eigen_growth(alpha,om,v,g*g,k)
        err=max(err,abs(ref-test))
        if ref>np.sqrt(alpha@(v*v))+1.e-9: raise ValueError("Frozen growth bound failed")
    if err>2.e-7:raise ValueError("Direct eigenvalues disagree with Schur root")
    return {"direct_stiffness_cases":75,"largest_absolute_growth_difference":err}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--history",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    meta=json.loads((args.history/"markers.json").read_text())
    charge=meta["charge_to_mass"]
    runs={}
    for name, ref in meta["runs"].items():
        runs[name]={}
        with np.load(args.history/(name+".npz")) as raw:
            trace={k:raw[k] for k in ("t","D","I","a","B")}
        for fraction,marker in ref["markers"].items():
            t=marker["time"]
            low=snapshot(trace,t,charge,8)
            row=snapshot(trace,t,charge,16)
            row["rate_order_refinement_max_absolute"]=max(abs(value-low["rates_for_comoving_k"][k]) for k,value in row["rates_for_comoving_k"].items())
            row["instantaneous_rate_ceiling_integral"]=growth_clock(trace,t,charge,256)
            row["interval_integration_absolute_change"]=abs(row["instantaneous_rate_ceiling_integral"]-growth_clock(trace,t,charge,128))
            if fraction=="0.01":
                row["injection_momentum_sensitivity"]={str(p):growth_clock(trace,t,charge,256,injection=p) for p in [-.32015621187164256,.32015621187164256]}
            runs[name][fraction]=row
    result={"model":"Frozen local mode, axial ballistic cohorts, birth momentum zero unless varied, inherited full-current normalization",
            "physical_k":"comoving k divided by the instantaneous radius", "controls":controls(),"runs":runs}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
