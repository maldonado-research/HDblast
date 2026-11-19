"""Linear trace-coupling widths at the two archived interior minima."""
import argparse
import json
from math import atanh, ceil, cos, cosh, exp, log, pi, sin, sqrt, tanh
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


def phase_factor(x):
    if x < 0:
        raise ValueError("A squared mass ratio cannot be negative")
    return sqrt(max(0.,1.-4*x))


def scalar_ratio(x, xi):
    return (1-6*xi+2*x)**2 * phase_factor(x)


def fermion_ratio(x):
    return 4*x * phase_factor(x)**3


def scalar_maximum(xi):
    # Include endpoints and the interior peak of [(1-6ξ)+2x]^2 sqrt(1-4x).
    choices=[0.,.25]
    root=(1+6*xi)/10
    if 0.<root<.25:
        choices.append(root)
    x=max(choices,key=lambda r:scalar_ratio(r,xi))
    return x,scalar_ratio(x,xi)


def check_spin_and_phase_space():
    """Independent slash trace and direct radial two-body delta-function Jacobian, m_sigma=1."""
    i2=np.eye(2)
    zero=np.zeros((2,2))
    pauli=[np.array([[0.,1.],[1.,0.]]),np.array([[0.,-1j],[1j,0.]]),np.diag([1.,-1.])]
    g0=np.block([[i2,zero],[zero,-i2]])
    gamma=np.array([np.block([[zero,p],[-p,zero]]) for p in pauli])
    eye=np.eye(4)
    nodes,weights=np.polynomial.legendre.leggauss(8)
    width_unit=1/(32*pi)
    worst=0.
    for x in (0.,0.01,.1,.2,.249):
        beta=phase_factor(x)
        p=beta/2
        e=.5
        # ∫d³p/(2π)^3 2E ∫d³k/(2π)^3 2E (2π)^4 δ³ δ(E1+E2−M)
        radial=(4*pi)/(2*pi)**2 * p**2/((2*e)**2) / (2*p/e)
        np.testing.assert_allclose(radial,beta/(8*pi),rtol=1e-14)
        spin_sum=0.
        for u,wu in zip(nodes,weights):
            for phi in np.linspace(0.,2*pi,5,endpoint=False):
                v=p*np.array([sqrt(1-u*u)*cos(phi),sqrt(1-u*u)*sin(phi),u])
                gdot=np.tensordot(v,gamma,axes=1)
                p1=e*g0-gdot
                p2=e*g0+gdot
                # Spin-summed squared amplitude in units α²/Mpl².
                spin_sum+=wu/10*x*np.trace((p1+sqrt(x)*eye)@(p2-sqrt(x)*eye)).real
        fermion=radial*spin_sum/2/width_unit
        worst=max(worst,abs(fermion-fermion_ratio(x)))
        for xi in (0.,1/6,.25):
            # Direct trace using p1·p2 and the improvement's total momentum.
            dot=e*e+p*p
            amp=2*dot+4*x-6*xi  # m_sigma=1 and mostly-minus invariant dot
            physical=radial*amp**2/4/width_unit  # includes identical-particle 1/2!
            worst=max(worst,abs(physical-scalar_ratio(x,xi)))
    if worst>8e-14:
        raise ValueError("Rest-frame phase-space amplitude check")
    return float(worst)


def warp_integral(y):
    def primitive(v):
        return v+2*log(cosh(v))/3+tanh(v)**2/6
    norm=primitive(y)
    # For u>70 and either branch the omitted exponentially falling tail is <1e-32.
    return quad(lambda u:exp((-2/3)*(primitive(y+u)-norm)),0.,70.,epsabs=3e-14,epsrel=3e-14)[0]


def geometry(phi):
    """Independent quadrature at the nominal minimum; f0 and location remain unrounded."""
    y=-atanh(phi)
    f0=2*warp_integral(0.)
    f=2*warp_integral(y)
    W=1-phi+phi**3/3
    D=1-W*f/3
    z=2*W*D/f+6*D**2/f**2
    # A=1/sqrt(f); dTheta/dy=+sqrt(z) on the tanh branch.
    alpha=D/(f*sqrt(z))
    step=1.e-4
    def log_A(u):
        return -.5*log(2*warp_integral(u))
    dA=(-log_A(y+2*step)+8*log_A(y+step)-8*log_A(y-step)+log_A(y-2*step))/(12*step)
    if not np.isclose(dA/sqrt(z),alpha,rtol=1e-9):
        raise ValueError("Independent logarithmic matter-metric derivative")
    b=4/f0-4/3+2/phi
    d2V=2*(1+b*phi)/phi**2
    r=f0*sqrt(3*d2V*(1-phi**2)**2/(f**2*z))
    return dict(f_at_min=f,alpha=alpha,m_over_Htop=r,Htop_over_Mpl_per_sqrt_t=1/(sqrt(3)*f0),
                scalar_width_per_Htop_over_Mpl_squared=alpha**2*r**3/(32*pi))


def run(ledger):
    data=json.loads(ledger.read_text())
    channels={}
    for name,xi in (("minimal_scalar",0.),("conformal_scalar",1/6.)):
        x,ratio=scalar_maximum(xi)
        check=minimize_scalar(lambda z:-scalar_ratio(z,xi),bounds=(0.,.25),method="bounded",
                              options={"xatol":1e-12})
        if abs(x-check.x)>3e-8 or abs(ratio+check.fun)>2e-14:
            raise ValueError("Scalar mass optimum")
        channels[name]=dict(xi=xi,massless_ratio=scalar_ratio(0.,xi),x_at_max=x,maximum_ratio=ratio)
    control=minimize_scalar(lambda x:-fermion_ratio(x),bounds=(0.,.25),method="bounded",
                            options={"xatol":1e-12})
    if abs(control.x-.1)>3e-8:
        raise ValueError("Dirac mass optimum")
    channels["dirac"]=dict(massless_ratio=fermion_ratio(0.),x_at_max=.1,maximum_ratio=fermion_ratio(.1))
    if any((scalar_ratio(0.,1/6),fermion_ratio(0.),scalar_ratio(.25,0.),fermion_ratio(.251))):
        raise ValueError("Conformal zero or closed threshold")
    results={}
    for key in ("phi_m=0.5","phi_m=-0.5"):
        row=data["B"][key]
        alpha=row["alpha"]
        r=row["m_over_H_top"]
        old=alpha**2*r**3
        c=old/(32*pi)
        exact=geometry(row["phi_m"])
        if abs(old/row["Gamma_over_Htop_per_HoverMpl_squared"]-1)>5e-14:
            raise ValueError("Historical ledger units")
        if abs(exact["scalar_width_per_Htop_over_Mpl_squared"]/c-1)>0.002:
            raise ValueError("Archive and nominal exact-geometry width")
        scales={}
        for detuning,point in row["scales"].items():
            h=point["H_top_over_Mpl"]
            ratio=c*h*h
            if abs(32*pi*ratio/point["Gamma_over_H_top"]-1)>5e-14:
                raise ValueError("Original benchmark scale")
            scales[detuning]=dict(Htop_over_Mpl=h,m_GeV=point["m_GeV"],
                                  Gamma_minimal_real_scalar_over_Htop=ratio,
                                  Gamma_minimal_real_scalar_GeV=ratio*point["H_top_GeV"],
                                  Gamma_best_colorless_Dirac_over_Htop=ratio*fermion_ratio(.1))
        results[key]=dict(alpha=alpha,f_m=row["f_m"],m_over_Htop=r,historical_scaling_coefficient=old,
                          scalar_width_per_Htop_over_Mpl_squared=c,
                          best_Dirac_width_per_Htop_over_Mpl_squared=c*fermion_ratio(.1),
                          scale_points=scales,exact_minimum=exact,
                          exact_vs_archived_scalar_width_fraction=exact["scalar_width_per_Htop_over_Mpl_squared"]/c-1)
    return dict(convention="Reduced Mpl=2.435e18 GeV; x=m_daughter^2/m_sigma^2; each width in rest frame",
                channels_in_units_of_one_massless_minimal_real_scalar=channels,
                classical_massless_gauge_ratio=0.,scalar_to_original_ledger_ratio=1/(32*pi),
                minimum_minimal_massless_species_to_match_original=ceil(32*pi),
                direct_trace_and_phase_space_max_abs_error=check_spin_and_phase_space(),
                branches=results)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(run(args.ledger),indent=2,allow_nan=False)+"\n")
