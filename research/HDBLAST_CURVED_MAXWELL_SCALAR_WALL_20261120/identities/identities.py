"""Curved Maxwell triplet, spherical Einstein equations and moving-frame checks."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


def derivatives(func, t, step):
    a, b, c, d, e = [func(t+i*step) for i in (-2,-1,0,1,2)]
    return (a-8*b+8*d-e)/(12*step), (-a+16*b-30*c+16*d-e)/(12*step**2)


def cartesian_stress(c, amp):
    # On a flat four-dimensional spatial section take three constant two-forms.
    # Each is (2 amp/sqrt(3)) times a sum of two orthogonal unit wedges.
    metric = np.diag([-1]+[c*c]*4)
    inv = np.linalg.inv(metric)
    stress = np.zeros((5,5))
    for i, j, k in [(1,2,3),(2,3,1),(3,1,2)]:
        form = np.zeros((5,5))
        form[4,i], form[j,k] = 1., -1.
        form = (form-form.T) * (2*amp/math.sqrt(3))
        raised = inv @ form @ inv
        stress += form @ inv @ form.T - .25*metric*np.sum(raised*form)
    return stress


def calculate():
    ell, kap, k = 1.7, .8, 1.15
    C = 2*kap*k*k/3
    def scale(t):
        return math.sqrt(ell*math.sqrt(C)*math.sin(2*t/ell))

    errors = {}
    def record(key, value):
        errors[key] = max(errors.get(key, 0.), float(np.max(np.abs(value))))

    angles=[.32,.71,1.24,1.57,1.98,2.59,2.82]
    count = 0
    for angle in angles:
        t = angle*ell/2
        c = scale(t)
        cd, cdd = derivatives(scale, t, ell*3e-4)
        H=cd/c
        rho = 4*k*k/c**4
        T = cartesian_stress(c,k)
        record("cartesian_triplet_stress",T-np.diag([rho,0,0,0,0]))
        # Independent Einstein tensor of ds²=-dt²+c² δij dXi dXj, all five indices.
        einstein = np.diag([6*H**2+6/ell**2]+[c*c*(-3*cdd/c-3*H**2-6/ell**2)]*4)
        record("cartesian_einstein_equations",einstein-kap*T)
        for chi in [.12,.33,.61]:
            a=c*chi
            cen=c*cd
            cenn=c*cd*cd+c*c*cdd
            ru,rv=(cen*chi-c)/math.sqrt(2), (cen*chi+c)/math.sqrt(2)
            ruu,rvv=(cenn*chi-2*cen)/2,(cenn*chi+2*cen)/2
            ruv=cenn*chi/2
            magnitude=c*c
            A=k*chi**2
            au,av=-math.sqrt(2)*k*chi,math.sqrt(2)*k*chi
            auv=-k
            m=C*chi**4
            F=-2*ru*rv/magnitude
            traceL=cd/math.sqrt(2)
            record("geometric_mass",a*a*(1+a*a/ell**2-F)-m)
            record("double_null_constraints",
                   [ruu-2*traceL*ru+kap*au*au/(3*a), rvv-2*traceL*rv+kap*av*av/(3*a)])
            record("characteristic_maxwell",2*a*auv+ru*av+rv*au+4*magnitude*A/a)
            record("characteristic_radius",ruv+magnitude*(a/ell**2+m/a**3-2*kap*A*A/(3*a**3)))
            record("angular_einstein",c*cdd/2 - (
                magnitude*(-.5/ell**2+1.5*m/a**4-5*kap*A*A/(3*a**4))-kap*au*av/(6*a**2)))
            for ucomp, other, grad, sign in [(ru,rv,au,-1),(rv,ru,av,1)]:
                computed=(2*kap/3)*(2*A*A*ucomp/a-a*grad**2*other/magnitude)
                exact=sign*4*C*chi**3/math.sqrt(2)
                record("null_mass_flux",computed-exact)
            for velocity in [-.83,-.2,0.,.45,.8]:
                gamma=1/math.sqrt(1-velocity*velocity)
                X=2*k*gamma*velocity/c**2
                Y=2*k*gamma/c**2
                B=2*k/c**2
                density=(X**2+Y**2+B**2)/2
                pressure=(X**2+Y**2-B**2)/2
                q=-X*Y
                speed=gamma*(H*a+velocity)
                normal=gamma*(1+velocity*H*a)
                record("boosted_time_mass_flux",
                       -(2*kap/3)*a**3*(pressure*speed+q*normal)-4*C*chi**3*gamma*velocity/c)
                record("boosted_normal_mass_flux",
                       (2*kap/3)*a**3*(density*normal+q*speed)-4*C*chi**3*gamma/c)
                record("boosted_dust_pressure",pressure-rho*gamma**2*velocity**2)
                record("pressure_inequality_violation",max(0.,abs(q)-B**2/2-pressure))
                count+=1

    # A vacuum Schwarzschild-AdS exterior provides an independent static control.
    for r in [1.1,1.7,2.4,4.]:
        mass=.56
        f=lambda rr: 1+rr**2/ell**2-mass/rr**2
        df, ddf=derivatives(f,r,3e-4*r)
        record("vacuum_base_curvature",df/2-(r/ell**2+mass/r**3))
        record("vacuum_angular_curvature",-ddf/4-(-.5/ell**2+1.5*mass/r**4))
    if max(errors.values())>3e-7:
        raise ArithmeticError(errors)
    spec=importlib.util.spec_from_file_location('inward_collar',Path(__file__).with_name('collar.py'))
    collar=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(collar)
    return dict(collar=collar.calculate(),kappa5_squared=kap,ads_length=ell,field_amplitude=k,
                sampled_times=len(angles),sampled_radii=3,boosted_checks=count,
                max_absolute_residuals=errors,units="benchmark coordinate units",
                exact_bulk="homogeneous spatial Maxwell triplet on flat FLRW spatial sections",
                comoving_sheet_kinetic_coefficient_over_radius_required=-1.)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,sort_keys=True)+"\n")
