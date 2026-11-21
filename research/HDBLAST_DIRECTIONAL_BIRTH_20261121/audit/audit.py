"""Independent hard-clipping angular quadrature, compact-isotropy and source-ledger checks."""
import argparse
from fractions import Fraction
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec

HERE=Path(__file__).resolve().parent


def rational_loss(c, d, rate):
    """Exact polynomial antiderivatives for loss moments 1,u,u²,u⁴."""
    c, d, rate = map(Fraction, (c, d, rate))
    points = {Fraction(-1), Fraction(1)}
    if d:
        for threshold in (-rate, 0, rate):
            crossing = (c-threshold)/d
            if -1 < crossing < 1:
                points.add(crossing)
    powers = (0, 1, 2, 4)
    moments = [Fraction(0) for _ in powers]
    ordered = sorted(points)
    for left, right in zip(ordered, ordered[1:]):
        value = c - d*(left+right)/2
        if abs(value) >= rate:
            continue
        sign = 1 if value >= 0 else -1
        intercept, slope = 1-sign*c/rate, sign*d/rate
        for j, power in enumerate(powers):
            moments[j] += (intercept*(right**(power+1)-left**(power+1))/(power+1)
                           +slope*(right**(power+2)-left**(power+2))/(power+2))/2
    m0, m1, m2, m4 = moments
    return (m0, m1, m2, (35*m4-30*m2+3*m0)/8)


def moments(c, d, rate):
    """Integrate the triangular loss profile using three Gauss points on each linear piece."""
    c, d, rate = np.broadcast_arrays(c, d, rate)
    center = np.divide(c, d, out=np.zeros_like(c, dtype=float), where=d != 0)
    width = np.divide(rate, np.abs(d), out=np.full_like(c, np.inf, dtype=float), where=d != 0)
    edges = [np.clip(bound, -1., 1.) for bound in (center-width, center, center+width)]
    ans = np.zeros((4, *c.shape))
    for left, right in zip(edges, edges[1:]):
        half = .5*(right-left)
        mid = .5*(right+left)
        for z, weight in ((-math.sqrt(.6), 5/9), (0., 8/9), (math.sqrt(.6), 5/9)):
            u = mid+z*half
            loss = .5*half*weight*np.maximum(0., 1.-np.abs(c-d*u)/rate)
            u2 = u*u
            for row, value in enumerate((1., u, u2, ((35*u2-30)*u2+3)/8)):
                ans[row] += loss*value
    return ans


def audit_history(file, rest, mu, action):
    history = np.memmap(file, dtype=np.float64, mode="r").reshape(-1, 10)
    coverage = 0.
    sums = np.zeros(4)
    peak = dict(loss=-1.)
    low, high = .00001, .99999
    for block in np.array_split(history, max(1, math.ceil(len(history)/32768))):
        times, dt, D0, D1, k0, k1, electric0, electric1, _, _ = block.T
        fraction = (D1-D0)/rest
        start = np.clip((low-D0/rest)/fraction, 0., 1.)
        end = np.clip((high-D0/rest)/fraction, 0., 1.)
        active = end > start
        k = k0 + (.5*(end+start))*(k1-k0)
        c = ((k1-k0)/dt)[active]
        t = (2*math.pi**2*action*(D1-D0)/(12*mu*dt*k**2))[active]
        d0 = (electric0+start*(electric1-electric0))[active]
        d1 = (electric0+end*(electric1-electric0))[active]
        w = ((end-start)*fraction)[active]
        coverage += float(np.sum(w))
        # For c>=t/2 the angular-loss extremum inside the mixed clipping region
        # occurs at |d|=sqrt(c²+2|c|t−t²). Also check the boundaries and d=0.
        smooth_peak = np.sqrt(np.maximum(0., c*c+2*np.abs(c)*t-t*t))
        for test in (np.zeros_like(c),smooth_peak,-smooth_peak,d0,d1):
            d_test=np.clip(test,np.minimum(d0,d1),np.maximum(d0,d1))
            frac=moments(c,d_test,t)[0]
            i=int(np.argmax(frac))
            if frac[i]>peak["loss"]:
                u=0. if d1[i]==d0[i] else (d_test[i]-d0[i])/(d1[i]-d0[i])
                s=start[active][i]+u*(end-start)[active][i]
                peak=dict(loss=float(frac[i]),time=float(times[active][i]+dt[active][i]*s),
                          shell_speed_over_requested=float(c[i]/t[i]),force_speed_magnitude_over_requested=float(abs(d_test[i])/t[i]))

        # On steep ramps, resolve all absolute-value branch endpoints explicitly
        # before adaptive integration. Most other steps change negligibly.
        danger = np.abs(d1-d0) > .05*np.maximum(np.minimum(np.abs(d0), np.abs(d1)), t)
        regular = ~danger
        for s, weight in ((.5*(1-math.sqrt(.6)), 5/18), (.5, 4/9), (.5*(1+math.sqrt(.6)), 5/18)):
            args = (c[regular], d0[regular]+s*(d1[regular]-d0[regular]), t[regular])
            sums += weight*np.sum(moments(*args)*w[regular], axis=1)
        for c1, x0, x1, target, weight in zip(c[danger],d0[danger],d1[danger],t[danger],w[danger]):
            cuts = sorted({(sign*value-x0)/(x1-x0)
                           for value in (0., c1, c1+target, c1-target) for sign in (-1,1)
                           if 0 < (sign*value-x0)/(x1-x0) < 1})
            def integrand(s):
                return moments(c1, x0+s*(x1-x0), target)
            integral, _ = quad_vec(integrand, 0., 1., points=cuts, epsabs=2e-10, epsrel=1e-8)
            sums += weight*integral
    loss, loss_u, loss_u2, loss_P4 = map(float, sums)
    accepted = coverage-loss
    if abs(coverage-(high-low)) > 2e-10 or not 0 < accepted <= coverage:
        raise ArithmeticError("Coverage or sign of angular deficit changed")
    return dict(coverage=coverage, loss_per_total_requested=loss,
                accepted_per_total_requested=accepted, signed_cosine_of_losses=loss_u,
                mean_accepted_axis_cosine=-loss_u/accepted,
                mean_accepted_cosine_squared=(coverage/3-loss_u2)/accepted,
                cubic_fourth_moment_of_accepted=-.4*loss_P4/accepted,
                peak_directional_snapshot=peak)


def discrete_checks():
    # Span flat caps, tangencies, full loss, clipped grazing bands and both signs.
    tests = [(c,d,t) for c in (-7, -1, 0, 1, 7)
             for d in (-20, -3, 0, 3, 20) for t in (1,3,12)]
    worst = 0.
    tensor_error = 0.
    dipole_sign_count = 0
    for c, d, target in tests:
        exact = np.array(list(map(float,rational_loss(c,d,target))))
        trial = moments(c,d,target)
        worst = max(worst, float(np.max(np.abs(trial-exact))))
        A0, A1, A2 = 1-exact[0], -exact[1], 1/3-exact[2]
        if c*d*A1 > 1e-13:
            raise ArithmeticError("Dipole sign inconsistent with the physical acceleration axis")
        dipole_sign_count += int(c*d*A1 < -1e-13)
        trace_matrix = np.zeros((3,3))
        for axis in np.eye(3):
            dyad = np.outer(axis,axis)
            trace_matrix += A2*dyad + .5*(A0-A2)*(np.eye(3)-dyad)
        tensor_error=max(tensor_error,float(np.max(np.abs(trace_matrix - A0*np.eye(3)))))
    if worst>3e-14 or tensor_error>3e-14 or dipole_sign_count<10:
        raise ArithmeticError("Angular quadrature or triad check failed")

    # Reconstruct the actual four-vector source for independent individual
    # signed charge classes, then compare against the stream variables.
    ledger_error=0.
    smallest_rest_heat=math.inf
    negative_mix_count=0
    for drift in (-.82,-.1,.1,.82):
        gamma=1/math.sqrt(1-drift**2)
        for c,d,target in tests:
            loss=np.array(rational_loss(c,d,target),dtype=float)
            acc, first=1-loss[0],-loss[1]
            R=1.4
            p=math.sqrt(R*R-1)*first
            E=R*acc
            q0=acc
            for sign in (-1,1):
                energy_prime=gamma*(E-sign*drift*(sign*p))
                gap=energy_prime-q0
                smallest_rest_heat=min(smallest_rest_heat,gap)
                mix=energy_prime-E
                negative_mix_count+=int(mix<-1e-14)
                ledger_error=max(ledger_error,abs(mix-((gamma-1)*E-gamma*drift*p)))
            D,S=1.2,-.69
            cold=math.hypot(D,S)
            expected=E-(D*q0+S*p)/cold
            formula=(R-D/cold)*q0-S*p/cold
            ledger_error=max(ledger_error,abs(expected-formula))
            if E*E+1e-14 < q0*q0+p*p:
                raise ArithmeticError("Birth source is spacelike")
    if ledger_error>2e-14 or smallest_rest_heat< -1e-15 or negative_mix_count==0:
        raise ArithmeticError("Source identities failed")
    return dict(rational_angular_cases=len(tests),max_error_in_loss_moments=worst,
                max_summed_triad_tensor_error=tensor_error,directed_dipole_sign_cases=dipole_sign_count,
                max_source_identity_error=ledger_error,
                minimum_rest_frame_birth_excitation=smallest_rest_heat,
                signed_mixing_negative_examples=negative_mix_count)


def header_checks():
    """Test the production angular header via independent numerical integration."""
    header=HERE.parent/"feedback/angular.hpp"
    driver=r'''#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "angular.hpp"
int main() {
    double c,d,t;
    std::cout << std::setprecision(17);
    while(std::cin >> c >> d >> t) {
        auto h=hard_moments(c,d,t), v=vacancy_moments(c,d,t);
        std::cout << h.accept << ' ' << h.first << ' ' << v.accept << ' ' << v.first << '\n';
    }
}'''
    cases=[(c,d,t) for c in (-7.,-1.,0.,1.,7.) for d in (-20.,-3.,0.,3.,20.) for t in (1.,3.,12.)]
    cases.extend([(1.,1e-16,1.), (1.,1e-7,1.), (1.,-1e-3,1.), (1.,-1.0000001,1.e-5),
                  (-.989999999999,1.,.01),(0.,1.e9,1.),(1.1,1e8,1.),
                  (1.,-1.e-4,1.e4),(1e-14,-.01,1.),(.01,10.,1.e-4)])
    with tempfile.TemporaryDirectory(prefix="hdblast-audit-") as tmp:
        src=Path(tmp)/"header.cpp"
        exe=src.with_suffix("")
        src.write_text(driver)
        subprocess.run([os.environ.get("CXX","g++"),"-std=c++17","-O2","-I",str(header.parent),str(src),"-o",str(exe)],check=True)
        raw=subprocess.check_output([str(exe)],input="".join("%r %r %r\n"%case for case in cases),text=True)
    max_error=np.zeros(4)
    for (c,d,t),row in zip(cases,raw.splitlines()):
        hard=rational_loss(c,d,t)
        cuts={0.}
        if d!=0:
            for n in range(-2,17):
                speed=t*10.**n
                cuts.update((c/d,(c-speed)/d,(c+speed)/d))
        def vacancy(u):
            speed=abs(c-d*u)
            return .5*speed/(speed+t)*np.array((1.,u))
        vac, _=quad_vec(vacancy,-1,1,points=sorted(v for v in cuts if abs(v)<1.),
                        epsabs=1e-13,epsrel=1e-13)
        expected=[float(1-hard[0]),float(-hard[1]),*vac]
        max_error=np.maximum(max_error,np.abs(np.fromstring(row,sep=" ")-expected))
    if max(max_error)>2e-12:
        raise ArithmeticError(("Actual header differs from test cases",max_error))
    return dict(production_header_cases=len(cases),max_absolute_errors=max_error.tolist())

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native",type=Path,required=True)
    ap.add_argument("--prior",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    opts=ap.parse_args()
    env=json.loads((opts.prior/"envelope/ENVELOPE.json").read_bytes())
    runs={}
    for n in ("mass_4e10","mass_4e12"):
        scan=env["screen"][n]
        action=scan["witness_action"] if n=="mass_4e10" else 5*scan["elastic_drag_action_floor"]
        runs[n]={}
        for name,cfl in (("native","0.08"),("half_step","0.04")):
            runs[n][name]=audit_history(opts.native/n/("forward."+cfl+".bin"),
                env["runs"][n]["end_rest"],env["parameters"][n]["daughter_mass"],action)
    results=dict(fixtures=discrete_checks(),header=header_checks(),histories=runs)
    opts.output.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
