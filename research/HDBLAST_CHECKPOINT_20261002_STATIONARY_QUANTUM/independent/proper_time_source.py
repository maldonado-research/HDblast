#!/usr/bin/env python3
"""Independent S4 proper-time determinant/source check; no producer imports."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp

ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def heat_coefficients(order):
    """Bracket coefficients of (4 pi s)^2 K_z(s)/Vol, in powers of z*s.

    Mellin/Euler–Maclaurin half-integer sum:
    sum n^p exp(-t*n*n) ~ Gamma((p+1)/2)/(2*t^((p+1)/2))
      + sum (-t)^k/k! zeta(-p-2*k,1/2).
    The n=1/2 term of n^3-n/4 vanishes identically.
    """
    base = [mp.mpf(0) for _ in range(order + 1)]
    base[0], base[1] = mp.mpf(1), -mp.mpf(1)/4
    for k in range(order-1):
        # zeta(-m,q) = -B_(m+1)(q)/(m+1), exact rational inputs.
        c = (-mp.bernpoly(4+2*k, mp.mpf(1)/2)/(4+2*k)
             + mp.bernpoly(2+2*k, mp.mpf(1)/2)/(4*(2+2*k)))
        base[k+2] = 2*(-1)**k*c/mp.factorial(k)
    return [mp.fsum(base[j]*(mp.mpf(9)/4)**(n-j)/mp.factorial(n-j)
                    for j in range(n+1)) for n in range(order+1)]


def uv_integrals(x, z, r, join, order):
    a = heat_coefficients(order)
    h = x-r
    A, B = 2*z-h, h*h/2-2*h*z+mp.mpf(29)*z*z/15
    Az, Bz = 2, -2*h+mp.mpf(58)*z/15
    result = [mp.mpf(0)]*3
    cancellations = []
    for n in range(order+1):
        d = mp.fsum(a[j]*z**j*(-x)**(n-j)/mp.factorial(n-j)
                    for j in range(n+1))
        dz = mp.fsum(j*a[j]*z**(j-1)*(-x)**(n-j)/mp.factorial(n-j)
                     for j in range(1,n+1))
        sn = (-r)**n/mp.factorial(n)
        sz = mp.mpf(0)
        qn = (-r)**n/mp.factorial(n)
        if n >= 1:
            sn += A*(-r)**(n-1)/mp.factorial(n-1)
            sz += Az*(-r)**(n-1)/mp.factorial(n-1)
            qn += (2*z-h)*(-r)**(n-1)/mp.factorial(n-1)
        if n >= 2:
            sn += B*(-r)**(n-2)/mp.factorial(n-2)
            sz += Bz*(-r)**(n-2)/mp.factorial(n-2)
        if n <= 2:
            cancellations.append(abs(d-sn))
        if n >= 3:
            result[0] -= (d-sn)*join**(n-2)/(32*mp.pi**2*(n-2))
            result[1] -= ((d-sn)-z*(dz-sz)/2)*join**(n-2)/(32*mp.pi**2*(n-2))
        if n >= 2:
            result[2] += (d-qn)*join**(n-1)/(16*mp.pi**2*(n-1))
    require(max(cancellations) < mp.mpf(10)**(-mp.mp.dps+5), "UV heat cancellation failed")
    return result


def gaussian_tail(t, N, power):
    """Bounds sum_(l=N)^infinity (l+2)^power exp(-t*l²)/3.

    N is beyond the maximum of the summand. Use first term plus integral.
    """
    require(N > mp.sqrt(power/(2*t)), "Tail monotonicity condition failed")
    integral = mp.fsum(mp.binomial(power,j)*2**(power-j)*
                      mp.gammainc((j+1)/2,t*N*N,mp.inf)/(2*t**((j+1)/2))
                      for j in range(power+1))
    return ((N+2)**power*mp.exp(-t*N*N)+integral)/3


def heat_sums(t, cutoff):
    k = mp.mpf(0)
    moment = mp.mpf(0)
    for ell in range(cutoff):
        lam = ell*(ell+3)
        degeneracy = (ell+1)*(ell+2)*(2*ell+3)//6
        value = degeneracy*mp.exp(-t*lam)
        k += value
        moment += lam*value
    return k, moment


def sources(x, z, r, settings):
    mp.mp.dps = settings["dps"]
    x,z,r = map(mp.mpf, (x,z,r))
    join, upper = mp.mpf(settings["join"])/r, mp.mpf(settings["upper"])/r
    order = settings["order"]
    h = x-r
    A, B = 2*z-h, h*h/2-2*h*z+mp.mpf(29)*z*z/15
    invvol = 3*z*z/(8*mp.pi**2)
    N = int(mp.ceil(mp.sqrt(mp.log(10)*(settings["dps"]+10)/(z*join))))+5
    spectral0 = gaussian_tail(z*join,N,3)
    spectral1 = gaussian_tail(z*join,N,5)
    # Uniform-in-s bounds integrated from join to upper.
    spectral_bounds = [invvol*spectral0*mp.log(upper/join)/2,
                       z*invvol*spectral1*(upper-join)/4,
                       invvol*spectral0*(upper-join)]

    cache = {}
    def integrands(s):
        key = str(s)
        if key in cache:
            return cache[key]
        k, moment = heat_sums(z*s,N)
        physical = mp.exp(-x*s)*invvol*k
        pref = mp.exp(-r*s)/(16*mp.pi**2*s*s)
        subtraction = pref*(1+A*s+B*s*s)
        w = -(physical-subtraction)/(2*s)
        q = physical-pref*(1+(2*z-h)*s)
        rho = (-z*mp.exp(-x*s)*invvol*moment +
               pref*(2+(2*z-2*h)*s+(h*h-2*z*h)*s*s)/s)/4
        cache[key] = (w,rho,q)
        return cache[key]

    nodes = [join]
    while nodes[-1] < 1/r:
        nodes.append(min(nodes[-1]*4,1/r))
    for value in [4/r,16/r,upper]:
        if value > nodes[-1]:
            nodes.append(value)
    uv = uv_integrals(x,z,r,join,order)
    # A vector cache means each spectral sum is shared across the three
    # distinct integrands, but no determinant values are shared with producer.
    bulk = [mp.quadgl(lambda s: integrands(s)[j], nodes) for j in range(3)]
    values = [uv[j]+bulk[j] for j in range(3)]

    ku, mu = heat_sums(z*upper,N)
    ku += gaussian_tail(z*upper,N,3)
    mu += gaussian_tail(z*upper,N,5)
    def exponential_moment(power):
        return r**(-power-1)*mp.gammainc(power+1,r*upper,mp.inf)
    ir_bounds = [
        invvol*ku*mp.e1(x*upper)/2 +
        (exponential_moment(-3)+abs(A)*exponential_moment(-2)+abs(B)*exponential_moment(-1))/(32*mp.pi**2),
        z*invvol*mu*mp.exp(-x*upper)/(4*x) +
        (2*exponential_moment(-3)+abs(2*z-2*h)*exponential_moment(-2)+abs(h*h-2*z*h)*exponential_moment(-1))/(64*mp.pi**2),
        invvol*ku*mp.exp(-x*upper)/x +
        (exponential_moment(-2)+abs(2*z-h)*exponential_moment(-1))/(16*mp.pi**2),
    ]
    anomaly = B/(16*mp.pi**2)
    trace = -4*values[1]+x*values[2]-anomaly
    wrong_source = -4*values[1]-x*values[2]-anomaly
    # Changing rho=W-z*W_z/2 to W+z*W_z/2 gives rho_wrong=2W-rho.
    wrong_curvature = -4*(2*values[0]-values[1])+x*values[2]-anomaly
    scale = mp.mpf("1e-9")+mp.mpf("1e-7")*max(abs(v) for v in values)
    return {"x":str(x),"z":str(z),"r":str(r),"W":mp.nstr(values[0],40),
            "rho":mp.nstr(values[1],40),"Q":mp.nstr(values[2],40),
            "p":mp.nstr(-values[1],40),"trace_residual":mp.nstr(trace,12),
            "negative_current_sign_residual":mp.nstr(wrong_source,12),
            "negative_metric_sign_residual":mp.nstr(wrong_curvature,12),
            "negative_controls_detected": bool(abs(wrong_source)>scale and abs(wrong_curvature)>scale),
            "spectral_tail_bounds":list(map(lambda v:mp.nstr(v,12),spectral_bounds)),
            "IR_tail_bounds":list(map(lambda v:mp.nstr(v,12),ir_bounds)),
            "spectral_cutoff":N,"integrand_nodes":len(cache),"settings":settings,
            "uv_remainder_certified":False,"arithmetic_interval_certified":False,
            "quadrature_error_certified":False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",default="INDEPENDENT_PROPER_TIME.json")
    args = parser.parse_args()
    start = time.time()
    prereg = ROOT/"PREREGISTRATION.md"
    settings = [{"dps":40,"join":".01","upper":"60","order":10},
                {"dps":55,"join":".005","upper":"80","order":12}]
    results = []
    for setting in settings:
        for x in [1,2]:
            for z in [".5","1"]:
                value = sources(x,z,1,setting)
                results.append(value)
                print(json.dumps(value),flush=True)
    gates = []
    for first,second in zip(results[:4],results[4:]):
        for name in ["W","rho","Q"]:
            delta=abs(mp.mpf(first[name])-mp.mpf(second[name]))
            tolerance=mp.mpf("1e-9")+mp.mpf("1e-7")*abs(mp.mpf(second[name]))
            gates.append({"x":first["x"],"z":first["z"],"observable":name,
                          "difference":mp.nstr(delta,12),"pass":bool(delta<=tolerance)})
    # Coefficient perturbation must leave an s² mismatch, hence a logarithmic
    # divergence in W. This does not rely on the numerically evaluated values.
    a=heat_coefficients(2)
    wrong_a2=mp.mpf(28)/15
    heat_negative=abs(a[2]-wrong_a2)>mp.mpf(".01")
    document={"status":"PASS" if all(g["pass"] for g in gates) and
              all(v["negative_controls_detected"] for v in results) and heat_negative else "FAIL",
              "preregistration_sha256":hashlib.sha256(prereg.read_bytes()).hexdigest(),
              "results":results,"refinement_gates":gates,
              "wrong_29_over_15_detected":bool(heat_negative),
              "total_error_enclosure":False,"elapsed_seconds":time.time()-start}
    (ROOT/args.output).write_text(json.dumps(document,indent=2)+"\n")
    require(document["status"]=="PASS","Independent proper-time review gate failed")


if __name__=="__main__":
    main()
