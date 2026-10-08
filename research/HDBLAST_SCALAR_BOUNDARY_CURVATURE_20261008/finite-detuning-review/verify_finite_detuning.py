#!/usr/bin/env python3
"""Independent diagnostic integration; no validated numerical bound is claimed.

Run with NumPy, SciPy and mpmath. Outputs JSON beside this source. Does not
import any attached or companion calculation. Dimensionless t=ky and r=k rho.
"""
import json
import math
from pathlib import Path
import platform
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import root
import mpmath as mp

mp.mp.dps = 80
MODELS = [
    dict(name="exact-w5-low", k=1., alpha=5., f0=1., dbar=.1, v=0., g2=0.),
    dict(name="exact-w5-medium", k=.2, alpha=5., f0=1., dbar=1., v=0., g2=0.),
    dict(name="exact-w5-high", k=1., alpha=5., f0=2., dbar=10., v=0., g2=0.),
    dict(name="asymmetric-w8", k=1., alpha=8., f0=1., dbar=.1, v=2., g2=.5),
    dict(name="near-threshold", k=.1, alpha=4.2, f0=.7, dbar=3., v=-1., g2=-.5),
    dict(name="large-detuning", k=.5, alpha=6., f0=1.3, dbar=30., v=4., g2=2.),
]

def potential(eta, model):
    alpha, v = model["alpha"], model["v"]
    w = 3 + alpha*eta**2/2 + v*eta**3/6
    wp = alpha*eta + v*eta**2/2
    wpp = alpha + v*eta
    u = wp**2/2 - 2*w**2/3
    up = wp*wpp - 4*w*wp/3
    upp = wpp**2 + wp*v - 4*(wp**2 + w*wpp)/3
    return w, wp, u, up, upp

def linear_data(model):
    alpha, dbar, f0 = map(mp.mpf, (model["alpha"], model["dbar"], model["f0"]))
    aa = 1 + dbar*f0/6
    tb = mp.acoth(aa)
    z = -mp.sinh(tb)**2
    a, b = 2-alpha/2, alpha/2
    f = mp.hyp2f1(a, b, mp.mpf(5)/2, z)
    ft = (a*b/(mp.mpf(5)/2))*mp.hyp2f1(a+1,b+1,mp.mpf(7)/2,z)*(-2*mp.sinh(tb)*mp.cosh(tb))
    d = ft/f
    dp = alpha*(alpha-4)-d*d-4*aa*d
    e = 4*aa*(d+alpha)-dp
    k, delta = mp.mpf(model["k"]), mp.mpf(model["k"])*dbar
    c2 = -delta**2*e/(48*(d+alpha)**2)
    # A second, floating-point route integrates only the logarithmic derivative.
    mass2 = float(alpha*(alpha-4))
    tstart = 1e-6
    p1 = mass2/5
    p3 = -(p1*p1+4*p1/3)/7
    ric = solve_ivp(lambda t, x: [mass2-x[0]**2-4*x[0]/math.tanh(t)],
                    (tstart,float(tb)), [p1*tstart+p3*tstart**3], method="DOP853",
                    rtol=3e-13,atol=1e-14,max_step=.01)
    if not ric.success:
        raise RuntimeError(ric.message)
    output = dict(tb=float(tb), F=float(f), D_over_k=float(d), Dprime_over_k2=float(dp),
                  E_over_k2=float(e), C2=float(c2),
                  riccati_abs_error=float(abs(mp.mpf(ric.y[0,-1])-d)),
                  horizon_amplitude_per_c=float(-dbar/(2*(d+alpha)*f)),
                  H0_squared=float(k*k*(aa*aa-1)))
    if model["alpha"] == 5:
        c2_explicit = -delta**2*(3+20*aa+1/aa**2)/(48*(1/aa+5)**2)
        output["exact_w5_D_error_80dps"] = str(abs(d-1/aa))
        output["exact_w5_C2_error_80dps"] = str(abs(c2-c2_explicit))
    return output

def integrate(model, eta_h, tb, refined):
    tstart = 5e-6 if refined else 1e-5
    _, _, u, up, upp = potential(eta_h, model)
    b2 = up/10
    a3 = -u/36
    b4 = b2*(upp-16*a3)/28
    a5 = u*u/4320-up*up/750
    y0 = [tstart+a3*tstart**3+a5*tstart**5,
          1+3*a3*tstart**2+5*a5*tstart**4,
          eta_h+b2*tstart**2+b4*tstart**4,
          2*b2*tstart+4*b4*tstart**3]
    def rhs(t, x):
        r, rp, eta, p = x
        _, _, u, up, _ = potential(eta, model)
        return [rp, -r*(p*p/4+u/6), p, up-4*rp*p/r]
    sol = solve_ivp(rhs, (tstart,tb), y0, method="DOP853", dense_output=True,
                    rtol=3e-13 if refined else 2e-11,
                    atol=[1e-14,1e-14,1e-17,1e-17],
                    max_step=.005 if refined else .015)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol

def nonlinear(model, linear, c, refined):
    ah = linear["horizon_amplitude_per_c"]
    def residual(z):
        scale, tb = z
        if tb <= 1e-4:
            return np.array([1e6,1e6])
        sol = integrate(model, c*ah*scale, tb, refined)
        r, rp, eta, p = sol.y[:,-1]
        w, wp, _, _, _ = potential(eta,model)
        g = eta + model["g2"]*eta*eta/2
        gp = 1+model["g2"]*eta
        sigma = 2*w+model["dbar"]*(model["f0"]+c*g)
        sigmap = 2*wp+model["dbar"]*c*gp
        return np.array([rp/r-sigma/6,(p+sigmap/2)/c])
    rt = root(residual,[1.,linear["tb"]],method="hybr",options={"xtol":1e-10})
    res = residual(rt.x)
    if np.max(np.abs(res)) > 1e-9:
        raise RuntimeError(f"Shooting did not converge: {model['name']} {c} {rt.message} {res}")
    scale, tb = rt.x
    sol = integrate(model,c*ah*scale,tb,refined)
    r,rp,eta,p = sol.y[:,-1]
    w,wp,u,_,_ = potential(eta,model)
    g = eta+model["g2"]*eta*eta/2
    gp = 1+model["g2"]*eta
    sigma = 2*w+model["dbar"]*(model["f0"]+c*g)
    sigmap = 2*wp+model["dbar"]*c*gp
    hsq = model["k"]**2/r**2
    shell_hsq = model["k"]**2*(sigma*sigma/36-sigmap*sigmap/48+u/6)
    # Evaluate the exact constraint identity relative to the base background
    # without subtracting two O(H0^2) floating-point numbers. This is separate
    # from the geometric 1/r^2 readout, whose cancellation floor is retained.
    dw = model["alpha"]*eta*eta/2+model["v"]*eta**3/6
    dsigma = 2*dw+model["dbar"]*c*g
    du = wp*wp/2-4*dw-2*dw*dw/3
    hshift_stable = model["k"]**2*((2*(6+model["dbar"]*model["f0"])*dsigma+dsigma*dsigma)/36-sigmap*sigmap/48+du/6)
    sample = sol.sol(np.linspace(1e-4,tb,201))
    rr,rrp,etas,pp = sample
    uu = potential(etas,model)[2]
    constraint = rrp**2-1-rr**2*(pp**2/12-uu/6)
    relconstraint = np.max(np.abs(constraint)/(1+rrp**2+rr**2*np.abs(pp**2/12-uu/6)))
    return dict(c=c, refined=refined, root_success=bool(rt.success),
                horizon_scale=float(scale), tb=float(tb), eta_boundary=float(eta),
                H_squared=float(hsq), shell_identity_abs_error=float(abs(hsq-shell_hsq)),
                H_squared_shift_stable_constraint=float(hshift_stable),
                junction_residual_max=float(max(abs(res[0]),abs(c*res[1]))),
                sampled_relative_constraint_max=float(relconstraint))

def main():
    output = {"status":"DIAGNOSTICS_ONLY_NOT_A_VALIDATED_BOUND",
              "scope":"Fixed positive detuning, independently varied weak scalar coupling; no claim of uniform detuning-zero continuation.",
              "runtime":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__},
              "models":[]}
    for model in MODELS:
        lin = linear_data(model)
        runs = [nonlinear(model,lin,c,refined) for refined in (False,True)
                for size in (.04,.02,.01) for c in (size,-size)]
        centered = []
        for size in (.04,.02,.01):
            pair = [q for q in runs if q["refined"] and abs(q["c"])==size]
            qval = (sum(q["H_squared"] for q in pair)-2*lin["H0_squared"])/(2*size*size)
            qstable = sum(q["H_squared_shift_stable_constraint"] for q in pair)/(2*size*size)
            centered.append(dict(c=size, centered_coefficient=qval,
                                 error=qval-lin["C2"],relative_error=abs(qval/lin["C2"]-1),
                                 stable_constraint_coefficient=qstable,
                                 stable_constraint_error=qstable-lin["C2"],
                                 stable_constraint_relative_error=abs(qstable/lin["C2"]-1)))
        output["models"].append(dict(parameters=model,linear=lin,runs=runs,centered=centered))
        print(model["name"], "C2",lin["C2"],"stable centered relative errors",[q["stable_constraint_relative_error"] for q in centered],flush=True)
    allruns = [r for m in output["models"] for r in m["runs"]]
    settingdiffs = [abs(m["runs"][i]["H_squared"]-m["runs"][i+6]["H_squared"])
                   for m in output["models"] for i in range(6)]
    output["summary"] = dict(nonlinear_configurations=len(allruns),
        max_junction_residual=max(q["junction_residual_max"] for q in allruns),
        max_sampled_relative_constraint=max(q["sampled_relative_constraint_max"] for q in allruns),
        max_shell_identity_abs_error=max(q["shell_identity_abs_error"] for q in allruns),
        max_H_squared_settings_difference=max(settingdiffs),
        all_centered_errors_decrease=all(abs(m["centered"][i+1]["error"])<abs(m["centered"][i]["error"]) for m in output["models"] for i in (0,1)),
        all_stable_centered_errors_decrease=all(abs(m["centered"][i+1]["stable_constraint_error"])<abs(m["centered"][i]["stable_constraint_error"]) for m in output["models"] for i in (0,1)),
        all_C2_negative=all(m["linear"]["C2"]<0 for m in output["models"]),
        max_linear_riccati_hypergeometric_difference=max(m["linear"]["riccati_abs_error"] for m in output["models"]))
    path = Path(__file__).with_name("FINITE_DETUNING_DIAGNOSTICS.json")
    path.write_text(json.dumps(output,indent=2)+"\n")
    print(json.dumps(output["summary"],indent=2))

if __name__ == "__main__":
    main()
