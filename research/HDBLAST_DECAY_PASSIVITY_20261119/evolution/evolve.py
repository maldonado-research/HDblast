"""Test alternative Gaussian loss targets on pinned reciprocal B1 wall trajectories."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline


WALL = Path("research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip")
WALL_HASH = "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"
SUMMARY = Path("research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/fine/quiet/rolling_summary.json")
SUMMARY_HASH = "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810"
CUT = 2.10
ARMS = ("quiet", "plus", "minus", "far")


def checked(path, pin):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=pin:
        raise ValueError(f"Unexpected source {path}")
    return raw


def read_csv(archive, key):
    with archive.open("fine/y1_b1/"+key+".csv") as stream:
        return np.genfromtxt(stream, delimiter=",", names=True, dtype=float)


def background(pre, post, reference):
    def column(key):
        return np.r_[pre[key],post[key]]
    clock = column("H0tau")
    if np.min(np.diff(clock))<=0 or not clock[0] < CUT < clock[-1]:
        raise ValueError("Invalid wall clock")
    scale_log = CubicHermiteSpline(clock, column("A_b")-reference, column("H_over_H0"))
    mass = CubicHermiteSpline(clock, 100*(column("phi_b")-.5), 100*column("v_over_H0"))
    def geometry(t):
        return np.exp(scale_log(t)), scale_log(t, 1), mass(t), mass(t, 1)
    return geometry, float(clock[0]), float(np.interp(CUT, post["H0tau"], post["rad_hat"]))


def evolve(geometry, start, mode, order, tol=1e-9, curves=None):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    k = 4*(nodes+1)
    measure = 4*weights*k**2/(2*np.pi**2)
    a0, h0, mu0, dmu0 = geometry(start)
    e0 = np.hypot(k/a0,mu0)
    b0 = a0**3*e0
    # The first two experiments retain the published preparation. The third
    # prepares the pure target of the sheared bath and measures its physical cost.
    d0 = h0 + (h0*mu0**2+mu0*dmu0)/(2*e0**2) if mode=="sheared_prepared" else np.zeros_like(e0)
    initial_energy = float(measure@(d0**2/(4*e0)))
    initial = np.r_[.5/b0,-d0/(2*e0),.5*b0 + a0**3*d0**2/(2*e0),np.zeros(8)]

    def values(t, state):
        a,h,mu,dmu = geometry(t)
        v = a**3
        p2 = (k/a)**2
        e2 = p2 + mu**2
        e = np.sqrt(e2)
        b = v*e
        gamma = mu**2/(8*np.pi*e)
        X,C,P = state[:3*order].reshape(3,order)
        x0 = .5/b
        p0 = b/2
        dx, dp = X-x0, P-p0
        U = .5*(dp/v+v*e2*dx)
        pressure = .5*(dp/v-v*(mu**2+p2/3)*dx)
        source = v*mu*dmu*dx
        # For the u=a*chi target, D is minus the real logarithmic slope of chi.
        D = (h + (h*mu**2+mu*dmu)/(2*e2)) if mode!="physical" else np.zeros_like(e)
        Cb = -D/(2*e)
        Pb = p0 + v*D**2/(2*e)
        Hlost = gamma*(U-D**2/(4*e))
        counted = gamma*(U+D*C+.5*v*D**2*X)
        # Calculate the work independently of subtracting the two reported powers.
        drive = gamma*(D*C+.5*v*D**2*X+D**2/(4*e))
        flux = np.array([measure@Hlost, measure@counted, measure@drive])
        return a,h,e2,v,gamma,X,C,P,x0,Cb,Pb,U,source,pressure,flux

    def rhs(t, state):
        a,h,e2,v,g,X,C,P,xb,cb,pb,_,source,pressure,flux = values(t,state)
        return np.r_[2*C/v-g*(X-xb), P/v-v*e2*X-g*(C-cb), -2*v*e2*C-g*(P-pb),
                     flux, a*flux, measure@source, -3*h*(measure@pressure)]

    def positive_crossing(t,state):
        return values(t,state)[-1][0]
    positive_crossing.direction = 1
    sol = solve_ivp(rhs, (start,CUT), initial, method="DOP853",rtol=tol,atol=tol*1e-3,
                    max_step=.001,dense_output=True,events=positive_crossing if mode!="physical" else None)
    if not sol.success:
        raise ValueError(sol.message)
    def diagnostics(t):
        state=sol.sol(t)
        a,h,e2,v,g,X,C,P,xb,cb,pb,U,source,pressure,flux=values(t,state)
        J=state[3*order:3*order+3]
        Z=state[3*order+3:3*order+6]
        S,A=state[-2:]
        nu=measure@U
        if min(X*P-C*C)<.25-2e-7 or min(U/np.sqrt(e2)) < -2e-8:
            raise ValueError("Unphysical evolved covariance")
        errors = [abs(nu-initial_energy+J[0]-S-A),abs(J[1]-J[0]-J[2]),abs(Z[1]-Z[0]-Z[2]),abs(flux[1]-flux[0]-flux[2])]
        if errors[0] > 3e-7*(1+abs(S)+abs(A)) or max(errors[1:])>1e-8:
            raise ValueError(("balance",errors))
        return dict(s=t,scale_factor=float(a),field_energy_comoving=nu,
                    # H_loss is signed energy removed from the physical field, count is
                    # the chosen positive oscillator quantum, drive=count-H_loss.
                    rate_comoving=dict(zip(("H_loss","count","drive"),map(float,flux))),
                    accumulated_comoving=dict(zip(("H_loss","count","drive"),map(float,J))),
                    redshifted_density=dict(zip(("H_loss","count","drive"),map(float,Z/a**4))),
                    radiation_a4=dict(zip(("H_loss","count","drive"),map(float,Z))),
                    mass_work=float(S),expansion_work=float(A),max_balance_defect=max(errors),
                    min_uncertainty=float(min(X*P-C*C)))
    z0, z1, z2, end = [diagnostics(t) for t in (start,start+.0001,start+.001,CUT)]
    d=dict(initial=z0,early_1e_4=z1,early_1e_3=z2,end=end)
    # d(a^4 rho_H)/ds = a L_H, so negative-to-positive power crossings
    # give the interior candidate minima. Also retain both interval endpoints.
    if mode!="physical":
        crossings=sol.t_events[0]
        if mode=="sheared" and (z0["rate_comoving"]["H_loss"]>=0 or z1["rate_comoving"]["H_loss"]>=0 or not len(crossings)):
            raise ValueError("Missing transient negative physical-power interval")
        extrema = [z0,end]+[diagnostics(float(t)) for t in crossings]
        if mode=="sheared":
            d["first_H_loss_negative_to_positive"]=extrema[2]
        d["min_signed_a4_radiation"]=min(extrema,key=lambda r:r["radiation_a4"]["H_loss"])
    else:
        d["min_signed_a4_radiation"]=z0
    count_index = 3*order + 4
    d["min_sampled_count_a4"] = float(sol.y[count_index].min())
    if d["min_sampled_count_a4"] < -1e-9:
        raise ValueError("Positive-quantum reservoir became negative")
    if curves is not None:
        ts = np.linspace(start, CUT, 401)
        states = sol.sol(ts)
        curves["s"] = ts
        for j,name in enumerate(("signed","count","drive")):
            curves[mode+"_"+name+"_a4"] = states[3*order+3+j]
    return d


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("--order",type=int,default=24)
    ap.add_argument("--tol",type=float,default=1e-9)
    ap.add_argument("--curves",type=Path,help="Also write the quiet proper-time reservoir curves")
    args=ap.parse_args()
    curves = {} if args.curves else None
    if args.order<8 or args.tol<=0 or args.tol>1e-6:
        raise ValueError("Expected at least eight nodes and a positive accurate tolerance")
    reference = json.loads(checked(args.repository/SUMMARY,SUMMARY_HASH))["matter"]["events"][0]["A_star"]
    out={"order":args.order,"relative_tolerance":args.tol,"cut":CUT,"K":8,"grid":"fine","arms":{}}
    with zipfile.ZipFile(io.BytesIO(checked(args.repository/WALL,WALL_HASH))) as archive:
        pre=read_csv(archive,"pre")
        for arm in ARMS:
            geom,start,native=background(pre,read_csv(archive,arm),reference)
            rec={model:evolve(geom,start,model,args.order,args.tol,curves if arm=="quiet" else None) for model in ("physical","sheared","sheared_prepared")}
            rec["archived_physical_density"]=native
            # Stored native radiation uses linear interpolation to the proper-time cut.
            if abs(rec["physical"]["end"]["redshifted_density"]["count"]-native)>3e-4:
                raise ValueError("Physical-target control changed")
            out["arms"][arm]=rec
    quiet=out["arms"]["quiet"]
    for arm,row in out["arms"].items():
        row["contrasts"]={m:{k:100*(v/quiet[m]["end"]["redshifted_density"][k]-1)
                    for k,v in row[m]["end"]["redshifted_density"].items() if k!="drive"}
                    for m in ("physical","sheared","sheared_prepared")}
    if curves is not None:
        with args.curves.open("w", newline="") as output:
            writer=csv.writer(output,lineterminator="\n")
            writer.writerow(curves)
            writer.writerows(zip(*curves.values(),strict=True))
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
