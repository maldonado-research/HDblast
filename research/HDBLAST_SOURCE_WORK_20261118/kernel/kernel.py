"""Retarded linear and exact finite-amplitude occupation differences for paired B1 histories."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, simpson, solve_ivp
from scipy.interpolate import CubicSpline


START, END, DS = .8, 2.10, .001
ARMS = {"quiet": ("classical", "quiet"), "near_plus": ("classical", "0.001"),
        "near_minus": ("classical", "-0.001"), "far_plus": ("far", "0.001"),
        "far_half": ("far", "0.0005")}


class History:
    def __init__(self, cap, grid, name):
        folder, branch = ARMS[name]
        base = cap / "classical" / "runs" / grid
        event = json.loads((base/"quiet/rolling_summary.json").read_bytes())["matter"]["events"][0]
        with np.load(base/"checkpoint.npz", allow_pickle=False) as z:
            pre = z["hist"]
        with np.load(cap/folder/"runs"/grid/branch/"rolling_timeseries.npz", allow_pickle=False) as z:
            obs = dict(zip(z["cols"], z["rec"].T, strict=True))
        self.cut_clock = lambda t: float(np.interp(t, obs["T"], obs["H0tau"]))
        self.quiet_crossing = event["s_star"]
        join = obs["H0tau"][0]
        early = pre[:,0] < join-1e-10
        t = np.r_[pre[early,0],obs["H0tau"]]
        self.phi = CubicSpline(t,np.r_[pre[early,1],obs["phi_b"]])
        self.h = CubicSpline(t,np.r_[pre[early,2],obs["H_over_H0"]])
        self.int_h = self.h.antiderivative()
        self.offset = obs["A_b"][0]-event["A_star"]-self.int_h(join)

    def pieces(self, s, k):
        # Last axis is k. The common time is each wall's dimensionless proper time.
        log_a = self.offset + self.int_h(s)
        mass = (100 * (self.phi(s)-.5))**2
        geometry = (np.expand_dims(np.exp(-2*log_a), -1)*k**2
                    - np.expand_dims(2.25*self.h(s)**2+1.5*self.h(s,1),-1))
        return np.expand_dims(mass, -1), geometry


def make_basis(V, dV, t):
    """Local order-two adiabatic analyzer and its exact derivative for V+lambda*dV."""
    ts = t+DS*np.arange(-2,3)
    values = V(ts)
    if np.any(values<=0):
        raise ValueError("Endpoint adiabatic frequency is nonpositive")
    o = np.sqrt(values)
    do = dV(ts)/(2*o)
    diff = lambda x:(x[2:]-x[:-2])/(2*DS)
    diff2 = lambda x:(x[2:]+x[:-2]-2*x[1:-1])/DS**2
    vel, acc = diff(o), diff2(o)
    dvel, dacc = diff(do), diff2(do)
    w0, dw0 = o[1:-1], do[1:-1]
    w = w0-acc/(4*w0**2)+3*vel**2/(8*w0**3)
    dw = (dw0-dacc/(4*w0**2)+acc*dw0/(2*w0**3)
          +3*vel*dvel/(4*w0**3)-9*vel**2*dw0/(8*w0**4))
    wp, dwp = diff(w)[0], diff(dw)[0]
    w, dw = w[1], dw[1]
    f = 1/np.sqrt(2*w)
    df = -.5*f*dw/w
    y = -1j*w - .5*wp/w
    dy = -1j*dw - .5*(dwp/w-wp*dw/w**2)
    return (f, f*y), (df, df*y+f*dy)


def det(v, z):
    return v[0]*z[1]-v[1]*z[0]


def evolve(V, seed, ts, backward=False):
    points = ts[::-1] if backward else ts
    n = len(seed[0])
    def rhs(t, z):
        return np.r_[z[n:],-V(t)*z[:n]]
    result = solve_ivp(rhs, (points[0],points[-1]), np.r_[*seed].astype(complex),
                       method="DOP853", t_eval=points, rtol=3e-11, atol=5e-13, max_step=.003)
    if not result.success:
        raise ValueError(result.message)
    fields = result.y.T[::-1] if backward else result.y.T
    return (fields[:,:n],fields[:,n:])


def point(state, idx):
    return (state[0][idx],state[1][idx])


def design_controls(ts, K0, qV, q_i, q_f, n0, count):
    # Equal positive areas, fixed quiet geometry and no change near the crossing or either endpoint.
    width = .08
    half = width/2
    def window(s, center):
        phase = np.clip((s-center+half)/width, 0, 1)
        return np.where((s > center-half) & (s < center+half), np.sin(np.pi*phase)**4, 0.)
    candidates = np.linspace(1.60,1.90,16)
    forecast = [float(simpson(K0*window(ts,c),x=ts)) for c in candidates]
    controls = {}
    for label, idx in (("lowest",np.argmin(forecast)),("highest",np.argmax(forecast))):
        center = float(candidates[idx])
        pred = forecast[idx]
        actual = {}
        for height in (1., .5):
            def perturb(s):
                return qV(s)+height*window(s,center)
            state=evolve(perturb,q_i,np.array([START,END]))
            actual[str(height)] = count(det(q_f,point(state,-1)))-n0
        if pred*actual["1.0"] <= 0:
            raise ValueError("Designed phase sign disagrees")
        controls[label] = {"center":center,"predicted_delta_N_per_height":pred,
                           "direct_delta_N":actual}
    if controls["lowest"]["direct_delta_N"]["1.0"]*controls["highest"]["direct_delta_N"]["1.0"]>=0:
        raise ValueError("Phase control failed to separate")
    return {"width":width,"area_per_unit_height":3*width/8,"candidates":candidates.tolist(),
            "trial_predictions":forecast,"controls":controls}


def run_mesh(cap, grid, order, num_intervals, finite_difference, export=None):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    k = 4*(nodes+1)
    measure = 4*weights*k*k/(2*math.pi**2)
    count = lambda b:float(np.dot(measure,abs(b)**2))
    quiet = History(cap,grid,"quiet")
    cuts = {"pre_front":quiet.cut_clock(1.20), "quiet_crossing":quiet.quiet_crossing,
            "before_near_compensation":quiet.cut_clock(1.25),"s_1.60":1.60,"s_1.75":1.75,
            "s_1.95":1.95,"readout":END}
    ts = np.unique(np.r_[np.linspace(START,END,num_intervals+1),list(cuts.values())])
    indices = {key:int(np.searchsorted(ts,cut)) for key,cut in cuts.items()}
    qV = lambda t:sum(quiet.pieces(t,k))
    zeros = lambda t:np.zeros_like(qV(t))
    q_i = make_basis(qV,zeros,START)[0]
    q_f = make_basis(qV,zeros,END)[0]
    q = evolve(qV,q_i,ts)
    adj = evolve(qV,q_f,ts,True)
    b0 = det(q_f,point(q,-1))
    if np.max(abs(det(point(adj,0),q_i)-b0))>2e-8:
        raise ValueError("Forward/backward symplectic check")
    n0 = count(b0)
    variance = lambda db:float(np.dot(measure,2*np.real(b0.conj()*db)))
    response_weight = -2*np.real(adj[0]*q[0]*b0.conj())
    K0 = response_weight @ measure
    K2 = response_weight @ (k**2*measure)
    if export is not None:
        # Smooth response on this interval. The table is sufficient for independent reuse.
        tab = np.linspace(START,END,2049)
        np.savetxt(export, np.stack([tab,np.interp(tab,ts,K0),np.interp(tab,ts,K2)],axis=1),
                   fmt="%.15g", delimiter=",", header="s,K0,K2", comments="")
    qM,qG = quiet.pieces(ts,k)
    arms = {}
    for name in ("near_plus","near_minus","far_plus","far_half"):
        arm=History(cap,grid,name)
        pV=lambda t:sum(arm.pieces(t,k))
        delta=lambda t:pV(t)-qV(t)
        p_i,_ = make_basis(pV,zeros,START)
        p_f,_ = make_basis(pV,zeros,END)
        _,di = make_basis(qV,delta,START)
        _,df = make_basis(qV,delta,END)
        p=evolve(pV,p_i,ts)
        b1=det(p_f,point(p,-1))
        pm,pg=arm.pieces(ts,k)
        arrays={}
        components={}
        for tag,base in (("linear",q),("exact",p)):
            integ={}
            for key,dv in (("mass",pm-qM),("geometry",pg-qG)):
                integrand=-adj[0]*base[0]*dv
                # Complex cumulative integration; each part uses the same Simpson rule.
                result=(cumulative_simpson(integrand.real,x=ts,axis=0,initial=0)
                        +1j*cumulative_simpson(integrand.imag,x=ts,axis=0,initial=0))
                integ[key]=result[-1]
                if tag=="linear":
                    arrays[key]={c:variance(result[i]) for c,i in indices.items()}
            if tag=="linear":
                a=det(point(adj,0),di)
                e=det(df,point(q,-1))
            else:
                a=det(point(adj,0),tuple(v-u for v,u in zip(p_i,q_i)))
                e=det(tuple(v-u for v,u in zip(p_f,q_f)),point(p,-1))
            db=integ["mass"]+integ["geometry"]+a+e
            components[tag]={"mass_evolution":variance(integ["mass"]),
                             "geometry_evolution":variance(integ["geometry"]),
                             "initial_basis":variance(a), "terminal_basis":variance(e),
                             "total_interference":variance(db)}
            if tag=="exact":
                components[tag]["quadratic_amplitude"]=count(db)
                components[tag]["max_amplitude_closure_error"]=float(max(abs(b1-b0-db)))
                components[tag]["reconstructed_delta_N"]=variance(db)+count(db)
        common_delta=count(det(q_f,point(p,-1)))-n0
        common_tangent=components["linear"]["total_interference"]-components["linear"]["terminal_basis"]
        exact=count(b1)-n0
        lin=components["linear"]["total_interference"]
        item={"delta_N":exact, "delta_percent":100*exact/n0,
              "linear_prediction_percent":100*lin/n0,"linear_to_finite_ratio":lin/exact,
              "common_quiet_analyzer":dict(exact_percent=100*common_delta/n0,
                                           linear_percent=100*common_tangent/n0),
              "components":components, "cumulative_linear_delta_N":arrays,
              "proper_clock_at_T1_255":arm.cut_clock(1.255)}
        if finite_difference:
            def homotopy(lam):
                V=lambda s:qV(s)+lam*delta(s)
                fake_i=make_basis(V,zeros,START)[0]
                fake_f=make_basis(V,zeros,END)[0]
                mode=evolve(V,fake_i,np.array([START,END]))
                return count(det(fake_f,point(mode,-1)))
            eps=.01
            step=(homotopy(eps)-homotopy(-eps))/(2*eps)
            item["centered_derivative_at_step_0.01"]=step
            if abs(step-lin)>2e-6:
                raise ValueError(f"Functional derivative check: {name}, {step}, {lin}")
        if components["exact"]["max_amplitude_closure_error"]>3e-8 or abs(exact-components["exact"]["reconstructed_delta_N"])>6e-8:
            raise ValueError("Finite Duhamel closure failed")
        # Two instrument kernels characterize all wave propagation for this clock, band and state.
        qh,ph=quiet.h(ts),arm.h(ts)
        dg=(np.exp(-2*(arm.offset+arm.int_h(ts)))-np.exp(-2*(quiet.offset+quiet.int_h(ts))))
        dU=(pm-qM)[:,0]-2.25*(ph**2-qh**2)-1.5*(arm.h(ts,1)-quiet.h(ts,1))
        kernel_integral=float(simpson(K0*dU+K2*dg,x=ts))
        item["two_kernel_bulk"]=kernel_integral
        if abs(kernel_integral-components["linear"]["mass_evolution"]-components["linear"]["geometry_evolution"])>1e-10:
            raise ValueError("Two-kernel projection identity")
        # Estimate error in the distributed 2049-point response table independently of the arm.
        tab=np.linspace(START,END,2049)
        tab0=np.interp(ts,tab,np.interp(tab,ts,K0)); tab2=np.interp(ts,tab,np.interp(tab,ts,K2))
        item["table_linear_interpolation_error"]=float(simpson(tab0*dU+tab2*dg,x=ts)-kernel_integral)
        arms[name]=item
    return {"momentum_order":order,"time_intervals":num_intervals,
            "cuts":cuts,"quiet_N":n0,"arms":arms,
            "designed_fixed_geometry_controls":design_controls(ts,K0,qV,q_i,q_f,n0,count)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capsule",type=Path,help="Path to the prior particle-response study")
    parser.add_argument("--write-kernel",type=Path,help="Write the fine background response weights as CSV")
    args=parser.parse_args()
    for name, digest in json.loads(Path(__file__).with_name("INPUTS.json").read_bytes()).items():
        if hashlib.sha256((args.capsule/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unexpected historical source: {name}")
    out={key:run_mesh(args.capsule,grid,order,time,fd,args.write_kernel if key=="fine" else None) for key,grid,order,time,fd in
         [("coarse","coarse",32,4096,False),("fine","fine",48,8192,True),
          ("fine_quadrature_check","fine",32,4096,False)]}
    print(json.dumps(out,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
