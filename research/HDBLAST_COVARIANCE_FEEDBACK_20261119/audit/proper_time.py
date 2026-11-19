"""Independent proper-time reconstruction of the scalar sources on four archived closed-field walls."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.interpolate import CubicHermiteSpline, CubicSpline

WALL = "research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip"
SUMMARY = "research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/{}/quiet/rolling_summary.json"
DIGESTS = {
    WALL: "7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15",
    SUMMARY.format("fine"): "abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810",
    SUMMARY.format("coarse"): "91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c",
}
ARMS = ("quiet", "plus", "minus", "far")
CUT = 2.10
LABELS = ("rho", "p", "j")


def checked(repo, name):
    data = (repo/name).read_bytes()
    if hashlib.sha256(data).hexdigest() != DIGESTS[name]:
        raise ValueError(name)
    return data


def table(z, grid, part):
    with z.open(f"{grid}/y0_b1/{part}.csv") as stream:
        return np.genfromtxt(stream, delimiter=",", names=True)


def interpolate(rows, reference):
    s = rows["H0tau"]
    loga = CubicHermiteSpline(s, rows["A_b"] - reference, rows["H_over_H0"], extrapolate=False)
    phi = CubicHermiteSpline(s, rows["phi_b"], rows["v_over_H0"], extrapolate=False)
    def wall(t):
        return np.exp(loga(t)), loga(t, 1), phi(t)-0.5, phi(t, 1)
    return wall


def evolve(q, wall, domain, initial):
    shape = initial.shape
    def rhs(s, z):
        a, h, d, _ = wall(s)
        field, velocity = z.reshape(shape)
        return np.stack((velocity, -3*h*velocity - ((q/a)**2 + 1e4*d**2)*field)).ravel()
    sol = solve_ivp(rhs, domain, initial.ravel(), method="DOP853", rtol=2e-12, atol=2e-14,
                    max_step=.002, dense_output=True)
    if not sol.success:
        raise ValueError(sol.message)
    def sample(s):
        return sol.sol(s).reshape((*shape, *np.shape(s)))
    return sample


def mode_matrices(q, w, wall, s0, sf):
    ai, _, d0, _ = wall(s0)
    ei = np.hypot(q/ai, 100*d0)
    initial_field = (2*ai**3*ei)**(-0.5)
    mode = evolve(q, wall, (s0, sf), np.stack((initial_field, -1j*ei*initial_field)))(sf)
    a, _, _, _ = wall(sf)
    field, momentum = mode[0], mode[1]*a**3
    x, p, c = abs(field)**2, abs(momentum)**2, (field.conj()*momentum).real
    impulse = -2*(w@c)/(w@x)
    states = {
        "retain": (x,p,c), "erase": (x,p,np.zeros_like(c)), "reflect": (x,p,-c),
        "shear": (x,p+2*impulse*c+impulse**2*x,c+impulse*x),
    }
    return states, float(impulse)


def observables(q, w, wall, s, x, p, c):
    # x,p,c may have node and time axes. All public array evaluations use both.
    a, _, d, _ = wall(s)
    k2 = (q[:,None]/a)**2
    m2 = 1e4*d*d
    e = np.sqrt(m2+k2)
    dx = x - 1/(2*a**3*e)
    dv = (p-a**3*e/2)/a**6
    return np.stack((w@((dv+e**2*dx)/2), w@((dv-(m2+k2/3)*dx)/2), 1e4*d*(w@dx)))


def ready(wall, fork, q, w, states):
    a, h, d, v = wall(fork)
    k2 = (q/a)**2
    mass2 = 1e4*d**2
    measure = w/a**3
    reference = states["erase"]
    output = {}
    for label, (x,p,c) in states.items():
        column = (x[:,None],p[:,None],c[:,None])
        src = observables(q, w, wall, np.array([fork]), *column)[:,0]
        dj = 2e4*d*(measure@c)
        dp = -2*(measure@((mass2+2*k2/3)*c))
        # The integrated P difference vanishes in all four states. This removes
        # its expansion-rate term from the source derivative difference.
        delta_rho = (w@(p-reference[1]))/(2*a**6)
        e = np.sqrt(mass2+k2)
        delta_n = (p-reference[1])/(2*a**3*e)
        output[label] = {
            "sources": dict(zip(LABELS,map(float,src))),
            "derivatives_relative_to_erase": {"j_prime":float(dj),"p_prime":float(dp),
                                              "rho_prime":0.,"rho_double_prime":float(v*dj-3*h*dp)},
            "min_determinant": float(np.min(x*p-c*c)), "energy_change":float(delta_rho),
            "max_abs_mode_occupation_change":float(max(abs(delta_n))),
            "comoving_number_change":float(w@delta_n),
        }
    return output


def replay_arm(q, w, wall, start, states, count):
    ai = wall(start)[0]
    f0 = np.stack((np.ones_like(q),np.zeros_like(q)))
    v0 = np.stack((np.zeros_like(q),np.full_like(q,1/ai**3)))
    propagator = evolve(q, wall, (start,CUT), np.stack((f0,v0)))
    dt = 1e-5
    times = np.r_[np.linspace(start,CUT,count),start+dt*np.arange(5)]
    a,h,d,phi_prime = wall(times)
    fg,DS = propagator(times)
    f,g = fg
    df,dg = a**3*DS
    det = float(np.max(abs(f*dg-g*df-1)))
    output, corners = {}, {}
    for label,(x,p,c) in states.items():
        x, p, c = x[:,None], p[:,None], c[:,None]
        xx = f*f*x + g*g*p + 2*f*g*c
        pp = df*df*x + dg*dg*p + 2*df*dg*c
        cc = f*df*x + g*dg*p + (f*dg+g*df)*c
        rho,pr,j = observables(q,w,wall,times,xx,pp,cc)
        before = slice(0,count)
        ws,wp = a**3*j*phi_prime,3*a**3*h*pr
        work = float(simpson(ws[before],x=times[before]))
        expansion = float(simpson(wp[before],x=times[before]))
        energy = a**3*rho
        error = float(energy[count-1]-energy[0]-work+expansion)
        endpoints = dict(zip(LABELS,map(float,(rho[count-1],pr[count-1],j[count-1]))))
        output[label] = {"sources_at_readout":endpoints,"scalar_work":work,"expansion_work":expansion,
                         "comoving_energy_change":float(energy[count-1]-energy[0]),
                         "relative_work_residual":abs(error)/max(1.,abs(work),abs(expansion))}
        corners[label] = np.stack((rho[count:],pr[count:],j[count:]))
    for label,row in output.items():
        diff = corners[label]-corners["erase"]
        slope = (-25*diff[:,0]+48*diff[:,1]-36*diff[:,2]+16*diff[:,3]-3*diff[:,4])/(12*dt)
        acceleration = (35*diff[:,0]-104*diff[:,1]+114*diff[:,2]-56*diff[:,3]+11*diff[:,4])/(12*dt**2)
        row["endpoint_differences_relative_to_erase"] = {
            "p_prime":float(slope[1]),"j_prime":float(slope[2]),"rho_prime":float(slope[0]),
            "rho_double_prime":float(acceleration[0]),
        }
    return dict(max_transfer_determinant_error=det,states=output)


def compute(repo,nodes,count):
    x, weights = np.polynomial.legendre.leggauss(nodes)
    q = 4*(1+x)
    w = 4*weights*q*q/(2*np.pi**2)
    grids = {}
    with zipfile.ZipFile(io.BytesIO(checked(repo,WALL))) as archive:
        for grid in ("fine","coarse"):
            ref = json.loads(checked(repo,SUMMARY.format(grid)))["matter"]["events"][0]["A_star"]
            pre = table(archive,grid,"pre")
            arms = {arm:table(archive,grid,arm) for arm in ARMS}
            geom = interpolate(np.r_[pre,arms["quiet"][:1]], ref)
            start, fork = pre["H0tau"][0], arms["quiet"]["H0tau"][0]
            states,kick = mode_matrices(q,w,geom,start,fork)
            initials = ready(geom,fork,q,w,states)
            native = {k:float(arms["quiet"][k+"_hat"][0]) for k in LABELS}
            record={"fork_s":float(fork),"shear":kick, "C_min_max":list(map(float,(min(states["retain"][2]),max(states["retain"][2])))),
                    "fork":initials, "native_fork_sources":native, "histories":{}}
            for arm,rows in arms.items():
                history = replay_arm(q,w,interpolate(rows,ref),fork,states,count)
                history["retained_native_sources_at_readout"] = {
                    k:float(CubicSpline(rows["H0tau"],rows[k+"_hat"],extrapolate=False)(CUT))
                    for k in LABELS
                }
                record["histories"][arm] = history
            grids[grid]=record
    return {"nodes":nodes,"samples":count,"cut":CUT,"grids":grids}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repository",type=Path)
    parser.add_argument("--nodes",type=int,default=24)
    parser.add_argument("--samples",type=int,default=2049)
    args=parser.parse_args()
    if args.nodes<4 or args.samples<65 or args.samples%2!=1:
        raise ValueError("Choose at least four Gaussian nodes and an odd work quadrature size >= 65")
    print(json.dumps(compute(args.repository,args.nodes,args.samples),indent=2,allow_nan=False))
