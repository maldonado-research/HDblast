"""Integrate the native five-dimensional B1 equations with the explicit first matter variation."""
import argparse
from functools import lru_cache
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import time
import zipfile

import numpy as np
from scipy.interpolate import CubicHermiteSpline

HERE = Path(__file__).resolve().parent
VENDOR = HERE / "vendor"
sys.path.insert(0, str(VENDOR))
import evolve_b1
import static_w
from quantum_shell import B, QuantumShell
from rolling_data import add_pulse


def archive(grid, arm):
    with zipfile.ZipFile(VENDOR / "references.zip") as saved:
        return np.genfromtxt(io.BytesIO(saved.read(f"{grid}_{arm}.csv")), delimiter=",", names=True)


def authenticate():
    for line in (VENDOR / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split()
        if hashlib.sha256((VENDOR / name).read_bytes()).hexdigest() != digest:
            raise ValueError(name)


class Source:
    """At fixed native T reconstruct the archived baseline clock and its explicit matter variation."""
    def __init__(self, path, reference, rb):
        self.clock = CubicHermiteSpline(reference["T"], reference["H0tau"], np.exp(reference["B_b"])/rb)
        raw = path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        with np.load(io.BytesIO(raw), allow_pickle=False) as d:
            keys = ("rho", "p", "j", "rad")
            v = np.stack([d[k] for k in keys], axis=1)
            dv = np.stack([d["d"+k] for k in keys], axis=1)
            self.support = (float(d["s"][0]), float(d["s"][-1]))
            self.spl = CubicHermiteSpline(d["s"], v, dv, extrapolate=False)

    @lru_cache(maxsize=6)
    def evaluate(self, T):
        s = float(self.clock(T))
        if s < self.support[0] - 1.e-10 or s > self.support[1] + 1.e-10:
            raise ValueError(f"Explicit matter history does not cover s={s}")
        s = np.clip(s, *self.support)
        return self.spl(s), self.spl(s, 1) * float(self.clock(T, 1))


class LiveMatter(QuantumShell):
    def initialize_source(self, strength, source):
        self.strength = strength
        self.source = source

    def totals(self, m, Ab, phib, T=None, **kwargs):
        base = super().totals(m, Ab, phib)
        if self.strength == 0:
            return base
        return base + self.mult*self.strength*self.source.evaluate(T)[0][:3]

    def totals_T(self, m, dm, Ab, phi, AT, phiT, eB, T=None):
        base = super().totals_T(m, dm, Ab, phi, AT, phiT, eB)
        if self.strength == 0:
            return base
        return base + self.mult*self.strength*self.source.evaluate(T)[1][:3]

    def rhs(self, m, Ab, phib, AT, phiT, eB, T=None):
        out = super().rhs(m, Ab, phib, AT, phiT, eB)
        if self.strength:
            _, p, j, _ = self.strength*self.source.evaluate(T)[0]
            vol = math.exp(3*(Ab-self.ref))
            out[1] += vol*j*phiT
            out[2] += 3*vol*AT*p
        return out


def run_case(grid, arm, strength, cutoff, scheme, tables, output, cache):
    authenticate()
    t0 = time.monotonic()
    reference = archive(grid, arm)
    ref = json.loads((VENDOR / "refs.json").read_bytes())[grid]
    fine = grid == "fine"
    # Only the static vacuum shell shooting cache depends on this path.
    evolve_b1.HERE = cache
    run = evolve_b1.Run(delta=.1, c=.5975949350280132, d=-3.106933495673783, dc=.01,
                        dzf=.0005 if fine else .001, dzc=.002 if fine else .004, L=16.,
                        kappa=10., project=None, cfl=.5, order=4, xc=1.9,
                        matter=lambda rb: LiveMatter(rb, ref, nodes=24, feedback=True, decay=0.),
                        log=lambda msg: None)
    sector = run.matter
    src = Source(tables/f"{grid}_{arm}_{cutoff}_{scheme}.npz", reference, run.rb) if strength else None
    sector.initialize_source(strength, src)
    with zipfile.ZipFile(VENDOR/"FORK_INPUTS.zip") as forks:
        with np.load(io.BytesIO(forks.read(grid+".npz")), allow_pickle=False) as cp:
            u, m, T, s = cp["bulk"], cp["matter"], float(cp["T"]), float(cp["tau"])
            run.mode = str(cp["mode"].item())
    if run.mode != "full":
        raise ValueError("Expected native full-field postfork state")
    if arm != "quiet":
        options = {"control": -.85, "spread": .25} if arm == "far" else {}
        u, match = add_pulse(run, T, u, m, -.001 if arm == "minus" else .001, options)
    else:
        match = {}
    den = run.rb
    dt = run.p["cfl"]*run.zp.min()
    # One record per native RK4 step. The RHS supplies its full, stage-one geometry.
    rows = []
    z = run.z
    ii = run.inner
    near = z[ii] > -1.
    while True:
        k1, r1, f = run.rhs(T, u, m)
        Ab, Bb, phi, AT, BT, phit = u[:, -1]
        eB = math.exp(Bb)
        h, v = AT*den/eB, phit*den/eB
        vol = math.exp(3*(Ab-ref))
        rho, p, j = np.asarray(sector.totals(m, Ab, phi, T=T))/sector.mult
        rhoT, pT, jT = np.asarray(sector.totals_T(m, r1, Ab, phi, AT, phit, eB, T=T))/sector.mult
        hc, mc, sh, sm = run.constraints(f, u[3], m, r1, k1[5,-1], T)
        hn, mn = hc/sh, mc/sm
        sig, sg1 = run.ten.s(phi), run.ten.s1(phi)
        junction = np.array([f["gA"], f["gB"], f["gF"]])
        junction_T = run.neumann_t(Ab,Bb,BT,AT,phi,phit,k1[5,-1],m,r1,T)
        explicit = strength*src.evaluate(T)[0] if strength else np.zeros(4)
        rad = m[0]+explicit[3]
        wy = (AT/eB)**2-(sig+sector.mult*rho)**2/36-(phit/eB)**2/12+(sg1+sector.mult*j)**2/48-static_w.U(phi)/6
        rows.append(dict(T=T,s=s,A=Ab,B=Bb,phi=phi,phit=phit,h=h,v=v,rho=rho,p=p,j=j,
                         work=m[1],pressure_work=m[2],E=rho*vol,sigma=sig,sigma1=sg1,inside_phi=u[2,-2:-7:-1].copy(),
                         Hshell=max(abs(hn[-6:])),Mshell=max(abs(mn[-6:])),Hnear=max(abs(hn[ii][near])),
                         Mnear=max(abs(mn[ii][near])),Hsigned=hn[-6:].copy(),Msigned=mn[-6:].copy(),
                         Hraw=hc[-1],Mraw=mc[-1],actual_grad=[f[k][-1] for k in ("A_Z","B_Z","phi_Z")],
                         junction=junction,junction_T=junction_T,rhoT=rhoT,pT=pT,jT=jT,rad=rad,
                         nonrad=rho-rad,explicit=explicit,Wy_over_H02=wy*den**2))
        if not np.isfinite(np.fromiter((rho,p,j,h,v), dtype=float)).all():
            raise ValueError("Run became nonfinite")
        if s >= 2.105:
            break
        if T >= 1.62:
            raise ValueError("Native source window exhausted")
        k2, r2, _ = run.rhs(T+dt/2, u+dt/2*k1, m+dt/2*r1)
        k3, r3, _ = run.rhs(T+dt/2, u+dt/2*k2, m+dt/2*r2)
        k4, r4, _ = run.rhs(T+dt, u+dt*k3, m+dt*r3)
        u = u + dt/6*(k1+2*k2+2*k3+k4)
        m = m + dt/6*(r1+2*r2+2*r3+r4)
        T += dt
        rn = run.ref_nodes(T, np.arange(run.nfix))
        for i, name in enumerate(("A","B","phi","A_T","B_T","phi_T")):
            u[i,:run.nfix] = rn[name]
        s += dt*.5*(eB+math.exp(u[1,-1]))/den
    output.mkdir(parents=True,exist_ok=True)
    tag = f"{grid}_{arm}_zero" if strength == 0 else f"{grid}_{arm}_{cutoff}_{scheme}_{'p' if strength>0 else 'm'}{abs(strength):g}"
    target = output/(tag+".npz")
    with tempfile.NamedTemporaryFile(dir=output, delete=False) as tf:
        np.savez_compressed(tf, **{key:np.array([row[key] for row in rows]) for key in rows[0]},
                            b=B,rb=den,Aref=ref,source_sha256=src.sha256 if strength else "",z_inside=z[-2:-7:-1],lam=strength,dt=dt,match=json.dumps(match))
    os.replace(tf.name, target)
    print(f"{tag}  steps={len(rows)-1}  seconds={time.monotonic()-t0:.2f}", flush=True)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("grid", choices=("coarse","fine"))
    ap.add_argument("arm", choices=("quiet","plus","minus","far"))
    ap.add_argument("strength", type=float)
    ap.add_argument("cutoff", type=int, nargs="?", default=120)
    ap.add_argument("scheme", nargs="?", choices=("bare","static"), default="static")
    ap.add_argument("--tables", type=Path, default=HERE.parent/"source"/"tables")
    ap.add_argument("--output", type=Path, default=HERE/"runs")
    ap.add_argument("--cache", type=Path, default=Path(tempfile.gettempdir())/"hdblast-live-static-cache")
    a=ap.parse_args()
    run_case(a.grid,a.arm,a.strength,a.cutoff,a.scheme,a.tables,a.output,a.cache)


if __name__=="__main__":
    main()
