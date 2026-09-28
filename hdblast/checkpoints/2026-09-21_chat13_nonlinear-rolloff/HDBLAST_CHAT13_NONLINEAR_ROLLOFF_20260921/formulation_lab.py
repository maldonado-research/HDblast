#!/usr/bin/env python3
"""Chat 13 - formulation laboratory: linearised semi-discrete operator of the 1+1 system about the static shell, its eigenvalues,
and the constraint content of the unstable eigenvectors, for different boundary treatments.  Small domain, coarse grid: a design tool."""
import sys, math, time, json
import numpy as np
import rolloff5d as R

def build(d=0.0, dz=2e-3, L=1.2, bc="neumann", ko_eps=0.05, damp=0.0, t_det=1e-3):
    Ip = 1.0357712571566784; c = 2/Ip - 4/3; ten = R.Tension(t_det, c, d)
    ph_h, y_b, _ = R.solve_shell(ten, guess=(math.log10(0.2207*t_det**1.8), 8.2282 + 0.9*math.log(1e-3/t_det))); N = int(round(L/dz)); z = -L + dz*np.arange(N + 1)
    S = R.static_on_grid(ph_h, y_b, z, hz=min(dz/4, 2.5e-5)); rho, Hc, phs, phz = S["rho"], S["Hc"], S["phi"], S["phiz"]
    rho2 = rho*rho; rb, pb, Hb = rho[-1], phs[-1], Hc[-1]
    Aszz = Hc**2 - 1 - phz**2/3
    def ext(f, g=None):                                  # ghost values: Neumann (g given) or cubic extrapolation (g None)
        if g is None:
            e1 = 4*f[-1] - 6*f[-2] + 4*f[-3] - f[-4]; e2 = 4*e1 - 6*f[-1] + 4*f[-2] - f[-3]
        else:
            e1 = f[-2] + 2*dz*g; e2 = f[-3] + 4*dz*g
        return np.concatenate(([0.0, 0.0], f, [e1, e2]))
    D1 = lambda F: (F[:-4] - 8*F[1:-3] + 8*F[3:-1] - F[4:])/(12*dz)
    D2 = lambda F: (-F[:-4] + 16*F[1:-3] - 30*F[2:-2] + 16*F[3:-1] - F[4:])/(12*dz*dz)
    def KO(f, g=None):
        F = ext(f, g); F = np.concatenate(([0.0], F, [2*F[-1] - F[-2]]))
        return +ko_eps*(F[:-6] - 6*F[1:-5] + 15*F[2:-4] - 20*F[3:-3] + 15*F[4:-2] - 6*F[5:-1] + F[6:])/(64*dz)
    def rhs(state):
        a, b, f, pa, pbv, pf = state
        eB = rb*math.exp(b[-1]); p = pb + f[-1]
        gA = (ten.s(p)*eB - ten.s(pb)*rb)/6; gF = -(ten.s1(p)*eB - ten.s1(pb)*rb)/2
        gB = gA if bc == "neumann" else None
        Fa, Fb, Ff = ext(a, gA), ext(b, gB), ext(f, gF)
        az, bz, fz = D1(Fa), D1(Fb), D1(Ff); azz, bzz, fzz = D2(Fa), D2(Fb), D2(Ff)
        dU, dU1, U0, U10 = R.dU_exact(phs, f); e2b = np.expm1(2*b)
        SU = rho2*(e2b*(U0 + dU) + dU); SU1 = rho2*(e2b*(U10 + dU1) + dU1)
        kin = 2*pa + pa*pa; grad = 2*Hc*az + az*az
        ra = azz - 3*kin + 3*grad + (2/3)*SU
        rb_ = bzz + 3*kin - 3*grad - 0.5*pf*pf + 0.5*(2*phz*fz + fz*fz) - (1/3)*SU
        rf = fzz - 3*(1 + pa)*pf + 3*((Hc + az)*(phz + fz) - Hc*phz) - SU1
        M = -3*D1(ext(pa, None)) - 3*(1 + pa)*(az - bz) + 3*(Hc + az)*pbv - pf*(phz + fz)
        out = [pa + KO(a, gA), pbv + KO(b, gB), pf + KO(f, gF), ra + KO(pa), rb_ + KO(pbv), rf + KO(pf)]
        if bc == "hamB":                                 # b at the shell from the Hamiltonian constraint (first-order ODE); pbv[-1] slaved
            G = (ten.s(p)*eB)/6; Gf = -(ten.s1(p)*eB)/2
            num = 2*SU[-1] - 6*kin[-1] + 6*(G*G - Hb*Hb) + 6*azz[-1] + pf[-1]**2 + (Gf*Gf - phz[-1]**2)
            Bt = num/(6*(1 + pa[-1]))
            out[1][-1] = Bt; out[4][-1] = (Bt - pbv[-1])/ (5*dz)      # relax pbv at the shell toward the constraint value
        if damp:
            out[3] = out[3] + damp*M/3.0*0                # placeholder for constraint-damping experiments
        for o in out: o[0:2] = 0.0
        return out, M
    return rhs, z, N

def spectrum(bc, d=0.0, dz=2e-3, L=1.2, t_det=1e-3):
    rhs, z, N = build(d=d, dz=dz, L=L, bc=bc, t_det=t_det); n = N + 1
    x0 = [np.zeros(n) for _ in range(6)]
    def flat(o): return np.concatenate(o)
    J = np.zeros((6*n, 6*n)); e = 1e-7
    for j in range(6*n):
        xp = [v.copy() for v in x0]; xm = [v.copy() for v in x0]
        xp[j // n][j % n] += e; xm[j // n][j % n] -= e
        J[:, j] = (flat(rhs(xp)[0]) - flat(rhs(xm)[0]))/(2*e)
    w, V = np.linalg.eig(J)
    order = np.argsort(-w.real); out = []
    for i in order[:8]:
        v = V[:, i]; st = [np.real(v[k*n:(k+1)*n]) for k in range(6)]
        sc = max(np.abs(st[2]).max(), np.abs(st[1]).max(), 1e-300)
        st = [s_/sc*1e-6 for s_ in st]
        M = rhs(st)[1]
        out.append(dict(eig=[float(w[i].real), float(w[i].imag)], constraint_ratio=float(np.abs(M[10:-6]).max()/1e-6),
                        where_f_peaks=float(z[np.argmax(np.abs(st[2]))])))
    return out

if __name__ == "__main__":
    res = {}
    for bc in ("neumann", "hamB"):
        t0 = time.time(); res[bc] = spectrum(bc)
        print("== bc = %s  (%.0f s)" % (bc, time.time() - t0))
        for r in res[bc]: print("   eig = %+9.3f %+9.3fi   |M|/|mode| = %.2e   f peaks at z = %.3f" % (r["eig"][0], r["eig"][1], r["constraint_ratio"], r["where_f_peaks"]))
    json.dump(res, open("runs/formulation_lab.json", "w"), indent=1)
