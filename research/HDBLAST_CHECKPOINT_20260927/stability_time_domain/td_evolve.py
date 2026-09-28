"""Small-amplitude time-domain evolutions of the registered 1+1 PDE system on a static background.

Initial data: time-symmetric (all velocities zero) and constraint-satisfying.  The scalar
disturbance is two smooth compact bumps in the bulk, f = eps*(b1 + beta*b2); a=b is obtained by
solving the discrete linearised Hamiltonian constraint with the solver's own finite-difference
operators and quintic-Hermite shell closure; beta is tuned so that a vanishes at the shell, which
(with f=0 near the shell) makes a vanish in a neighbourhood of the shell and removes the
junction-compatibility (corner) obstruction.  The momentum constraint then holds identically
(a=b, zero velocities).

Evolution: classical RK4 with dt=0.4*hmin (as in the registered solver), using either the unchanged
nonlinear RHS (default) or the verified linear matrix (--linear).
Recorded: shell scalar f_b (gauge invariant under the residual shell-preserving conformal gauge),
shell Hubble perturbation pa_b-b_b (also gauge invariant), constraint diagnostics.
"""
import argparse, json, math, time
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
import td_core as T


def bump(z, c, w):
    x = (z - c)/w
    out = np.zeros_like(z); m = np.abs(x) < 1
    out[m] = np.exp(-1/(1 - x[m]**2) + 1)
    return out


def hamiltonian_operator(s):
    """Sparse matrices (Ka, Kf) with H_lin = Ka@a + Kf@f for a=b, zero velocities, on all n nodes."""
    n = s.n; I = sparse.eye(n, format='csr')
    E = sparse.csr_matrix(([1.], ([0], [n - 1])), shape=(1, n))
    # boundary data ga = rb*(s0*b_N + s10*f_N)/6 with b=a
    ga_a = s.rb*s.s0/6*E; ga_f = s.rb*s.s10/6*E
    gf_a = -s.rb*s.s10/2*E; gf_f = -s.rb*s.s20/2*E
    G1 = sparse.csr_matrix(s.g1[:, None]); G2 = sparse.csr_matrix(s.g2[:, None])
    Dz_a = s.D1 + G1@ga_a; Dzz_a = s.D2 + G2@ga_a          # derivative of a (and of b=a)
    Dz_a_from_f = G1@ga_f; Dzz_a_from_f = G2@ga_f
    Dz_f = s.D1 + G1@gf_f; Dz_f_from_a = G1@gf_a
    hc = sparse.diags(s.hc); pz = sparse.diags(s.phz); r2 = sparse.diags(s.rho2)
    U0 = sparse.diags(s.pot[0]); U1 = sparse.diags(s.pot[1])
    # H = -2 rho^2 (2 b U + U' f) - 18 hc az + 6 hc bz - 6 azz - 2 phz fz   (a=b, velocities 0)
    Ka = -4*r2@U0 - 12*hc@Dz_a - 6*Dzz_a - 2*pz@Dz_f_from_a
    Kf = -2*r2@U1 - 12*hc@Dz_a_from_f - 6*Dzz_a_from_f - 2*pz@Dz_f
    return Ka.tocsr(), Kf.tocsr()


def two_bump_data(s, z1, w1, z2, w2):
    n = s.n; Ka, Kf = hamiltonian_operator(s)
    free = np.arange(2, n)
    K = Ka[free][:, free].tocsc()
    sols = []
    for (zc, w) in [(z1, w1), (z2, w2)]:
        f = bump(s.z, zc, w); f[:2] = 0
        a = np.zeros(n); a[free] = spsolve(K, -(Kf@f)[free])
        sols.append((a, f))
    beta = -sols[0][0][-1]/sols[1][0][-1]
    a = sols[0][0] + beta*sols[1][0]; f = sols[0][1] + beta*sols[1][1]
    v = np.zeros((6, n)); v[0] = a; v[1] = a; v[2] = f
    scale = np.max(np.abs(f)); v /= scale
    H, M, Hs, Ms = T.linear_constraints(s, v)
    mask = (s.z > -.8*s.L) & (np.arange(n) > 5)
    info = dict(beta=float(beta), a_shell=float(v[0, -1]), f_shell=float(v[2, -1]),
                max_abs_a=float(np.max(np.abs(v[0]))), max_abs_a_near_shell=float(np.max(np.abs(v[0][s.z > -0.005]))),
                initial_H_rel_l2=float(np.linalg.norm(H[mask])/np.linalg.norm(Hs[mask])),
                initial_M_rel_l2=float(np.linalg.norm(M[mask])/max(np.linalg.norm(Ms[mask]), 1e-300)),
                initial_H_max_abs=float(np.max(np.abs(H[mask]))))
    return v, info


def scalar_energy(s, v, wq, near):
    """Weighted scalar energy E=sum dz rho^3 (pf^2 + fz^2 + rho^2 U'' f^2), whole domain and z>-0.1.
    For a scalar decoupled from the metric (the +1 branch, where the background gradient is ~1e-4)
    the substitution X=e^{3t/2} f removes the damping, so a purely continuum (Re lambda=-3/2)
    spectrum makes E decay like e^{-3t} up to oscillating cross terms.  Positive only where U''>0."""
    f, pf = v[2], v[5]
    gf = -s.rb*(s.s10*v[1, -1] + s.s20*f[-1])/2
    fz = s.D1@f + s.g1*gf
    dens = s.rho**3*(pf*pf + fz*fz + s.rho2*s.pot[2]*f*f)
    return float(np.sum(wq*dens)), float(np.sum((wq*dens)[near]))


def compensated_energy(s, v, wq, t):
    """E_X for X=e^{3t/2} f.  If the scalar decouples from the metric, X_tt = rho^-3 (rho^3 X_z)_z
    - (V - 9/4) X with V=rho^2 U'', and the Robin junction X_z = -rho_b (sigma''/2) X gives the conserved
    E_X = 1/2 sum dz rho^3 [X_t^2 + X_z^2 + (V - 9/4) X^2] + rho_b^4 sigma''/4 X_b^2 (Dirichlet far end).
    Constancy of E_X is equivalent to every excited scalar mode having Re(lambda) = -3/2."""
    f, pf = v[2], v[5]
    gf = -s.rb*(s.s10*v[1, -1] + s.s20*f[-1])/2
    fz = s.D1@f + s.g1*gf
    bulk = .5*np.sum(wq*s.rho**3*((pf + 1.5*f)**2 + fz*fz + (s.rho2*s.pot[2] - 2.25)*f*f))
    return float(np.exp(3*t)*(bulk + s.rb**4*s.s20/4*f[-1]**2))


def metric_norm(s, v, wq):
    return float(np.sqrt(np.sum(wq*(v[0]**2 + v[1]**2 + v[3]**2 + v[4]**2))))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bg', choices=['plus', 'original'], required=True)
    p.add_argument('--hmin', type=float, default=2e-4); p.add_argument('--L', type=float, default=6.)
    p.add_argument('--stretch', type=float, default=.05); p.add_argument('--tdet', type=float, default=1e-3)
    p.add_argument('--tf', type=float, default=4.); p.add_argument('--eps', type=float, default=1e-6)
    p.add_argument('--z1', type=float, default=-.03); p.add_argument('--w1', type=float, default=.02)
    p.add_argument('--z2', type=float, default=-.12); p.add_argument('--w2', type=float, default=.04)
    p.add_argument('--linear', action='store_true'); p.add_argument('--sample', type=float, default=.002)
    p.add_argument('--diag-every', type=float, default=.05); p.add_argument('--s20-scale', type=float, default=1.)
    p.add_argument('--max-f', type=float, default=1e-2)
    p.add_argument('--phi-poly', action='store_true')
    p.add_argument('--out', required=True)
    a = p.parse_args(); Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    bg = T.plus_background(a.tdet) if a.bg == 'plus' else T.original_background(a.tdet)
    s = T.GeneralSolver(bg, tdet=a.tdet, hmin=a.hmin, L=a.L, stretch=a.stretch, precise_eta=not a.phi_poly)
    if a.s20_scale != 1.: s.s20 *= a.s20_scale
    v0, info = two_bump_data(s, a.z1, a.w1, a.z2, a.w2)
    v = a.eps*v0
    if a.linear:
        M = s.linear_matrix(); rhs = lambda x: (M@x.ravel()).reshape(6, s.n)
    else:
        rhs = s.rhs
    dt = .4*s.zp.min(); steps = math.ceil(a.tf/dt); dt = a.tf/steps
    every = max(1, round(a.sample/dt)); dev = max(every, round(a.diag_every/dt)//every*every)
    mask = (s.z > -.8*s.L) & (np.arange(s.n) > 5)
    ts, fb, hb, ab, fmax = [], [], [], [], []; diag = []; t0 = time.monotonic(); reason = 'tf'
    wq = s.zp.copy(); wq[0] *= .5; wq[-1] *= .5; wq[:2] = 0   # trapezoid weights in the index coordinate
    near = s.z > -0.1; E, En, MN, EX = [], [], [], []
    for i in range(steps + 1):
        if i % every == 0 or i == steps:
            t = i*dt; ts.append(t); fb.append(v[2, -1]); hb.append(v[3, -1] - v[1, -1]); ab.append(v[0, -1])
            fmax.append(float(np.max(np.abs(v[2]))))
            e1, e2 = scalar_energy(s, v, wq, near); E.append(e1); En.append(e2); MN.append(metric_norm(s, v, wq)); EX.append(compensated_energy(s, v, wq, t))
            if i % dev == 0 or i == steps:
                d, H, Mm, cc = s.diagnostics(v, 0.)
                Hl, Ml, Hs, Ms = (T.linear_constraints if a.linear else T.nonlinear_constraints)(s, v)
                diag.append(dict(t=t, registered_H_max=d['hamiltonian_max'], registered_M_max=d['momentum_max'],
                                 H_rel_l2=float(np.linalg.norm(Hl[mask])/max(np.linalg.norm(Hs[mask]), 1e-300)),
                                 M_rel_l2=float(np.linalg.norm(Ml[mask])/max(np.linalg.norm(Ms[mask]), 1e-300)),
                                 H_rel_max=float(np.max(np.abs(Hl[mask]))/max(np.max(Hs[mask]), 1e-300)),
                                 M_rel_max=float(np.max(np.abs(Ml[mask]))/max(np.max(Ms[mask]), 1e-300)),
                                 max_abs_f=float(np.max(np.abs(v[2]))), max_abs_a=float(np.max(np.abs(v[0]))),
                                 max_abs_b=float(np.max(np.abs(v[1])))))
            if not np.isfinite(v).all() or np.max(np.abs(v[2])) > a.max_f: reason = 'guard'; break
        if i == steps: break
        k1 = rhs(v); k2 = rhs(v + .5*dt*k1); k3 = rhs(v + .5*dt*k2); k4 = rhs(v + dt*k3)
        v = v + dt*(k1 + 2*k2 + 2*k3 + k4)/6
    out = dict(script='td_evolve.py', args=vars(a), background=s.meta, nodes=s.n, dt=dt, pot_mode=s.pot_mode, stop_reason=reason,
               initial_data=info, runtime_seconds=time.monotonic() - t0, diagnostics=diag)
    np.savez_compressed(a.out + '.npz', t=np.array(ts), f_b=np.array(fb), h_b=np.array(hb), a_b=np.array(ab),
                        fmax=np.array(fmax), E=np.array(E), E_near=np.array(En), E_X=np.array(EX), metric_norm=np.array(MN), z=s.z, final=v, initial=v0)
    T.dump(a.out + '.json', out)
    print(json.dumps(dict(out=a.out, nodes=s.n, steps=steps, runtime=out['runtime_seconds'], reason=reason, init=info)), flush=True)


if __name__ == '__main__':
    main()
