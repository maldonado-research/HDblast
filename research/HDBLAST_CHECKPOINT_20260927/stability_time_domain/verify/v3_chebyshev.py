"""Audit V3: independent discretisation of the continuum linear problem.
Multi-domain Chebyshev collocation (not the audited stretched-grid finite differences / quintic Hermite
closure) of the linearised 1+1 system whose coefficients were re-derived symbolically in v1:
  a_tt = a'' + 6hc a' + (4/3)rho^2 U b + (2/3)rho^2 U' f - 6 a_t
  b_tt = b'' - 6hc a' + phz f' - (2/3)rho^2 U b - (1/3)rho^2 U' f + 6 a_t
  f_tt = f'' + 3hc f' + 3phz a' - 2rho^2 U' b - rho^2 U'' f - 3 f_t
with the linearised junctions at z=0:  a' = b' = rho_b(s0 b + s10 f)/6,  f' = -rho_b(s10 b + s20 f)/2,
Dirichlet a=b=f=0 at z=-L, C^1 matching at subdomain interfaces.  Eigenproblem A x = lam B x
(x=(q,p), singular B on boundary rows).  Backgrounds: +1 branch from the independent v2 solver;
original shell from the frozen registered solver (hash-checked).
Also: decoupled-scalar Sturm-Liouville check and the Liouville-transform bound.
Writes verify/v3_chebyshev.json."""
import json, math, sys, time, hashlib, importlib.util
from pathlib import Path
import numpy as np, scipy.linalg as sl
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import v2_background as V2

FROZEN = HERE.parent/'frozen_input'/'registered_solver.py'
assert hashlib.sha256(FROZEN.read_bytes()).hexdigest().startswith('4147239c40b63c2d')
spec = importlib.util.spec_from_file_location('rs', FROZEN); RS = importlib.util.module_from_spec(spec); spec.loader.exec_module(RS)
C = V2.C
UC = np.array([-1/6, 4/3, -5/3, -4/9, 17/18, 0., -2/27])


def Uphi(p, k):
    c = UC.copy()
    for _ in range(k): c = np.arange(1, len(c))*c[1:]
    return np.polynomial.polynomial.polyval(p, c)


def cheb(N):
    x = np.cos(np.pi*np.arange(N + 1)/N); c = np.hstack([2, np.ones(N - 1), 2])*(-1)**np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T; dX = X - X.T
    D = np.outer(c, 1/c)/(dX + np.eye(N + 1)); D -= np.diag(D.sum(1))
    return x[::-1], D[::-1, ::-1]      # increasing nodes


def background(kind, d=1e-3):
    if kind == 'plus':
        roots, _ = V2.solve_plus(d); lx = roots[0]
        samp, info = V2.sampler(d, lx)
        def bg(z):
            s = samp(z); e = s['eta']
            return dict(rho=s['rho'], hc=s['hc'], phz=s['phiz'], U=V2.Ud(e), U1=V2.Ud(e, 1), U2=V2.Ud(e, 2)), 1 + e[-1] if np.ndim(e) else 1 + e
        pb = 1 + info['eta_b']; rb = info['rho_b']
        def bgf(z):
            s = samp(z); e = s['eta']
            return dict(rho=s['rho'], hc=s['hc'], phz=s['phiz'], U=V2.Ud(e), U1=V2.Ud(e, 1), U2=V2.Ud(e, 2))
        return bgf, pb, rb, dict(eta_h=-10**lx, **info)
    samp, meta = RS.shell_background(d)
    def bgf(z):
        s = samp(np.asarray(z)); p = s['phi']
        return dict(rho=s['rho'], hc=s['Hc'], phz=s['phiz'], U=Uphi(p, 0), U1=Uphi(p, 1), U2=Uphi(p, 2))
    return bgf, meta['phi_b'], meta['rho_b'], meta


def build(bgf, pb, rb, d, breaks, N, s20_scale=1.0):
    x, Dx = cheb(N); K = len(breaks) - 1; M = K*(N + 1)
    z = np.concatenate([breaks[k] + (x + 1)*(breaks[k + 1] - breaks[k])/2 for k in range(K)])
    D1 = sl.block_diag(*[Dx*2/(breaks[k + 1] - breaks[k]) for k in range(K)]); D2 = D1@D1
    B = bgf(z); r2 = B['rho']**2
    s0 = 2*(1 - pb + pb**3/3) + d*(1 + C*pb); s10 = 2*(pb*pb - 1) + d*C; s20 = 4*pb*s20_scale
    I = np.eye(M); Z = np.zeros((M, M)); dg = np.diag
    hc, phz = dg(B['hc']), dg(B['phz'])
    Kt = np.block([[D2 + 6*hc@D1, (4/3)*dg(r2*B['U']), (2/3)*dg(r2*B['U1'])],
                   [-6*hc@D1, D2 - (2/3)*dg(r2*B['U']), phz@D1 - (1/3)*dg(r2*B['U1'])],
                   [3*phz@D1, -2*dg(r2*B['U1']), D2 + 3*hc@D1 - dg(r2*B['U2'])]])
    Ct = np.block([[-6*I, Z, Z], [6*I, Z, Z], [Z, Z, -3*I]])
    n3 = 3*M
    A = np.block([[np.zeros((n3, n3)), np.eye(n3)], [Kt, Ct]]); Bm = np.eye(2*n3)
    def setrow(field, node, coeffs):   # coeffs: dict field -> row vector over M nodes (acting on q)
        row = n3 + field*M + node
        A[row, :] = 0; Bm[row, :] = 0
        for fld, vec in coeffs.items(): A[row, fld*M:(fld + 1)*M] += vec
    e = lambda i: np.eye(M)[i]
    for fld in range(3):
        setrow(fld, 0, {fld: e(0)})                                    # Dirichlet at z=-L
        for k in range(K - 1):
            iL = k*(N + 1) + N; iR = (k + 1)*(N + 1)
            setrow(fld, iL, {fld: e(iL) - e(iR)})                     # continuity
            setrow(fld, iR, {fld: D1[iL] - D1[iR]})                   # C1 matching
    last = M - 1
    setrow(0, last, {0: D1[last], 1: -rb*s0/6*e(last), 2: -rb*s10/6*e(last)})
    setrow(1, last, {1: D1[last] - rb*s0/6*e(last), 2: -rb*s10/6*e(last)})
    setrow(2, last, {2: D1[last] + rb*s20/2*e(last), 1: rb*s10/2*e(last)})
    return A, Bm, z, B, dict(s0=s0, s10=s10, s20=s20, M=M)


def analyse(A, Bm, z, B, M, sel_min=-1.5 + 1e-6):
    w, V = sl.eig(A, Bm)
    fin = np.isfinite(w) & (np.abs(w) < 1e5)
    rows = []
    for k in np.nonzero(fin)[0]:
        lam = w[k]
        if lam.imag < -1e-9: continue
        q = V[:3*M, k].reshape(3, M); sc = np.max(np.abs(q))
        fb = abs(q[2, -1])/sc
        with np.errstate(all='ignore'):
            ch, sh = np.cosh(lam*z), np.sinh(lam*z)
            g = np.array([ch + B['hc']*sh, lam*ch + B['hc']*sh, B['phz']*sh])
            al = np.vdot(g.ravel(), q.ravel())/np.vdot(g.ravel(), g.ravel())
            gres = float(np.linalg.norm(q.ravel() - al*g.ravel())/np.linalg.norm(q.ravel()))
        if not np.isfinite(gres): gres = 1.0
        far = float(np.max(np.abs(q[:, z < -0.9*abs(z[0])]))/sc)
        rows.append(dict(real=float(lam.real), imag=float(lam.imag), shell_phi_rel=float(fb), gauge_res=gres, far_weight=far))
    rows.sort(key=lambda r: -r['real'])
    above = [r for r in rows if r['real'] > sel_min]
    line = [r for r in rows if abs(r['imag']) > 1e-6 and abs(r['imag']) < 2000]
    return rows, above, line


def decoupled_scalar(bgf, pb, rb, d, breaks, N, s20_scale=1.0):
    """Scalar-only Sturm-Liouville problem  -(f''+3hc f') + rho^2 U'' f = mu2 f,  f'(0) = -rho_b s20 f(0)/2,
    f(-L)=0; and the Liouville potential Veff = 9/4 + rho^2 (U'' - 5U/8 - 3 phi_y^2/16) (exact identity
    using the static equations), plus boundary coefficient beta = rho_b (s20/2 - s0/4)."""
    x, Dx = cheb(N); K = len(breaks) - 1; M = K*(N + 1)
    z = np.concatenate([breaks[k] + (x + 1)*(breaks[k + 1] - breaks[k])/2 for k in range(K)])
    D1 = sl.block_diag(*[Dx*2/(breaks[k + 1] - breaks[k]) for k in range(K)]); D2 = D1@D1
    B = bgf(z); r2 = B['rho']**2
    s0 = 2*(1 - pb + pb**3/3) + d*(1 + C*pb); s20 = 4*pb*s20_scale
    L = -(D2 + np.diag(3*B['hc'])@D1) + np.diag(r2*B['U2']); Bm = np.eye(M)
    def setrow(i, vec): L[i] = vec; Bm[i] = 0
    setrow(0, np.eye(M)[0])
    for k in range(K - 1):
        iL = k*(N + 1) + N; iR = (k + 1)*(N + 1)
        setrow(iL, np.eye(M)[iL] - np.eye(M)[iR]); setrow(iR, D1[iL] - D1[iR])
    setrow(M - 1, D1[M - 1] + rb*s20/2*np.eye(M)[M - 1])
    mu2 = sl.eigvals(L, Bm); mu2 = np.sort(mu2[np.isfinite(mu2)].real)
    py = B['phz']/B['rho']
    Veff = 9/4 + r2*(B['U2'] - 5*B['U']/8 - 3*py*py/16)
    # direct Liouville potential from hc for comparison (identity check)
    hcz = D1@B['hc']
    Veff_direct = r2*B['U2'] + 1.5*hcz + 2.25*B['hc']**2
    beta = rb*(s20/2 - s0/4)
    return dict(lowest_mu2=mu2[:4].tolist(), min_Veff_minus_9_4=float(np.min(Veff) - 9/4),
                Veff_identity_max_rel_err=float(np.max(np.abs(Veff - Veff_direct)/np.abs(Veff))), beta=float(beta))


NS = [int(x) for x in (sys.argv[1] if len(sys.argv) > 1 else '24,32').split(',')]


def main():
    t0 = time.time(); out = {}
    brk_plus = [-6, -4, -2.5, -1.5, -0.8, -0.4, -0.2, -0.1, -0.05, -0.025, -0.012, -0.005, 0.]
    brk_orig = brk_plus
    GN = 1.6571936312453648
    for kind in ['original', 'plus']:
        bgf, pb, rb, meta = background(kind)
        res = {}
        for N in NS:
            A, Bm, z, B, info = build(bgf, pb, rb, 1e-3, brk_plus if kind == 'plus' else brk_orig, N)
            rows, above, line = analyse(A, Bm, z, B, info['M'])
            ent = dict(M=info['M'], top_above_minus1p5=above[:12],
                       n_line_modes=len(line), max_dev_from_line=float(max(abs(r['real'] + 1.5) for r in line)) if line else None,
                       real_axis_above_minus3=[r for r in rows if abs(r['imag']) < 1e-6 and r['real'] > -3][:12])
            if kind == 'original':
                ent['top_real'] = rows[0]['real']; ent['diff_vs_GN'] = rows[0]['real'] - GN
            else:
                ent['shell_supported_above_line'] = [r for r in above if r['shell_phi_rel'] > 1e-6 and abs(r['imag']) < 1000]
            res['N%d' % N] = ent
            print(kind, N, [(round(r['real'], 7), round(r['imag'], 4), '%.1e' % r['shell_phi_rel'], '%.1e' % r['gauge_res']) for r in above[:8]], flush=True)
        if kind == 'plus':
            flip = {}
            for N in NS:
                A, Bm, z, B, info = build(bgf, pb, rb, 1e-3, brk_plus, N, s20_scale=-1)
                rows, above, _ = analyse(A, Bm, z, B, info['M'])
                ra = sorted([r['real'] for r in rows if abs(r['imag']) < 1e-6])
                top = max((r for r in rows if abs(r['imag']) < 1e-6), key=lambda r: r['real'])
                flip['N%d' % N] = dict(top_real_axis=top, min_real_axis=ra[0], pair_sum=top['real'] + ra[0])
                print('flip', N, top['real'], ra[0], flush=True)
            res['control_s20_flipped'] = flip
            dec = {}
            for N in [32, 40]:
                dec['N%d' % N] = decoupled_scalar(bgf, pb, rb, 1e-3, brk_plus, N)
                dec['flip_N%d' % N] = decoupled_scalar(bgf, pb, rb, 1e-3, brk_plus, N, -1)
            for key in ['N40', 'flip_N40', 'N32', 'flip_N32']:
                m = dec[key]['lowest_mu2'][0]
                dec[key]['lambda_from_lowest_mu2'] = (-1.5 + math.sqrt(2.25 - m)) if m < 2.25 else None
            res['decoupled_scalar'] = dec
            print('decoupled', json.dumps(dec)[:600], flush=True)
        else:
            dec = {}
            for N in [32, 40]:
                dec['N%d' % N] = decoupled_scalar(bgf, pb, rb, 1e-3, brk_orig, N)
            res['decoupled_scalar'] = dec
        res['background'] = dict(phi_b=float(pb), rho_b=float(rb))
        out[kind] = res
        out['N_values'] = NS
        (HERE/'v3_chebyshev.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    out['runtime_seconds'] = time.time() - t0
    (HERE/'v3_chebyshev.json').write_text(json.dumps(out, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
