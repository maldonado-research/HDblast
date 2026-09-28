"""Core utilities for the time-domain / discrete-operator stability study (Method B).

The PDE system, gauge, stretched grid, quintic-Hermite shell closure and the nonlinear
right-hand side are those of the frozen registered solver (copied unchanged into
frozen_input/registered_solver.py; its SHA-256 is checked on import).  This module only
(i) supplies a second static background, the regular-cone "+1 branch", on the same grid,
(ii) builds the linear operator on either background and verifies it against a complex-step
Jacobian of the unchanged nonlinear RHS, and (iii) classifies eigenvectors into
constraint-violating, pure-gauge and physical modes.

Conventions (inherited): metric e^{2B}(-dt^2+dz^2)+e^{2A}dx_3^2, shell at z=0, bulk z<0,
background A=t+ln rho(z), B=ln rho(z); perturbations a=A-t-ln rho, b=B-ln rho, f=phi-phi_bg.
Coordinate time t is measured in units of the static shell's own Hubble time 1/H (H=1/rho_b).
Floating-point study; nothing here is an interval certificate.
"""
from pathlib import Path
import hashlib, importlib.util, json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy import sparse

HERE = Path(__file__).resolve().parent
FROZEN = HERE / 'frozen_input' / 'registered_solver.py'
FROZEN_SHA = '4147239c40b63c2de89f5091e732e19f09c4644a0c08ae5735202a548cdcd6b6'
if hashlib.sha256(FROZEN.read_bytes()).hexdigest() != FROZEN_SHA:
    raise RuntimeError('frozen registered_solver.py hash mismatch')
_spec = importlib.util.spec_from_file_location('registered_solver', FROZEN)
RS = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(RS)

C_REG = RS.C  # 2/1.0357712571566784-4/3 = 0.5975949350280132


# ---------------------------------------------------------------- +1 branch background
def pot_eta(e):
    """U and derivatives written in eta=phi-1 (exact polynomial expansion of U about phi=1)."""
    w = 1/3 + e*e + e**3/3; wp = e*(2 + e); p = 1 + e
    u = .5*wp*wp - (2/3)*w*w
    up = wp*(2*p - (4/3)*w)
    upp = 4*p*p + 2*wp - (4/3)*(wp*wp + 2*p*w)
    return u, up, upp


def plus_background(tdet=1e-3, c=C_REG, rtol=3e-13, y0=5e-5, guess=None):
    """Regular-cone +1 static branch, re-solved here with a second-order (rho, rho_y, eta, eta_y, z)
    system (no first-integral square root).  Unknowns: log10(-eta_h), y_b.  Both junctions
    rho_y/rho = sigma/6 and eta_y = -sigma'/2 are solved by root finding."""
    sig = lambda e: 2*(1/3 + e*e + e**3/3) + tdet*(1 + c*(1 + e))
    sig1 = lambda e: 2*e*(2 + e) + tdet*c

    def integrate(x, dense=False):
        eh = -10**x[0]; yb = x[1]; u, up, upp = pot_eta(eh)
        A_ = -u/36; B_ = u*u/4320 - up*up/750; Cc = up/10; D_ = up*(upp/280 + u/630)
        ini = [y0 + A_*y0**3 + B_*y0**5, 1 + 3*A_*y0**2 + 5*B_*y0**4,
               eh + Cc*y0*y0 + D_*y0**4, 2*Cc*y0 + 4*D_*y0**3, 0.0]
        def rhs(y, v):
            r, ry, e, p, _ = v; u, up, _ = pot_eta(e)
            return [ry, -r*(p*p/4 + u/6), p, up - 4*ry*p/r, 1/r]
        sol = solve_ivp(rhs, (y0, yb), ini, method='DOP853', rtol=rtol,
                        atol=[1e-14, 1e-14, 1e-44, 1e-44, 1e-14], max_step=.05, dense_output=dense)
        if not sol.success: raise RuntimeError(sol.message)
        return sol

    def residual(x):
        r, ry, e, p, _ = integrate(x).y[:, -1]
        return np.array([ry/r - sig(e)/6, p + sig1(e)/2])/tdet

    if guess is None:
        arch = json.loads((HERE/'frozen_input'/'PLUS_BRANCH_RESULTS.json').read_text())
        row = min(arch['rows'], key=lambda q: abs(q['delta'] - tdet))['refined']
        guess = [math.log10(-row['eta_h']), row['y_b']]
    fit = root(residual, guess, tol=1e-12, options={'eps': 1e-9})
    res = residual(fit.x)*tdet
    if np.max(np.abs(res)) > 1e-11: raise RuntimeError('plus-branch junction root failed %r' % res)
    sol = integrate(fit.x, True); r, ry, e, p, z_end = sol.y[:, -1]
    u, _, _ = pot_eta(e)
    ham = ry*ry - 1 - r*r*(p*p/12 - u/6)

    def sample(z):
        """Background on conformal nodes z<=0 (shell at z=0).  Newton inversion of z(y)."""
        targets = np.asarray(z, float) + z_end
        yy = np.linspace(y0, fit.x[1], 200001); zz = sol.sol(yy)[4]
        if targets.min() < zz[0]: raise ValueError('grid extends beyond the cone table')
        ys = np.interp(targets, zz, yy)
        for _ in range(6):
            v = sol.sol(ys); ys = np.clip(ys - (v[4] - targets)*v[0], y0, fit.x[1])
        v = sol.sol(ys); rho, rhoy, eta, py = v[:4]
        err = float(np.max(np.abs(v[4] - targets)))
        uu = pot_eta(eta)[0]
        hc_first = np.sqrt(1 + rho*rho*(py*py/12 - uu/6))
        return dict(rho=rho, phi=1 + eta, eta=eta, phiz=rho*py, Hc=rhoy, Hc_first_integral=hc_first,
                    lnrho=np.log(rho), z_inversion_error=err)

    meta = dict(branch='plus', tdet=tdet, c=c, eta_h=-10**fit.x[0], y_b=float(fit.x[1]), rho_b=float(r),
                phi_b=float(1 + e), eta_b=float(e), H2=float(1/r**2), junction_residual=res.tolist(),
                hamiltonian_first_integral_at_shell=float(ham), z_end=float(z_end), rtol=rtol, y0=y0)
    return sample, meta


def original_background(tdet=1e-3, dc=0.):
    """Original unstable registered shell, from the frozen solver (unchanged)."""
    sample, meta = RS.shell_background(tdet, dc)
    meta = dict(meta); meta['branch'] = 'original'
    return sample, meta


# ---------------------------------------------------------------- solver on any background
UETA = np.array([-2/27, 0., 14/9, 50/27, -1/6, -4/9, -2/27])  # U(1+eta) coefficients


def eta_derivatives(eta):
    """[U, U', ..., U^(6)] at phi=1+eta, evaluated as polynomials in eta."""
    out = []; c = UETA.copy()
    for j in range(7):
        out.append(np.polynomial.polynomial.polyval(eta, c)); c = np.arange(1, len(c))*c[1:]
    return out


class GeneralSolver(RS.Solver):
    """Registered Solver with an externally supplied static background.  The RHS, junction
    data, linear matrix and diagnostics are inherited unchanged."""
    def __init__(self, bg, tdet=1e-3, hmin=2e-4, L=6, stretch=.05, ko=0., c=C_REG, precise_eta=True):
        sample, meta = bg
        self.z, self.zp, self.zpp = RS.grid(hmin, L, stretch); self.n = len(self.z)
        self.D1, self.D2, self.g1, self.g2, self.ghost = RS.matrices(self.zp, self.zpp)
        self.meta = meta; self.bg = sample(self.z)
        self.rho = self.bg['rho']; self.hc = self.bg['Hc']; self.phz = self.bg['phiz']; self.phi = self.bg['phi']
        self.rb = self.rho[-1]; self.pb = self.phi[-1]; self.rho2 = self.rho**2
        self.tdet = tdet; self.ko = ko; self.hmin = hmin; self.L = L; self.stretch = stretch
        self.pot = RS.derivatives(self.phi)
        self.pot_mode = 'phi-polynomial (registered)'
        if meta.get('branch') == 'plus' and precise_eta:
            # U and its derivatives about phi=1 evaluated as polynomials in eta=phi-1 (exact
            # coefficients), avoiding cancellation when |eta| is below ~1e-12 in the deep bulk.
            self.pot = eta_derivatives(self.bg['eta']); self.pot_mode = 'eta-polynomial (precise)'
        if c != RS.C:
            raise ValueError('the inherited bc() uses the module constant C; only c=C_REG supported')
        self.s0 = 2*RS.W(self.pb) + tdet*(1 + RS.C*self.pb)
        self.s10 = 2*(self.pb*self.pb - 1) + tdet*RS.C; self.s20 = 4*self.pb
        self.state0 = np.zeros((6, self.n)); self.seed_meta = None
        # Static junction residuals on the sampled grid (background consistency).
        self.junction_on_grid = [float(self.hc[-1] - self.rb*self.s0/6), float(self.phz[-1] + self.rb*self.s10/2)]


def free_index(n):
    """DOFs that are not frozen (first two nodes of each of the six blocks are held at zero)."""
    keep = np.ones(6*n, bool)
    for blk in range(6): keep[blk*n:blk*n + 2] = False
    return np.nonzero(keep)[0]


def verify_matrix(s, seed=4801):
    M = s.linear_matrix()
    rng = np.random.default_rng(seed); v = rng.normal(size=(6, s.n)); v[:, :2] = 0
    num = s.rhs(1j*1e-30*v).imag/1e-30
    return M, float(np.max(np.abs(M@v.ravel() - num.ravel()))/(1 + np.max(np.abs(num))))


# ---------------------------------------------------------------- linear diagnostics
def linear_constraints(s, v):
    """Linearised Hamiltonian and momentum constraints (same expressions as the frozen
    diagnostics(), linearised by hand) and the sum of absolute values of their terms, used
    as a scale.  v: (6,n) real or complex."""
    a, b, f, pa, pb, pf = v
    ga = s.rb*(s.s0*b[-1] + s.s10*f[-1])/6
    gf = -s.rb*(s.s10*b[-1] + s.s20*f[-1])/2
    gpa = s.rb*(s.s0*pb[-1] + s.s10*pf[-1])/6
    az = s.D1@a + s.g1*ga; bz = s.D1@b + s.g1*ga; fz = s.D1@f + s.g1*gf
    azz = s.D2@a + s.g2*ga; paz = s.D1@pa + s.g1*gpa
    one = np.array([25/12, -4, 3, -4/3, 1/4])/s.zp[-1]
    paz = paz.copy(); paz[-1] = np.dot(one, pa[-1:-6:-1])
    su = s.rho2*(2*b*s.pot[0] + s.pot[1]*f)
    Mterms = [-3*paz, -3*(az - bz), 3*s.hc*pb, -pf*s.phz]
    Hterms = [-2*su, 12*pa, 6*pb, -18*s.hc*az, 6*s.hc*bz, -6*azz, -2*s.phz*fz]
    M = sum(Mterms); H = sum(Hterms)
    Ms = sum(np.abs(t) for t in Mterms); Hs = sum(np.abs(t) for t in Hterms)
    return H, M, Hs, Ms


def nonlinear_constraints(s, state):
    """Full nonlinear Hamiltonian/momentum constraint residuals (the expressions of the frozen
    diagnostics(), unchanged) together with the sum of the absolute values of their perturbative
    terms, as a scale.  The background-only pieces cancel analytically in these expressions."""
    a, b, f, pa, pb, pf = state
    ga, gf, gpa, gpf = s.bc(b[-1], f[-1], pb[-1], pf[-1])
    az = s.D1@a + s.g1*ga; bz = s.D1@b + s.g1*ga; fz = s.D1@f + s.g1*gf
    azz = s.D2@a + s.g2*ga; paz = s.D1@pa + s.g1*gpa
    one = np.array([25/12, -4, 3, -4/3, 1/4])/s.zp[-1]
    paz = paz.copy(); paz[-1] = np.dot(one, pa[-1:-6:-1])
    du = sum(s.pot[j]*f**j/math.factorial(j) for j in range(1, 7))
    su = s.rho2*(np.expm1(2*b)*(s.pot[0] + du) + du)
    Mterms = [-3*paz, -3*(1 + pa)*(az - bz), 3*(s.hc + az)*pb, -pf*(s.phz + fz)]
    Hterms = [-2*su, 12*pa + 6*pa*pa, 6*(1 + pa)*pb, -18*s.hc*az, 6*s.hc*bz, -12*az*az + 6*az*bz, -6*azz,
              -pf*pf, -2*s.phz*fz - fz*fz]
    return sum(Hterms), sum(Mterms), sum(np.abs(t) for t in Hterms), sum(np.abs(t) for t in Mterms)


def gauge_template(s, lam):
    """Residual conformal gauge mode preserving z=0: xi^t=e^{lam t}cosh(lam z), xi^z=e^{lam t}sinh(lam z).
    Returns the (6,n) perturbation it generates on the static background."""
    z = s.z; ch = np.cosh(lam*z); sh = np.sinh(lam*z)
    a = ch + s.hc*sh; b = lam*ch + s.hc*sh; f = s.phz*sh
    return np.array([a, b, f, lam*a, lam*b, lam*f])


def classify(s, lam, vec, mask_nodes=6, zcut=0.8):
    """Diagnostics of one eigenpair.  vec: flat (6n) complex eigenvector."""
    v = vec.reshape(6, s.n)
    mask = (s.z > -zcut*s.L) & (np.arange(s.n) > mask_nodes)
    H, M, Hs, Ms = linear_constraints(s, v)
    rH = float(np.max(np.abs(H[mask]))/max(np.max(Hs[mask]), 1e-300))
    rM = float(np.max(np.abs(M[mask]))/max(np.max(Ms[mask]), 1e-300))
    # Weighted-l2 relative constraint residual
    rH2 = float(np.linalg.norm(H[mask])/max(np.linalg.norm(Hs[mask]), 1e-300))
    rM2 = float(np.linalg.norm(M[mask])/max(np.linalg.norm(Ms[mask]), 1e-300))
    scale = np.max(np.abs(v[:3]))
    fb = v[2, -1]/scale; hb = (v[3, -1] - v[1, -1])/scale
    # Least-squares projection on the gauge template with the same eigenvalue.
    with np.errstate(all='ignore'):
        g = gauge_template(s, lam).astype(complex)
        gm = g[:, mask].ravel(); vm = v[:, mask].ravel()
        alpha = np.vdot(gm, vm)/np.vdot(gm, gm)
        gauge_res = float(np.linalg.norm(vm - alpha*gm)/np.linalg.norm(vm))
    if not np.isfinite(gauge_res): gauge_res = 1.0  # template overflow (|lam| L >> 700): not a gauge match
    return dict(real=float(lam.real), imag=float(lam.imag), constraint_H_max_rel=rH, constraint_M_max_rel=rM,
                constraint_H_l2_rel=rH2, constraint_M_l2_rel=rM2,
                shell_phi_rel=float(abs(fb)), shell_hubble_rel=float(abs(hb)),
                gauge_template_residual=gauge_res,
                far_weight=float(np.max(np.abs(v[:3, :s.n//10]))/scale))


def label(d, ctol=1e-3, gtol=1e-3, stol=1e-6):
    if max(d['constraint_H_l2_rel'], d['constraint_M_l2_rel']) > ctol: return 'constraint-violating'
    if d['gauge_template_residual'] < gtol or (d['shell_phi_rel'] < stol and d['shell_hubble_rel'] < stol):
        return 'gauge'
    return 'physical-candidate'


def dump(path, obj):
    def conv(x):
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, complex): return [x.real, x.imag]
        raise TypeError(type(x))
    Path(path).write_text(json.dumps(obj, indent=1, default=conv) + '\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
