# Scalar-sector spectrum of the registered dS shell - independent Gaussian-normal (ADM radial split) derivation

16 September 2026. Blind derivation: `RESULT.json` was written before any orchestrator / effective-theory / literature
file of Chat 9 was opened (comparison in sec. 9 was added afterwards). Floating point (numpy RK4), not an interval certificate.

**Result.** For the registered shell (t = 1e-3, c = c_star) there is exactly one bound state below the continuum
(mu^2 < 9/4, scanned -400 < mu^2 < 2.2499):

    mu^2 = m^2/H^2 = -7.71787163 (+- 1e-8),   m^2 = -1.242037e-3,   growth e^{1.657 H tau} on the shell.

It has POSITIVE norm (tachyon with healthy kinetic term, not a ghost); the norm is the manifestly positive bulk integral
N = 2 int rho^2 X^2 dy + 12 int rho^2 Psi^2 dy (sec. 7). The static dS shell is linearly unstable.

## 1. Set-up and radial split

Gaussian-normal coordinates adapted to the shell: `ds^2 = dy^2 + g_mn(y,x) dx^m dx^n`, `K_mn = (1/2) d_y g_mn`.
With `Gamma^y_mn = -K_mn`, `Gamma^m_yn = K^m_n` one finds (derived by hand, checked on the background):

    R_yy      = -d_y K - K^m_n K^n_m                          = phi_y^2 + (2/3) U
    R_ym      = D_n K^n_m - D_m K                             = phi_y d_m phi
    R^m_n(5)  = R^m_n[g] - d_y K^m_n - K K^m_n                = g^{ml} d_l phi d_n phi + (2/3) U delta^m_n
    2 G_yy    = K^2 - K^m_n K^n_m - R[g]                      = phi_y^2 - (d phi)_4^2 - 2U          (Hamiltonian constraint)
    phi_yy + K phi_y + Box_g phi = U_phi .

Background `g = rho^2 gamma`, `K^m_n = H delta^m_n`, `H = rho'/rho`: these reproduce the three registered background equations
(e.g. `R^m_n`: `3 - 3 rho'^2 - rho rho'' = (2/3) U rho^2`, true by `rho''/rho = -phi'^2/4 - U/6`).

Perturbation (scalar harmonic `Box_gamma Y = mu^2 Y`):

    g_mn = rho^2 [ (1 + 2 psi Y) gamma_mn + 2 E nabla_m nabla_n Y ],   phi = phi_0 + chi Y,   delta g_yy = delta g_ym = 0 .

To first order `delta K^m_n = psi' Y delta^m_n + E' nabla^m nabla_n Y`, `delta K = (4 psi' + mu^2 E') Y`.
4D curvature: the E-part is a diffeo of the Einstein space gamma, so `delta R^m_n = 0`; the psi-part is a conformal rescaling,
`delta R^m_n[g] = rho^-2 [ -6 psi delta^m_n Y - 2 psi nabla^m nabla_n Y - psi delta^m_n Box Y ]`, `delta R[g] = -(24 + 6 mu^2) psi Y / rho^2`.
The connection variation drops out of the momentum constraint because the background `K^m_n` is pure trace, and
`nabla_n nabla^n nabla_m E' - nabla_m Box E' = R_m^l nabla_l E' = 3 nabla_m E'`.

## 2. Linearized equations (Q := E', P := psi')

    (MC)   E' = psi' + phi' chi / 3                                                         [coefficient of nabla_m Y]
    (HC)   3H (4 psi' + mu^2 E') + 3 (4 + mu^2) psi / rho^2 = phi' chi' - U_phi chi
    (EV-E) E'' + 4H E' = -2 psi / rho^2                                                     [coefficient of nabla^m nabla_n Y]
    (EV-p) psi'' + 8H psi' + mu^2 H E' + (6 + mu^2) psi / rho^2 + (2/3) U_phi chi = 0       [coefficient of delta^m_n]
    (SC)   chi'' + 4H chi' + (4 psi' + mu^2 E') phi' + (mu^2/rho^2 - U_phiphi) chi = 0
    (YY)   -(4 psi' + mu^2 E')' - 2H (4 psi' + mu^2 E') = 2 phi' chi' + (2/3) U_phi chi .

Closed system actually integrated: state `(chi, chi', Q, psi)` with `psi' = Q - phi' chi/3` (MC), (EV-E) and (SC).
(HC), (EV-p), (YY) are NOT used in the flow; they are monitored and stay satisfied to 1e-9 ... 1e-13 (normalised) from cone to
shell for all 3595 scanned mu^2 (`run_spectrum.py`) - a numerical Bianchi-identity check of the whole set.
For `mu^2 = -4` (`nabla_m nabla_n Y = -gamma_mn Y`) psi and E are degenerate (only psi - E is defined); (HC) then forces
`chi = C phi'`, i.e. the regular bulk solution is pure gauge (eps = C), so mu^2 = -4 is not a physical eigenvalue unless
`phi'' + sigma'' phi'/2 = 0` at the shell. The integrated flow has no singularity at mu^2 = -4 (no division by 4 + mu^2).

## 3. Residual gauge of Gaussian-normal gauge

`xi^y = eps(x)`, `xi^m = -F(y) nabla^m eps + nabla^m L(x)`, `F' = 1/rho^2` preserve `delta g_yA = 0`:

    delta psi = H eps,   delta E = -F eps + L,   delta chi = phi' eps      (so  delta Q = -eps/rho^2,  delta chi' = phi'' eps).

L is a pure 4D diffeo (E itself never appears, only Q). eps moves the y = const surfaces: in SHELL-adapted GN gauge eps is fixed,
but then the fields are singular at the cone (`psi ~ eps/y`). I therefore integrate in CONE-adapted GN gauge (all fields regular,
shell displaced) and transform at the shell with an unknown constant eps: `f_s = f + delta_eps f`.
Gauge invariants: `X = chi + phi' rho^2 Q`, `Psi = psi + H rho^2 Q`.

## 4. Cone behaviour and continuum

Near y = 0: `rho ~ y`, `phi' ~ a y` (a = U_phi(phi_h)/5), (SC) -> `s(s+3) + mu^2 = 0`, `s_pm = -3/2 +- sqrt(9/4 - mu^2)`.
Weight of mu^2 in (SC) is rho^2, so normalisability needs `2s + 3 > 0`: only `s_+`. For mu^2 > 9/4 the exponents are complex
(continuum). Frobenius data used (`k^2 = -U_h/6`):

    chi = y^s (1 + c2 y^2),  c2 = (U'' + mu^2 k^2/3 - 4 s k^2/3)/((s+2)(s+5) + mu^2),
    Q = 2a y^{s+1} / (3 (s+3)(s+4)),   psi = -(a/3)(s+5) y^{s+2} / ((s+3)(s+4))

(the apparent pole of psi at mu^2 = -4 cancels because 4 + mu^2 = -(s+4)(s-1)). The third local solution `Q ~ 1/y^2` is the
eps gauge mode. Grid: geometric in y from y0 to 0.5 (RK4-stable against the 4H ~ 4/y stiffness), then uniform.
WARNING found on the way: a uniform grid with 4 H dy > 2.8 at the cone produced a SPURIOUS second root near mu^2 = 2.177;
it disappears on the stable grid and is not physical.

## 5. Junction conditions and the gauge-invariant mismatch

Israel (Z2, bulk on y < y_b): `K^m_n = sigma_t(phi)/6 delta^m_n`, scalar: `phi' = -sigma_t'/2`. Perturbed, in shell gauge:

    (J1) E_s' = 0                         (trace-free part)
    (J3) psi_s' = sigma_t' chi_s / 6      (trace part)   - identical to (J1) by (MC) and the background phi' = -sigma_t'/2
    (J2) chi_s' = -(sigma_t''/2) chi_s .

(J1) fixes `eps = rho_b^2 Q(y_b)`; inserting into (J2):

    M(mu^2) = chi' + (sigma_t''/2) chi + (phi'' + sigma_t'' phi'/2) rho^2 E'   at y_b
            = X' + (sigma_t''/2) X + 2 phi' Psi .

Under the residual gauge `delta M = phi'' e + (sigma''/2) phi' e - (phi'' + sigma'' phi'/2) e = 0` identically.
Numerical confirmation (`run_gauge_test.py`): contaminating the cone data with a gauge mode 1000x larger than chi changes M by
< 4e-12 relative; pure gauge data integrate to the analytic gauge mode at the shell to 1.5e-8 for mu^2 in {-7.7, -4.3, -1, 0.7, 2}
(this tests MC + EV-E + SC for mu^2 != 0, complementary to T1) and give M = O(1e-11).

## 6. Spectrum

`run_spectrum.py`, `run_roots.py`: M has a single sign change on -400 < mu^2 < 2.2499.

| grid n (uniform part), y0 | mu^2 |
|---|---|
| 1000, 2e-3 | -7.7178720572 |
| 2000, 2e-3 | -7.7178716515 |
| 4000, 2e-3 | -7.7178716252 |
| 4000, 5e-4 | -7.7178716252 |
| 4000, 8e-3 | -7.7178716255 |

The eigenfunctions X and Psi have no nodes (ground state). `phi_b'' + sigma'' phi_b'/2 =: B phi_b'` with B = -2.68e-4 < 0.

## 7. Sign of the norm (three routes, all agree)

**Route A/B - on-shell quadratic action.** For a quadratic action `S2 = (1/2) <f, D f>`, if f solves every linearised equation
except one junction condition, S2 equals one half of (that junction's equation) x (its conjugate field) on the shell.
First variations (two Z2 copies, outward normal +y, action (1/2)R - (1/2)(d phi)^2 - U - sigma delta(shell)):
`delta S = int_shell sqrt(-g) { -2 (phi' + sigma'/2) delta phi - (1/2)[2(K^mn - K g^mn) + sigma g^mn] delta g_mn }`.
* A (eps from J1, scalar junction left open):  `S2 = (1/2) int sqrt(-gamma) q^2 * [-2 rho_b^4 chi_s M]`.
* B (eps from J2, metric junction left open):  `S2 = (1/2) int sqrt(-gamma) q^2 * [6 (mu^2 + 4) rho_b^4 psi_s E_s']`
  (the E_s-dependent pieces cancel after integration by parts on dS_4 using `int (nabla nabla Y)^2 = mu^2(mu^2+3) int Y^2`: 4D gauge invariance).
Near the root `S2 = (1/2) N int q (Box - mu_0^2) q`, N = d/dmu^2 of the bracket. A healthy field has N > 0
(test-field check: for a probe scalar with Robin condition, `-2 rho_b^4 chi dM/dmu^2 = 2 int rho^2 chi^2 > 0`).
Both reductions restrict the same 2x2 quadratic form in (amplitude, eps) to different lines through its null vector, so their
derivatives at the root must coincide. Numerically (cone normalisation chi -> y^s):

    N_A = 5.33072493e15,   N_B = 5.33072493e15   (agree to 1e-9 at n = 4000)  > 0.

Canonical normalisation `Z = N / (rho_b^2 chi_s^2) = 0.896383` (kinetic coefficient of delta phi_b in physical units).

**Route C - positive bulk integral (proved).** From (MC), (HC), (EV-E) the gauge invariants obey the first-order system

    Psi' = -2H Psi - (phi'/3) X,     X' = (phi''/phi') X + (3 lam/(rho^2 phi') - 2 phi') Psi,    lam = mu^2 + 4,

and `M = B X + 3 lam Psi/(rho_b^2 phi_b')`, `B = phi''/phi' + sigma''/2`. Two identities follow (cone terms vanish for s_+):
`[ (rho^2/phi') (Psi dX/dlam - X dPsi/dlam) ]' = 3 Psi^2/phi'^2` and `[rho^2 Psi X/phi']' = -rho^2 X^2/3 - 2 rho^2 Psi^2 + 3 lam Psi^2/phi'^2`. They give, at a root,

    N_A = 12 lam N_FK = 12 Q_FK,
    Q_FK = int (rho^2 X^2/6 + rho^2 Psi^2) dy > 0,     N_FK = int 3 Psi^2/(2 phi'^2) dy + 3 Psi_b^2/(2 phi_b'^2 B).

Hence **N = 2 int rho^2 X^2 dy + 12 int rho^2 Psi^2 dy > 0 for every bound state: there are no scalar ghosts**, and an eigenvalue
with mu^2 < -4 requires the Frolov-Kofman weight N_FK < 0, i.e. B < 0 (true here). The negative FK "norm" is not a ghost signal:
the physical Klein-Gordon norm is lam x N_FK. Numerics (`run_norm_bulk.py`, `run_fk_identity.py`, trapezoid, n = 2000 -> 4000):
`N_bulk/N_A - 1 = 4.9e-6 -> 1.3e-6`, `Q_FK/(lam N_FK) - 1 = 4.9e-6 -> 1.2e-6` (second-order convergence).
Status of the proof: rigorous given the first-order system and the on-shell-action definition of N; the statement
"S2 = (1/2) f.EOM with these boundary terms" is standard but the boundary term at the cone was only argued (flux ~ y^{2s+3} -> 0).

## 8. Validation

* **T1** (`run_T1.py`): finite difference in phi_h of two exact backgrounds = exact mu^2 = 0, Y = const, E = 0 GN perturbation.
  For constant Y, (MC) (a 4-gradient) and (EV-E) (coefficient of nabla nabla Y) are vacuous; (HC), (SC), (EV-p), (YY) must hold.
  At alpha = 0.3: residuals 1e-8 ... 1e-11 (finite-difference limited). At the registered alpha = 8.8e-7 the same holds to
  1e-9 for y > 6; for y < 3 the test is void because |delta ln rho| < 1e-10 is below double-precision resolution (SC still passes to 1e-7).
* **Bianchi/constraint propagation** and **gauge-mode** tests: secs. 2, 5.
* **T2** (`run_tscan.py`, c = c_star, n = 1500):

| t | mu^2 | m^2 = mu^2 H^2 | Z |
|---|---|---|---|
| 1e-2 | -7.7005615 | -1.24059e-2 | 0.896300 |
| 3e-3 | -7.7140236 | -3.72515e-3 | 0.896365 |
| 1e-3 | -7.7178717 | -1.24204e-3 | 0.896383 |
| 3e-4 | -7.7192191 | -3.72645e-4 | 0.896390 |
| 1e-4 | -7.7196057 | -1.24218e-4 | 0.896391 |
| 3e-5 | -7.7197495 | -3.72659e-5 | 0.896392 |
| 1e-5 | -7.71982 (precision-limited) | -1.24220e-5 | 0.896385 |

  mu^2(t) = -7.719796 + 1.924 t: mu^2 tends to a constant, so `m^2 = -1.2422 t -> 0` linearly with finite positive Z:
  the mode is continuously connected to the massless brane-position modulus of the flat BPS wall.

## 9. Comparison with the orchestrator (read only after RESULT.json was on disk)

* Orchestrator longitudinal-gauge master equation: mu^2 = -7.7178716187 (n = 12000); mine -7.7178716252. Difference 6.5e-9. AGREE.
  My gauge-invariant first-order system (sec. 7) reproduces their constraint (C) and boundary condition (BC) exactly
  (`Psi` = their psi, `X` = their chi, `xi = -2 psi` <-> `(rho^2 Q)' = -2 Psi`).
* Their t-scan: t = 1e-2: -7.700561458 vs -7.700561461; 3e-3: -7.7140236 both; 1e-4: -7.7196034 vs my -7.7196057 (my n = 1500 run is
  less converged). t -> 0 intercept -7.719796 vs HJ closed form -7.7197959. AGREE.
* Norm: the HJ effective theory gives Jordan-frame Z = c I_+, f = 2 I_+, f' = c I_+, so the kinetic coefficient of delta phi_b after
  removing the mixing with the conformal mode is `Z + (3/2) f'^2/f = c I_+ (1 + 3c/4) = 0.8963924`; my 5D value tends to 0.896391-0.896392. AGREE
  (this identification of the combination is mine, made after reading; the number 0.89638 was on disk before).
* `literature/fk_quadratic_form_check_output.json` reports `Q/(lam N) = 0.737` (identity not closing). With my eigenfunction the
  identity closes to 1.2e-6 and is proved above, so that script's Q or N evaluation is likely in error (units/normalisation); its
  qualitative conclusion (N_FK < 0 because B < 0, boundary term dominant: my N_FK,bulk/N_FK,bdy = -8.2e-5) AGREES.

## 10. Prior art (from memory; not re-verified online in this session)

Frolov & Kofman, hep-th/0309002 (scalar perturbations of inflating branes with bulk scalar; continuum at 9H^2/4, tachyonic radion,
self-adjoint form with eigenvalue-dependent b.c.); Gen & Sasaki, gr-qc/0011078 (radion on dS brane, m^2 = -4H^2);
Garriga & Vilenkin, PRD 44 (1991) 1007 (wall fluctuation m^2 = -4H^2 for a 3-brane); Garriga & Tanaka, hep-th/9911055 and
Charmousis, Gregory & Rubakov, hep-th/9912160 (brane bending / radion in Gaussian-normal gauge); DeWolfe, Freedman, Gubser &
Karch, hep-th/9909134 (flat-wall scalar sector). The method is standard; new here are the registered-model numbers, the explicit
GN-gauge mismatch function and the identity N_phys = 12 lam N_FK = 12 Q_FK > 0 (I do not know whether the latter is in the literature).

## 11. Files

`gn_lib.py` (equations, integrator, mismatch), `run_spectrum.py` (scan + constraint monitors, `scan.npy`), `run_roots.py`
(`roots_convergence.json`), `run_T1.py` (`T1_result.json`), `run_tscan.py` (`tscan.json`, `tscan.log`), `run_gauge_test.py`
(`gauge_test.json`), `run_norm_bulk.py` (`norm_bulk.json`), `run_fk_identity.py` (`fk_identity.json`), `RESULT.json`.
