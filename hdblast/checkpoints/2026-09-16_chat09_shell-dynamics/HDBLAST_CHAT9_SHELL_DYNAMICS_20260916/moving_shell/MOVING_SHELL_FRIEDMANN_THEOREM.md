# Moving-shell Friedmann theorem in the frozen SO(4,1) HDBLAST bulk

Date: 2026-09-16. Folder: `HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/moving_shell/`.
Status: **exact identities proved by hand and machine-checked (exact rationals + floating trajectory checks); physical consequences clearly graded below; no novelty claim for the core identity** (it is the 00-component of the Shiromizu-Maeda-Sasaki / Maeda-Wands projected equations specialised to a conformally flat bulk, see Section 7).

Files: `verify_moving_shell.py` (run `python3 verify_moving_shell.py`, 4 s), `MOVING_SHELL_CHECKS.json` (all results), `MOVING_SHELL_VERIFY_LOG.txt` (printed log of the final run). `run_log.txt` and `verify_run_output.txt` are stale logs of earlier intermediate runs that the sandbox would not let me overwrite or delete; ignore them.

Grading used throughout: **[P]** proved (hand derivation, with exact-rational sampled control of the algebra), **[N]** numerically supported on the registered background, **[C]** conjecture / interpretation.

---

## 1. Setting, sides and signs

Model: `R_AB = d_A phi d_B phi + (2/3) U g_AB`, `Box phi = U_phi`, `kappa_5 = 1`, signature `(-++++)`.

Two charts for the same SO(4,1)-symmetric bulk:

* **y-chart**: `ds^2 = dy^2 + rho(y)^2 [ -d eta^2 + cosh^2(eta) dOmega_3^2 ]`, `phi = phi(y)`, with
  `rho'^2 = 1 + rho^2 (phi'^2/12 - U/6)`, `(rho'/rho)' = -1/rho^2 - phi'^2/3`, `phi'' + 4 (rho'/rho) phi' = U_phi`.
* **conformal chart**: `g = Omega(sigma)^2 eta`, `sigma = R^2 - T^2`, `w = Omega_sigma/Omega`, `Hs = Phi_sigma`, with `w_sigma = w^2 - Hs^2/3` and `4w + sigma w_sigma + 3 sigma w^2 + Omega^2 U/6 = 0`.

Dictionary (checked exactly, control E5, and in floating point, B2): `sqrt(sigma) = b`, `T = b sinh eta`, `R = b cosh eta`, `rho = Omega b`, `dy = Omega db`, `rho' = 1 + 2 sigma w`, `phi' = 2 b Hs/Omega`, `w = (rho' - 1)/(2 sigma)`, `ln sigma = 2 ln y + 2 int_0^y (1/rho - 1/y') dy'` (normalisation `Omega(0) = 1`). With this normalisation the registered shell is at `sigma_b = 12.4474`, `Omega_b = 22.343` [N] (matches the project's "sigma ~ 12.45").

**Shell kinematics.** S^3-symmetric timelike shell, proper time tau, scale factor `a = Omega R = rho cosh eta`.
y-chart: `ydot = sinh psi`, `rho etadot = cosh psi` (psi = rapidity of the shell relative to the static leaves y = const). Unit tangent `u = (sinh psi, cosh psi/rho)`, unit normal pointing toward increasing y: `n = (cosh psi, sinh psi/rho)`. Define

    alpha = rho'/rho,  beta = tanh(eta)/rho,
    X = n.grad ln a = alpha cosh psi + beta sinh psi        ( = K^theta_theta for the normal n )
    H = u.grad ln a = alpha sinh psi + beta cosh psi        ( = adot/a )

Conformal chart: `u = (Tdot, Rdot)`, `n = eps (Rdot, Tdot)`, `Omega^2 (Tdot^2 - Rdot^2) = 1`,
`K^theta_theta = n^A d_A ln(Omega r) = eps [ Tdot/R + 2 w (R Tdot - T Rdot) ]`, `H = Rdot/R + w sigmadot`, `sigmadot = 2 (R Rdot - T Tdot)`. The two charts give identical X and H (E5 exact; trajectory check to 2e-15).

**Which side is kept, and eps.** Junction convention of `BRANE_REHEATING_BRIDGE.md` eq. (2) (= Maartens-Koyama): directed normal n from M- to M+, `[K_mu nu] = -(S_mu nu - S h_mu nu/3)`, `S^mu_nu = diag(-lam - rho_m, (-lam + p_m) x 3)`, `[n phi] = lam' + J`, Z2: `K+ = -K-`. Therefore, on the **kept side with the normal pointing out of the kept bulk**,

    K^theta_theta = +(lam + rho_m)/6,     K^tau_tau = +(lam - 2 rho_m - 3 p_m)/6,     n.grad phi = -(lam' + J)/2        (E8 exact)

The registered configuration keeps the interior `y < y_b` (the side containing the cone), so the outward normal is `+d_y` (eps = +1 evaluated with the interior geometry) and the static conditions are exactly the registered ones `rho'/rho = sigma_t/6`, `phi' = -sigma_t'/2`. The bridge memo's form `eps[...] = -(lam + rho_m)/6` is the same equation written for `K+` (directed normal entering the mirror copy, i.e. eps = -1 in interior coordinates). Keeping the exterior instead would need `1 + 2 sigma w < 0` for positive tension, which does not occur on the registered background (`rho' > 0`). **The Friedmann identity below is quadratic in K and independent of eps**; eps only matters for the first-order system of Section 4.

---

## 2. Theorem 1 (moving-shell Friedmann identity) [P]

> For ANY timelike S^3-symmetric trajectory in the frozen SO(4,1) bulk, if the theta-theta Israel condition `X = (lam + rho_m)/6` holds at an instant, then at that instant
>
>     H^2 + 1/a^2 = (lam + rho_m)^2/36 - 4 w (1 + sigma w)/Omega^2                                  (2.1)
>                 = (lam + rho_m)^2/36 + U(phi_b)/6 - (grad Phi)^2_b/12,     (grad Phi)^2 = 4 sigma Hs^2/Omega^2 = phi'(y_b)^2.   (2.2)

This confirms the orchestrator's claim, including all coefficients.

**Proof A (conformal chart).** With `X = Tdot/R + 2w(R Tdot - T Rdot)`, `Y := H = Rdot/R + 2w(R Rdot - T Tdot)`:

    X^2 - Y^2 = (Tdot^2 - Rdot^2)/R^2 + (4w/R)[Tdot(R Tdot - T Rdot) - Rdot(R Rdot - T Tdot)] + 4w^2[(R Tdot - T Rdot)^2 - (R Rdot - T Tdot)^2]
              = (Tdot^2 - Rdot^2) [ 1/R^2 + 4w + 4 sigma w^2 ],

because the first bracket equals `R (Tdot^2 - Rdot^2)` and the second `(R^2 - T^2)(Tdot^2 - Rdot^2)`. Using `Tdot^2 - Rdot^2 = Omega^-2` and `a = Omega R` gives (2.1). Then `w_sigma = w^2 - Hs^2/3` in the registered constraint gives `4w(1 + sigma w) = sigma Hs^2/3 - Omega^2 U/6`, hence (2.2); `(grad Phi)^2 = Omega^-2 eta^{AB} (2 x_A Hs)(2 x_B Hs) = 4 sigma Hs^2/Omega^2`. (Controls E1, E1b, E2, E10.)

**Proof B (covariant, independent).** For any S^3-symmetric 5D metric `q_ij dx^i dx^j + a(x)^2 dOmega_3^2` and any shell with orthonormal (u, n) in the 2D orbit space, `K^theta_theta = n.grad ln a` and `H = u.grad ln a`, so

    X^2 - H^2 = (grad a)^2 / a^2      =>      H^2 + 1/a^2 = X^2 + [1 - (grad a)^2]/a^2.                      (2.3)

In the y-chart `a = rho cosh eta`, `(grad a)^2 = rho'^2 cosh^2 eta - sinh^2 eta`, so `[1 - (grad a)^2]/a^2 = (1 - rho'^2)/rho^2 =: F(y)`, and the Hamiltonian constraint gives `F = U/6 - phi'^2/12`. Equivalently `X^2 - H^2 = alpha^2 - beta^2` (E4). Because of the SO(4,1) symmetry the "5D mass function" `F` depends on the invariant y only, not on a and eta separately. Chart equivalence `(1 - rho'^2)/rho^2 = -4w(1 + sigma w)/Omega^2` is control E5 and B2 (2.8e-16).

**Remarks.**
1. (2.3) is the standard generalised-Birkhoff statement; what is specific here is that F is a function of the dS-invariant position, not of the scale factor. Hence, unlike the AdS-Schwarzschild case (F = -1/l^2 + mu/a^4), **(2.2) is not a closed equation for a(tau)**: one also needs y_b(tau) (Section 4).
2. The same computation with flat or open slicings of the dS_4 leaves (`a = rho e^eta`, `a = rho sinh eta`) gives `H^2 + k/a^2 = (lam + rho_m)^2/36 + F(y_b)`, k = 0, -1, with the same F [P, by hand; not machine-checked].
3. Only the theta-theta junction is used. The identity is therefore purely kinematical + Hamiltonian constraint; the tau-tau and scalar junctions enter the dynamics (Section 4).

## 3. Corollaries [P]

**3.1 Static shell = M462 (7.1).** psi = 0, eta-independent: `H^2 + 1/a^2 = 1/rho_b^2 = h`, `phi' = -sigma_t'/2`, so `h = sigma_t^2/36 + U/6 - sigma_t'^2/48`. With `sigma_t = 2W + t(1 + c phi)`, `U = W_phi^2/2 - 2W^2/3`:
`sigma_t^2/36 = W^2/9 + tW(1+c phi)/9 + t^2(1+c phi)^2/36`, `U/6 = W_phi^2/12 - W^2/9`, `sigma_t'^2/48 = W_phi^2/12 + t c W_phi/12 + t^2 c^2/48`; the sum is exactly M462 (7.1), and (7.2) for the cubic W (controls E3, E3b; on the registered float solution `theorem rhs - h = -5e-14`, `(7.1) - h = -4e-13`, B1).

**3.2 BPS no-force.** At t = 0 on the flat BPS flow (`sigma = 2W`, `phi'^2 = W_phi^2`): `W^2/9 + U/6 - W_phi^2/12 = 0` identically (E7). So `Lambda_eff = 0` at every position: the whole effective vacuum energy of the registered shell is O(t).

**3.3 Defect form.** With `Z = rho'/rho - W/3`, `Fd = phi' + W_phi` (this orientation: BPS flow is `phi' = -W_phi`, `rho'/rho = +W/3`): `F + W^2/9 = -Fd (Fd - 2 W_phi)/12`, hence
`Lambda_eff/3 = t W (1 + c phi)/9 + t^2 (1 + c phi)^2/36 - Fd (Fd - 2W_phi)/12`.

**3.4 The static shell is a stationary point of Lambda_eff(y).** For `Lambda_eff(y) := 3[sigma_t(phi(y))^2/36 + F(y)]`, using `F' = 2 alpha phi'^2/3`:
`d(Lambda_eff/3)/dy = (2 phi'/3) [ sigma_t sigma_t'/12 + alpha phi' ]`, which vanishes when both static junction conditions hold. On the registered background the stationary point is a **maximum** [N] (table below). Do not read this as an instability proof: Lambda_eff(y) is not an effective potential for a consistent motion (Section 4.3); stability is decided by the perturbation problem.

---

## 4. Theorem 2: the closed moving-shell system in the frozen bulk

### 4.1 Equations [P]
State `(y, eta, psi, rho_m)` plus the brane data `lam(phi)`, equation of state `p_m`, coupling `J`:

    (K)  ydot = sinh psi,   etadot = cosh psi / rho(y)
    (A)  psidot + alpha cosh psi = (lam - 2 rho_m - 3 p_m)/6            [tau-tau Israel; n.a = psidot + alpha cosh psi]
    (E)  rho_m' + 3H (rho_m + p_m) = J phidot_b,   phidot_b = sinh(psi) phi'(y)
    (C)  C := alpha cosh psi + beta sinh psi - (lam + rho_m)/6 = 0       [theta-theta Israel]
    (S)  S := cosh(psi) phi'(y) + (lam'(phi) + J)/2 = 0                  [scalar junction]

On C = 0, (A) is equivalent to the **rapidity law** (E12)

    psidot = tanh(eta) sinh(psi)/rho - (rho_m + p_m)/2.                                               (4.1)

So any matter with `rho_m + p_m > 0` decelerates outward motion and drives the shell toward the cone (toward smaller warp factor) on the time scale `2/(rho_m + p_m)`; this is the dynamical version of the known statement that the static umbilic leaf carries `rho_m + p_m = 0`.

**Off-shell Codazzi law** (derived by hand using `alpha' = -1/rho^2 - phi'^2/3`; control E6 exact at 400 random points; float check B8 to 2e-13 against a drift of 2.8e-5):

    Cdot = -(phidot_b/3) S - H C.                                                                       (4.2)

Consequences: (i) if (A), (E), (S) hold then C = 0 is preserved - the system (K, A, E, S) with C = 0 initial data is consistent and implies Theorem 1 for all tau; (ii) conversely if (K, A, E, C) hold and `phidot_b != 0` then (S) follows - this is the Codazzi/energy identity `rho_m' + 3H(rho_m + p_m) + lam-dot = [n phi] phidot_b` of the bridge memo.

### 4.2 Over-determination for a generic prescribed sigma_t(phi) [P]
With `lam(phi)`, `p_m(rho_m)` and `J(rho_m, phi)` all prescribed, (K, A, E) is a closed ODE flow on the 4D state space, but C = 0 **and** S = 0 must both be preserved. (4.2) propagates C once S = 0; S = 0 must then propagate by itself, which for `sinh psi != 0`, J = 0 requires the additional pointwise relation

    lam''(phi_b) = -2 psidot - 2 cosh(psi) (U_phi - 4 alpha phi')/phi'      along the trajectory.      (4.3)

For a generic function lam (in particular the registered `sigma_t = 2W + t(1 + c phi)`) this fails: C = 0, S = 0, Sdot = 0 cut the state space down to isolated points, the only persistent solution being the static shell. **A moving shell with a generic prescribed tension is incompatible with the frozen SO(4,1) bulk; the bulk must respond** (the shell radiates scalar/gravitational disturbances and the bulk is only SO(4)-symmetric: a 1+1 PDE problem in (y, eta)). Demonstration B8 [N]: registered sigma_t, J = rho_m = 0, kick psi_0 = 0.05, evolved with (K, A): `max|C| = 4.2e-4`, `max|S| = 1.2e-3`, growing exactly as (4.2) predicts.

### 4.3 The two consistent "reconstructed" families in the frozen bulk [P + N]
* **Reconstructed tension.** Treat `lam(tau)` as a state variable with `lam-dot = -(2 cosh psi phi' + J) sinh psi phi'` (this is (S) times phidot_b). Then (K, A, E, lam-dot) is a closed 5D ODE system, C = 0 is conserved, and `lam_rec(phi) := lam(tau(phi_b))` is a legitimate tension function **on any segment where phi_b(tau) is monotone**. These are exact solutions of the full bulk + junction system for the tension lam_rec. Runs B3-B7 [N]: |C| <= 8e-10 (vacuum shells 3e-14), Friedmann residual with an independent finite-difference Hubble rate <= 3e-7 on a scale of 10 (vacuum shells 5e-13), conformal-chart formulae for normalisation/K/H reproduced to 2e-14.
  Caveat: the matter runs B5-B7 turn around (phidot_b changes sign), and lam_rec is then double-valued in phi (flag `phi_b_monotone = false`); across a turning point these are consistency checks of the identities, not solutions of a single theory with one lam(phi).
* **Reconstructed coupling.** Keep any lam(phi) (e.g. the registered sigma_t) and define `J := -lam'(phi_b) - 2 cosh(psi) phi'(y_b)`. Then (K, A, E) close, C = 0 is conserved (run B11 [N]: dust shell on the registered tension, |C| <= 6e-16, J in [-0.091, -0.030]). But J is then a function of the state, not derived from a matter Lagrangian. For a Lagrangian conformal coupling `J = xi'(phi)(rho_m - 3 p_m)` the scalar junction becomes one more constraint; this is the known restriction "equation of state tied to the coupling" of Langlois & Rodriguez-Martinez for moving branes in static dilatonic bulks (hep-th/0106245; review hep-th/0209261 eqs. (104)-(106)).

### 4.4 Theorem 3: vacuum shells are time-shifted hyperboloids, with closed-form tension [P + N]
If `rho_m + p_m = 0` then `K^tau_tau = K^theta_theta`: the shell is totally umbilic. Umbilicity is conformally invariant and the bulk is conformally flat, so the shell is the image of a totally umbilic Minkowski hypersurface; S^3 symmetry and timelike character leave exactly

    R^2 - (T - T_0)^2 = b^2       (two parameters T_0, b; T_0 = 0 is the static leaf),

with unit normal `n^A = (x - x_0)^A/(Omega b)` and `(x - x_0).x = (sigma + b^2 + T_0^2)/2`, hence

    lam_{b,T0}(sigma)/6 = [ 1 + w(sigma) (sigma + b^2 + T_0^2) ] / (Omega(sigma) b),     lam' Hs = d lam/d sigma  (uses only w_sigma = w^2 - Hs^2/3),
    n.grad Phi = Hs (sigma + b^2 + T_0^2)/(Omega b).

[N]: along the integrated vacuum runs B3, B4 the quantities `T_0 = T - R Rdot/Tdot` and `b^2` are constant to 1e-13 and 8e-12 and the closed-form tension matches to 3e-14. T-translation is not a symmetry of the bulk (Omega depends on sigma), so these are genuinely different, non-de Sitter brane cosmologies: `H^2 + 1/a^2 = lam_{b,T0}^2/36 + F(sigma(tau))` with a time-dependent right-hand side. The registered sigma_t(phi(sigma)) coincides with `lam_{b,T0}` only for `T_0 = 0`, `b^2 = sigma_b`, and there only to first order in (sigma - sigma_b); this is Section 4.2 again.

---

## 5. Lambda_eff along the registered background [N]

`Lambda_eff(phi_b) := 3 [ sigma_t(phi_b)^2/36 + U(phi_b)/6 - (grad Phi)^2_b/12 ]` is the value of `3(H^2 + 1/a^2)` that an empty (`rho_m = 0`) shell with the registered tension would record **at the instant it passes** the position y, given the theta-theta condition. It is a local, instantaneous quantity; by Section 4.2 only y = y_b is a consistent permanent location in the frozen bulk. `3h = 4.82790e-4`.

| y | phi | rho | sigma (conf.) | sigma_t | sigma_t^2/36 | F = U/6 - phi'^2/12 | Lambda_eff | Lambda_eff/3h | sigma_t/6 - rho'/rho | phi' + sigma_t'/2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | -1.00000 | 0.020 | 0.000 | 3.33374 | 0.308716 | -0.308642 | 2.2357e-04 | 0.463 | -4.9e+01 | 2.97e-04 |
| 1.0 | -1.00000 | 1.052 | 0.951 | 3.33374 | 0.308716 | -0.308642 | 2.2357e-04 | 0.463 | -5.5e-01 | 2.98e-04 |
| 3.0 | -0.99994 | 4.595 | 6.033 | 3.33374 | 0.308717 | -0.308642 | 2.2359e-04 | 0.463 | -4.1e-02 | 2.95e-04 |
| 5.0 | -0.99684 | 14.419 | 10.103 | 3.33372 | 0.308713 | -0.308638 | 2.2469e-04 | 0.465 | -4.2e-03 | 2.78e-04 |
| 6.0 | -0.97701 | 25.195 | 11.236 | 3.33270 | 0.308525 | -0.308447 | 2.3229e-04 | 0.481 | -1.3e-03 | 2.48e-04 |
| 7.0 | -0.84195 | 43.859 | 11.941 | 3.28650 | 0.300030 | -0.299934 | 2.8921e-04 | 0.599 | -3.9e-04 | 1.81e-04 |
| 7.5 | -0.62183 | 57.315 | 12.182 | 3.08399 | 0.264195 | -0.264070 | 3.7520e-04 | 0.777 | -1.7e-04 | 1.28e-04 |
| 8.0 | -0.22424 | 72.431 | 12.373 | 2.44183 | 0.165626 | -0.165469 | 4.6941e-04 | 0.972 | -4.2e-05 | 5.23e-05 |
| **8.2282 (shell)** | 0.00002 | 78.828 | 12.447 | 2.00095 | 0.111217 | -0.111056 | **4.8279e-04** | 1.000 | -8e-14 | -2e-12 |
| 8.5 | 0.26523 | 85.276 | 12.530 | 1.48314 | 0.061103 | -0.060947 | 4.6691e-04 | 0.967 | 3.7e-05 | -9.0e-05 |
| 9.0 | 0.64773 | 93.771 | 12.670 | 0.88710 | 0.021860 | -0.021725 | 4.0501e-04 | 0.839 | 7.2e-05 | -3.9e-04 |
| 9.5 | 0.85382 | 100.034 | 12.802 | 0.70883 | 0.013957 | -0.013837 | 3.5779e-04 | 0.741 | 8.2e-05 | -1.0e-03 |

Reading the table:
* Two O(0.01 - 0.3) terms cancel to O(t) everywhere (BPS no-force, 3.2): Lambda_eff is between 2.2e-4 and 4.8e-4 across the whole bulk, i.e. `Lambda_eff/3h` between 0.46 and 1.
* Deep AdS side (phi -> -1): `Lambda_eff -> 3[(20/3) t (1 - c) + t^2 (1 - c)^2]/36 = 2.2357e-4` exactly as tabulated (F -> -25/81 = -1/l^2).
* Lambda_eff is maximal at the registered shell (3.4). The last two columns show how badly a static shell would violate the two junctions elsewhere: the Israel mismatch is negative for y < y_b, so `cosh psi = sigma_t/(6 alpha) < 1` has no solution at eta = 0 - an empty registered-tension shell cannot even instantaneously sit at rest inside y_b.

### 5.1 The rho^2 term, high- and low-energy regimes
Expanding Theorem 1: `H^2 + 1/a^2 = Lambda_eff(phi_b)/3 + (sigma_t(phi_b)/18) rho_m + rho_m^2/36`.
* **High energy**, `rho_m >> 2 sigma_t` (here sigma_t ~ 2-3.3, so rho_m >> 4-7 in kappa_5 = 1 units): `H ~ rho_m/6`, the usual Binetruy-Deffayet-Langlois behaviour. [P] as an identity.
* **Low energy**, `rho_m << 2 sigma_t`: the linear term has the standard coefficient `8 pi G_F/3 = sigma_t(phi_b)/18`, i.e. `8 pi G_F = sigma_t(phi_b)/6` (kappa_5^4 lam/6 of Maeda-Wands eq. (23)); it is position dependent: 0.556 at the AdS end, 0.3335 at the registered shell, 0.118 at y = 9.5.
* **But `G_F` is not the Newton constant that governs the response of H^2 to brane energy**, because phi_b and F respond too. Exact statement for vacuum-type energy:

**Theorem 4 (Planck-mass sum rule) [P + N].** With `Z = rho'/rho - W/3`, `Fd = phi' + W_phi`, the bulk equations give `(rho^4 Z)' = rho^2 + rho^4 (2 Z^2 - Fd^2/6)` (E11 exact). Integrating from the regular cone to a static shell with `sigma_t = 2W + dlam(phi)`, where `Z_b = dlam(phi_b)/6`:

    dlam(phi_b)/6 = h I_kept + int_0^{y_b} (rho/rho_b)^4 (2 Z^2 - Fd^2/6) dy,      I_kept = int_0^{y_b} (rho/rho_b)^2 dy.

On the registered solution: lhs 1.6666894e-4, `h I_kept` = 1.6638013e-4, quadratic term 2.888e-7, residual 1.9e-15 (B10). Hence `h = dlam(phi_b)/(6 I_kept) [1 - 1.7e-3]`, i.e. `H^2 = (8 pi G_N/3) dlam` with

    8 pi G_N = 1/(2 I_kept) = 0.48362      (zero-mode value; flat-wall I_+ gives 0.48273),     versus     8 pi G_F = sigma_t/6 = 0.33349.

The 45% difference is exactly compensated in (2.2) by the shift of `F(phi_b)`: M462 (7.1) at phi_b ~ 0 reads `h/t = 1/9 + c/12`, and `c = c_star` makes this `1/(6 I_+)`. An exploratory scan [N, not in the script; float Newton solves with dv = 1e-3] for c in {0.2, 0.3, 0.4, 0.5, c_star, 0.8, 1.0} gave phi_b from -0.58 to +0.24 and `h/t` equal to `(1 + c phi_b)/(6 I_kept)` within 0.18% in every case, while `sigma_t/18` varied from 0.168 to 0.085. (It also shows that c_star is the value for which the detuned shell sits at the wall centre phi_b ~ 0; for other c the O(t) detuning selects an O(1)-different position along the BPS modulus.)
* [C] For genuine matter (`rho_m + p_m != 0`) there is no frozen-bulk solution (Section 4), so the low-energy Newton constant must come from the back-reacted problem; the natural expectation, consistent with Theorem 4, is `8 pi G_N = 1/(2 I_kept)` for the massless-graviton exchange plus a Brans-Dicke-like correction from the bulk scalar/radion sector. This has not been derived here.
* [N] In the frozen bulk with reconstructed tension, matter shells obey (4.1): runs B5-B7 (radiation, dust with J, stiff) reverse within tau ~ 1-2 and plunge toward the cone with |psi| > 2.5 by tau = 1.7-4.0 while `rho_m` grows (contraction, H down to -3.3). This is a statement about the frozen-bulk toy system only.

---

## 6. What is and is not established

Proved [P]: Theorem 1 in both forms and both charts, signs/sides as in Section 1; static reduction to M462 (7.1)/(7.2); BPS no-force; stationarity 3.4; rapidity law (4.1); off-shell Codazzi law (4.2); over-determination criterion (4.3); hyperboloid theorem and closed-form vacuum tension; sum-rule integrand identity. "Proved" means a complete hand derivation given above, with every algebraic step sampled in exact rational arithmetic at 400 random points (14 controls, zero failures). No proof assistant or external referee has checked it.

Numerically supported [N]: everything evaluated on the registered background (floating RK4, step 2e-4 bulk and 1e-3 shell, cubic Hermite interpolation; not interval arithmetic): table of Lambda_eff, maximum at the shell, sigma_b = 12.447, sum-rule residual 2e-15, trajectory checks, c-scan.

Not established / not claimed: any statement about the true moving-shell dynamics with the registered tension (requires bulk back-reaction); stability of the static shell; the physical 4D Newton constant for ordinary matter; any cosmological fit; any novelty or priority.

---

## 7. Prior art (honest statement)

The core identity is **not new**. Verified references (arXiv abstract pages fetched 2026-09-16; full texts only where stated):
* P. Binetruy, C. Deffayet, D. Langlois, "Non-conventional cosmology from a brane-universe", hep-th/9905012; with U. Ellwanger, "Brane cosmological evolution in a bulk with cosmological constant", hep-th/9910219 - the `rho^2/36` Friedmann law and its first integral.
* P. Kraus, "Dynamics of anti-de Sitter domain walls", hep-th/9910149; D. Ida, "Brane-world cosmology", gr-qc/9912002 - moving wall in a static (Schwarzschild-)AdS bulk; this is the method of Proof B.
* P. Bowcock, C. Charmousis, R. Gregory, "General brane cosmologies and their global spacetime structure", hep-th/0007177 - generalised Birkhoff theorem, `H^2 + k/a^2 = rho^2/36 + F(a)`.
* K. Maeda, D. Wands, "Dilaton-gravity on the brane", hep-th/0008188 (HTML full text consulted) - **closest general formula**: projected equations with a bulk scalar, eq. (22) `Lambda_4 = (1/2)[Lambda_5 + lam^2/6 - (d lam/d phi)^2/8]` at kappa_5 = 1 (their Lambda_5 is, to my reading, the bulk potential term; the placement of kappa_5 factors was read through an automated page summariser and should be re-checked against the paper before quoting with units) and eq. (23) `8 pi G_N = kappa_5^4 lam/6`. Control E9 checks that our static `3[lam^2/36 + U/6 - lam'^2/48]` is identically their Lambda_4; the moving-shell term `+phidot_b^2/12` contained in `-(grad Phi)^2/12 = -(n phi)^2/12 + phidot_b^2/12` is the 00-component of their `(2 kappa_5^2/3) T-hat_mu nu(phi)` term with `T-hat_mu nu = D_mu phi D_nu phi - (5/8) g_mu nu (D phi)^2` (checked by hand: `T-hat_00 = (3/8) phidot^2`, so `3(H^2 + k/a^2)` receives `(2/3)(3/8) phidot^2 = phidot^2/4`, i.e. `H^2` receives `phidot^2/12`), and their undetermined Weyl term `E_mu nu` vanishes here because the registered bulk is conformally flat. So Theorem 1 = Maeda-Wands 00-equation with `E_mu nu = 0`, written in closed form.
* A. Mennim, R. Battye, "Cosmological expansion on a dilatonic brane-world", hep-th/0008192 - same covariant projection with non-minimal matter coupling (abstract only consulted).
* D. Langlois, M. Rodriguez-Martinez, "Brane cosmology with a bulk scalar field", hep-th/0106245 (abstract only; PDF could not be parsed here), and D. Langlois, "Brane cosmology: an introduction", hep-th/0209261, Sec. 8.1 eqs. (101)-(106) (HTML consulted): moving Z2 brane in a static dilatonic bulk, three junction conditions, generalised Friedmann equation, non-conservation `rho-dot + 3H(rho + p) = (1 - 3w) xi' rho phi-dot`, and the resulting **constraint tying the equation of state to the coupling** - the direct precedent for Section 4.2/4.3.
* H. Chamblin, H. Reall, "Dynamic dilatonic domain walls", hep-th/9903225 - walls moving in static dilaton bulks with Liouville potentials (tension exponent must match the bulk: an instance of the "reconstructed tension" restriction).
* S. C. Davis, "Cosmological brane world solutions with bulk scalar fields", hep-th/0106271; E. Flanagan, S.-H. Tye, I. Wasserman, "Brane world models with bulk scalar fields", hep-th/0110070 - superpotential-generated bulk solutions with cosmological branes (abstracts only).
* Sum rule (Theorem 4): of the same type as the consistency conditions of G. Gibbons, R. Kallosh, A. Linde, "Brane world sum rules", hep-th/0011225 (abstract verified: integrated Einstein-equation constraints relating brane tensions and bulk scalar gradient energy), combined with the standard zero-mode normalisation `M_Pl^2 = int e^{2A}`; **I did not find or verify a reference containing this exact dS-shell identity**; treat its attribution as open.

What may be specific to HDBLAST (modest, application-level): the identification `F = U/6 - (grad Phi)^2/12` as a function of the SO(4,1)-invariant position for the registered cone-to-shell background; the closed form of all vacuum moving shells as shifted hyperboloids with explicit reconstructed tension; the off-shell law (4.2); and the numerical Lambda_eff profile with its maximum at the certified shell. I found no paper stating the hyperboloid family for a conformally flat scalar bulk, but the argument is elementary and I make no priority claim.
