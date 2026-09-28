# The O(t) coefficient of the shell tachyon mass: analytic derivation

16 September 2026. Status: analytic derivation (curvature expansion to second order, equivalent to a four-derivative
Hamilton-Jacobi expansion restricted to de Sitter slices plus linear perturbations), reduced to **two explicit
quadratures over the BPS wall**, validated against full 5D cone-to-shell numerics to ~5e-7 relative.
Floating point only; no interval certificate.

Result for the registered shell (c = c_star, linear detuning):

    mu^2(t) = mu0 + s1 t + O(t^2),   mu0 = -4(3c^2-4c+8)/(c(3c+4)) = -7.7197959184,
    s1 = 1.9243896400   (analytic, digits stable under step halving of the quadratures)
    5D numerics:  s1 = 1.9243896(5)  (least-squares extrapolation, t in [2e-3, 1.2e-2])

## 1. Set-up: the bulk scalar as radial coordinate, expansion in X = 1/rho^2

Orientation: y grows from the cone/throat (phi = -1) to the shell, so the BPS flow is phi' = -W_phi, A' = +W/3
(phi = tanh(y - const)). Let P(phi) = phi', H(phi) = rho'/rho, X(phi) = 1/rho^2 (primes below = d/dphi). Exact equations:

    H^2 = X + P^2/12 - U/6,      P H' = -X - P^2/3   (=> P P' + 4HP = U_phi),      X' = -2HX/P.

In the wall region X = O(t). Expand  H = W/3 + X z1 + X^2 z2 + ...,  P = -W_phi + X p1 + X^2 p2 + ... .

* Order X:   W_phi z1' - (2/3) W z1 = -1, IR-regular solution z1 = I(phi) = e^{-2A} int_{-inf}^{y} e^{2A} dy
  (the HJ function, f = 2I; I(0) = I_plus = a).    p1 = 6(1 - 2WI/3)/W_phi = -6 I'.
* Order X^2: all z2-dependent terms cancel and the source vanishes identically:  **z2' = 0**.
  The constant is the coefficient of the bulk black-hole ("dark radiation") term; regularity of the cone
  (exact AdS_- throat, H^2 = 1/l^2 + X) fixes

      z2 = -l_-^3/8 = -729/1000   (l_- = 9/5),      p2 = (6/W_phi)[ -2W z2/3 - I^2 + p1^2/12 ].

  `check_background_expansion.py`: on the full 5D background (t = 1e-3) the residuals H - W/3 - XI - X^2 z2 and
  phi' + W_phi - X p1 - X^2 p2 are O(X^3) with O(1) coefficients all along the wall (in the throat the H residual
  tends to l^5/16 X^3 = 1.18 X^3 as it must).

Junctions with sigma_t = 2W + t V(phi) (registered: V = 1 + c phi):

    X I + X^2 z2 = t V/6,       X p1 + X^2 p2 = -t V'/2       at phi_b.

Leading order: V'/V = 2I'/I (the HJ stationarity condition, giving c_star at phi_b = 0) and X_b = x1 t, x1 = V/(6I).
Next order: phi_b = phi0 + beta t,

    beta = - x1 (p2/I - p1 z2/I^2) / [ (p1/I)' + 3 (V'/V)' ] .

Registered case: beta = [a^3(4-3c^2) - 4 zeta] / [2 a^3 (3c^2-4c+8)] = 0.0227845923 (zeta = 729/1000, a = 6/(3c+4));
5D value (phi_b - 0)/t -> 0.02279, and phi_b - beta t = 5.46e-3 t^2 smoothly over t in [1e-3, 1.6e-2].

## 2. Perturbations: exact mass formula and a Riccati equation

From the orchestrator's longitudinal-gauge system (Codazzi + Gauss) one finds the first-order pair (d/dy)

    psi_y = -2H psi - P chi/3,      chi_y = [ (3 mu^2+12) X/P - 2P ] psi + (phi_yy/P) chi,

so the scalar junction chi_y + 2 phi_y psi + (sigma_t''/2) chi = 0 is **exactly**

    mu^2 = -4 - E_b / (3 X_b s_b),    E = phi_yy + P sigma_t''/2 = P (F' + t V''/2),  F = P + W_phi,   s = psi/chi.     (*)

(the "-4" is the Frolov-Kofman shift; E vanishes identically on the BPS wall - this is why mu^2 is O(1) although the
shell condition is degenerate at t = 0). s obeys the Riccati equation (d/dphi)

    s' = [2 - (3mu^2+12) X/P^2] s^2 - ((P' + 2H)/P) s - 1/3 .

* Order 0: the IR-regular flat-wall zero mode psi_C = 1 - 2WI/3, chi_C = 2 W_phi I gives **s0 = -I'/(2I)**.
  With (*) at first order this reproduces mu0 (= the HJ two-derivative result) including general V''.
* Order X: s = s0 + X s1, and the linearised Riccati equation has integrating factor I^2 W_phi:

      (I^2 W_phi s1)' = -(3mu^2+12) I'^2/(4 W_phi) + I^2 I' - 6 (W_phiphi - W) I I'^2 / W_phi .

  IR-regular solution (lower limit at the throat, where I^3 -> l^3/8 = -z2):

      s1 = [ (I^3 + z2)/3 - 6 K2 ] / (I^2 W_phi)  +  (3 mu^2 + 12) K1 / (4 I^2 W_phi),
      K1(phi_b) = int_{-inf}^{y_b} I_phi^2 dy,      K2(phi_b) = int_{-inf}^{y_b} (W - W_phiphi) I I_phi^2 dy,

  along phi = tanh y, with I_phi = (2WI/3 - 1)/W_phi. At phi_b = 0:
  K1 = 0.029376124556, K2 = 0.050017165134 (step-halving stable to 1e-13).
  Check: 5D psi/chi at the shell, t = 1e-3: -0.14936484; s0(phi_b) + X s1 = -0.1493648.

Why the expansion is analytic in t: in the throat the master equation has exponents e^{4y} (regular series) and
e^{-2y/l}; 4l + 2 = 9.2 is not an even integer, so there are no logarithms and the non-analytic admixture at the
wall is O(h^{4.6}). The scalar profile's own response mode is O(t^{5.6}). Hence mu^2(t) is a power series up to these.

## 3. The slope

With T0 = G0 - 3 W_phi V'' I/V, T1 = G1 + 3V''(p1 I - W_phi z2)/V, G0 = -W_phi(k0 p1 + p1'),
G1 = p1(k0 p1 + p1') - W_phi(k1 p1 + 2 k0 p2 + p2'), k0 = 2W/(3W_phi), k1 = k0(3I/W + p1/W_phi):

    mu^2 = -4 - [T0 + X T1] / (3 [s0 + X s1])   at phi_b,   mu0(phi) = -4 - T0/(3 s0),
    s1_slope = beta mu0'(phi0) - x1 [ T1/(3 s0) - T0 s1/(3 s0^2) ]_{phi0, mu^2 = mu0} .

Registered shell (phi0 = 0, V = 1 + c phi, a = 6/(3c+4), zeta = 729/1000, Q = 3c^2-4c+8):

    beta      = [a^3(4-3c^2) - 4 zeta]/(2 a^3 Q)
    mu0'      = 8ac/3 - 56a/3 + (32a/9)(2-c)/c^2
    G0 = 4a(c-1),   G1 = a^2(6c^2+12c-8) + 28 zeta/3,   s0 = -c/4
    S  = -[(a^3 - zeta)/3 - 6 K2]/a^2 - (3 mu0 + 12) K1/(4a^2)
    s1 = beta mu0' + (2/(9ac)) G1 + (8/(9 a c^2)) G0 S
       = (-0.0734055) + (2.9488104) + (-0.9510152) = 1.9243896400

(three pieces: shift of the equilibrium position; curvature correction of the junction/background; curvature
correction of the mode profile, the only place where the new quadratures K1, K2 enter).
`slope_closed_form.py` checks this closed form against the general-phi0 routine (agreement 1e-10).

## 4. Validation against the full 5D problem (`validate_5d.py`)

Full cone-to-shell background Newton solve + shooting of the master equation (Richardson in step number), several t,
polynomial extrapolation of (mu^2 - mu0)/t:

| case | c | mu0 | s1 analytic | s1 from 5D |
|---|---|---|---|---|
| phi0 = 0 (registered) | 0.5975949 | -7.7197959 | 1.92438964 | 1.9243896(5) |
| phi0 = -0.3 | 0.3514241 | -6.9527741 | 1.47634094 | 1.4763407 |
| phi0 = -0.6 | 0.1909850 | -6.6431420 | 1.31445779 | 1.3144588 |
| phi0 = +0.3 | 1.1679402 | -8.8835813 | 3.11390860 | 3.1139080 |
| phi0 = 0, d = 1 | 0.5975949 | -0.7868640 | -2.72392556 | -2.72393(1) |

5D mu^2 values at t <= 1e-3 carry noise of a few 1e-9 (phi_h + 1 ~ t^1.8 loses relative precision in double
arithmetic), which limits the extrapolation to ~1e-6; t = 1e-3 alone was excluded from the high-accuracy fit.

## 5. By-products

* `resummed_predictor.py`: solving the two truncated junction conditions for (phi_b, X_b) at finite t and using (*)
  with s0 + X s1 gives mu^2 with an error 1.03 t^2 (registered), 0.10 t^2 (d = 1), ~0.6 t^2 (tuned d).
* Near the slow-roll tuning d -> d0 the plain t-series has a small radius of convergence (beta ~ 1/mu0; s1 = -50.85 at
  d = 1.105924): at t = 1e-3 the 5D mass is mu^2 = -0.08988, not -0.0525. More strongly (`slowroll_floor.py`): the
  O(t) displacement term turns d = d0 into an avoided saddle-node; hilltop and minimum coexist and the hilltop mass
  obeys mu^2 <= -2.31 sqrt(t) for any d (t = 1e-3: -0.07309 at d = 1.1135; 5D check -0.0730905). With a purely
  quadratic detuning, mu^2 = -0.0525 therefore needs t <~ 5e-4 (or a cubic term in the detuning). Numerically
  supported (predictor at two t, one 5D point); not proved.

## 6. Status and prior art

Derived: sections 1-3 (hand algebra, each step checked numerically). Numerically supported: the IR-regularity
prescription (lower limits at the throat, z2 = -l^3/8), confirmed by the 1e-6-level 5D agreement for five cases.
Not attempted: closed forms for K1, K2 (they are nested integrals of e^{2A} on the tanh wall; no reduction found);
the O(t^2) coefficient (5D fit: about -0.094 for the registered shell); interval certification.
Prior art (from memory, not re-verified online in this session): Hamilton-Jacobi/holographic RG derivative expansion
(de Boer-Verlinde-Verlinde hep-th/9912012), braneworld gradient expansion to next order (Kanno-Soda hep-th/0207029),
de Sitter brane radion mass and the -4 shift (Frolov-Kofman hep-th/0309002), curved-wall first-order/fake-superpotential
formalism (Skenderis-Townsend hep-th/0602260). The method is standard; the explicit second-order result for the
registered model, formula (*), z2' = 0 and the K1, K2 reduction are new to this project, external novelty not claimed.
