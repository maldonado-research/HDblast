# Adversarial mathematical verification of R1-R4 (math agent, 16 Sep 2026)

Status: floating-point + exact-rational checks and hand algebra. No interval certificate. Attempted refutation FAILED on every line;
two minor defects in the orchestrator material were found (section G).

## A. Bulk perturbation equations (R1) - brute force
`einstein_linear_check.py` builds the full perturbed 5D metric/scalar (flat-sliced unit dS4, harmonic Y = F(tau) cos(k x), generic mu2),
computes R_AB - d phi d phi - (2/3)U g and Box phi - U_phi numerically (complex-step in eps, 4th-order stencils on exact 2-jets)
and tests the claimed system  xi = -2 psi,  psi' = -2(rho'/rho)psi - phi' chi/3,  chi' = 3(mu2+4)psi/(rho^2 phi') - 2phi' psi + (phi''/phi')chi.
All 15 Einstein components + scalar equation vanish to <= 1.4e-9 (terms O(1-10)) at 20 random points, registered W and generic U.
Sabotage of any coefficient by 2.5-5 % gives residuals 1e-4 ... 2e-2.
Hand algebra: differentiating the Codazzi equation and using the chi' equation gives
psi'' = [2/rho^2 + (4/3)phi'^2 - (mu2+4)/rho^2 + 4(rho'/rho)phi''/phi'] psi - 2(rho'/rho - phi''/phi') psi'  (the 2/rho^2 comes from -2(rho'/rho)'),
which is exactly (M_y) with the coefficient (2+mu2)/rho^2. Cone indicial equation: alpha(alpha-1) = 2 - mu2  => alpha = 1/2 +- sqrt(9/4 - mu2). Threshold 9/4 confirmed.
Horizon regularity: alpha = s + 2 with s = -3/2 + sqrt(9/4 - mu2), so psi ~ (y e^tau)^s y^2, chi ~ (y e^tau)^s: regular on the future cone; the other root diverges.

## B. Gauge / boundary conditions
* Longitudinal gauge is completely fixed for mu2 != -4 (eps_par = 0 from the traceless mu-nu part, then eps^y = 0 from the y-mu part), so any
  non-zero solution is physical. At mu2 = -4 (harmonics with nabla nabla Y = -gamma Y) a residual gauge family exists:
  rho^2 e'' + 4 rho rho' e' + 2e = 0, psi_g = -(rho'/rho)rho^2 e' - e, chi_g = -phi' rho^2 e'. `toy_known_answers.py` (A) feeds it through the
  orchestrator's own `full_rhs`: agreement 1.6e-13 (known-answer test of M_y); mu2 = -3.9 instead of -4 changes psi by 6e-2.
* Traceless Israel condition: induced metric perturbation is pure trace, brane stress is tension x metric, so (nabla nabla - gamma Box/4) zeta = 0: zeta = 0 for mu2 != -4. Correct.
* Trace Israel: delta k = psi' - (rho'/rho) xi = sigma' chi/6 = -phi' chi/3 is identical to the Codazzi constraint at the shell. Correct (no new information).
  The tau-tau / ij components are both delta K^mu_nu = delta k delta^mu_nu (umbilic), so nothing else remains.
* Scalar junction: (1 - xi)(phi' + chi') = -sigma'(phi+chi)/2  =>  chi' + 2 phi' psi + sigma'' chi/2 = 0. Correct. Using the Gauss equation it is
  equivalent to the Frolov-Kofman-type form  3(mu2+4) psi/(rho^2 phi') + (phi'' + sigma'' phi'/2) chi/phi' = 0.
* Pure AdS5 + dS brane: (C) gives psi = C/rho^2 and (H) gives (mu2+4)psi = 0: no scalar mode except the cone-singular mu2 = -4 gauge/radion mode. Known answer reproduced.
The mu2 = -7.72 mode is physical: it is far from -4, has chi_b != 0, and its pole residue matches the 4D modulus (section C).

## C. Norm: tachyon, not ghost
1. Sturm-Liouville form (new here): P = rho^2 psi, p = 1/(rho^2 phi'^2), w = 1/(rho^4 phi'^2), lam = mu2+4:
   -(pP')' + 2P/(3 rho^2) = lam w P,  shell: pP' = (lam/beta)P, beta = rho^4 phi'(phi'' + sigma''phi'/2)|_b.
   Identity  int[pP'^2 + 2P^2/(3rho^2)] = lam (int wP^2 + P_b^2/beta).  LHS > 0, so a mode with mu2 < -4 needs beta < 0 and has negative w-norm N,
   but lam*N > 0 for EVERY mode. Numerically beta = -1.03e4 (phi''+sigma''phi'/2 = -2.68e-4 = O(t)), identity satisfied to 1e-8.
   The "(mu2+4)" sign flip is a property of the variable psi, not of the physical norm.
2. Physical test: couple a source j(phi - phi_b0) on the shell (no stress perturbation at linear order). Junction becomes B = j/2, so
   G(mu2) = chi_b/j = chi_b/(2B). A healthy 4D scalar has Res G = -1/(kinetic coefficient) < 0. Numerics: Res = -6932.17;
   EFT prediction -rho_b^2/(f Z_E) = -6932.10 (f = 2 I_plus, Z_E = c/2 + 3c^2/8). Ratio 1.00001. Sign healthy and the magnitude
   independently confirms the HJ kinetic normalisation.
Not done: a full second-order 5D action.

## D. Hamilton-Jacobi algebra (R3)
Own derivation of the order-2 constraint: (2/3)W[Phi R - 3 Box Phi + M(dphi)^2/2] - W'(Phi' R - M'(dphi)^2/2 - M Box phi) - R/2 + (dphi)^2/2 = 0,
giving the same three equations as the orchestrator. `hj_exact_check.py` (python Fractions, 200 random points, registered and random cubic W,
free integration constant): consistency identity, Z = -W f'/W', f'' relation, mu2 = 3V''/(V Z_E) - 4 - Delta: 0 failures. At phi_b = 0 as identities in I:
c = 2/I - 4/3, Z = cI, Z_E = c/2 + 3c^2/8 = 3/(2I^2) - 1/I (the form transcribed from Brax-van de Bruck-Davis-Rhodes), mu2 closed form: 0 failures.
I could confirm that hep-th/0209158 and hep-th/0309002 exist and their abstracts; I could NOT read the equations of hep-th/0209158 (no PDF reader), so
the statement "equals the published moduli-space metric" rests on the literature agent's transcription.

## E. Toy with exact answer
Test scalar on dS-sliced AdS5 with Robin brane (Langlois-Sasaki-type): exact 2F1 solution (checked against ODE to 1e-9) vs the orchestrator's shooting
logic: roots agree to all 8 printed digits (1.20852285, -4.25644740; none/none for the other two cases; massless Neumann zero mode -1e-14).

## F. Independent numerics (R2, R4)  `independent_spectrum.py`
Own background (start 2e-3, graded steps, own Newton) reproduces phi_h+1 = 8.785515e-7, y_b = 8.228201550, rho_b = 78.828177.
First-order (psi, chi/phi') system: single root mu2 = -7.717871623 (N = 20000 and 40000; orchestrator -7.717871627). One sign change on [-60, 2.2].
Argument principle with ANALYTIC normalisation B (y0/y_b)^alpha: winding +1 on both rectangles; psi_b and chi_b have winding 0.
t-scan: mu2 - mu2_EFT = 1.9234e-2, 5.7723e-3, 1.9243e-3, 5.7731e-4, 1.9246e-4 at t = 1e-2 ... 1e-4; fit -7.71979591 + 1.9244 t - 0.09 t^2 (EFT -7.71979592).

## G. Defects found (minor)
* complex_plane_winding.py divides by psi_b, so it counts zeros of B minus zeros of psi_b; harmless here because psi_b has winding 0 (checked), but the
  script did not check it. On the large rectangle the phase step between samples reaches ~2 rad (< pi, but coarse).
* T_SCAN_ORCHESTRATOR.json loses accuracy at small t: t = 1e-4 value -7.719603854 vs -7.719603462 here (4e-7), and its t = 3e-5 point implies slope 1.95;
  the quoted slope 1.9243 is right, the claim "agreement ~1e-7" holds with the present data (fit intercept differs from EFT by 1e-8).
