# Oscillation / comparison theorem for the scalar sector of a Z2 dS shell (Task B, Chat 10, 18 Sep 2026)

Status legend: **[P]** proved here from the stated equations and hypotheses; **[S]** standard theorem cited, not re-proved;
**[N]** floating-point evidence only (scripts in this folder); nothing in [N] is a proof.

Input taken from Chat 9 (not re-derived here): the gauge-invariant first-order system and shell condition

    Psi' = -2H Psi - (phi'/3) X,   X' = (phi''/phi') X + (3 lam/(rho^2 phi') - 2 phi') Psi,   H = rho'/rho,  lam = mu^2 + 4,
    m := B X + 3 lam Psi/(rho^2 phi') = 0  at y_b,     B = phi''/phi' + sigma_t''/2 ,

(in longitudinal gauge Psi = psi, X = chi, and m = chi' + 2 phi' psi + sigma_t'' chi/2). "Bound state" = nontrivial solution
with the cone behaviour psi ~ y^alpha, alpha = 1/2 + sqrt(9/4 - mu^2), satisfying m = 0.

## 0. Hypotheses

* **(H0)** (rho, phi) solves the background equations on (0, y_b] with rho > 0, and is real-analytic at the regular cone:
  rho = y + O(y^3) (odd), phi = phi_h + a y^2/2 + O(y^4) (even), **a = U_phi(phi_h)/5 != 0**. [S: analyticity of the regular
  solution of an analytic ODE system at a regular singular point (Briot-Bouquet); not re-proved. See sec. 7, flag F1.]
* **(H1)** phi' != 0 on (0, y_b].
  *Sufficient condition [P]:* phi_h > -1 and phi(y) < phi_0 on [0, y_b], phi_0 = 0.40878... the root of 30 x - 12 - 4 x^3.
  Proof: U_phi = (phi^2 - 1) g(phi), g = (10/3)phi - 4/3 - (4/9)phi^3, g' > 0 on (-1,1), g(phi_0) = 0, so U_phi > 0 on (-1, phi_0).
  The scalar equation is (rho^4 phi')' = rho^4 U_phi with rho^4 phi' -> 0 at the cone. Let y_1 = sup{y <= y_b : phi' > 0 on (0,y)} (> 0 since a > 0).
  On (0, y_1] phi is increasing, so phi in [phi_h, phi_0) and rho^4 phi'(y_1) = int_0^{y_1} rho^4 U_phi > 0; hence y_1 = y_b. QED.
  [N] all shells used below have min phi' = 1.5e-8 > 0 (at the first grid point) and |phi_b| < 3e-4.
* **(H2)** B != 0 (B = 0 is the degenerate case in which mu^2 = -4 is on the verge of becoming physical).

## 1. Sturm-Liouville form and dictionary [P]

Put **P = rho^2 Psi**, p = 1/(rho^2 phi'^2), w = 1/(rho^4 phi'^2), q = 2/(3 rho^2) (all > 0 on (0,y_b] by H0-H1). Then

    P' = rho^2 (Psi' + 2H Psi) = -(rho^2 phi'/3) X    =>   X = -3 phi' (p P'),
    (X/phi')' = (3 lam/(rho^2 phi'^2) - 2) Psi        =>   -(p P')' + q P = lam w P                                  (SL)
    m = -3 B phi' pP' + 3 lam P/(rho^4 phi') = -(3/(rho_b^4 phi_b')) D(lam),   D := beta (pP')(y_b) - lam P(y_b),   beta := rho_b^4 phi_b'^2 B.

So the shell condition is **p P' = (lam/beta) P at y_b** (eigenparameter in the boundary condition), sign(beta) = sign(B), and with

    F(lam) := (pP')(y_b)/P(y_b) = -X_b/(3 phi_b' rho_b^2 Psi_b),     G(lam) := F(lam) - lam/beta,
    m/(B phi_b' Psi_b) = -3 rho_b^2 G(lam).                                                                           (1.1)

The physical (Klein-Gordon/on-shell-action) norm of Chat 9 is N_phys = 12 Q_FK with
Q_FK = int(rho^2 X^2/6 + rho^2 Psi^2) = (3/2) int (p P'^2 + q P^2) =: (3/2) E[P], i.e. **the Dirichlet form of (SL)**, and
N_FK = (3/2) [P,P], where [P,Q] := int_0^{y_b} w P Qbar dy + P_b Qbar_b/beta.

**Remark 1.2 [P].** Using the junction conditions (H = sigma_t/6, phi' = -sigma_t'/2) and the constraint, B and rho_b are
functions of phi_b alone:  B = -2U_phi/sigma_t' - (2/3)sigma_t + sigma_t''/2,  1/rho_b^2 = sigma_t^2/36 - sigma_t'^2/48 + U/6 (at phi_b).
So an enclosure of phi_b gives enclosures of B, rho_b, phi_b', beta by rational arithmetic. ([N] check: `B_closed_form_check.py`, agreement 1e-9.)
B is a cancellation of O(1) terms (1.3339 - 1.3341 + 0.0007 for d = 8/5), but the O(1) parts of dB/dphi_b cancel at phi_b = 0
(-U_phiphi(0) + 4/3 + 2 = 0): in exact rational arithmetic for t = 1/1000, c = 5975949350280/10^13, d = 8/5: B(phi_b = 0) = 5.318490e-4 and dB/dphi_b = -3.27e-4 at phi_b = -5.2e-5,
so B = B(0) + O(t phi_b) > 0 for any |phi_b| << 1; numerically B = -2.682e-4 + t d/2 for all shells of sec. 6, B = 0 at d = 0.536. (Evaluate B in mean-value form, or with a
phi_b enclosure narrower than ~1e-6, to avoid interval dependency blow-up of the O(1) cancellation.)

## 2. The cone endpoint [P, given H0]

Write (SL) as P'' + (p'/p) P' + ((lam w - q)/p) P = 0. By H0, p'/p = -2H - 2phi''/phi' = -4/y + b(y), (lam w - q)/p = lam/rho^2 - (2/3)phi'^2
= lam/y^2 + c(y,lam) with b, c analytic at y = 0 (c linear in lam). y = 0 is a regular singular point with indicial equation

    gamma^2 - 5 gamma + lam = 0,   gamma_pm = 5/2 +- nu,   nu := sqrt(25/4 - lam) = sqrt(9/4 - mu^2)      (alpha = gamma - 2).

**Lemma 1 (Frobenius).** For Re nu > -1/2 there is a solution P_+(y,nu) = y^{5/2+nu} h(y,nu), h(0,nu) = 1, h analytic in (y,nu)
for |y| < R; it extends to (0,y_b] as a solution analytic in nu. For nu != 0 a second solution is P_- = y^{5/2-nu} h_-(y) + K P_+ ln y
(K = 0 unless 2nu is a positive integer); for nu = 0, P_- = P_+ ln y + y^{5/2} k(y).
Proof: the Frobenius recursion has denominators I(gamma_+ + n) = n(n + 2nu), |n(n+2nu)| >= n^2 for Re nu >= 0 (and != 0 for Re nu > -1/2);
the standard majorant argument gives convergence uniform on compact nu-sets, hence analyticity in nu; continuation to (0,y_b] is
analytic dependence of a regular ODE on a parameter. QED. For real lam < 25/4 we normalise P(y,lam) := P_+ (positive near y = 0; this is psi = y^alpha(1+o(1))).

**Lemma 2 (normalisability; limit point).** Let lam in C. A nontrivial solution has int_0 w|P|^2 < infinity iff E[P] < infinity near 0
iff lam not in [25/4, infinity) and P is a multiple of P_+ (with Re nu > 0). Hence (SL) is in the limit-point case at y = 0 (no boundary condition
at the cone may or need be imposed), and there is no normalisable solution at all for real lam >= 25/4 (mu^2 >= 9/4).
Proof: w ~ 1/(a^2 y^6), p ~ 1/(a^2 y^4). For P = c_+P_+ + c_-P_-, Re nu > 0, c_- != 0: |P|^2 w ~ |c_-|^2 y^{-1-2Re nu}/a^2 and p|P'|^2 ~ |gamma_-|^2 |c_-|^2 y^{-1-2Re nu}/a^2 likewise (for gamma_- = 0, i.e. lam = 0, use q|P|^2 ~ 2|c_-|^2/(3y^2) instead): divergent;
for c_- = 0: ~ y^{-1+2Re nu}: convergent. For nu = 0: |P|^2 w ~ |c_+ + c_- ln y|^2/(a^2 y): divergent. For nu = i kappa, kappa > 0:
|P|^2 w = |c_+ e^{i kappa s} + c_- e^{-i kappa s}|^2 (1+O(y))/(a^2 y), s = ln y, and int ds of the bracket over a length S is (|c_+|^2+|c_-|^2) S + O(1): divergent.
Limit point: at lam = 0 the solution P_- ~ 1 is not in L^2(w). QED. So the three notions "psi ~ y^alpha", "finite w-norm", "finite physical norm Q_FK" coincide.

**Lemma 3 (cone boundary terms vanish; closes the Chat 9 gap).** For lam in C \ [25/4, infinity) and P = P_+, as y -> 0:

    (a) P (pP')              = (gamma_+/a^2) y^{2 nu} (1 + O(y^2))   -> 0,
    (b) p (P dP'/dlam - P' dP/dlam) = -(1/(2 nu a^2)) y^{2 nu} (1 + O(y^2))      -> 0,
    (c) p (P_1' P_2bar - P_1 P_2bar') = O(y^{nu_1 + nubar_2})        -> 0   for two such solutions (lam_1, lam_2).

Proof: (a), (c): insert P = y^{gamma}h, P' = y^{gamma-1}(gamma h + y h'), p = y^{-4}(1+O(y^2))/a^2. (b): d/dlam = -(1/(2nu)) d/dnu and
dP/dlam / P = -(1/(2nu)) (ln y + d_nu ln h), so p(P P_lam' - P' P_lam) = p P^2 (P_lam/P)' = -(1/(2nu)) pP^2 (1/y + d_y d_nu ln h), with d_y d_nu ln h = O(y). QED.
In Chat 9 variables (a) is -rho^2 Psi X/(3 phi') and (b) is -(rho^2/(3phi'))(Psi X_lam - X Psi_lam): exactly the two cone terms
that were "only argued (flux ~ y^{2s+3})" in GN_DERIVATION sec. 7; 2s + 3 = 2 nu. The decay is a power y^{2 nu}, nu = sqrt(9/4 - mu^2) > 0: it
vanishes for every mode below the continuum threshold and fails exactly at/above it. [N] `osc_validate.py`: measured slopes 6.3094/6.3090 (2nu = 6.3087),
3.0014/3.0004 (3), 1.7340/1.7326 (1.7321) and prefactors gamma_+/a^2, -1/(2 nu a^2) reproduced to 1e-4 at y = 0.01.

## 3. Theorem 1 (reality, positivity, simplicity) [P]

Let B != 0 and let P be a bound state for some lam in C (normalisable at the cone, pP' = (lam/beta)P at y_b). Then

1. **E[P] = lam [P,P]**, E[P] = int_0^{y_b} (p|P'|^2 + q|P|^2) dy > 0.
2. lam is real, lam != 0, lam < 25/4, and **sign [P,P] = sign lam**. N_phys = 12 Q_FK = 18 E[P] > 0: no ghosts, for either sign of B.
3. **If B > 0 then lam > 0, i.e. mu^2 > -4**, for every scalar bound state.
4. Eigenfunctions to different eigenvalues are [.,.]-orthogonal. If B < 0 there is at most one eigenvalue with lam < 0.
5. Every eigenvalue is geometrically simple, and a simple zero of D(lam): P_b D'(lam) = -beta [P,P] = -beta E[P]/lam != 0.

Proof. (1) Multiply (SL) by Pbar, integrate over (eps, y_b): int_eps (p|P'|^2 + q|P|^2) - [pP'Pbar]_eps^{y_b} = lam int_eps w|P|^2. At y_b,
pP'Pbar = (lam/beta)|P_b|^2; at eps the term -> 0 by Lemma 2 + Lemma 3(a). (2) E > 0 is real, so Im: (Im lam)[P,P] = 0 and Re: (Re lam)[P,P] = E > 0 force
[P,P] != 0, Im lam = 0, lam != 0, sign[P,P] = sign lam; lam < 25/4 by Lemma 2. (3) beta > 0 makes [P,P] > 0. (4) Green's identity with Lemma 3(c):
(lam_1 - lam_2)[P_1,P_2] = 0. If lam_1 != lam_2 were both < 0, [.,.] would be negative definite on span{P_1,P_2}; but that span contains
a nonzero element with vanishing shell value, on which [.,.] = int w|P|^2 > 0. Contradiction. (5) Geometric: Lemma 2 (P_+ is unique up to scale). Differentiating (SL) in lam:
[p(P_lam' P - P' P_lam)]' = -wP^2; integrate with Lemma 3(b):

    p (P_lam' P - P' P_lam)(y) = - int_0^y w P^2 dy        (lam < 25/4, all y in (0,y_b]).                                 (3.1)

At a root pP' = lam P/beta, so P_b D' = beta pP_lam' P - P^2 - lam P P_lam = beta p(P_lam'P - P'P_lam) - P_b^2 = -beta [P,P]. P_b != 0 at a root
(else P_b = P'_b = 0). QED.

Remarks. (i) (3.1) gives for all real lam < 25/4 with P_b != 0:  **dF/dlam = -int w P^2/P_b^2 < 0** and **dG/dlam = -[P,P](lam)/P_b^2**. (ii) For beta > 0 the
problem is the standard right-definite one in the Hilbert space L^2(w) (+) C with weight 1/beta on the second component (operator (P,P_b) -> (w^{-1}(-(pP')'+qP), beta pP'(y_b)));
with Lemma 2 (limit point) it is self-adjoint and > 0 [S: Walter 1973, Fulton 1977]. For beta < 0 it is a Pontryagin space with one negative square.

## 4. Pruefer angle and the counting theorem for B > 0

For real lam <= 25/4 define theta(y,lam) by P = r sin theta, pP' = r cos theta, continuous in y, theta(0+) = 0 (possible since P, pP' > 0 near 0 and
tan theta = P/(pP') ~ a^2 y^5/gamma_+ -> 0). Then theta' = cos^2 theta/p + (lam w - q) sin^2 theta.

**Lemma 4 [P].** (a) theta = k pi is crossed only upwards (theta' = 1/p > 0 there): #zeros of P(.,lam) in (0,y_b] = floor(theta_b/pi), theta_b := theta(y_b,lam);
the number in the open interval (0,y_b) is Z(lam) := ceil(theta_b/pi) - 1. (b) theta_b is continuous on lam <= 25/4 and strictly increasing on lam < 25/4:
d theta/d lam = int_0^y w P^2/r^2 > 0 (from (3.1); continuity at 25/4 from analyticity of P_+ in nu at nu = 0; the branch is fixed uniformly because theta in (0,pi/2) for small y, locally uniformly in lam).
(c) For lam <= 0: P > 0 and pP' > 0 on (0,y_b], so theta_b in (0,pi/2), F(lam) > 0.
Proof of (c): pP'(0+) >= 0 (gamma_+ >= 5). While pP' > 0, P increases and stays > 0, so (pP')' = (q - lam w)P > 0 and pP' increases: it can never reach 0. QED.
(d) Consequently the Dirichlet eigenvalues (P_b = 0) are the solutions of theta_b = k pi, k >= 1; they are > 0, and #{Dirichlet eigenvalues < lam*} = Z(lam*).

**Theorem 2 [P] (B > 0).** Let beta > 0, arccot: R -> (0,pi) decreasing, Theta(lam) := theta_b(lam) - arccot(lam/beta). Theta is continuous and strictly
increasing on lam <= 25/4, Theta(lam) in (-pi,0) for lam <= 0, and lam is an eigenvalue iff Theta(lam) in pi Z, lam < 25/4. Hence for lam* < 25/4

    N(lam*) := #{bound states with lam <= lam*} = floor(Theta(lam*)/pi) + 1 = Z(lam*) + [ G(lam*) <= 0 ],                          (4.1)
    N_tot := #{bound states, mu^2 < 9/4}     = max(0, ceil(Theta(25/4)/pi))  = Z(25/4)  + [ G(25/4) <  0 ],                          (4.2)

where [.] = 1 if true, 0 otherwise, Z = number of zeros of the cone-regular solution in the open interval (0,y_b), and if P_b(lam*) = 0 the bracket is read as 1
(G = -infinity). The k-th eigenfunction (k = 0,1,...) has Theta = k pi and exactly k zeros in (0,y_b); eigenvalues interlace with Dirichlet eigenvalues.
Proof. cot theta_b = F, so the shell condition is cot theta_b = lam/beta, i.e. theta_b = arccot(lam/beta) mod pi. Monotonicity: Lemma 4(b) and beta > 0.
For lam <= 0: theta_b in (0,pi/2), arccot in [pi/2,pi). Count of k >= 0 with k pi <= Theta(lam*) gives the first equality. Write theta_b = Z pi + th, th in (0,pi]
(th = pi iff P_b = 0), cot th = F: th - arccot(lam*/beta) in (-pi,pi) is >= 0 iff F <= lam*/beta. For (4.2) eigenvalues need lam < 25/4 strictly: count k pi < Theta(25/4). QED.

**Corollary (finite checkable criterion) [P].** Let B > 0, mu*^2 <= 9/4, lam* = mu*^2 + 4, and let psi be the cone-regular solution at mu*^2 normalised psi = +y^alpha(1+o(1)). Then

    there is NO scalar bound state with mu^2 <= mu*^2 (for mu*^2 = 9/4: none at all, none embedded in mu^2 >= 9/4 either)
        <=>   (C1) psi > 0 on (0, y_b]      and      (C2) G(lam*) > 0,
    (C2)  <=>  beta (pP')_b - lam* P_b > 0  <=>  m(mu*^2)/(B phi_b' psi_b) < 0  <=>  **phi_b' m(mu*^2) < 0**   (given C1, B > 0),
    explicitly  B chi_b/(phi_b' psi_b) + 3 lam*/(rho_b^2 phi_b'^2) < 0,   m = B chi + 3 lam psi/(rho^2 phi') = chi' + 2 phi' psi + sigma_t'' chi/2.

**The sign is s = -sign(phi_b'); for the HDBLAST shells phi' > 0, so the criterion is: psi has no zero and m(mu*^2) < 0** (with psi > 0; equivalently the Chat 9
`shoot_vec` value is negative). It is an "iff", it is monotone (if it holds at mu*^2 it holds at every smaller value), both conditions are open (robust under
interval enclosures), and m, Psi, X are gauge invariant (Chat 9 GN sec. 5). For mu*^2 = 9/4 equality G = 0 (threshold resonance) also gives no bound state.
For linear stability of the shell alone (no tachyon) mu*^2 = 0 suffices: all bound states then have mu^2 > 0; the continuum starts at 9/4.
Riccati version (no zero <=> F(y) := pP'/P stays finite): F' = q - lam w - F^2/p on (0,y_b], F = gamma_+/(a^2 y^5) (1 + O(y^2)) -> +infinity at the cone, criterion F(y_b) > lam*/beta.

## 5. Theorem 3 [P] (B < 0: exactly one extra mode, below lam = 0)

Let beta < 0. Then (a) there is exactly one eigenvalue lam_- < 0 (mu^2 < -4); it is simple, nodeless (Lemma 4c), has [P,P] < 0 (negative Frolov-Kofman weight) but
positive physical norm 18 E[P]; (b) for lam* <= 0: N(lam*) = [G(lam*) >= 0]; (c) for 0 < lam* < 25/4:

    N(lam*) = 1 + Z(lam*) + [ G(lam*) <= 0 ],        N_tot = 1 + Z(25/4) + [ G(25/4) < 0 ]   (the last assuming G(25/4) != 0).                (5.1)

Unified with Theorem 2:  **N = [B<0] + Z + [ m/(B phi_b' psi_b) >= 0 ]**  (by (1.1)).
Proof. G is C^1 between the poles of F (the Dirichlet eigenvalues, all > 0), and at any root G' = -[P,P]/P_b^2 = -E[P]/(lam P_b^2): roots with lam > 0 are
down-crossings, roots with lam < 0 are up-crossings; so each pole-free interval of fixed sign of lam contains at most one root. On (-infinity,0): G(0) = F(0) > 0, and
G -> -infinity as lam -> -infinity: on J = [y_b - delta, y_b] the Riccati variable F(y) = pP'/P > 0 obeys F' = q - lam w - F^2/p <= K - F^2/p_M, K = max_J(q + |lam| w),
p_M = max_J p; g = c coth(c (y - y_b + delta)/p_M), c = sqrt(K p_M), solves g' = K - g^2/p_M with g(y_b - delta) = +infinity, and u = F - g satisfies u' <= -(F+g)u/p_M with u < 0
initially, so u < 0: F(y_b) <= c coth(c delta/p_M) = O(|lam|^{1/2}), whereas -lam/beta = -|lam|/|beta|. So there is exactly one negative root; (b) follows from the sign pattern of G.
For lam > 0: F jumps from -infinity to +infinity across each pole and F_lam < 0 nearby, so each interval between consecutive Dirichlet eigenvalues, and (0, lam^D_1) (where G(0) > 0), contains
exactly one root; the last, incomplete interval contains one iff G(lam*) <= 0. Lemma 4(d) converts the number of poles below lam* into Z. QED.
Continuity through B = 0: the extra root is lam ~ beta F(0) -> 0 from either side, i.e. mu^2 crosses -4 exactly when B changes sign
([N]: B = 0 at d = 0.536; tachyon at -4.95 for d = 0.4 and -2.87 for d = 0.7; EFT: mu^2 = -4 at d = 0.5365).

## 6. Consistency with the known counts [N] (`osc_shells.py`, `osc_validate.py`, `osc_criterion_S85.py`; float RK4, t = 1e-3, c = 0.5975949350280)

| d | B | roots of m on -60 < mu^2 < 9/4 | Z(9/4) | G(25/4) | m/psi_b at 9/4 | N_tot from (4.2)/(5.1) |
|---|---|---|---|---|---|---|
| 0   | -2.6818e-4 | -7.7178716 | 0 | +9.63e-4 | +4.81e-3 | 1 + 0 + 0 = 1 |
| 0.4 | -6.819e-5  | -4.9452698 | 0 | +2.73e-3 | +3.47e-3 | 1 |
| 0.7 | +8.179e-5  | -2.8662898 | 0 | -1.62e-3 | +2.47e-3 | 0 + 1 = 1 |
| 1.0 | +2.3169e-4 | -0.7895812 | 0 | -3.40e-4 | +1.47e-3 | 1 |
| 1.3 | +3.8192e-4 | +1.2949033 | 0 | -6.48e-5 | +4.61e-4 | 1 |
| 1.4 | +4.3189e-4 | +1.9871330 | 0 | -1.58e-5 | +1.27e-4 | 1 |
| **1.6** | **+5.3187e-4** | none | 0 | **+5.47e-5** | **-5.42e-4** | **0** |
| 2.0 | +7.3185e-4 | none | 0 | +1.38e-4 | -1.88e-3 | 0 |

On all 373 grid values of mu*^2 per shell the theorem count N(lam*) equals the number of roots of m below mu*^2 (0 mismatches); F and theta_b are monotone as proved;
dF/dlam = -int wP^2/P_b^2 to 1.3e-9; E = lam[P,P] at every root to 1e-8 with sign[P,P] = sign lam (negative only for d = 0, 0.4); identity (1.1) to 1e-13.
B < 0 iff mu^2 < -4, as Theorem 1.3 / Theorem 3 require. For d = 8/5 (float, guidance for the certificate):

| mu*^2 | zeros | F | lam*/beta | G | G/F | m/psi_b |
|---|---|---|---|---|---|---|
| 0    | 0 | 3.59215e-4 | 1.94882e-4 | +1.6433e-4 | 0.457 | -1.629e-3 |
| 2    | 0 | 3.59199e-4 | 2.92323e-4 | +6.688e-5  | 0.186 | -6.63e-4 |
| 2.2  | 0 | 3.59198e-4 | 3.02067e-4 | +5.713e-5  | 0.159 | -5.66e-4 |
| 9/4  | 0 | 3.59197e-4 | 3.04503e-4 | +5.469e-5  | 0.152 | -5.42e-4 |

F is almost independent of lam (3.5925e-4 ... 3.5920e-4 on -4 <= mu^2 <= 9/4), so the eigenvalue is lam ~ beta F and a bound state exists iff beta F(25/4) > 25/4, i.e.
B > 4.51e-4 ... equivalently d < ~1.44 [N], in agreement with the Chat 9 scan (last bound state at d = 1.43) and with the EFT value 3(d - d_0)/Z_E = 9/4 at d = 1.438.
Margin for the certificate: the relative enclosure width of F and of beta (i.e. of B, see Remark 1.2) must be well below 15 % at mu*^2 = 9/4 (46 % at mu*^2 = 0).

## 7. Flags: what is NOT proved here

* **F1.** (H0) analyticity of the background at the cone is cited, not proved. Everything in secs. 2-5 uses it only through Lemma 1. (A C^2 version via a Volterra
  contraction, as in Chat 9 `cert_core.py`, would also do; the exponents and Lemma 3 are unchanged.)
* **F2.** For a given shell, (H1) phi' > 0 and the sign of B are verified only in floating point here; they must come from the interval background enclosure
  (sufficient condition in sec. 0; B from phi_b by Remark 1.2).
* **F3.** The first-order system and the shell condition are taken from Chat 9 (two gauges + 1e-9 numerical check against the linearised Einstein equations); not re-derived.
  At lam = 0 (mu^2 = -4) longitudinal gauge is incomplete; in the (SL) problem lam = 0 is never an eigenvalue when B != 0 (Theorem 1.2).
* **F4.** "No bound state => linear stability of the scalar sector" additionally uses self-adjointness/completeness in L^2(w) (+) C (B > 0) [S], and the continuum mu^2 >= 9/4
  being non-tachyonic. Tensor and vector sectors are outside this task.
* **F5.** N_tot for B < 0 assumes G(25/4) != 0 (non-generic equality excluded); for B > 0 no such assumption is needed.
* **F6.** All numbers in sec. 6 are floating point. The statement "S_8/5 has no scalar bound state" becomes a theorem only when (C1), (C2) at some mu*^2 (ideally 9/4; 0 suffices for
  stability), together with F2, are verified in interval arithmetic.
