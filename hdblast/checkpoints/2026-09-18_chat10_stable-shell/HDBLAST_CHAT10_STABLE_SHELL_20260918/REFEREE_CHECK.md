# REFEREE CHECK - Chat 10 stable shell S_8/5 (Task D, 18 Sep 2026)

Scope: Tasks A (existence/), B (oscillation_theorem/), C (stability_certificate/). Nothing in those folders was modified;
replays were run on copies in referee/ (existence_copy, stability_copy). Logs: referee/replay_existence_log.txt,
referee/replay_stability_log.txt, referee/float_spotcheck.py + float_spotcheck_output.json.

## 1. Replays (run by the referee)

| replay | result | wall time |
|---|---|---|
| existence/REPLAY.py | PASS: 29/29 gates, certificate identical except runtime; negative control (box shifted 4 radii) rejected at "Gt_1 < 0 on p_lo", Gt_1 in [2.3665e-30, 2.3666e-30] | 13.6 s (certify alone 7 s) |
| stability_certificate/REPLAY.py | PASS: 48/48 gates, certificate identical except runtime; negative control (d = 13/10 in B only) rejected at "b < 0, mu*^2 = 2", b = +1.137155e-4 | 4 min 46 s (certify alone 143 s) |

Also checked: source_sha256 in both delivered certificates match the current source files; cert_core.py is byte-identical to the
Chat 9 core and identical in both folders; exist_lib.py identical in both folders; C's inputs/EXISTENCE_CERTIFICATE_S85.json is
byte-identical to A's certificate. (A's run_log.txt is stale, as A disclosed; replay_log_final.txt is the valid log.)

## 2. Attack on the logic

### 2A. Existence (Task A) - no defect found
* Miranda hypotheses are verified on the WHOLE boundary: Gt_2 on both y-faces for all p at once (Taylor model in u over the full
  p-range), Gt_1 on both p-faces for all y in [ya,yb] (state at exact u = -1/+1, Taylor polynomial in tau composed with
  tau = (1+v)/2, Lagrange remainder over a verified a-priori box). Correct component on the correct pair of faces. det C != 0 is an
  exact rational gate, so Gt = 0 => G = 0. C, p_center, y_center enter as exact Fractions parsed from strings; the float entries of
  CENTER_GUESS.json are never read into the proof path. It is a Miranda test, not Krawczyk: existence only, no uniqueness (as stated).
* No float in the proof path: checked certify_existence.py and exist_lib.py line by line; floats only in print-outs and detail
  strings. iv_sqrt uses integer isqrt with outward +1.
* A-priori boxes: the verified inclusion x0 + [0,h] f(B) subset B is on the AUGMENTED polynomial system with G carried as an
  extra variable (G' = U2 - G^2 + 4HG, s' = s(G-4H)). I checked that this is legitimate: along the augmented flow
  d/dy (U1(phi) - G s) = 0, so G = U1/s persists, and s = s0 exp(int(G-4H)) cannot vanish. Hence s > 0 on (0,y*] really follows
  (cone gate sigma = s/y > 0 on (0,1/20], then every verified step). The gate label "s > 0 on every a-priori box" only tests the
  last tube; the statement is nevertheless true for the reason just given. The v2 inflation heuristic changes only the search.
* Cone start: exact rational self-map/contraction gates, uniform in |p+1| <= 1e-6; Cauchy-estimate tails. Inherited analytic
  premise: fixed point = the regular cone solution (uniqueness of bounded solutions). Accepted as standard; not machine-checked.
* Not re-audited line by line: the Chat 9 core (TM.cauchy truncation bookkeeping, IV rounding). I read recip, taylor,
  vector_field, apriori_box and step and found nothing wrong; the core is unchanged since the Chat 9 tachyon certificate.
* Continuity of (p,y) -> state on BOX (needed by Miranda) is standard continuous dependence; fine.

### 2B. Oscillation theorem (Task B) - proof complete modulo (H0); sign convention consistent
* Singular endpoint: handled properly. Lemma 2 shows the non-regular Frobenius branch is non-normalisable for every nu with
  Re nu > 0 (|P_-|^2 w ~ y^{-1-2Re nu}), so the cone is limit-point and no boundary condition is hidden there; Lemma 3 gives the
  vanishing of the cone terms as exact powers y^{2nu}. Lemma 3(b) degenerates at nu = 0 (prefactor 1/(2nu)); B uses it only for
  lam < 25/4 and takes lam = 25/4 by continuity of theta_b (analyticity of P_+ in nu at nu = 0). That step is stated in one line,
  not written out; I believe it, but it is the thinnest point of B.
* Eigenparameter-dependent boundary condition: handled by the modified Pruefer function Theta = theta_b - arccot(lam/beta),
  strictly increasing for beta > 0, Theta in (-pi,0) for lam <= 0. The count floor(Theta/pi)+1 = Z + [G <= 0] follows. Correct.
  Reality of lam for complex candidates: Theorem 1 uses only Lemma 2 + 3(a) and B != 0; correct for complex lam.
* Sign convention against known counts (B's own table, float): d = 0 (B<0, G>0): 1+0+0 = 1; d = 0.7, 1.0, 1.3, 1.4 (B>0, G<0): 1;
  d = 1.6, 2.0 (B>0, G>0): 0. m/psi_b = -3 rho_b^2 B phi_b' G gives m > 0 for d <= 1.4 and m < 0 for d >= 1.6, as tabulated.
  C's b = lam w + B R satisfies m = (3 psi/s) b, so C's "b < 0" is exactly B's (C2). Consistent.
* Remaining soft spots: (H0) cone analyticity cited (C's disc-algebra contraction does supply holomorphy on |y| < 1/2 and
  a = U_phi(phi_h)/5 > 0 follows from the gated sigma > 0); Frobenius majorant argument "standard"; Theorem 3 (B < 0) is not
  needed for S_8/5; self-adjointness/completeness cited only. osc_shells.py exceeds the 6-minute guideline (cached).
* IMPORTANT: C does not use B's counting theorem. C needs from B only: Lemma 1-2 (bound state <=> regular branch, none for
  real mu^2 >= 9/4) and Theorem 1 (mu^2 real when B > 0). Those parts are the solid ones.

### 2C. Stability certificate (Task C) - covers every background in BOX and every real mu^2 < 9/4
Re-derived by hand (all agree with the script):
* X' and the Riccati equation R' = R^2 + (2G-6H)R + lam w - (2/3)s^2 from the master equation in the task statement;
  shell condition = (3 psi/s)(lam w + B R); B = G - 4H + 2 phi + t d/2; cone limit yR -> -5/2 - nu.
* y r' = (r-r_+)(r-r_-) + delta1 r + delta0 with delta1 = 2y sigma'/sigma + 2y beta, delta0 = lam(y^2 w - 1) - (2/3)y^2 s^2,
  and the gate function M(k); Horner/tail bounds for beta/y, sigma, sigma'/y (needs s_2 = 0 exactly: gated; justified by parity).
* Lemma C incl. nu = 0. I tested it against the obvious counterexample: at nu = 0 the logarithmic solution also has r -> r_+,
  but r - r_+ = -1/ln y makes E ~ |ln y| unbounded, so the argument correctly does NOT apply to it; it applies only to the
  log-free solution (r - r_+ = O(y^2)), which is the normalisable-limit branch. OK.
* Lemma M (a),(b): correct; R cannot reach -infinity forwards, D' = (...)D + (lam2-lam1)w forbids a first zero of D, and
  r_1(0+) < r_2(0+) strictly since nu_1 > nu_2 >= 0. One upper Riccati solution at lam* = 25/4 therefore bounds ALL real
  mu^2 < 9/4 (including mu^2 <= -4); final inequality uses B > 0 and w > 0, both gated on the tube.
* Coverage: BOX rationals are reproduced exactly from A's certificate (gated); same P_TM; the final tube covers all p and all
  y in [ya,yb]; R components are inside verified a-priori boxes at every step (finiteness => psi > 0). Floats appear only in the
  choice of k_lo, k_hi, which are then exact rationals and gated. Covered.
* Weak points (none fatal): (i) Lemmas P, C, M are pencil-and-paper; (ii) the perturbation equations and junction condition
  are INHERITED and verified only in floating point (1e-9) plus agreement with Frolov-Kofman - this is the main non-rigorous
  premise of the whole chain; (iii) only REAL mu^2 is covered by the computation; complex mu^2 is excluded by B's Theorem 1
  with the certified B > 0; (iv) the negative control alters d only inside B (not a physical d = 1.3 shell) - it shows
  non-vacuity only; (v) if one distrusts the nu = 0 barrier, the fallback is "no bound state with mu^2 <= 11/5".

## 3. Independent float spot-check (guidance only; referee/float_spotcheck.py, 16 s)
Own RK4 loop (y0 = 2e-3, n = 20000, psi rescaled each step, psi zero counter), 495 values of mu^2 in [-200, 2.2499]:

| shell | B | psi zeros | sign changes of b | roots | b(0) | b(2.2499) |
|---|---|---|---|---|---|---|
| S_8/5 at Task A centre (phi_h, y_b) | +5.318661e-4 | 0 | 0 | none | -5.428004e-4 | -1.80673e-4 |
| d = 1.6 (Task B float shell) | +5.318661e-4 | 0 | 0 | none | -5.428004e-4 | -1.80673e-4 |
| d = 1.3 (Task B float shell) | +3.819207e-4 | 0 | 1 | mu^2 = +1.29490 | -2.084e-4 | +1.537e-4 |

b < 0 on the whole grid for d = 1.6 (max = the threshold end value). Agrees with C's certified b(0) = -5.4280037e-4 and
b(9/4) = -1.8065687e-4 and with the known count (exactly one state, mu^2 = 1.2949, at d = 1.3). Junction residual with Task A's
certified centre in the float background: 2e-14.

## 4. Strongest SAFE theorem statement

**Theorem (computer-assisted; exact rational parameters t = 1/1000, c = 5975949350280/10^13, d = 8/5).**
Consider the regular-cone SO(4,1) solutions of the 5D Einstein-scalar background system with W = 1 - phi + phi^3/3 and the
Z2 junction conditions rho'/rho = sigma_t/6, phi' = -sigma_t'/2, sigma_t = 2W + t(1 + c phi + d phi^2/2).
(a) There exists a solution (phi_h, y_b) of both junction conditions with phi_h = -1 + 8.78437001888956254016e-7 (+- 2^-100),
y_b = 8.22819172487514216035 (+- 2^-84), rho_b in [78.8290595254422237385095, ...096], phi_b = -5.2155053891198384e-5,
phi' > 0, rho'/rho > 0 on (0, y_b], and B = phi''/phi' + sigma_t''/2 in [5.31866099143642445554e-4, ...561e-4] > 0.
(b) For every such solution inside the certified box, and ASSUMING the longitudinal-gauge scalar master system and scalar
junction condition of Chat 9 (= Frolov-Kofman), the scalar perturbation problem has no normalisable mode with
mu^2 in C \ [9/4, infinity): no tachyon, no zero mode, no massive bound state; psi of the regular solution is nodeless and the
junction mismatch b = (mu^2+4)/rho_b^2 + B R_b is <= -1.8065e-4 < 0 for all real mu^2 < 9/4 (margin >= 18 % of lam w).
Without the threshold (nu = 0) barrier the same is certified for mu^2 <= 11/5.
Proof ingredients: outward-rounded 2^-384 integer arithmetic (Poincare-Miranda for (a); cone barrier + validated Riccati
integration + monotone comparison for (b)); hand lemmas P, C, M; Task B Lemmas 1-2 and Theorem 1 (normalisable = regular
branch; mu^2 real since B > 0; nothing normalisable on mu^2 >= 9/4).

Do NOT claim more than this. In particular "S_8/5 is stable" should be worded "S_8/5 has no scalar bound state below the
continuum threshold 9/4 H^2; the scalar sector is mode-stable".

## 5. What remains unproved
1. Uniqueness of the junction root: not in BOX (no Krawczyk/interval Newton; float det J = 0.0256 only), not globally, and not
   that y* is the first junction point along the flow. (b) holds for every root in BOX, so this does not weaken (b).
2. The perturbation equations and the scalar junction condition themselves: derived in two gauges and float-checked to 1e-9,
   never verified symbolically/exactly. This is the largest remaining gap; an exact (rational/symbolic) verification that the
   master system implies the linearised Einstein + Israel equations would close it.
3. Hand lemmas P, C, M and B's Lemma 1-3 / Theorem 1 are not machine-checked; (H0) analyticity and the identification of the
   contraction fixed point with the regular solution are cited as standard.
4. Completeness / self-adjointness in L^2(w) (+) C (needed to pass from "no unstable mode" to linear stability of arbitrary
   scalar data) is cited, not proved. Continuum mu^2 >= 9/4 is non-tachyonic but its decay properties are not analysed.
5. Tensor and vector sectors, and the mu^2 = -4 (lam = 0) sector where longitudinal gauge degenerates (B shows lam = 0 is not an
   eigenvalue of the SL problem, but the gauge-complete statement is not re-derived).
6. Nonlinear stability; stability against SO(4,1)-breaking background deformations; what selects d (d = 8/5 is an input).
7. Parameter identity: the exact c used differs from the Chat 9 c_star enclosure beyond the 13th digit, so S_8/5 is "a" nearby
   exact-parameter shell. The window "stable for d >= 1.44", the threshold d ~ 1.1135, and all other d values remain float/EFT.
8. B's counting theorem at lam = 25/4 (continuity step) and Theorem 3 are proved on paper only; not needed for (b).
