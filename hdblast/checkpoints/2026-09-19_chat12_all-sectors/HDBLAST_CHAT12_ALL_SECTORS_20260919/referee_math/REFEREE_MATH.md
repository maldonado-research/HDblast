# Referee report (mathematics) - HDBLAST Chat 12, "all sectors" for the shell S_8/5

19 Sep 2026. Skeptical AI referee (not external peer review). Everything I wrote is in this folder; nothing else was modified
(the three scripts were re-run on copies in `rerun/`, because they write their JSON into the working directory).

Legend: **[P]** proved (exact symbolic identity or short analytic argument given here, under the stated hypotheses);
**[N]** floating-point evidence only; **[C]** cited/standard, not re-proved; **[open]** not established.

## Verdict

I tried to refute the three new sector analyses and did not succeed. All symbolic checks reproduce, 22 deliberately wrong
variants are rejected, and the logical gaps I found could be closed (items G1-G4 below) rather than turned into counterexamples.
Two of the arguments as written in the scripts are weaker than they need to be and should be replaced by the stronger versions
given here: (i) the tensor equation was checked for one polarisation only, (ii) the vector mode `Bv = C/rho^2` was excluded by
"cone regularity", which is heuristic because in the Lorentzian section y = 0 is a horizon, not a point. Neither changes the conclusion.

The combined claim is supported in this form (see sec. 5 for the exact wording): **no normalisable unstable or bound mode in the
tensor, vector and special-harmonic sectors; together with the Chat 10 scalar certificate, S_8/5 is mode-stable in every sector.**
It is conditional on the same hand lemmas as Chat 10, and it is mode stability, not linear or non-linear stability.

## 1. Reproduction and negative controls

* `rerun/`: the three scripts pass (7 + 7 + 8 checks, 1-2 s each, SymPy 1.14.0); the JSON files are byte-identical to the archived ones.
* `negative_controls.py` (text mutation of copies, `NEGATIVE_CONTROLS_RESULT.json`): 22 mutants, all rejected by the expected check:
  tensor `4H->3H`, `m2->2m2`, `9/16->1/2`, `U/8->U/6`, `Q` coefficient `3/2->1`, weight `rho^{3/2}->rho^2`, U-range `3/20->1/5` (U has a root near 0.16),
  U-bound range `[-1,0]->[-1,1/10]`; vector `2H->3H`, `(Box+3)->(Box+2)`, `(Box+3)->(Box+4)`, gauge rule `rho^2->rho`, non-transverse harmonic;
  special: Hessian sign, wrong harmonic `e^{-tau}`, `mu2=-4 -> -3`, Codazzi `3->2`, gauge ODE `4->3`, `sigma''/2 -> sigma''`, `B -> B+1`,
  background junction `-2phi' -> -phi'`, `6H -> 4H`.
  (Two of my own first-draft mutants were faulty, not the scripts: one replaced text inside a comment, one changed W to a potential that is also negative
  on the interval, so it was not a violation. Both were fixed/removed.)

## 2. Tensor sector

**G1 (gap closed) - all polarisations.** T1 checks a single helicity-2 component in flat slicing. `general_tensor_check.py` takes a *general* symmetric
gamma-traceless `TT_{mu nu}(tau,x)` (nine arbitrary functions; not assumed transverse, not assumed an eigenfunction) and proves, for all 15 components,

    EQ_{mu nu} = -(rho^2/2)(h'' + 4H h') TT - (h/2)(Box - 2)TT + h nabla_(mu D_nu),   EQ_{y nu} = (1/2) h' D_nu,   EQ_yy = 0,   D_nu = nabla^mu TT_{mu nu}.   [P]

So every transverse-traceless harmonic with `(Box - 2)TT = m2 TT` (all helicities, including the helicity-1 and helicity-0 parts of massive
gravitons) obeys `h'' + 4Hh' + m2 h/rho^2 = 0`, and the script's `m2` is the Fierz-Pauli mass (G4 check). Controls `4H->3H`, `Box-2 -> Box-3` fail.

**Shell condition [P, by hand].** For a Z2 shell with action `-int sigma_t(phi) sqrt(-q)` the junction is `K_{mu nu} = (sigma_t/6) q_{mu nu}`. In the TT sector
the shell cannot bend (bending is a scalar), `delta phi = 0`, and `K_{mu nu} = (1/2) d_y g_{mu nu}`:
`rho rho' h TT + (rho^2/2) h' TT = (sigma_t/6) rho^2 h TT`, i.e. `h'(y_b) = 0` by the background junction. Correct. (Assumes no anisotropic stress on the shell.)

**SUSY-partner argument - complete? Yes, with these details [P given (H): rho>0, U<=0 along the bulk, background analytic at the cone].**
`V1 - 9/4 = rho^2(-5U/8 - 3phi'^2/16)`, `V2 - 9/4 = rho^2(9phi'^2/16 - U/8)`; both are `O(e^{2z})` as `z -> -infinity`, so (Levinson) solutions behave as
`e^{+-kappa z}`, `kappa = sqrt(9/4 - m2)`. For `m2 < 9/4` an L^2 eigenfunction is `u ~ e^{kappa z}`, with `u_z` likewise; the endpoint `z = -infinity` is limit point
(bounded potential), so no condition is imposed there. Then (a) `m2 int u^2 = int (Qu)^2 - [u Qu]` with `Qu(z_b) = 0` and exponential decay: `m2 >= 0`, and `m2 = 0`
iff `Qu = 0` iff `u = c rho^{3/2}` (h = const): the zero mode is simple, **a second zero mode is excluded** (the other m2 = 0 solution `h' = C/rho^4` violates
`h'(y_b)=0` and is not normalisable). (b) For `0 < m2 < 9/4`: `v = Qu` is not identically 0, `v(z_b) = 0`, `v ~ (kappa - 3/2)e^{kappa z}`, `QQ^+ v = m2 v`, boundary terms
vanish, hence `m2 int v^2 = int v_z^2 + V2 v^2 > (9/4) int v^2`: contradiction. Note that V1 itself is **not** above 9/4: next to the shell
`V1 - 9/4 ~ -0.083 rho^2 ~ -5e2` [N]; a naive "potential above threshold" argument would be wrong, the partner is genuinely needed.

**Independent elementary proof (new) [P given rho' >= 1, i.e. U <= phi'^2/2].** `R = rho h'/h` obeys `R_z = -3 rho' R - m2 - R^2`, `R(-inf) = s_+ = (-3+sqrt(9-4m2))/2`.
For `0 < m2 < 9/4`: `s_+ in (-3/2, 0)`; at `R = 0`, `R_z = -m2 < 0`; at `R = -3/2`, `R_z >= 9/4 - m2 > 0`. So R is trapped in (-3/2, 0): h has no node and
`h'(y_b) != 0`. For `m2 < 0`: `s_+ > 0` and `R_z = -m2 > 0` at `R = 0`, so `R > 0`. For `m2 = 0`, `R = 0`. Same conclusion without SUSY and with a weaker hypothesis on U.

**Numerical spot check [N]** (`tensor_numeric_spotcheck.py`, numpy; float background reproduces `phi_b = -5.2157e-5`, `rho_b = 78.82906`, `B = 5.320e-4`, junction
residuals < 1e-7): `R(z_b)` is negative for all 9 sampled `m2 in (0, 9/4)` (approx `-m2/(3 rho'_b)`), positive for `m2 < 0`; finite-volume spectrum of
`-(rho^4 h')' = m2 rho^2 h` with Neumann at the shell: eigenvalues `{<1e-9, 2.365, 2.702, ...}` for y_min = 1e-4 and `{<1e-9, 2.453, ...}` for y_min = 1e-3,
i.e. the zero mode and a box-quantised continuum `9/4 + (pi/L_z)^2` that approaches 9/4 from above. No eigenvalue in (0, 9/4). (My first attempt on a uniform
z-grid was under-resolved - the wall next to the shell has z-thickness ~0.013 - and gave a spurious negative eigenvalue; the y-grid version resolves it. Reported for honesty.)

Not needed for stability but stated in the claim: `m2 >= 9/4` contains no L^2 solutions (oscillatory at the cone), so "only the massless graviton" is right. Any `m2 >= 0`
would be stable anyway.

## 3. Vector sector

* Decomposition `delta g_{y mu} = Bv V_mu`, `delta g_{mu nu} = 2 rho^2 Fv nabla_(mu V_nu)`, `delta g_yy = delta phi = 0` is complete for transverse V (no scalar can be built
  linearly from V). `Fv = 0` is always reachable (`zeta = -Fv`), and since `v^y = 0` the shell is not moved. [P]
* **G2 (gap closed) - the exclusion should not rest on cone regularity.** In the Lorentzian section y = 0 is a horizon; normalisable tensor modes themselves diverge there
  like `y^{s_+}`, so "singular at the cone" is not by itself a criterion. `vector_special_extra_checks.py` proves in a *general* gauge: the bulk equations involve only the gauge
  invariant `sigma_V = Bv - rho^2 Fv'` (X1a-c, including the identity `(Box-2)[nabla_(mu V_nu)] = nabla_(mu [(Box+3)V]_nu)`), and the perturbed Israel condition is

      delta( K_{mu nu} - (sigma_t/6) g_{mu nu} ) = - sigma_V nabla_(mu V_nu)      at the shell   (X2).   [P]

  Hence `sigma_V(y_b) = 0`, so `C = 0` in `sigma_V = C/rho^2` for every non-Killing V. The cone argument is then a second, independent reason
  (`delta g_{y-hat mu-hat} ~ C/y^3`, more singular than any normalisable mode).
* Killing V: `nabla_(mu V_nu) = 0`, both bulk equations and the junction are empty, `Bv` is arbitrary and equals `rho^2 zeta'`: pure gauge (a y-dependent isometry). Correct.
* Non-Killing solutions of `(Box+3)V = 0` do exist on Lorentzian dS_4 (k != 0), so the case is not vacuous. For them `nabla_(mu V_nu)` is TT with `m2 = 0` (X3): the would-be
  vector mode is exactly the second, non-normalisable `m2 = 0` tensor solution `h' ~ 1/rho^4`, excluded in sec. 2 by `h'(y_b) = 0`. The two sector analyses agree.
* Consistency with KK gravitons: their helicity-1 components live inside the TT tensors of the 4D-covariant decomposition and are covered by G1. "No vector modes" agrees with
  Frolov-Kofman, arXiv:hep-th/0209133 (abstract: massless scalar and vector projections of bulk gravitons are absent; abstract checked online today).

## 4. Special harmonics (Hess Y = -gamma Y)

* **G3 (gap closed) - "every bulk solution is pure gauge" without gauge fixing or solution counting.** The script's route (gauge `dy dx = 0`, `xi + 2psi = 0`, then a 2-dim
  solution space of the master equation) is reachable and correct (the gauge ODE is always solvable on (0,y_b]; injectivity: X5 gives `psi_g = 0 => rho^2 phi'^2 b/(3rho'^2) = 0`),
  but it leans on the Chat 11 master equation at a value where `xi = -2psi` is a choice, not a consequence. Direct proof (X4a/b): for arbitrary `(xi, Bs, psi, chi)` put
  `a = chi/phi'`, `b = Ha - psi`; subtracting the Lie derivative leaves a remainder `(xi_r, Bs_r, 0, 0)`, and the linearised Einstein equations for it are
  `xi_r = Bs_r/(rho rho')`, `c1 Bs_r = 0`, `c2 Bs_r = 0` with `2c2 - c1 = 3 rho^2 phi'^2`. So `Bs_r = xi_r = 0` wherever `phi' != 0`, `rho' != 0` (both certified on (0,y_b]). [P]
  Regular vs singular branch is then irrelevant for the junction: whatever the generating vector, the configuration is the unperturbed bulk with the shell at
  `y_b + eps(zeta - a(y_b))Y`. (The singular branch `psi ~ y^{-2}` is in addition non-normalisable.) "Both gauge modes needed" is covered by the same statement.
* Displaced shell [P; re-derived by hand]: `delta K^mu_nu = zeta Y (H' + 1/rho^2) delta = -(phi'^2/3) zeta Y delta = (sigma_t' phi'/6) zeta Y delta`; scalar junction mismatch
  `zeta Y phi'(phi''/phi' + sigma_t''/2)`. Normal tilt enters `n.dphi` only at `O(eps^2)` because `d_mu phi0 = 0`. Correct. So `B != 0` and `phi'_b != 0` give no mode, for either sign of B.
* **G4 (clarification, important for wording).** The partition is by `Hess_TF Y = 0`, *not* by `mu2 = -4`. On Lorentzian dS_4 there are infinitely many harmonics with
  `Box Y = -4Y` and `Hess_TF Y != 0` (k != 0); they belong to the generic scalar sector. I checked that the Chat 10 certificate does cover them: it proves `b = lam w + B R < 0`
  for every real `mu2 < 9/4` by monotone comparison, which includes `lam = 0`. X6 shows the generic mismatch on the gauge profile at `mu2 = -4` is `phi' B a(y_b)`: the same
  obstruction as S3b. For these Y, `Hess_TF Y` is a TT tensor with `m2 = 2`, and the residual-gauge ODE `rho^2 b'' + 4 rho rho' b' + 2b = 0` *is* the tensor equation at `m2 = 2`:
  the scalar/tensor split is not unique there (nor at `m2 = 0`, sec. 3), but each overlapping configuration is excluded by both analyses, consistently.
* Missing item (minor): **Y = const** (mu2 = 0, `grad Y = 0`). Neither Codazzi nor `xi = -2psi` is forced, so no sector covers it. It is a static, dS-invariant deformation, i.e. a
  neighbouring background solution; it exists iff the background root is degenerate. It is not an instability and not normalisable on dS_4, but it is the linearised form of the
  "uniqueness of the root" question, which remains **[open]** (non-singularity of the Poincare-Miranda Jacobian would settle local uniqueness; I did not check it).

## 5. Does the combined claim follow?

Yes, as: *"For S_8/5 no sector of the 4D-covariant decomposition contains a normalisable mode with m2 < 9/4 other than the massless graviton: scalar (Chat 10 certificate, all real
mu2 < 9/4, including non-special mu2 = -4 harmonics), special harmonics (none, B != 0), vector (none; junction forces sigma_V = 0), tensor (zero mode only)."*

Still not proved:
1. Hand lemmas of Chat 10 (cone barrier, monotone comparison, counting theorem, Frobenius analyticity H0) - not machine-checked. The tensor/vector/special results use only
   `phi' > 0`, `H > 0`, `phi in (-1, 0)` (so U <= -1/6) and `B > 0` from the existence certificate, plus the exact Sturm counts T4.
2. Completeness / self-adjointness: passing from "no unstable mode" to "every finite-energy perturbation stays bounded" needs the spectral theorem for each radial operator
   (tensor: standard, limit point at the horizon, regular Robin endpoint [C]; scalar: eigenparameter in the boundary condition, cited in Chat 10) and decay estimates for the
   continuum. Mode stability is not linear stability.
3. Horizon regularity: y = 0 is a horizon of the Lorentzian geometry. "Normalisable in z" is the standard criterion (Garriga-Sasaki arXiv:hep-th/9912118, Langlois-Maartens-Wands
   arXiv:hep-th/0006007, Frolov-Kofman arXiv:hep-th/0209133; abstracts checked today: mass gap 3H/2 above the zero mode) but its equivalence with regularity on the future
   horizon for growing modes was not examined here.
4. Starting points of the junction analysis (first-order expansion of `K = sigma_t/6`, `n.dphi = -sigma_t'/2`; Z2 symmetry imposed on perturbations - Z2-odd perturbations
   are excluded by assumption, not by proof).
5. Uniqueness of the background root / the Y = const deformation; non-linear stability; quantum (tunnelling) stability.

## 6. The Chat 9 registered shell (phi_b = +2.28e-5, B = -2.68e-4)

* Tensor: identical argument; needs `U <= 0` on `[phi_h, phi_b]`, supplied by T4a (`U < 0` on [-1, 3/20]); float: `U_max = -0.1666`, `R(z_b) < 0` on (0, 9/4), spectrum
  `{~0, 2.365, ...}` [N]. The Riccati proof needs only `rho' >= 1`. Conditional on the archived M462 enclosure and `phi' > 0` (H1 of the oscillation theorem), not on anything new.
* Vector: background-independent; none.
* Special harmonics: `B = -2.68e-4 != 0`, so none. The instability of the registered shell is therefore confined to the generic scalar sector (`mu2 = -7.72`); the five special
  harmonics do not add a mode. (B is float here, not interval-certified for this shell, but its sign is what Chat 9/10 already rely on.)

## 7. Prior art

The tensor result (zero mode + gap 9H^2/4, via supersymmetric quantum mechanics) and the absence of vector modes are standard for de Sitter branes (references above; Gen-Sasaki
arXiv:gr-qc/0011078 for the negative radion mass-squared - abstract checked, the coefficient -4H^2 is from memory, unverified). Nothing here is novel physics; the value is that it
is verified exactly for this model and shell. Frolov-Kofman "Can inflating braneworlds be stabilized?" is relevant to the scalar sector; I could not confirm its arXiv number today.

## Files

`negative_controls.py/.log`, `NEGATIVE_CONTROLS_RESULT.json`; `general_tensor_check.py`, `GENERAL_TENSOR_RESULT.json` (4 s); `vector_special_extra_checks.py`,
`VECTOR_SPECIAL_EXTRA_RESULT.json`; `tensor_numeric_spotcheck.py`, `TENSOR_NUMERIC_SPOTCHECK.json` (system python3 + numpy, ~1 min); `rerun/` (copies + logs); `negctl/` (mutants).
SymPy scripts: run from this folder with the venv python (they import the local copy of `lin_gr.py`).
