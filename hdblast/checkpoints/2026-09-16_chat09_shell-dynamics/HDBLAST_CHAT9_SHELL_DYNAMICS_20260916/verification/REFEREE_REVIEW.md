# Referee review - HDBLAST Chat 9 "shell dynamics" package

16 September 2026. Adversarial physics-and-claims review of the seven track reports (independent-derivation,
verify-math, literature, tensor-calibration, blast-dynamics, slope-analytic, certified-instability) and of the
orchestrator results R1-R5. Written as a skeptical referee for a Phys. Rev. D submission.
Everything I ran myself is floating point (numpy), not certified. New files from this review:
`verification/referee/referee_stable_shell_test.py`, `verification/referee/REFEREE_STABLE_SHELL_TEST.json`.

## 0. Verdict in one paragraph

The linear-stability result is solid and unusually well cross-checked: the registered static Z2 dS shell
(t = 1e-3, c = c_star, LINEAR detuning) has exactly one scalar bound state, mu^2 = m^2/H^2 = -7.71787, with
positive norm (a tachyon, not a ghost), a healthy tensor sector, and an O(t)-accurate two-derivative effective
description. It is a model-specific instance of the known Frolov-Kofman instability of de Sitter branes with a bulk
scalar, in a "soft brane potential" regime (mu^2 < -4) that FK's quantitative analysis did not cover.
What is NOT supported: that this roll-off is a viable "blast"/hot-big-bang mechanism. Inside the project's own EFT the
roll toward phi_b -> +1 ends in another eternally inflating dS brane (no exit, no reheating, produced-particle energy
fraction <= 1e-8 to 1e-5), the roll toward the throat heads to negative V_E / brane-frame collapse with unknown 5D fate,
only ~3 e-folds occur, and the slow-roll variants give n_s <= 0.93 unless a slope is tuned to ~1e-5 and still never
end. Moreover (new referee test, section 2.1) the tachyon is a property of the chosen linear detuning: adding a
brane-tension curvature of relative size ~1e-3 (d >= 1.3) removes it in the full 5D problem. Recommended decision
for a paper: publishable as a careful stability/EFT analysis of a specific thick-wall dS-brane model with honest negative
cosmological conclusions; not publishable with any "origin of the Big Bang" framing.

## 1. Cross-track consistency

### 1.1 Numbers that agree (checked by me from the reports / JSON files, some re-run)

| quantity | tracks | values | status |
|---|---|---|---|
| mu^2 at t=1e-3 | orchestrator, GN, verify-math, literature(FK eqs), slope, my re-run, certificate | -7.7178716187 / -7.7178716252 / -7.717871623 / -7.7178716231 / -7.7178716 / -7.7178716186 (n=6000) / certified bracket (-7.71788,-7.71786) | consistent to < 1e-8; four different formulations (longitudinal master eq., Gaussian-normal 4-variable flow, first-order FK system, Riccati) |
| t -> 0 intercept | HJ closed form, MSA (Brax et al.), GN t-scan, verify-math fit | -7.719795918 (all) | consistent; I re-derived -4(3c^2-4c+8)/(c(3c+4)) numerically |
| slope | R4 1.9243, verify-math 1.9244, slope-analytic 1.92438964 (analytic) vs 1.92438965 (5D) | consistent | the certified bracket excludes the t->0 value at t=1e-3, as it must |
| kinetic normalisation | GN 5D norm Z = 0.896383 -> 0.896392; HJ f Z_E = c I_+ (1+3c/4) = 0.896392; verify-math pole residue -6932.17 vs -rho_b^2/(f Z_E) = -6932.10 | consistent (I checked 6213.88/0.896392 = 6932.1) | three independent determinations of the sign AND magnitude of the kinetic term |
| FK quadratic form | GN: Q/(lam N) = 1 + 1.2e-6; literature (corrected file on disk): 1.000000013; B rho_b^2 = -1.66646 both; bulk/boundary = -8.2e-5 both | consistent | the 0.737 value quoted by the GN track refers to a superseded first attempt; the file on disk is the corrected one. Disagreement is RESOLVED. |
| B = phi''/phi' + sigma''/2 at shell | GN -2.68e-4, verify-math -2.681e-4, literature -2.6818e-4, certificate [-2.681831e-4,-2.681824e-4] | consistent | |
| growth | s = -3/2 + sqrt(9/4+7.71787) = 1.65719; e-fold time 0.6034/H; m^2 = -1.24204e-3 | consistent in all tracks | |
| Planck mass | tensor: M_4^2 = 2 I_kept = 2.06773, 8 pi G = 0.483621; EFT f = 2I; instanton -8 pi^2 M_4^2/H^2 = -24 pi^2 f^2/V | mutually consistent conventions | |
| landscape | hilltop 0.23303 t, throat value 0.12420 t = (1-c)(5/9)^2 t, Theta(phi_b=-1) = 0.54791 | orchestrator, literature(MSA), tensor, blast | consistent |
| calibration | H L_0 = sqrt(h) = 0.0126858; case A M_5 = 1.31e17 GeV (I re-derived from M_Pl^2 = 2 I_kept M_5^3 L_0); rescued-inflation A_s: 2H^2/(pi^2 M_Pl^2 r) = 2.3e-9 with H = 1.2e11 GeV, r = 2.1e-7 | self-consistent | |

Exact-identity script `exact_checks/verify_exact_identities.py`: re-run, ALL PASS (0.2 s).

### 1.2 Inconsistencies and defects found

1. **`REPLAY.py` on disk does not run (SyntaxError).** In stage 5 an inline `#` comment swallows the closing
   arguments of the `gate(...)` call (`... < 1e-4   # exact identity; ... , dict(...))`). `ast.parse` fails at line 46.
   `REPLAY_RESULT.json` (all_pass: true) therefore was not produced by the file as it now stands (same-minute
   timestamps: the comment was evidently added after the run). Must be fixed and re-run before the package is cited
   as "replayable". (I did not modify it - outside my write scope.)
2. **n_s claim in `effective_theory_orchestrator/HJ_EFFECTIVE_THEORY_DERIVATION.md` line 83-84 is wrong twice.**
   (a) "mu^2 = -0.0525 i.e. n_s ~ 0.965 at d = 1.105924" uses n_s = 1 + 2 eta and ignores V'''/V = +7.7; the tensor
   track's exact-background integration gives n_s = 0.910/0.918/0.924 (N = 50/55/60), bounded by 1 - 4/N.
   (b) mu^2 = -0.0525 is a t -> 0 statement; at the registered t = 1e-3 the 5D value at that d is -0.08988, and no
   purely quadratic detuning gives a hilltop with |mu^2| < 0.0731 (floor ~ 2.31 sqrt t; one 5D point supports it).
   R5's "slow roll needs d within ~1% of d0" must be replaced by the tensor/slope tracks' statements.
3. **"Field space ENDS at Theta = +0.548" (R5) is contradicted by the blast track** (Z(phi_b=-1) = 9/23 finite, EFT
   continues onto the -coth branch). I re-derived Z(-1) = W I'' with I'' = 54/230, i.e. 9/23: correct. Both the
   tensor track and the literature track repeated the "ends" wording; only the blast track tested it. Conversely the
   "+1" side, described in R5 as a "long gentle slope", also ends at finite distance (Theta = -3.4653): I checked the
   mechanism (homogeneous solution I - 9/2 ~ (1-phi)^{1/9}, so Z ~ (1-phi)^{-17/9}, distance ~ (1-phi)^{1/18},
   V_E - V_end ~ DeltaTheta^18, f = 9 - O(DeltaTheta^2)): consistent with the blast track's numbers. Physically the
   end point is the brane receding forever into AdS_+ (infinite y, finite Theta), not a wall of field space.
4. **V_E(phi_b -> +1) = 0.01976 t (R5) should read (1+c)/81 t = 0.019723 t** (two tracks agree; R5's value is a
   finite-distance table entry).
5. **Instanton identity numerics differ between files**: `SHELL_INSTANTON_ENTROPY_IDENTITY.md` quotes ratio
   1.0000006 at t=1e-3, `REPLAY_RESULT.json` 1.0000198 (different quadrature, 1e-4 tolerance). Harmless, but the
   identity is exact and the replay gate is loose because of a 1e11 cancellation; say so in the gate text.
6. **Negative-mode statements conflict.** The literature report (sec. 6) calls the single homogeneous tachyon "the
   natural candidate" for the Coleman negative mode; the orchestrator's instanton note correctly observes that
   mu^2 = -7.72 < -4 makes l(l+3)+mu^2 negative for l = 0 AND l = 1 (1+5 modes). The second statement is the right
   one at EFT level: the static shell is the analogue of a Hawking-Moss saddle with |V''| > 4H^2, which is known NOT to
   be a single-negative-mode bounce (a CDL-type O(4) solution is then expected instead). Any "creation from nothing of
   the shell universe" reading of the entropy identity is therefore unsupported. The 5D l = 1 count is genuinely open
   because l = 1 are the conformal-Killing harmonics where the longitudinal-gauge argument degenerates.
7. **Minor code/documentation defects already flagged by tracks and confirmed as still present**:
   `complex_plane_winding.py` normalises by psi_b without checking its winding (verify-math checked: 0);
   `T_SCAN_ORCHESTRATOR.json` loses accuracy for t <= 1e-4 (4e-7); `eft_landscape.py` comment has the sign of Theta
   reversed; `gn_lib` produces a spurious root near 2.177 on a uniform grid; `independent_spectrum.py` runs 8.5 min;
   stale `run_log*.txt` in `certified_instability`; `log_slowroll_first.txt` is a transcription, not a run log;
   `VALIDATE_cstar_hi.json` and `TENSOR_AMPLITUDE_TSCAN.json` do not exist.
8. **Citations.** Only the literature track (and partly verify-math and tensor) verified references online. The
   GN, blast, slope and certificate tracks cite from memory. Two arXiv numbers in the task list were wrong and were
   corrected by the literature track (BraneCode is hep-th/0309001; gr-qc/0212114 is Himemoto-Tanaka). Kanno-Soda is
   cited as hep-th/0303203 by one track and hep-th/0207029 by another (these are different papers of the same
   programme; check which is meant). No reference may go into a manuscript from the "from memory" lists without being
   opened.

### 1.3 Independence audit (how many genuinely independent confirmations?)

* Linearised equations: TWO independent derivations (orchestrator longitudinal gauge; GN Gaussian-normal ADM split)
  plus ONE brute-force numerical linearisation of the full 5D Einstein-scalar system with sabotage tests (verify-math).
  The literature track's FK system and the certificate's Riccati system are algebraic rewritings of R1 - not
  independent of it. The certificate is conditional on R1 + the archived M462 enclosures.
* Eigenvalue numerics: effectively three independent codes (orchestrator/slope share `hdblast_background.py`;
  verify-math and GN have their own background solvers).
* Norm sign: three routes (GN on-shell-action reductions; Sturm-Liouville identity lam*N = Q > 0; source-response
  residue matched to the EFT). None is a full second-order 5D action computed from scratch with all boundary terms;
  the cone boundary term is argued, not proved. I regard the residue test (sign and magnitude to 1e-5) as decisive at
  referee level, but the manuscript should say "established by ... ; a complete quadratic action with Gibbons-Hawking
  terms was not derived".
* Uniqueness of the bound state: real-axis scans (to -400), argument principle on finite rectangles, and the
  analytic "at most one mode below -4 + real spectrum" argument. Nothing excludes bound states with
  mu^2 < -400 numerically; the analytic argument does (given completeness of the FK system). Acceptable.

## 2. Physics assessment

### 2.1 Is "the static registered dS shell is unstable within about one Hubble time" justified?  YES, with scope limits.

* It is a linear statement about an SO(4,1)-invariant background with data on the complete t=0 Cauchy slice (the
  4-ball 0 <= y <= y_b doubled). The growing mode behaves as (y e^{tau})^{s} near the cone with the SAME exponent
  s = 1.657 in space and time, i.e. it is a function of the regular null coordinate, finite on the future cone but only
  C^1 there. This is the standard situation for dS-brane/bubble-wall bound states (Garriga-Vilenkin, Garriga-Sasaki),
  but the manuscript should state the regularity class explicitly and cite the project's cone-transparency note.
* Time scale: e-folding time 0.60/H; with a quantum seed the EFT gives 3.4 e-folds (k=0) / 2.6 (closed bounce)
  before roll-off, growing only logarithmically as t decreases (9.6 e-folds at t = 1e-12). "About one Hubble time" is
  the e-folding time, "a few Hubble times" is the lifetime; use the second phrase.
* **New referee test: the instability is a property of the registered LINEAR detuning, not of dS shells in this
  model.** Full 5D solves (machinery of `slope_analytic/validate_5d.py`, t = 1e-3, c = c_star, scan -60 < mu^2 < 2.2):

  | d (sigma_t = 2W + t(1 + c phi + d phi^2/2)) | 5D bound states | EFT prediction mu0 + s1 t |
  |---|---|---|
  | 0 (registered) | -7.7178716 | -7.7178715 |
  | 1.0 | -0.7895812 | -0.7895879 |
  | 1.3 | **+1.2949033** (no tachyon) | +1.2949052 |
  | 1.4 | **+1.9871331** (no tachyon) | +1.9871335 |
  | 1.6 | none below 9/4 (mode has merged with the continuum; no tachyon) | +3.3724 |

  So a brane-tension curvature sigma'' of absolute size ~1.3e-3 (in units where sigma ~ 2) already gives a linearly
  stable static dS shell with a light (0 < m^2 < 9H^2/4) or absent modulus. This agrees with Frolov-Kofman's message
  (stabilisation needs enough brane-potential curvature) and with the EFT formula mu^2 = 3(d-d0)/Z_E, now checked in 5D
  on the stable side to ~2e-6. Consequences: (i) "the dS shell is unstable" must always carry the qualifier "for the
  registered linear detuning"; (ii) whether the universe "blasts" is an input (choice of sigma_t at O(t)), not an output;
  (iii) the same family supplies a STABLE dS shell, which is arguably the more useful object for the project and
  should be registered and certified (gate G3).
* Sectors not covered: the vector sector is not discussed in any report (expected to be empty for Z2 + scalar bulk,
  but it has to be stated and referenced); the five mu^2 = -4, l = 1 harmonics were treated analytically only (pure
  gauge / translations) - accepted; the scalar continuum above 9/4 was not analysed for resonances (tensor continuum
  was).

### 2.2 Is "the roll-off is a concrete dynamical candidate for the blast" justified?  Only in a weak sense.

What is established (within the two-derivative EFT, leading order in t):
* + side: fast roll, max eps_H = 0.44 (acceleration never stops), all potential energy lost to Hubble friction, end
  state = RS2-like dS brane receding into AdS_+ with H_J = 0.606 H_top (I checked H_J^2 = (1+c)t/27 from the
  tension excess over 2W(1) = 2/3). No radiation, no reheating, produced spectator energy fraction <= 2.4e-8 (minimal)
  to 5e-5 (g = 30H). This is a transition between two de Sitter states, not a bang.
* throat side: EFT continues through phi_b = -1, V_E changes sign at phi_b = -1/c, brane-frame turnaround and
  collapse; this relies on the cubic W being unbounded below and on a branch (phi = -coth y) outside the registered
  solution's field range.
* e-folds: 3-10 in total; horizon/flatness NOT solved; closed shells with a V = 0 minimum recollapse after 4-6.5/H.
* slow-roll variants: n_s <= 0.93 for c = c_star (excluded at > 8 sigma by Planck 0.9649 +- 0.0042; the more recent
  ACT DR6 combination, n_s ~ 0.974 +- 0.003 - quoted from memory, to be verified - makes it worse); the rescued
  inflection variant needs c tuned to ~2e-5 (comparable to the uncomputed O(t) 5D corrections), has r ~ 2e-7 and
  running -2.6e-3 (both allowed) and still has no exit.

What is NOT known (and must be said in any write-up):
1. **5D nonlinear fate.** No 5D time-dependent solution exists in the project. Once the shell moves, SO(4,1) is
   broken, the bulk becomes dynamical, and the moving-shell theorem already shows the frozen-bulk ansatz is
   over-determined. The only published nonlinear study (BraneCode, two branes, quadratic potentials) finds
   "relax to lower H or collide/singularity"; the EFT + side is qualitatively the first option. The throat side could
   end in a bulk singularity, horizon formation, or brane-cone collision; nothing here decides.
2. **Non-local / KK sector.** At LINEAR order the 5D analysis is exact and includes the continuum; the O(t) agreement
   (slope 1.9244 derived analytically and confirmed) shows that the local EFT's error on the hilltop mass is controlled
   by t ~ (H l)^2, NOT by |mu^2| ~ O(1) or by the gap 3H/2 being comparable to the growth rate. That is a genuine and
   useful result. But the gap being ~H means there is no hierarchy protecting the NONLINEAR roll: energy transfer to
   the continuum (dark-radiation term; the slope track's integration constant z2 is exactly the bulk mass term and
   was fixed by cone regularity of the STATIC solution) is outside the two-derivative EFT and is O(1) uncertain in a
   fast roll with Theta'^2 ~ H^2.
3. **Validity of the two-derivative truncation during fast roll.** Tested only at quadratic order around the
   hilltop. Four-derivative HJ coefficients were not computed; their regularity at W_phi = 0 (phi_b = -1) is argued,
   not shown. The blast track's claim that the leaf-rapidity blow-up near the throat is a chart artefact is
   plausible (AdS boost isometry) but unproved; both readings must be kept open.
4. **Initial conditions.** The static shell is a measure-zero unstable configuration; nothing selects it. The
   instanton cannot be invoked for that purpose (section 1.2 item 6). The "H/2pi" seed for a tachyon with
   |m^2| = 7.7 H^2 > 9H^2/4 is heuristic (no dS-invariant state exists for it).
5. **Inhomogeneity.** With |m| = 2.8 H all harmonics grow at the same asymptotic rate; different Hubble patches roll to
   different sides. The homogeneous trajectories are therefore not the expected outcome; walls between (i)- and
   (ii)-patches, and their 5D meaning, are unexplored.

Conclusion: "candidate" is tolerable only as "the instability identifies the direction(s) in which the registered
static configuration would evolve; in the local EFT neither direction produces a radiation-dominated universe".

### 2.3 Is the comparison with Frolov-Kofman fair?

Mostly, thanks to the literature track, with three corrections of emphasis:
* R1 is FK's system equation-for-equation in another chart; R3 is dBVV / Brax-van de Bruck-Davis-Rhodes /
  Kanno-Soda. The method is prior art and must be presented as such.
* R3's wording "the -4 is the Frolov-Kofman universal shift" is an identification by value. FK's quantitative
  results use the rigid limit delta phi|_brane = 0, which forces m^2 > -4H^2; the registered mode lies BELOW -4 because
  the brane potential is soft (B < 0), a regime FK's Rayleigh bound does not cover. Fair wording: "consistent with FK's
  general conclusion that inflating branes with a bulk scalar are tachyonic unless the brane potential is stiff
  enough; our example sits in the opposite, soft, limit and violates the rigid-limit bound m^2 > -4H^2, as it may".
  Section 2.1's table is the constructive counterpart: stiffening sigma_t by O(t) stabilises, as FK anticipate.
* FK equation numbers were read through a summarising fetch of ar5iv, one fetch hallucinated and was discarded.
  A human must open hep-th/0309002 and hep-th/0209158 and confirm (5), (8a,b), (14), (16), (20), (26) and BBDR (25)
  before they are quoted. Likewise hep-th/0312087 (Du-Wang-Abdalla-Su) and the Kobayashi-Tanaka /
  Minamitsuji-Himemoto-Sasaki papers must be read in full before any statement that the single-brane thick-wall
  scalar spectrum is new. Banerjee-SenGupta 1705.05015 claims dS-brane self-stabilisation and should be addressed.

### 2.4 Calibration and bounds

* Table arithmetic is self-consistent (checked case A). The calibration is a unit conversion, not a prediction: H
  is free. The statement "cases A-C satisfy the sub-mm bound by 17 orders of magnitude" is true but empty for the
  same reason; the "sub-mm limit" row applies a present-day laboratory bound to the unstable t = 1e-3 shell and is only
  illustrative (the track says so). The relevant length for the 1/r^3 correction of this asymmetric thick wall is an
  unknown O(1); I would drop the 38.6 micron number (it is a Yukawa |alpha| = 1 bound) or label it clearly.
* r < 0.036 is irrelevant for the registered shell (no observable inflationary stage) and trivially satisfied by the
  rescued variant (r ~ 2e-7).
* F^2 - 1 = 0.18 % tensor enhancement: fine, and the 10 % deficit relative to Langlois-Maartens-Wands is a thick-wall
  effect with a conjectured log coefficient (one t-point breaks the trend; precision-limited).
* "A_s fixes t/(M_5 L_0)^3 = 2.9e-14, not t": correct; R5's "scale t from A_s" is ill-posed as stated.
* Reheating numbers (Gamma, T_rh) refer to modified, non-registered potentials whose hilltops are even more
  tachyonic (mu^2 = -30 to -207) and carry unknown O(1)-O(8 pi) prefactors; they must not be quoted as model results.

## 3. Claims

### 3.1 Safe to publish (exact wording matters)

S1. "In the 5D Einstein-scalar model with W = 1 - phi + phi^3/3 and a Z2 shell of tension 2W + t(1 + c_star phi),
    t = 1e-3, the SO(4,1)-symmetric static de Sitter shell has a normalisable scalar perturbation with
    -7.71788 < m^2/H^2 < -7.71786 (interval-arithmetic certificate, conditional on the archived M462 background
    enclosure and on the linearised equations, which were independently derived in two gauges and checked by numerical
    linearisation of the full field equations to 1e-9)."
S2. "Floating-point evidence (three independent codes, argument-principle count on Re mu^2 in [-400,2], |Im| <= 300,
    plus an analytic argument giving real spectrum and at most one mode below -4) indicates it is the only scalar bound
    state below the continuum threshold 9H^2/4."
S3. "The mode has positive Klein-Gordon norm (it is a tachyon, not a ghost): lam*N_FK = Q_FK > 0 identically, and the 5D
    source-response residue equals the effective-theory value -rho_b^2/(f Z_E) to 1e-5. A complete second-order 5D
    action including all boundary terms was not derived."
S4. "The static shell is therefore linearly unstable with e-folding time 0.60 H^-1."
S5. "The tensor sector contains one massless graviton (M_4^2 = 2 I_kept = 2.0677), a continuum above 9H^2/4 and no
    other bound state or tachyon (Q^+Q factorisation; numerically monotone shooting function)."
S6. "A two-derivative Hamilton-Jacobi effective action (f = 2I, Z = -W f'/W_phi, V = sigma_t - 2W) reproduces the 5D
    mass: mu^2 = -4(3c^2-4c+8)/(c(3c+4)) + 1.92439 t + O(t^2), the O(t) coefficient being derived analytically (two
    quadratures) and confirmed by 5D numerics to ~1e-6 in five detuning configurations; c_star is the stationarity
    condition of V_E at phi_b = 0." The method is that of de Boer-Verlinde-Verlinde and the moduli-space
    approximation of Brax et al.; only the model-specific results are new to the project.
S7. "The instability is specific to the soft (linear) detuning: adding t d phi^2/2 with d >= 1.3 yields, in the full
    5D linear problem at t = 1e-3, a static dS shell with no tachyon (floating point, real-axis scan -60 < mu^2 < 2.2)."
    [new, this review; needs an independent re-run and a wider scan before publication]
S8. "Within the two-derivative effective theory the registered shell rolls off in a few Hubble times; toward
    phi_b -> +1 it relaxes, without ending acceleration and without reheating, to a de Sitter brane in AdS_+ with
    H = 0.606 H_top; toward the AdS_- throat the effective potential becomes negative and the 5D fate is unknown."
S9. "Hilltop slow-roll variants within the quadratic detuning family give n_s <= 1 - 4/N (<= 0.93), are excluded by
    Planck, and have no graceful exit; an inflection variant can match n_s only with a ~1e-5 tuning of the slope and
    still does not end."
S10. "The Euclidean action of the O(5) continuation equals minus the 4D de Sitter entropy computed with the zero-mode
    Planck mass (exact, elementary identity; prior art Garriga-Sasaki, Hawking-Maldacena-Strominger)."
S11. "These results are an instance of the Frolov-Kofman instability of inflating branes, here in the soft-potential
    regime mu^2 < -4 outside FK's rigid-limit window; we found no prior computation for a single Z2 shell with a regular
    cone on a backreacted BPS wall in a bounded literature search (~30 papers)."

### 3.2 Overclaims - do not make these statements

O1. Any statement that this "explains", "models" or "is the origin of" the Big Bang, or that the roll-off "is the
    blast" in an observational sense. The EFT end state on the + side is de Sitter again; there is no radiation era,
    no reheating, ~3 e-folds, and no viable spectrum.
O2. "de Sitter shells in the model are unstable" without "for the registered linear detuning" (refuted by S7).
O3. "The instability is proved." Only the existence of the eigenvalue is certified, conditionally; uniqueness and
    norm sign are numerical/analytic-with-assumptions; the background enclosure is inherited.
O4. "First computation of ...", "new instability", "novel mechanism". The instability type (FK 2003), the EFT (dBVV
    1999, BBDR 2002, Kanno-Soda), the -4 (Gen-Sasaki, Garriga-Vilenkin), the 9/4 threshold, the SUSY-QM tensor
    factorisation, the entropy identity and the hilltop n_s bound are all prior art. Novelty of N_phys = 12(mu^2+4)N_FK
    and of "at most one mode below -4" is unknown; present them as remarks.
O5. "The -4 is the Frolov-Kofman universal shift" as a derivational claim; "mu^2 = m_0^2 - 4" with FK's m_0^2 > 0 is
    false here.
O6. "n_s ~ 0.965 at d = 1.105924" and "slow roll needs d within 1% of d0" (both refuted inside the project).
O7. "The field space ends at Theta = 0.548" and "long gentle slope toward +1" (both superseded).
O8. "The shell universe is created from nothing with probability exp(+-S_dS)" - the saddle has 1+5 negative modes at
    EFT level; not a bounce.
O9. "The two-derivative EFT is valid during the roll" - validated only to quadratic order about the hilltop.
    Likewise "the rapidity divergence is a chart artefact" (argued only).
O10. "Particle production reheats the universe" - fractions <= 1e-8..1e-5 into a dS end state; and the Chat-8 sech^2
    ceiling is outside its validity domain here (H tau_p ~ 0.8), usable only as an order-of-magnitude comparator.
O11. Quoting the calibration table as predictions of M_5, L_0 or of sub-mm signals; quoting T_rh ~ 6e15 GeV (it belongs
    to a non-registered illustrative potential with uncomputed prefactors).
O12. "Fully certified" for anything other than S1. "Replayable" until REPLAY.py is repaired and re-run.

## 4. Recommended next gates (in order of decisiveness)

G1. **5D nonlinear evolution of the shell** (the only thing that can decide the "blast" question). Minimal version:
    SO(4)-symmetric (closed FRW on the brane) 1+1 characteristic code for the Einstein-scalar bulk with a Z2 boundary,
    initial data = registered static solution + tachyon eigenfunction, both signs. Deliverables: end state on the + side
    (compare H_J -> 0.606 H_top and the dark-radiation term), fate on the throat side (singularity / horizon / cone
    collision), energy radiated into the bulk. A cheaper precursor: second-order perturbation theory (cubic coupling of
    the tachyon to the continuum) to estimate the leak rate.
G2. **Four-derivative HJ coefficients** along the wall, in particular their behaviour at phi_b -> -1 and their size
    relative to the two-derivative terms on the fast-roll trajectory (Theta'^2 ~ H^2). Decides O9 and the
    chart-artefact question without a PDE code.
G3. **Register and certify a STABLE dS shell** (e.g. d = 1.4: 5D mu^2 = +1.98713, or d = 1.6: no bound state): interval
    certificate that the Riccati mismatch has no zero on mu^2 < 0 (needs a uniform-in-mu^2 argument: monotonicity of b
    plus one sign, or the Rayleigh identity with B > 0 giving lam > 0 analytically). This gives the project a
    perturbatively stable dS brane and turns S7 from a float result into a theorem. Also scan larger d and the
    scalar continuum for resonances.
G4. **Negative-mode count of the Euclidean O(5) saddle in 5D**, including the l = 1 conformal-Killing sector, and a
    search for the CDL-type O(4)-symmetric solution that should exist when |mu^2| > 4. Decides whether any
    creation/tunnelling interpretation survives.
G5. **Repair REPLAY.py**, re-run, and add to it: the GN norm, the residue test, the d = 1.4 stable-shell test, and the
    corrected fk_quadratic_form check. Fix the documentation errors listed in 1.2 (items 2-4, 7).
G6. **Human literature pass**: open hep-th/0309002, hep-th/0309001, hep-th/0209158, hep-th/0312087,
    hep-th/0311197, gr-qc/0303108, 1705.05015 and confirm every quoted equation number and claim; verify all
    "from memory" citations of the GN, blast, slope and certificate tracks; verify the ACT DR6 n_s value before use.
G7. **Close the small analytic gaps**: vanishing of the cone boundary term in the on-shell-action norm; explicit
    statement of the regularity class of the growing mode on the future cone; vector sector; make h a second
    Taylor-model variable to tighten the certificate; O(t^2) coefficient (-0.094) analytically if cheap.
G8. Only after G1-G2: inhomogeneous roll-off (patch structure, walls between + and throat domains) and any
    exit/reheating construction. Without an exit mechanism that survives G1 the cosmological programme for this
    shell is negative, and that should be reported as such.
