# Prior art for the stability and effective theory of the registered HDBLAST de Sitter shell

16 September 2026. Literature agent report (second pass; the first agent was interrupted, its scripts are reused).

**How the sources were read.** No PDF text tool is available on this machine, so every paper below was read through the
ar5iv HTML rendering of the arXiv source (or the arXiv abstract page) via a summarising fetch tool. Equation numbers quoted
here were returned by that tool; the Frolov-Kofman numbers (5), (8a), (8b), (14), (16), (17), (20), (22), (26) were obtained
consistently in three separate reads and are, in addition, *numerically* confirmed to reproduce the project's eigenvalue.
Other equation numbers were read once and carry a transcription risk. Nothing is cited that was not opened. Statements about papers are paraphrases, not verbatim quotations. Two arXiv
numbers in the task list needed correction:

* BraneCode is **hep-th/0309001** (not hep-ph/0309001, which is an unrelated pentaquark paper).
* gr-qc/0212114 is Himemoto-Tanaka, *"Braneworld reheating in the bulk inflaton model"* (not a stability paper). The
  backreaction papers of that group are gr-qc/0112027 and gr-qc/0303108 (below).

Status vocabulary: **[known]** = in the literature; **[reproduced]** = re-derived/re-computed here against the literature
form; **[specific]** = we found no source for it, and it depends on the registered model; **[not verified]**.

---

## 1. Frolov and Kofman, "Can inflating braneworlds be stabilized?", hep-th/0309002

**What is derived.** 5D Einstein + bulk scalar, two Z2 branes with potentials U_i(phi), de Sitter slicing. Scalar
perturbations in generalized longitudinal gauge, reduced to a self-adjoint Sturm-Liouville problem; Rayleigh bound on the
lowest eigenvalue; conclusion that the radion of inflating branes is typically tachyonic, m^2 = -4H^2 + m_0^2(H).

**Their equations (FK numbering).**

| FK | content |
|---|---|
| (1) | ds^2 = a(w)^2 [dw^2 - dt^2 + e^{2Ht} dx^2] (conformal bulk coordinate w) |
| (2) | S = M_5^3 int sqrt(-g){R - (grad phi)^2 - 2V(phi)} - 2 M_5^3 sum int sqrt(-q){[K] + U(phi)} |
| (3a-c) | phi'' + 3(a'/a)phi' - a^2 V' = 0;  a''/a = 2a'^2/a^2 - H^2 - phi'^2/3;  6(a'^2/a^2 - H^2) = phi'^2/2 - a^2 V |
| (4),(5) | ds^2 = a^2[(1+2 Phi)dw^2 + (1+2 Psi) ds_4^2],  Psi = -Phi/2 |
| (6),(7) | separation on dS_4 harmonics, Box_4 Q_m = m^2 Q_m |
| (8a) | (a^2 Phi)' = (2/3) a^2 phi' dphi |
| (8b) | ((a/phi') dphi)' = (1 - (3/2)(m^2+4H^2)/phi'^2) a Phi |
| (9) | u'' + (m^2 + 4H^2 - V_eff)u = 0, V_eff = z''/z + (2/3)phi'^2, z = ((2/3) a phi'^2)^(-1/2) |
| (13) | background junctions a'/a^2 = -/+ U/6,  phi'/a = +/- U'/2 |
| (14) | (dphi' - phi' Phi) = +/- (1/2) U'' a dphi at the branes (brane displacement xi set to 0: FK argue that, without matter on the brane, oscillatory displacement modes are not excited) |
| (16) | rigid limit U'' -> infinity: dphi = 0 |
| (17a,b) | with Y = a^2 Phi: -(g Y')' + f Y = lambda g Y, Y'(w_+/-) = 0, f = 1/a, g = z^2, lambda = m^2 + 4H^2 |
| (20) | m^2 <= -4H^2 + (2/3) int(dw/a) / int(dw/(a phi'^2))  (Rayleigh, trial function Y = 1, rigid b.c.) |
| (21) | late-time growth exp[(sqrt(9/4 + |m^2|/H^2) - 3/2) H t] |
| (22) | 65 e-folds of inflation need m^2 >~ -H^2/20 |
| (26) | m^2 = -4H^2 + m_0^2(H) |

**Dictionary to the registered chart.** FK's action is 2 M_5^3 times ours, so V_FK = U (our bulk potential),
U_FK = sigma_t (our tension), phi identical, and (13) is exactly our junction pair rho'/rho = sigma_t/6, phi' = -sigma_t'/2
for the side y < y_b. Further: a = rho, dw = dz = dy/rho, H = 1 (unit dS_4), m^2 = mu^2, Phi = xi, Psi = psi, dphi = chi.

**Which FK equations reproduce R1.**

* R1 "xi = -2 psi" **is** FK (5).
* R1 Codazzi relation chi = -3(psi' + 2(rho'/rho)psi)/phi' **is** FK (8a) (insert Phi = -2 psi, convert w -> y).
* R1 master equation (M_y) is FK (8a)+(8b) with chi eliminated (equivalently FK (9)/(17a) in another variable). The first
  literature agent integrated FK (8a),(8b) directly (script `prior_art_crosscheck.py`) and the orchestrator integrated
  (M_y); both give mu^2 = -7.7178716 to 1e-8 (`prior_art_crosscheck_wide_output.json`), so the equivalence is numerically
  confirmed, not just asserted. **[reproduced]**
* R1 shell condition chi' + 2 phi' psi + (sigma_t''/2) chi = 0 **is** FK (14) (lower sign, a d/dw = d/dy, Phi = -2 psi).
* R1 "no brane bending in longitudinal gauge" is FK's statement below (13)/(14) that xi can be set to zero.
* The continuum threshold 9H^2/4 and the normalisable cone exponent alpha = 1/2 + sqrt(9/4 - mu^2) are the standard dS-brane
  gap (Garriga-Sasaki hep-th/9912118 for tensors; Himemoto-Sasaki gr-qc/0010035, Kobayashi-Koyama-Soda hep-th/0009160 and
  Langlois-Sasaki hep-th/0302069 for a bulk scalar). In FK variables Y ~ y^p with p = 2 + alpha.

**What is *not* in FK and is different in the registered shell.**

1. FK have two branes; the registered solution has one Z2 shell and a regular cone/horizon. The inner boundary condition
   is the Himemoto-Sasaki-type regularity condition, not a second brane. This combination (FK system + regular cone) is
   also the set-up of Du-Wang-Abdalla-Su hep-th/0312087 (section 3), but only for a perturbative background.
2. FK's quantitative results (17b), (20) assume the **rigid** limit (16). The registered shell is in the opposite, **soft**
   regime: sigma_t'' = 2W'' = 4 phi_b ~ 1e-4. Then (14) is an *eigenvalue-dependent* boundary condition,
   `Y_y B = lambda Y / rho^2` at y_b, with `B = phi''/phi' + sigma_t''/2` (B = -2.68e-4 at t = 1e-3, B rho_b^2 = -1.666).
3. With the rigid condition, FK (17) gives lambda > 0, i.e. **-4H^2 < m^2** always. The registered tachyon has
   mu^2 = -7.72 < -4: it lies *outside* the window accessible in FK's rigid analysis. FK's bound (20) cannot be applied:
   the trial function Y = 1 is not admissible (int dw/a = int dy/rho^2 diverges at the cone, and the b.c. is not Neumann).
4. New structural remark obtained while mapping FK (script `fk_quadratic_form_check.py`). Multiplying FK (17a) by Y and
   using the soft b.c. gives the identity

       Q[Y] = lambda N[Y],   Q = int [3Y_y^2/(2 rho^2 phi'^2) + Y^2/rho^2] dy > 0,
       N = int 3Y^2/(2 rho^4 phi'^2) dy + 3Y_b^2/(2 rho_b^4 phi'_b^2 B).

   Numerically Q/(lambda N) = 1.00000001 on the eigenmode (N_bulk = 1.2e-8, N_boundary = -1.449e-4 in units Y_b = 1).
   Consequences (elementary linear algebra, valid if FK (8),(14) are the complete linear system): (i) all eigenvalues are
   real, because Q > 0 forbids N = 0; (ii) a mode with mu^2 < -4 requires N < 0, hence B < 0; (iii) N has a single negative
   direction (the rank-one boundary weight), and eigenvectors are N-orthogonal, so there is **at most one** mode with
   mu^2 < -4. This is an analytic counterpart of the argument-principle count in R2. It does **not** exclude bound states in
   -4 < mu^2 < 9/4 (the numerics find none), and the sign of N is **not** the sign of the physical kinetic norm (in 4D
   cosmology the Bardeen-potential norm and the Mukhanov-Sasaki norm differ by a factor of the eigenvalue; the analogous
   conjecture here is that the physical norm is proportional to lambda N = Q > 0, i.e. tachyon not ghost, consistent with
   Z_E > 0 in R3; this is **not proved**).

**FK's physical conclusions relevant to us (paraphrased).** For inflating branes the radion mass squared is typically
negative, giving a strong tachyonic instability; stabilised inflation needs both terms of (20) comparable (large bulk
gradients); FK expect the unstable configuration to restructure violently into the other static configuration, the one
with lower brane curvature.

## 2. BraneCode: Martin, Felder, Frolov, Peloso, Kofman, hep-th/0309001

Full nonlinear 1+1 evolution, metric ds^2 = e^{2B(t,y)}(-dt^2 + dy^2) + e^{2A(t,y)}dx^2 (their (3)), **two** orbifold
branes fixed at y = 0, 1, V = m^2 phi^2/2 + Lambda, U_i = M_i(phi - sigma_i)^2/2 + lambda_i (their (2)). Static dS
configurations come in pairs (one unstable, one stable, in the parameter regimes they study); the one with larger H is unstable through
the tachyonic radion m^2 = -4H^2 + m_0^2(H) (their (30)) and, in their simulations, re-configures violently to the second static configuration
with lower 4D curvature; when no second static solution exists the branes collide with a universal Kasner-like asymptotic
(their (40)), independent of the potentials.

Relevance: this is the only nonlinear study of the FK instability we found. It is for two branes with quadratic
potentials and flat spatial sections. **No nonlinear evolution of a single Z2 shell with a regular cone and a
superpotential bulk was found.** The BraneCode dichotomy (relax to a lower-H static state / run away) is the natural
template for the registered landscape R5: toward phi_b -> +1 the modulus rolls to V_E -> 0.0197 t (lower H), toward
phi_b -> -1 the field space ends. Whether the 5D evolution follows the 4D EFT is **open**.

## 3. Bulk-inflaton models (single dS brane, bulk horizon)

* **Himemoto-Sasaki gr-qc/0010035.** Test scalar V = V_0 + m^2 phi^2/2 (m^2 < 0) on fixed AdS_5 with a dS brane; separation
  phi = psi(t)u(r); regularity u(0) = 0 at the centre/horizon and Neumann at the brane; slow roll for |m^2| << H^2. No metric
  backreaction. Prior art for the *inner boundary condition* and the dS harmonic decomposition used in R1.
* **Himemoto-Tanaka-Sasaki gr-qc/0112027**, **Minamitsuji-Himemoto-Sasaki gr-qc/0303108.** Backreaction to O(phi^2);
  late-time behaviour dominated by the bound-state pole; effective 4D mass m_eff^2 = m^2/2 for H l, |m| l << 1; KK
  continuum above 9H^2/4 decays as a^(-3/2); corrections O(H^2 l^2), anisotropic part only at O(H^4 l^4).
  **Himemoto-Tanaka gr-qc/0212114** adds brane dissipation (reheating) to the same framework.
* **Kobayashi-Koyama-Soda hep-th/0009160.** Quantum fluctuations of a bulk inflaton: massless zero mode plus continuum
  m > 3H/2; zero mode dominates (> 30x).
* **Langlois-Sasaki hep-th/0302069.** Test bulk scalar V = M^2 phi^2/2 with brane coupling sigma = sigma_0 + (alpha/l)phi^2:
  a bound state exists only between two critical couplings alpha_zm < alpha < alpha_bs (the reader returned the window
  0 <= m_4^2 <= (3/2)H^2 for its mass; that inequality was read once and is **not verified**), quasi-normal (complex-mass) modes outside it, and an effective 4D potential
  M_eff^2 = M^2/2 + 2alpha/(l l_0) - alpha^2/(2 l^2). This is the closest analogue of "brane-potential curvature shifts the
  bound-state mass" (our quadratic-detuning statement d vs d_0 in R5), but without backreaction and around phi = 0 of a
  Z2-symmetric potential rather than on a BPS wall.
* **Koyama-Takahashi hep-th/0301165, hep-th/0307073.** Exactly solvable dilatonic model: bulk Lambda(phi) ~ (Delta/8 +
  delta) lambda_0^2 e^{-2 sqrt2 b phi}, tension ~ e^{-sqrt2 b phi}; power-law inflation, full backreaction, master variable
  Box_5 omega_c = 0 with decoupled junction omega_c' = 0; spectrum: zero mode + gap m >= -H/(Delta+2); 4D limit is
  Brans-Dicke with omega_BD = 1/(2b^2). No tachyon there: the exponential (scaling) structure makes the modulus a flat
  direction of a power-law attractor. **Kobayashi-Tanaka hep-th/0311197** explain these as shadows of higher-dimensional
  vacuum solutions and note that massive scalar modes do not separate simply.
* **Du-Wang-Abdalla-Su hep-th/0312087.** Single dS brane, bulk inflaton, generalized longitudinal gauge (FK's gauge);
  conclude that for the dS brane the radion mass squared is not positive. Background is perturbative (slowly
  varying scalar on AdS_5), brane tension not coupled to the scalar. This is the closest published analogue of R1's
  *set-up* (FK equations + single brane + horizon regularity); it does not contain a backreacted domain-wall background,
  a field-dependent tension, or a numerical eigenvalue comparable to R2.

Earlier, scalar-free results for the "-4": **Gen-Sasaki gr-qc/0011078** (two-brane RS with dS branes: radion
m^2 = -n K, i.e. -4H^2 for n = 4; and: a *single* positive-tension dS brane in pure AdS has no radion) and **Chacko-Fox
hep-th/0102023** (radion m^2 negative for dS branes, positive for AdS branes; 4D effective theory). The registered shell has
a modulus only because the bulk scalar makes Z != 0; R3's remark "pure AdS gives Z = 0" is the Gen-Sasaki statement.

## 4. Effective actions

* **Brax, van de Bruck, Davis, Rhodes hep-th/0209158** (and the review **Brax-van de Bruck hep-th/0303095**, eqs
  (58)-(65), (140)-(145)). Action S = (1/2k^2) int [R - (3/4)((d psi)^2 + U)], BPS: a'/a = -U_B/4, psi' = dU_B/dpsi,
  U = (dU_B)^2 - U_B^2. Moduli-space approximation (their (25); review (140),(141)):
  S_MSA = int sqrt(-g)[ f R + (3/4k^2) a^2(phi)U_B(phi)(d phi)^2 - (3/4k^2) a^2(sigma)U_B(sigma)(d sigma)^2 ],
  f = (1/k^2) int_phi^sigma a^2 dz; detuned tensions give potentials a^4 V (their (43)-(46) for the exponential case).
  Explicit Einstein-frame results are given only for U_B = 4k e^{alpha psi}.
  Dictionary: phi_ours = (sqrt3/2) psi, U_B = 4W/3, so (3/4)a^2 U_B = a^2 W, f = int a^2.
  Single-shell specialisation used by the first agent (second brane removed to the AdS_- throat, a -> 0):
  G_zz = (3/2)a^4/f^2 - a^2 W/f = (a^2/f^2) int_z^inf W_phi^2 f dz' >= 0, V_E = a^4 V/(4 f^2), mu^2 = 3 (ln V_E)_zz/G_zz.
  At phi_b = 0: stationarity gives c = 2/I_+ - 4/3 (= registered c_star) and mu^2 = **-7.719795918**
  (`prior_art_crosscheck_output.json`), identical to R3's HJ closed form to all printed digits. The positivity identity for
  G_zz (second form) was derived by the first agent; we did not find it in the papers read. **[method known; single-brane
  general-W specialisation and numbers specific]**
* **Palma-Davis hep-th/0406091.** Systematic low-energy expansion around BPS configurations for general U_B; zeroth order is
  the bi-scalar-tensor MSA theory; detuning potentials v_1, v_2 (V_i = U_B + v_i); validity k_5^2|U_B| >> |T|. Supports the use
  of MSA at O(t) for a *general* superpotential; it also claims that U_B alone can stabilise the moduli in some two-brane
  cases.
* **Kanno-Soda hep-th/0303203** ("Low energy effective action for dilatonic braneworld"). Gradient expansion with
  U = -6/(k^2 l^2) + V(phi), sigma = sigma_0 + sigma~(phi): single brane, first order: Einstein-scalar gravity with
  V_eff = sigma~/l + V/2 (the "1/2" is Himemoto-Sasaki's m^2/2) plus a non-conserved dark radiation chi_mu nu coupled to the
  scalar (their (39)); second order (68) gives f = 1 + l^2 k^2 V/12, a modified kinetic function and an S_CFT remainder.
  Their expansion is around AdS_5 with a *slowly varying* scalar; the registered wall has O(1) variation of phi across one
  AdS length, so their truncated formulas do not apply, but the structure "local two-derivative action + non-local
  CFT/dark-radiation piece" is the same as R3.
* **de Boer-Verlinde-Verlinde hep-th/9912012.** L = V + R + (1/2)G_IJ d phi^I d phi^J; HJ constraint (13); ansatz (14),(15)
  S = S_loc + Gamma, S_loc = int sqrt g [U + Phi R + (1/2) M_IJ d phi d phi]; (16) V = U^2/3 - (1/2) dU.dU; (21)
  beta^I = (6/U) G^IJ d_J U; (23) beta^K d_K Phi = 2 Phi + 6/U; (24) an equation for M_IJ.
  Dictionary (checked, `dbvv_transport_check.py`): phi_d = sqrt2 phi, U_d = -2W, V_d = -2U, Phi_d = 2Phi = I. Then (16) is
  U = W_phi^2/2 - (2/3)W^2 identically (residual 9e-16) and (23) **is** R3's transport equation W_phi I' - (2/3) W I = -1
  (residual 3.5e-10 on the wall). So R3 orders 0 and 2(R) are dBVV (16),(23) verbatim. dBVV (24) could not be matched
  because its index structure was not reliably readable; R3's M = 2W Phi'/W_phi comes from the Box-phi coefficient.
  **[method and transport equations known; choice of the IR-regular solution I(phi_b) on a two-vacuum wall, identification
  f = 2I with the shell's Planck mass, and Z = -W f'/W_phi: specific, though implied by MSA]**
* The equality "HJ closed form = MSA value" (R3 vs first agent) shows the two prior-art routes agree for this model; that
  agreement is expected on general grounds (both are the two-derivative truncation of the same on-shell action) and is
  reproduced here numerically, not proved in general.

## 5. Thick dS walls and first-order formalisms

* **DeWolfe-Freedman-Gubser-Karch hep-th/9909134.** S = int[-R/4 + (d phi)^2/2 - V], V = (1/8)W_phi^2 - (1/3)W^2,
  phi' = W_phi/2, A' = -W/3; curved slices: A' = -(W/3) gamma, phi' = W_phi/(2 gamma), gamma = sqrt(1 + 9 Lambda_4 e^{-2A}/W^2);
  thin-brane jumps Delta A' = -(2/3) lambda, Delta phi' = (1/2) lambda_phi; SUSY-QM factorisation for TT modes. Origin of the
  "fake superpotential" and of the flat-slice limit of R1's master equation (as SCALAR_SECTOR_DERIVATION.md says).
  Our W normalisation differs (ours: phi' = W_phi, U = W_phi^2/2 - (2/3)W^2).
* **Sasakura hep-th/0201130, hep-th/0203032.** Exact dS thick wall (Jacobi elliptic functions) for V = a + b cos(sqrt(2/3)phi)
  with the gamma-deformed first-order system (his (7)-(10)); scalar fluctuation potential V_e = V_t + 6H^2 a^4/beta^2 (his
  (45),(46)) so the *smooth* dS wall is stable against scalar modes. Contrast: a regular smooth dS wall with no thin shell
  has no free modulus; the registered instability is tied to the thin shell with a nearly flat (soft) potential. The
  registered detuned bulk is not gamma-BPS with the registered W; it is integrated as a second-order system.

## 6. Brane-world creation

**Garriga-Sasaki hep-th/9912118.** dS-brane instanton = two AdS_5 balls glued on an S^4; action reduces to the 4D dS
instanton for H l << 1 (e.g. the black-cigar/Nariai-brane pair-creation rate reduces to the 4D value exp(-pi/(3 G H^2))); tensor spectrum: zero mode + continuum from m = 3H/2. Scalar
modes are not treated. The registered solution (regular first cone at y = 0, shell at y_b) is the Lorentzian continuation of
the scalar-dressed version of this instanton. Standard Coleman-type reasoning (not checked in any source read here) says an
instanton with one negative mode describes a decay/creation process and more than one signals a non-dominant saddle; the
single tachyonic homogeneous mode found in R2 is the natural candidate for that negative mode. **[conjectural link]**

## 7. Dark-bubble cosmology

Banerjee-Danielsson-Dibitetto-Giri-Schillo 1807.01570 (PRL 121, 261301): 4D cosmology on the wall of a bubble mediating
AdS_5 -> AdS_5 decay; the brane is *not* Z2 (inside/outside differ), gravity is induced through junction conditions with
non-normalisable bulk modes ("The dark bubbleography", JHEP 02 (2024) 102). Developments read (abstracts only):
Danielsson-Panizo 2311.14589 (experimental tests); Basile-Borys-Masias 2507.03748 (PRD 113, 026009: equivalence-principle
problem for the proton); Danielsson-Giri 2511.21362 (gravity weakens below microns; an inflation-like phase from radiation
alone); Danielsson-Giri 2606.20942 (dark dimension / fat graviton); 2606.16547 (self-gravitating EM waves).
We found **no linear stability analysis of a scalar-dressed dark-bubble wall** comparable to R1-R2. The registered model
differs structurally (Z2 shell, one-sided regular cone, RS-type localisation), so dark-bubble results do not transfer; the
shared idea is only "dS on a shell between/around AdS_5 regions".

## 8. 2025-2026 items

* Karmakar-SenGupta 2605.18403, "Tachyonic (in)stability in RS braneworld scenarios" (May 2026): radion as scalar-tensor
  field, Damour-Esposito-Farese-type instability; states conditions for radion stability on dS_4/AdS_4 branes (section 4.5;
  the text of that section could not be retrieved, so its formulas are **not verified** here). Two-brane RS/Goldberger-Wise.
* Banerjee-SenGupta 1705.05015 claim that the vacuum energy of a dS 3-brane by itself generates a modulus potential with a metastable minimum;
  this is in tension with FK/Gen-Sasaki and with R2-R3 and should not be cited as support for stability.
* Hassfeld-Hebecker-Schiller 2505.07934 ("Localized gravity, de Sitter, and the horizon criterion"): ETW branes with a
  brane scalar or a bulk modulus with a brane-localized potential; refined horizon criterion. Anastasi-Angius-Huertas-
  Uranga-Wang 2501.03310: swampland constraints for localized gravity ("relative quantum gravity"). Neither computes a
  fluctuation spectrum. Observation (ours, heuristic): R5's statement that every natural detuning gives |mu^2| >= 6.6, i.e.
  |V_E''|/V_E = |mu^2|/3 >= 2.2 in Planck units, is of the type allowed by the refined dS conjecture (Ooguri-Palti-Shiu-Vafa
  1810.05506: min V'' <= -c' V); the slow-roll window d ~ d_0 +/- 1% is the braneworld version of the eta problem. We found no
  2025-2026 paper deriving an eta-problem for a BPS-wall shell modulus.

---

## 9. Agreement with the orchestrator's numbers

| item | orchestrator | prior-art route | status |
|---|---|---|---|
| R2 mu^2(t=1e-3) | -7.7178716 | FK (8a),(8b),(14) + cone regularity: -7.71787162 (dv, v0 stable to 5e-9) | agree |
| R2 uniqueness | one zero by winding | FK scan mu^2 in [-40, 2.2]: one sign change; plus sec. 1 item 4: at most one mode below -4, spectrum real | agree |
| R3 mu^2(t->0) | -7.719795918 | MSA 3(ln V_E)_zz/G_zz = -7.719795918 | agree |
| R3 c_star | 2/I_+ - 4/3 | MSA stationarity, same | agree |
| R4 t-dependence | -7.7197959 + 1.9243 t | FK roots at t = 1e-3, 3e-4, 1e-4: slopes 1.924, 1.924 | agree |
| R5 d_0 | 1.1134966 | -(ln V_E)_zz = 1.1134966 | agree |
| R5 Theta end | 0.548 | MSA canonical distance 0.54791 | agree |
| R5 plateaus V_E/t | 0.12420, 0.01976 | (1-c)(5/9)^2 = 0.124199; (1+c)(1/9)^2 = **0.019723** | second differs by 0.2% |

The 0.01976 vs 0.019723 difference is small; the MSA grid value at z = -12 is 0.0222 (the approach is slow, ~e^{-2|z|/9}),
so the orchestrator figure is probably a finite-distance value rather than the asymptote. Worth a one-line check.

## 10. PRIOR-ART BOUNDARY

**R1 (linear scalar sector).** Entirely known in method and equations: gauge, xi = -2psi, the constraint, the bulk
system, the no-bending statement and the scalar junction are FK (4),(5),(8a),(8b),(14); the cone regularity and the 9/4
threshold are Himemoto-Sasaki / Garriga-Sasaki / Langlois-Sasaki; the single-brane + longitudinal-gauge combination
appears in Du et al. hep-th/0312087. The second-order master equation (M_y) for psi and its Schrodinger potential are
re-packagings. Not found in the literature: the explicit treatment of the soft regime (sigma_t'' ~ 0) in which (14) becomes
eigenvalue-dependent, and the resulting indefinite-form argument (sec. 1 item 4). That argument is elementary and should be
described as a remark, not a theorem, until the second-order action is checked.

**R2 (one tachyon, mu^2 = -7.7179).** Known qualitatively: inflating branes with a bulk scalar generically have a
tachyonic radion with |m^2| ~ 4H^2 (FK, BraneCode, Gen-Sasaki, Chacko-Fox, Du et al.). Specific: the number; the fact that it
lies *below* -4H^2 (impossible under FK's rigid b.c.); the single-shell + regular-cone backreacted background; the
complex-plane count. No published computation for a superpotential domain wall with a linearly detuned Z2 shell was found.

**R3 (HJ effective action).** Method known twice over: dBVV (16),(23) are literally R3's order-0 and transport equations;
Brax et al. (25) is the same two-derivative action in MSA language; Kanno-Soda and Palma-Davis give the systematic
expansions and the non-local remainder. Specific to this model: the single-shell reduction with the IR-regular I(phi_b),
the closed form mu^2 = -4(3c^2-4c+8)/(c(3c+4)), the explanation of c_star as modulus stationarity, the positivity identity
for G_zz, and the split "-4 - Delta". The interpretation of the "-4" as the FK/Gen-Sasaki universal shift is an
identification by value and structure; a derivation that the EFT's -4 *is* FK's -4H^2 term (e.g. from the conformal
coupling f R with R = 12H^2) was not located in the papers read. **[plausible, not proved here]**

**R4 (5D vs EFT at O(t)).** A quantitative test of the moduli-space/HJ truncation against the exact 5D linear spectrum on
a dS brane, with linear convergence in the detuning. We found no such test in the literature read (Palma-Davis and
Kanno-Soda discuss validity parametrically; BraneCode tests nonlinear dynamics, not the MSA mass). Treat as specific to
this project, with the caveat that absence of evidence after a bounded search is not proof of novelty.

**R5 (landscape, fast-roll hilltops, tuning d ~ d_0).** Known in spirit: FK (20)-(22) already state that sufficient
inflation needs a finely balanced m_0^2(H) ~ 4H^2; BraneCode shows the relax-or-run-away dichotomy; Langlois-Sasaki show how
a quadratic brane coupling shifts the bound-state mass. Specific: all numbers, the finite field-space distance to the AdS_-
throat, and the statement that natural detuning families give |mu^2| in [6.6, 38]. The eta-problem/swampland reading is
an interpretation, not a derived result.

**Not found at all (open):** nonlinear 5D evolution of the registered single shell; a second-order (symplectic) action
fixing the sign of the tachyon's norm; any dark-bubble or 2025-2026 swampland paper that computes this spectrum.
