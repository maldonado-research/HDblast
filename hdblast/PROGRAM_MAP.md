# HDBLAST program map

> **About this map.** Compiled 27 September 2026 from the author's research archive. Paths are relative to this `hdblast/` folder. Paths under `../new-files/` (`DB3` = `../new-files/D-Blast 3`) refer to the author's private research archive and are not part of this public repository; the published checkpoints are in [`checkpoints/`](checkpoints/). The recomputation script and its output are in [`../research/HDBLAST_CHECKPOINT_20260927/program_map/`](../research/HDBLAST_CHECKPOINT_20260927/program_map/).


Ricardo Maldonado's higher-dimensional blast (HDBLAST) research program. Map compiled 27 September 2026 from local files. This is a navigation and status document. It adds no new physics. Every number below is quoted from a named source. The numbers that can be recomputed cheaply were recomputed by `verify_program_map.py` (results in `program_map_checks.json`, section 10).

**Status labels used here**

| Label | Meaning |
|---|---|
| EXACT | Exact symbolic or rational identity, verified by a saved script (for example SymPy or exact rational arithmetic) |
| CERTIFIED | Computer-assisted interval proof, conditional on the premises its source lists |
| NUMERICAL | Floating-point result with resolution or independent-integration evidence; no error certificate |
| SCREENING | Early observational comparison using approximate likelihoods or inputs of uncertain provenance |
| CONDITIONAL | Holds only under stated extra assumptions |
| NEGATIVE | A tested route that fails |
| WITHDRAWN / DEMOTED | Publicly corrected |
| OPEN | Not established |

Path shorthand: **DB3** = `../new-files/D-Blast 3`. **LATEST** = `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922`. Section 8 gives the full paths.

---

## 1. The hypothesis in plain words

HDBLAST proposes that our universe is a four-dimensional "brane" (a shell) inside a five-dimensional spacetime. It asks whether a violent gravitational event in the fifth dimension (the "blast") could have transferred energy onto the brane and started the hot Big Bang. The project's own public guide puts the question this way: *"Could a powerful event in an extra direction of space have transferred energy into our Universe and helped create its hot beginning?"* It calls this a speculative, testable hypothesis, not a discovery (DB3/`HDBLAST_GENERAL_PUBLIC_GUIDE_V20.md`).

The program has run in three eras:

1. **Observational phenomenology, August 2025 to May 2026.** A smooth broken power law (SBPL) in the stochastic gravitational-wave background, with a "knee" in the pulsar-timing (PTA) band, was taken as the blast's proposed imprint. It was tested against public PTA products with kill-switch criteria.
2. **The PHYS-M mathematical series, June to early September 2026.** This era produced exact and certified mathematics for a registered five-dimensional Einstein–scalar model: bulk collisions, cone and tail analysis, the registered shell root, and a certified linear response.
3. **Shell dynamics, 16–22 September 2026 (Chats 9–14 and two further packages).** This era asked what the registered shell universe actually does: its stability, its nonlinear roll-off, the endpoint, and a matter extension.

**Key structural gap (OPEN):** the PTA-knee phenomenology of era 1 has **not** been derived from the registered 5D model of eras 2–3. The 8 Sept 2026 handoff says so directly: *"The SBPL is a phenomenological imprint. The optional brane-world mapping is illustrative, not a complete derivation of the blast mechanism."* (DB3/`VECTOR STORE/HDBLAST VECTOR STORE/DBlast 2 vector store/HDBLAST_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md`, section A2.6 item 7). The 22 Sept checkpoint also defers observational predictions until a radiation era has been demonstrated (LATEST/`00_READ_FIRST.md`, section 5).

## 2. Registered model definition

The sources are the M462R1 root certificate (DB3/`untitled folder 120/M462_MAIN_REPORT_RESEND.md`), Chat 9 section 1, Chat 11, the frozen solver (LATEST/`frozen/registered_solver.py`) and the action in LATEST/`matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md` section 1.

- **Field content.** Five-dimensional Einstein gravity with one canonical bulk scalar φ. The shell (brane) has tension σ(φ). The bulk is a mirror-symmetric (Z2) double copy with a single shell. The signature is (−++++). Units: M=κ=v=W₀=1 in M462; κ₅=1 in the numerical packages.
- **Action** (doubled bulk, one shell action):
  S = Σ_{a=1,2} { (1/2κ₅²)∫_{M_a}√−g [R − (∇φ)² − 2U(φ)] + (1/κ₅²)∫_Σ √−h K_a } + ∫_Σ √−h [−σ(φ)/κ₅² + L_m].
  In the registered (vacuum) model L_m = 0. The bulk field equations are R_AB − ∂_Aφ∂_Bφ − (2/3)U g_AB = 0 plus the scalar equation (Chat 11).
- **Superpotential and potential.** W(φ) = 1 − φ + φ³/3 and U = ½W_φ² − (2/3)W². The equivalent form is U = ½(φ²−1)² − (2/3)W². There are two AdS vacua, at φ = −1 (W = 5/3) and φ = +1 (W = 1/3). With zero detuning the background is the flat wall φ = −tanh y.
- **Tension.** σ(φ) = 2W(φ) + δ(1 + cφ). The detuning δ is written t in Chats 9–14.
  - Registered values: δ = t = 10⁻³ and c = c★ = 2/I₊ − 4/3 = 0.5975949350280132, with I₊ = 1.0357712571566784 (the `IP` constant in the frozen solver; recomputed in section 10).
  - Chat 9 showed that c★ is exactly the condition for a static shell at φ_b = 0 in the four-dimensional effective theory.
- **Junction conditions** (normal pointing from the bulk toward the shell).
  - Vacuum: K^μ_ν = (σ/6)δ^μ_ν, and the scalar condition n·∂φ = −σ′(φ)/2 (Chat 11; the solver residuals are `rpy/rho − sigma/6` and `py + sigma1/2`).
  - With shell matter (proposed 22 Sept extension, not part of the registered model): nA = (σ+κ₅²ρ)/6, nB = (σ−κ₅²(2ρ+3p))/6, nφ = −(σ′+κ₅²j)/2, and ρ̇ + 3H(ρ+p) = jφ̇.
- **Coordinates.**
  - Static shells: ds² = dy² + ρ(y)²γ, with γ the unit de Sitter metric. The regular "cone" is at y = 0 and the shell at y_b. The shell Hubble rate is H = 1/ρ_b.
  - Evolution chart: ds² = e^{2B(t,z)}(−dt²+dz²) + e^{2A(t,z)}dx₃², with the shell at z = 0 and the bulk at z < 0. The static solution has A = t + ln ρ and B = ln ρ.
- **Registered (unstable) shell.**
  - Certified root in M462R1, with local-box uniqueness. The cone lies near φ = −1 (φ_h = −1 + t^{9/5}δ in M462 notation).
  - The de Sitter radius is ρ₀ = 78.82817714224423 (quoted in LATEST/`static_branch/INDEPENDENT_BRANCH_REVIEW.md`), so H₀ = 1/ρ₀.
- **Variants that define different models.**
  - (i) Curved tension σ = 2W + t(1 + cφ + dφ²/2). The stable shell S₈⁄₅ has d = 8/5 (Chats 9–12).
  - (ii) Shell matter χ with a φ-dependent mass (22 Sept proposal).
  - The sources warn explicitly that adding matter, higher-curvature terms, extra fields or a timelike extra dimension makes a new branch of the investigation (DB3/`untitled folder 152/00_READ_FIRST.md`).
- **New static "+1 branch"** (22 Sept 2026, NUMERICAL; section 4). φ_b = 0.9999159473169134, ρ_b = 129.9247628497, H² = 5.92401479433×10⁻⁵ and H/H₀ = 0.606721732. The cone displacement is η_h = φ_h − 1 = −1.3367×10⁻²³.
  - Expansions: φ_b = 1 − (9c/64)δ + O(δ²), and H² = δ(1+c)/27 + δ²[(1+c)²/36 − c²/384] + O(δ³).

## 3. Dated chronology

"Rec." means a Zenodo record number found in local files. The HDBLAST Zenodo concept (all-versions) DOI is **10.5281/zenodo.17088132**. Version records 0–20 come from the 29 Aug 2026 provenance audit table (DB3/`VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md`); the verifier re-parses them into `program_map_checks.json`. Later records are located in the files named in section 7. Version labels mix internal, release and Zenodo numbering (audit item MC-004). For example, "HDBLAST v27/v28" are internal labels on concept versions 10–11.

### Era 1: observational phenomenology (SBPL knee, kill switches)

| Date | Checkpoint | Rec. | Content (as described in its source) | Label |
|---|---|---|---|---|
| 2025-08-13 → 08-21 | Brane-world two-link repro packs (separate concept 16866883) | 16868085, 16866884, 16868562, 16887277, 16896080, 16907982 | "one λ links a GW spectral break and ΔN_eff" | SCREENING; historical companion |
| 2025-08-22, 08-24 | HDBC v2.3 / v3.9 (concept 16929972) | 16929973, 16937520 | Reproducibility bundles (PISC/PLI overlays, DECIGO, sub-mm bounds, ΔN_eff) | historical companion |
| 2025-09-06 | Cross-PTA/LISA/CMB pipeline (concept 17069899) | 17069900 | Pipeline companion | historical companion |
| 2025-09-09 → 11-22 | HDBLAST concept versions 0–5: PTA→LISA meta-release, GEN5, v5.8.1 mini-kits, vNOW | 17088133, 17136572, 17211812, 17335636, 17547897, 17683486 | Phenomenology and toolkits | SCREENING |
| 2025-12-12 | vLASTMILE+ canonical snapshot | 17918095 | Canonical SBPL: f_k = 3.16×10⁻⁸ Hz, Ω_k = 8.0×10⁻⁹, α₁ = +3, α₂ = −2, Δ = 2 | SCREENING |
| 2025-12-17 → 12-29 | REALONLY v22 / GrandStatus; vNEXTLEVEL+ covariance frontier; vKILLSWITCH++ | 17968738, 18048513, 18089007 | ΔBIC, curvature, PCI, shared-knee BF, AR1 covariance stress | SCREENING (provenance caveats, section 5) |
| 2026-01-03, 01-06 | Internal "v27" kernel-family thresholds; "v28" kernel-free (Fréchet–Hoeffding) certificate on NANOGrav 15-yr KDE free spectra (HD) | 18157426, 18158842 | N = 10 bins: all four stationary kernels PASS; kernel-free n_eff lower bound 1.6707 | SCREENING (covariance robustness only) |
| 2026-01-13, 01-22 | RS2 curvature killswitch, cross-band closure, LISA-band knee; Theory Contract v3 (knee := f50, K2 filter) | 18236838, 18341882 | f50 = x50·f_scale, x50(K2) = 1.339139 | CONDITIONAL test harness |
| 2026-01-28 | v15.0 shape fingerprint + PTA KDE pipeline scan | 18406745 | Identity-inversion fingerprint | SCREENING |
| 2026-02-27 | v7 Bessel-K fingerprint + knee-scaling tests, explicit kill-switch criteria | 18807647 | Frozen gate (section 6) | registered test |
| 2026-03-06/07 | v7plus frozen gate; go/no-go memo (local only) | none | Three real-data-style pilots **FAIL / FAIL** under the strict gate and the stability sweep | NEGATIVE (method-vs-null undecided) |
| 2026-04-10 | Status snapshot v8 | 19499967 | Repeated IPTA_DR2B_TOP5 files judged non-authoritative; no authoritative export in hand | status |
| 2026-04-24 | V50/V51 proof-lane, sidecar falsification kit (label 2026.04.24-v18) | 19747409 | Methods | status |
| 2026-05-20 | V180 / Zenodo v19: D481 internal certificate, "No Proof Upgrade" | 20313949 | Methods | status |

### Era 2: PHYS-M mathematical series and certified linear response

The summaries below come from checkpoint titles and status lines. They were not re-audited for this map.

| Date | Checkpoint (folder under DB3) | Rec. | Content | Label |
|---|---|---|---|---|
| 2026-06-25 | PHYS-M51 covariant bulk–brane tensor transfer (earliest PHYS-M file name located) | none | Theory only | EXACT/theory |
| 2026-07-10 | M285–M312 handoffs (folders 64–68) | none | M305 derivative-redundancy audit; M306 planar-collision TT no-go; M307 KK energy theorem; **M309–M312 non-collinear 3D trigger and "same-event" TT energy ledger** | M309/M311 later WITHDRAWN/DEMOTED |
| ~2026-07-10 → 07-14 (inferred from the bracketing dated handoffs) | M313–M400 certificates (folders 69–88) | none | 4+1D numerical-readiness gates | stopped fail-closed at gate A10 |
| **2026-07-14** | **Zenodo v20 correction** (folder 89 ZENODO) | **21367209** | **Withdraws** the M309 radiative TT certification (the wave vector was absent from the declared quadratic stress). **Demotes** M311's F_KK^E = 1.33765% to an imposed-window diagnostic. A01–A09 true; A10–A12 false. | WITHDRAWN (partial); record stays published |
| mid-July → 2026-08-11 | M408–M453D development lanes (folders 90–111; dated files 07-17 to 08-11) | none | Development/synthetic-only lanes | development |
| 2026-08-16 → 08-17 | M455D program; M456R2; M457R1 tuned-shell flatness; M459R1–M461R1 formal response; **M462R1 certified registered shell root** (folders 112–120) | none | Fixes the registered shell: t = 10⁻³, c★ = 2/I₊ − 4/3, root box with local uniqueness | CERTIFIED |
| 2026-08-17 → 08-18 | M463–M467; **M468R1 → Zenodo v21** (folders 121–126) | **21998599** (concept 17088132) | Uniform TT and principal-L²(H⁴) scalar wave-packet quiescence on the certified M466 4+1D open-FRW kinetic tail. Isotropic attraction fails; logarithmic memory persists; close tensor precedent W. Li, arXiv:2401.08437 | CERTIFIED (computer-assisted), sector-limited |
| 2026-08-18 → 08-25 | M471R1 first-cone 16-field normal form; M472 characteristic Mellin-shift obstruction; M473 TT Goursat-to-terminal channel; M474 preflight (folders 127–129) | none | Cone-to-tail transfer program | EXACT pieces; end-to-end HOLD |
| 2026-08-27 → 09-02 | M482–M488G adjoint Green identity, center-safe adjoint, null-moment obstruction, operator closure, causal excision (folders 130–135) | none | Registered linear-response infrastructure | EXACT/CERTIFIED pieces |
| 2026-09-03 | **M489G-A → v22** (folder 136) | **22285737** | Riemann preconditioning; registers two sources and two detectors; response HOLD | CERTIFIED bounds |
| 2026-09-03 | M489G-B (folder 137): nonzero response entries, determinant unresolved | not separately published; its sealed package is bundled in M489G-C | Rank ≥ 1 | CERTIFIED |
| 2026-09-03 | **M489G-C → v23** (folder 138 ZENODO) | **22287013** | Rank two: det M ∈ [2.19870426928051, 2.46151588838361] | CERTIFIED (conditional on inherited premises) |
| 2026-09-04 | Five-field forward vs adjoint numerics (folder 138); continuum components (folder 139) | none | Numerical det: adjoint 2.330107384842 vs forward 2.330107321657 (2.7×10⁻⁸ apart); scalar reference det enclosed to width < 3×10⁻³⁶ | NUMERICAL / CERTIFIED |
| 2026-09-04 (published 09-05) | **Conditional global adjoint certification → v24** (folder 140) | **22347452** | det M ∈ [2.12278154373842557755059, 2.53743323308364125340774] over 176 characteristic cells | CERTIFIED, CONDITIONAL on inherited model, background-root membership, Green identities and stability bounds |
| 2026-09-05 | Forward error-control components (folder 141) | none | Startup contraction κ ≤ 1.2045×10⁻⁸; forward propagation constants | CERTIFIED components |
| 2026-09-06 | **Chat 8** (folders 142–144): fixed startup and scalar reduction; scalar closure; **independent global forward response** det R ∈ [2.3290685, 2.3311463] (inside folder 144 `inputs/PREVIOUS_GLOBAL_FORWARD.zip`); information recovery; quantum-production benchmark | none | See section 4 | CERTIFIED (linear, coefficient-normalized) plus NEGATIVE energy test |
| 2026-09-06 | Chat 7 full transcript/handoff (folder "142 new chat") | none | Chronology of M489G → v24 publication | record |
| 2026-09-09 | v24 APS author-review package (folder "145 APS JOURNAL") | none | Clean-copy replay PASS. **Not submitted**; an earlier PRD email was an editorial inquiry only | record |

### Era 3: shell dynamics (September 2026)

| Date | Checkpoint (folder under DB3) | Rec. | Main result | Label |
|---|---|---|---|---|
| 2026-09-16 (final ZIP 09-17) | **Chat 9** shell dynamics (146) | none | Registered shell has exactly one tachyonic scalar bound state, −7.71788 < μ² < −7.71786, growth e^{1.657Hτ}. Closed-form 4D effective theory; stable shells for d > d₀ = 1.1135. 4D roll-off: no radiation era; n_s ≤ 0.93. | CERTIFIED (tachyon) / NUMERICAL / NEGATIVE |
| 2026-09-18 | **Chat 10** stable shell (147) | none | S₈⁄₅ exists (Poincaré–Miranda) and has no scalar bound state below 9H²/4. Two-branch structure, 12 predicted shells found in 5D. | CERTIFIED (conditional) / NUMERICAL |
| 2026-09-18 | **Chat 11** symbolic verification (148) | none | All 15 linearized Einstein components, the scalar equation, the junction algebra and the certificate ODEs; 7 negative controls rejected | EXACT |
| 2026-09-19 | **Chat 12** all sectors (149) | none | S₈⁄₅ is mode-stable in the tensor, vector and special-harmonic sectors. The registered shell's only instability is the Chat 9 scalar mode. | EXACT + referee |
| 2026-09-21 | **Chat 13** nonlinear 5D roll-off at t = 0.1 and 0.03 (150) | none | Two fates: relaxation toward an empty Randall–Sundrum de Sitter brane, or expansion reversal and collapse. No radiation branch. | NUMERICAL / NEGATIVE |
| 2026-09-22 | **Chat 14** roll-off at the registered t = 10⁻³ (151) | none | Growth rate 1.65714–1.65719 vs certified 1.65719. Same two fates. +1 side H_J/H₀ = 0.63–0.64 at the chart freeze (exact RS 0.6067). Throat turnaround at H₀τ ≈ 5.894–5.895, φ_b = −1.942. | NUMERICAL / NEGATIVE |
| 2026-09-22 | **Registered-shell controls** (152; ZIP only, plus inputs) | none | New solver finest eigenvalue 1.657193663767 vs prior 1.657193631245. Exact constraint-transport identity. Single-bump initial-data obstruction; balanced two-bump candidates. | NUMERICAL / EXACT |
| 2026-09-22 | **Scalar-profile branch and matter** (153 = LATEST) | **22922928** (published 2026-09-23) | The constant φ = 1 endpoint fails the scalar junction. Corrected nonconstant static +1 branch. Balanced-disturbance evolutions show a constraint plateau. Matter extension with signed energy exchange. | NUMERICAL / EXACT / NEGATIVE / CONDITIONAL |

About record 22922928: its title ("HDBLAST: Scalar Junction Consistency and a Corrected Static de Sitter Branch") and its date (23 Sept 2026) come from God-Plays-Dice site data (DB3/`GPD-site-deploy/gpd5-articles.json`). That data calls it the newest HD-Blast record. **Its version label, its concept membership and its uploaded file list were not verified** (zenodo.org is not reachable from this machine, and a web search on 27 Sept 2026 returned no index entry). No local evidence was found that Chats 9–13 or the folder-152 package were deposited as separate Zenodo records. The LATEST package embeds the Chat 14 archive (`source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip`), and its `REFERENCE_INPUT.json` identifies the folder-152 ZIP by hash.

Provenance note: the folder-152 and LATEST packages record a workspace named `CHAT8_CONTINUATION_20260922` (LATEST/`frozen/PROVENANCE.json`). They audit and extend the Chat 9–14 folders; they are not "Chat 15".

## 4. Claims ledger

### 4a. Established within stated scope

| Claim | Label | Source |
|---|---|---|
| Registered shell root exists, with local-box uniqueness (t = 10⁻³) | CERTIFIED | M462R1 (folder 120) |
| Registered shell has a tachyonic scalar mode, −7.71788 < μ² < −7.71786, growth rate 1.65719 H. No unstable tensor, vector or special-harmonic mode exists. Uniqueness of the scalar bound state is numerical only. | CERTIFIED (conditional on the M462 enclosure and hand lemmas); equations EXACT (Chat 11); other sectors EXACT (Chat 12) | Chats 9, 11, 12 |
| Closed-form 4D effective action: μ²(t→0) = −4(3c²−4c+8)/(c(3c+4)) = −7.719796, plus 1.9243896 t; agrees with 5D to about 10⁻⁶ | EXACT derivation + NUMERICAL agreement | Chat 9 |
| Stable shell S₈⁄₅ (d = 8/5) exists and is mode-stable in every sector; the constant harmonic (uniqueness) is open | CERTIFIED + EXACT, conditional on cited self-adjointness/completeness | Chats 10, 12 |
| Linear response of the registered two-source/two-detector experiment has rank two (det > 0) | CERTIFIED, CONDITIONAL on inherited premises; coefficient-normalized, not physical energy/noise | v23, v24, Chat 8 forward |
| Two-sector linear kinetic-tail quiescence on the M466 background | CERTIFIED, sector-limited; tensor ingredients have prior art (W. Li 2024) | v21 |
| Constraint-transport identity (∂t ∓ ∂z){e^{3A}(C_H ± 2C_M)} = 0 (continuum, no dissipation) | EXACT | folder 152 |
| Single-bump initial-data obstruction: (ρ²D)_y = ερ⁴φ_y b/6; D_b ≠ 0 breaks the second time derivatives of the junctions | EXACT | folder 152 |
| Constant φ = 1 cannot satisfy the scalar junction when δc ≠ 0 | EXACT | LATEST |
| Matter-extension junctions and energy exchange ρ̇ + 3H(ρ+p) = jφ̇ (47 exact checks, 4 wrong-formula controls) | EXACT, as a proposed extension | LATEST/matter |

### 4b. Numerical results (floating point, with convergence evidence)

- Chats 13–14: the nonlinear roll-off has two fates at t = 0.1, 0.03 and 10⁻³. The turnaround (φ_b = −1.9412 to −1.9417, H₀τ = 5.8944–5.8947) is converged across three grids.
- Growth-rate calibration: 1.657193663767 (eigenvalue) and 1.657194158956 (evolution fit) against the prior 1.657193631245 (folder 152).
- New static +1 branch: six detunings; 18 independent second-order integrations; junction residuals of about 3.6×10⁻¹⁶ (metric) and 1.6×10⁻¹⁴ (scalar). These are residuals, not error bounds.
- Two-branch shell structure: 12 shells predicted by the 4D formula and found in 5D (Chat 10).
- Forward and adjoint numerical determinants agree to 2.7×10⁻⁸ (folder 138).

### 4c. Negative results

- **No hot, radiation-dominated branch** in the homogeneous, classical, single-scalar, matter-free roll-off at the registered parameters. The +1 side approaches an empty de Sitter brane; the throat side collapses (Chats 9, 13, 14). *This mechanism in this model is excluded as the sole cause of a hot Big Bang*, as stated in Chat 13.
- Shell-modulus inflation is not viable: n_s ≤ 1 − 4/N ≈ 0.93, with no exit (Chat 9, 4D). An earlier note saying "n_s ≈ 0.965" was wrong and was corrected.
- The constant-σ brane is excluded as a source of ordinary radiation under minimal junction assumptions (Chat 8, folder 142).
- A single positive sech² quantum-production pulse cannot meet the conditional heating target under the stated cutoff assumptions (at least 1,067 independent spectators would be needed) (Chat 8, folder 144).
- The v7plus Bessel-K PTA pilots **FAIL / FAIL** under the frozen strict gate (7 Mar 2026). The memo leaves open whether this reflects the observable mapping, the derivative estimation, or a true null.
- The balanced-disturbance evolutions have a Hamiltonian-constraint plateau: 0.00754 → 0.00655 on the last refinement at t = 0.5. They are not a validated late-time solver (LATEST).

### 4d. Withdrawn, demoted or corrected

| Item | Correction | Where recorded |
|---|---|---|
| PHYS-M309 radiative (Fourier-resolved) TT certification: the July 2026 gravitational-radiation claim | **Withdrawn** in v20, 14 Jul 2026: the TT matrix was evaluated at a wave vector absent from the declared quadratic stress | DB3/`untitled folder 89 ZENODO/HDBLAST_ZENODO_V20_CORRECTION_AND_CLAIM_BOUNDARY.md` |
| PHYS-M311 F_KK^E = 1.33765% "same-event" KK energy fraction | **Demoted** to an imposed-window/profile diagnostic | same |
| All pre-v20 HDBLAST Zenodo summaries | Marked `affected_by_later_claim_correction`; must not resurrect the M309 certification | Aug 29 audit |
| Chat 13 drafts | The plateau had been called "settled"; a quoted 0.603 was a script bug (correct value 0.638); the reversal had been attributed to U unbounded below (the true cause is sub-balanced tension); the t = 10⁻³ failure had been misdiagnosed | Chat 13 section 5 |
| Chat 14 drafts | Plateau 0.61 → 0.63–0.64 (the far grid under-resolved the receding wall); turnaround −1.953 → −1.942 (interpolated); "converged to −10H₀" was too strong | Chat 14 section 4 |
| Late-time Randall–Sundrum comparison with constant φ = 1 | Replaced by the nonconstant scalar-profile branch; the H change is only about −7.87 ppm | LATEST |
| Chats 9–11 wording "μ² = −4 modes" | Precise criterion: a vanishing trace-free Hessian of the harmonic | Chat 12 |

### 4e. Open

1. Stability, uniqueness, rigorous existence and dynamical attraction of the new static +1 branch. The branch's cone sits near φ = +1, while the original shell's cone sits near φ = −1; connecting the two needs the departing wall and the global spacetime.
2. The asymptotic +1 plateau beyond the chart freeze (0.63–0.64 at the freeze vs the exact 0.6067), the nature of the final singularity, and the origin of the late constraint feature at z ≈ −0.9.
3. Constraint-controlled finite-amplitude evolution, and a benchmarked proper-clock slicing for late times.
4. A matter sector with chosen physical scales and couplings; renormalized particle production with backreaction through all three junctions; thermalization; a sustained radiation era; the residual vacuum energy.
5. The primordial perturbation spectrum and any observable of the registered model. No observational prediction of the 5D model has been derived.
6. A derivation linking the registered 5D model to the PTA-knee phenomenology.
7. Physical normalization, noise model and observability of the certified linear response (coefficient norms only; condition number about 1.98×10⁵).
8. External novelty and priority. No source claims them, and all reviews are internal (no external peer review).

## 5. Early observational screening numbers (for context only)

The Dec 2025 layer is quoted in the 8 Sept 2026 handoff (section A2.5). It reports:

- ΔN_eff = 7.09×10⁻⁴ (PASS against 0.06/0.3 guardrails);
- Ω_gw(3 mHz) = 1.32×10⁻²³ and Ω_gw(25 Hz) = 2.28×10⁻³⁵;
- ΔBIC(SMBHB − SBPL) = 26.2 in a binned-Gaussian screen;
- shared-knee log₁₀BF = 1.07, with shared median f_k = 3.31×10⁻⁸ Hz.

The same source lists its limitations:

- the likelihoods are screening only;
- input provenance is ambiguous ("REALONLY" is an author label);
- a covariance loophole exists;
- the highest-frequency bin carries strong leverage;
- the ringdown inputs are placeholders;
- flexible SMBHB models were not tested.

These are **not** confirmations, and they are not derived from the registered 5D model.

## 6. Registered predictions and falsification tests

### 6a. Observational (era 1): the PTA "knee"

| Item | Registered content | Status |
|---|---|---|
| Canonical SBPL (vLASTMILE+, 12 Dec 2025) | f_k = 3.16×10⁻⁸ Hz, Ω_k = 8.0×10⁻⁹, α₁ = +3 (rises as f² below the knee), α₂ = −2 (falls as f⁻³ above it), Δ = 2 | SCREENING |
| Falsifiers listed for the SBPL | Full-likelihood PTA analyses prefer an unbroken power law or a flexible SMBHB turnover; the high-frequency downturn disappears with better noise models; the true bin covariance erases the preference; cross-PTA knee posteriors become inconsistent; the amplitude violates ΔN_eff; LISA/LVK see a large background incompatible with the steep tail | Registered; none decisively executed with full likelihoods |
| Knee convention (Theory Contract v3, 22 Jan 2026) | knee := f50 with F(f50) = 0.5. RS2 filter F(x) = ((x²K₂(x))/2)², so x50 = 1.339139 and ln x50 = 0.292027 (recomputed, section 10). Cross-band: f50_LISA(ℓ,ν) = x50(ν)·√(K_RS2/ℓ), ξ = f50_PTA/f50_LISA, y₀/ℓ = ln(1/ξ) | CONDITIONAL test harness |
| Kernel-free covariance certificate (6 Jan 2026, rec. 18158842) | NANOGrav 15-yr KDE free spectra (HD), N = 10: ρ_adj,max = 0.514 < thresholds (AR1 0.659, SE 0.767, Matérn-3/2 0.728, Matérn-5/2 0.744) | PASS as a covariance-robustness screen only |
| Bessel-K frozen gate (v7, 27 Feb 2026; v7plus revD/revE) | PASS needs all four: ΔBIC > 0 for Bessel-K vs logistic and SBPL; fingerprint residual RMS and ν-flatness inside a preregistered knee window; κ50_band_med < 4.0; stability-sweep PASS fraction ≥ 0.6 | **FAIL / FAIL** on three pilots (7 Mar 2026) |
| Empirical route after April 2026 | `v431 / WAIT_FOR_RETURN`: no authoritative IPTA_DR2B_TOP5 HD export, no official return data used, no detection claim (v8; v20 "empirical firewall") | OPEN |

### 6b. Model-internal registered predictions (eras 2–3)

| Prediction | Test performed | Outcome |
|---|---|---|
| Tachyon growth rate 1.65719 H (certified μ²) | Recovered by two evolution/spectral codes at t = 10⁻³: 1.65714–1.65719 (Chat 14 v2 code), 1.657193663767 (folder-152 replacement solver) | PASSED (NUMERICAL calibration; the target was known in advance) |
| 4D effective theory predicts the throat turnaround at φ_b = −1.942 and the crossing H_J = 0.797 H₀ | 5D: −1.9412 to −1.9417 and 0.798 | PASSED (NUMERICAL) |
| 4D predicts the +1 endpoint H = 0.606 H₀; the exact RS brane gives 0.6067 | 5D: 0.63–0.64 and still falling at the chart freeze | OPEN (consistent with relaxation, not demonstrated) |
| 4D predicts a second shell branch and the stability window (d₀ = 1.1135; merger near d = 1.438) | 12 shells found in 5D; window located at d ∈ (1.11, 1.12) and merger between 1.43 and 1.44 | PASSED (NUMERICAL) |
| A hot universe from the registered shell's roll-off | Full 5D nonlinear evolution | **FAILED** (no radiation branch) |
| The proposed late-time endpoint is constant φ = 1 | Scalar junction check | **FAILED**; replaced by the scalar-profile branch |

### 6c. What would count against or for the hypothesis next (from the 22 Sept checkpoint, section 5)

The next discriminating steps are the perturbation spectrum of the +1 branch; constraint-controlled late-time evolution; an explicit matter coupling with renormalized production fed back through all three junctions; and a demonstrated sustained radiation era **before** any observational prediction.

## 7. Zenodo records located (HDBLAST line)

| # | Record | Date | Label / content | Evidence file |
|---|---|---|---|---|
| 0–18 | 17088133 … 20313949 | 2025-09-09 → 2026-05-20 | See era 1 table | Aug 29 audit |
| 19 | 21367209 | 2026-07-14 | v20 correction with partial withdrawal | folder 89 ZENODO metadata; Aug 29 audit |
| 20 | 21998599 | 2026-08-18 | v21, M468R1 | Aug 29 audit (the folder 126 kit predates the record and says no record was yet created) |
| – | 22285737 | 2026-09-03 | v22, M489G-A | DB3/`untitled folder 137/M489GB_ZENODO_RECOMMENDATION_AND_PLAIN_LANGUAGE_SUMMARY.md` |
| – | 22287013 | 2026-09-03 | v23, M489G-C | Chat 7 transcript ("Your HDBLAST v23 (M489G-C) is now published!") |
| – | 22347452 | 2026-09-05 | v24, conditional full-domain certification | Chat 7 transcript; APS package replay |
| – | 22922928 | 2026-09-23 | Scalar junction consistency and corrected static de Sitter branch | GPD site data only (see section 3) |

Companion concepts: 16866883 (brane-world two-link), 16929972 (HDBC), 17069899 (pipeline).

Other record numbers found in D-Blast 3 are either the author's other programs or external datasets. Examples: 22319172 TOE-N00AK-r1, 22399940 DMDE, 22400470 AntiMatter, 22400638 HRF (horizon-response), 22398093 (a reserved draft, not verified published), and 10344086 (the NANOGrav 15-yr KDE dataset).

## 8. Index of key files and scripts

### Latest checkpoint (22 Sept 2026 → record 22922928)

- Overview and guides: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/00_READ_FIRST.md` (identical copy in DB3/untitled folder 153), with `REPRODUCE.md`, `FINAL_CLAIM_REVIEW.md` and `RESULTS_SUMMARY.md` alongside.
- Package verifier: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/verify_package.py`
- Frozen solver and seeds: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/frozen/registered_solver.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/frozen/balanced_constraint_seed.py`
- Static branch: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/static_branch/solve_plus_branch.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/static_branch/independent_branch_checks.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/static_branch/PLUS_BRANCH_RESULTS.json`
- Matter: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/verify_matter_extension.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/analyze_crossing_eligibility.py`
- Evolution: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/evolution/evolve_balanced.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/evolution_review/EVOLUTION_RESIDUAL_REVIEW.md`
- Chat 14 audit: `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/source_audit/recompute_chat14.py`, `checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/source_audit/CHAT14_RECOMPUTED_RESULTS.json`

### Shell dynamics (Chats 9–14; folder 152)

- Chat 9: `../new-files/D-Blast 3/untitled folder 146/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/REPLAY.py`; `../new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/certified_instability/REPLAY.py`; `../new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/background/hdblast_background.py`; `../new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/effective_theory_orchestrator/HJ_EFFECTIVE_THEORY_DERIVATION.md`; `../new-files/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/blast_dynamics/BLAST_DYNAMICS_REPORT.md`
- Chat 10: `../new-files/D-Blast 3/untitled folder 147/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 147/HDBLAST_CHAT10_STABLE_SHELL_20260918/existence/REPLAY.py`; `../new-files/D-Blast 3/untitled folder 147/HDBLAST_CHAT10_STABLE_SHELL_20260918/stability_certificate/REPLAY.py`; `../new-files/D-Blast 3/untitled folder 147/HDBLAST_CHAT10_STABLE_SHELL_20260918/oscillation_theorem/OSCILLATION_THEOREM.md`
- Chat 11: `../new-files/D-Blast 3/untitled folder 148/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 148/HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/REPLAY.py`; `../new-files/D-Blast 3/untitled folder 148/HDBLAST_CHAT11_SYMBOLIC_VERIFICATION_20260918/negative_controls.py`
- Chat 12: `../new-files/D-Blast 3/untitled folder 149/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 149/HDBLAST_CHAT12_ALL_SECTORS_20260919/REPLAY.py`
- Chat 13: `../new-files/D-Blast 3/untitled folder 150/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/rolloff5d.py`; `../new-files/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/derive_evolution_equations.py`
- Chat 14: `../new-files/D-Blast 3/untitled folder 151/00_READ_FIRST.md`; `../new-files/D-Blast 3/untitled folder 151/HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/rolloff5d_v2.py`; `../new-files/D-Blast 3/untitled folder 151/HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/analyse_v2.py`; `../new-files/D-Blast 3/untitled folder 151/HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/referee_numerics/ref_analyse.py`
- Registered-shell controls: `../new-files/D-Blast 3/untitled folder 152/00_READ_FIRST.md` and `../new-files/D-Blast 3/untitled folder 152/HDBLAST_REGISTERED_SHELL_CONTROLS_20260922.zip`. The ZIP contains `solver/registered_solver.py`, `numerics_review/verify_constraint_transport.py`, `initial_data/INITIAL_DATA_AND_CORNER_OBSTRUCTION.md`, `literature/PROPER_CLOCK_PROOF.md` and `verify_package.py`.

### Linear-response and M-series

- `../new-files/D-Blast 3/untitled folder 120/M462_MAIN_REPORT_RESEND.md` (registered shell root)
- `../new-files/D-Blast 3/untitled folder 126 ZENODO/HDBLAST_ZENODO_V21_METADATA.txt` (v21)
- `../new-files/D-Blast 3/untitled folder 136/00_READ_FIRST_PHYS_M489GA_20260903.md` (v22 content)
- `../new-files/D-Blast 3/untitled folder 138 ZENODO/PHYS_M489GC_CERTIFIED_RANK_TWO_REGISTERED_LINEAR_RESPONSE.md` (v23 content)
- `../new-files/D-Blast 3/untitled folder 140/HDBLAST_CONDITIONAL_GLOBAL_ADJOINT_RESPONSE_CERTIFICATION_20260904.zip` (v24 archive)
- `../new-files/D-Blast 3/untitled folder 143/00_READ_FIRST.md` and `../new-files/D-Blast 3/untitled folder 144/00_READ_FIRST.md` (Chat 8)
- `../new-files/D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/inputs/PREVIOUS_GLOBAL_FORWARD.zip` (global forward certificate)
- `../new-files/D-Blast 3/untitled folder 142 new chat/HDBLAST_CHAT7_COMPLETE_TRANSCRIPT_AND_PROGRESS_20260906.md`
- `../new-files/D-Blast 3/untitled folder 145 APS JOURNAL/HDBLAST-v24-APS-review-package/HDBLAST-v24-scientific-verification.md`
- v20 withdrawal: `../new-files/D-Blast 3/untitled folder 89 ZENODO/HDBLAST_ZENODO_V20_CORRECTION_AND_CLAIM_BOUNDARY.md`; M309–M312 origin: `../new-files/D-Blast 3/untitled folder 68/HDBLAST_NEW_CHAT_HANDOFF_PHYS_M285_M312_20260710.md`

### Observational / PTA knee

- `../new-files/D-Blast 3/ZENODO-PUBLICATIONS/Zenodo-18158842-jan-06-2026/HDBLAST_RESULTS_NANOGrav15yr_HD_20260106/README.md` and `../new-files/D-Blast 3/ZENODO-PUBLICATIONS/Zenodo-18158842-jan-06-2026/HDBLAST_free_spectrum_knee_estimator_20260106.py`
- `../new-files/D-Blast 3/ZENODO-PUBLICATIONS/Zenodo-18157426-jan-05-2026/HDBLAST_ZENODO_UPLOAD_BUNDLE_REALONLY_v27_20260105_vKERNELFAMILY_RELEASE1/HDBLAST_nanograv_kdegrid_killswitch_v3_20260103.py`
- `../new-files/D-Blast 3/ZENODO-PUBLICATIONS/ZENODO-18341882-Version-01-22-26/ZENODO_DESCRIPTION_HDBLAST_DBLAST_20260122.md`
- `../new-files/D-Blast 3/DBLAST FILES PURGE/Dblast 181/HDBLAST_fscale_vs_f50_clarification_20260122.md`
- `../new-files/D-Blast 3/D-BLAST-PROGRESS-3-FOLDERS/untitled folder 89/HDBLAST_v7plus_Authoritative_LastMile_20260306/README_CURRENT_STATE_20260306.md`
- `../new-files/D-Blast 3/DBLAST FILES PURGE/Dblast 571/HDBLAST_GO_NO_GO_MEMO_20260307.md`

### Provenance and registries

- `../new-files/D-Blast 3/VECTOR STORE/DARK MATTER DARK ENERGY VECTOR STORE/GPD-ARCHITECT-EXPLORE-ZENODO-AUDIT-2026-08-29.md` (Zenodo chain, with claim boundaries)
- `../new-files/D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 2 vector store/HDBLAST_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md` (SBPL equations and screening numbers)
- `../new-files/D-Blast 3/GPD-site-deploy/gpd5-articles.json` (public-site summaries mentioning 22922928)
- `../new-files/D-Blast 3/HDBLAST_GENERAL_PUBLIC_GUIDE_V20.md`

## 9. Known inconsistencies in the record

- **Version labels.** Internal labels ("v27", "v28", "REALONLY v22/v24/v26", "V180"), release labels and Zenodo version numbers overlap. Use immutable record IDs. Package 18341882 embeds the concept DOI as if it were its own (audit item LC-004).
- **Two assistant lines.** "Chat 7/Chat 8" (6 Sept) and the continuation workspace `CHAT8_CONTINUATION_20260922` belong to one line of work; Chats 9–14 belong to another. Both work on the same registered model.
- **Publication evidence for v22–v24 and 22922928** comes from transcripts and site data, not from a live API check made for this map.
- **Folder 153** holds only `00_READ_FIRST.md` and `source_audit/` beside its ZIP. The full extraction is under `new-files/latest-work`.
- **The Chat 13 ZIP** omits five intermediate files listed in its own manifest; they are intact in the extracted folder (folder 152 audit).
- **Two slightly different ρ₀ values** give the RS benchmark 0.6067265 (78.82817714 in the static-branch review; 78.82819001 in a Chat 14 referee run). The difference is below 2×10⁻⁷ in relative terms and does not affect any conclusion.

## 10. Verification of this map

`verify_program_map.py` writes `program_map_checks.json`. It recomputes:

- c from I₊;
- the Chat 9 closed-form μ² and its O(t) slope against the certified interval;
- the growth rate from μ² (control: s = √(−μ²) is rejected);
- the static-branch identities: H²ρ_b² = 1; the metric-only H² equals the exact RS formula; the expansions, tested by residual slopes of 2.00 and 3.00 (controls: coefficients 9c/32 and −c²/192 are detected);
- H/H₀ = ρ₀/ρ_b = 0.606721732 and the −7.87 ppm correction;
- the Chat 14 turnaround and +1 values from the 22 Sept recompute JSON;
- the ordering and containment of the linear-response determinant intervals;
- x50 of the K₂ filter (control: a K₁-type filter gives 0.761).

It also:

- checks 34 date and phrase anchors in the source files;
- rebuilds the 21-version Zenodo chain;
- locates the four later records (control: the fabricated record number 22922929 is not found);
- confirms that every absolute path quoted in this file exists (control: a fabricated path is reported missing).

Result on 27 Sept 2026: see `README.md` in this folder.
