# HDBLAST literature sweep 1: theory

Research round of 27 September 2026, written 28 September 2026. **Work in progress.** This is an annotated bibliography with a short synthesis. It adds no new physics. Nothing here is a claim of the HDBLAST programme until it appears in a dated checkpoint or a Zenodo record.

## Question

What published work (2022–2026, plus foundational papers) bears on the HDBLAST hypothesis that a five-dimensional gravitational event sparked our Big Bang? Specifically: which results support, constrain or challenge the registered model and its current status? That status is: a de Sitter shell in a Z2-doubled Einstein–scalar bulk; a certified tachyonic roll-off with two fates (relaxation to an empty Randall–Sundrum de Sitter brane, or reversal and collapse); no radiation era; and the corrected static "+1 branch" of 22 Sept 2026.

## How the sweep was done, and what "read" means here

- **111 distinct web searches** were run on 28 Sept 2026. All are listed in the appendix. They cover:
  - dark bubbles (de Sitter from decaying AdS);
  - brane-world creation and catalysed nucleation;
  - ekpyrotic, cyclic and colliding-bubble cosmologies;
  - Randall–Sundrum and Karch–Randall de Sitter branes;
  - stability of de Sitter-sliced walls;
  - Einstein–scalar de Sitter-sliced flows;
  - holographic and brane reheating;
  - swampland constraints;
  - 2025–2026 work on "a 5D event creating the universe".
- **No paper was read in full.** Pages could not be opened on arxiv.org, zenodo.org, link.springer.com, uu.diva-portal.org or www.uu.se; a fetch attempt returned `EGRESS_BLOCKED`. Every finding below paraphrases a **search-result snippet**. Snippets are machine-written summaries and can be wrong. Where two snippets conflicted, or a summary could not be attributed to one paper, the entry carries a caveat.
  - Authors are named only when a snippet named them.
  - A few years and journal details came from general knowledge; each such case is flagged.
  - Every URL cited appeared in the result list of the query or queries named in the entry.
- Each entry has an **implication** label:
  - `supports`: consistent with, or a precedent for, an HDBLAST mechanism or result;
  - `constrains`: sets a bound or requirement HDBLAST must meet;
  - `challenges`: a result that could undercut an HDBLAST ingredient;
  - `method`: a technique or benchmark HDBLAST could use;
  - `context`: background or prior art.
- One small calculation was added (`theory_checks/rs_thin_brane_crosscheck.py`). It places the registered model inside the standard literature formulas and is a consistency check, not a new result.

Totals: **95 annotated sources.** By implication: supports 7, constrains 8, challenges 4, method 41, context 35. By period: 2022–2026: 31, 2010–2021: 23, before 2010: 41.

**Relation to the programme's own earlier literature record.**
- The Chat 9 package (16 Sept 2026) contains `literature/PRIOR_ART_SHELL_STABILITY.md`. That report read papers through ar5iv at the time. It already covers: Frolov–Kofman; BraneCode; the bulk-inflaton papers; Koyama–Takahashi; DeWolfe–Freedman–Gubser–Karch; Garriga–Sasaki; Gen–Sasaki; and the dark-bubble papers up to June 2026.
- The Chat 9 moving-shell analysis already uses the Kraus/Ida moving-wall method and the Maeda–Wands projected equations.
- `theory_checks/prior_mention_scan.py` searched the programme files for each source's arXiv number and title: 98 Markdown/text files in D-Blast 3 folders 146–152, the DBlast 4 vector-store handoffs and the 22 Sept package. Result: **26 of 95 sources were already mentioned there, and 69 were not found.** Each entry below says which.
- The material that is new to the programme is concentrated in themes B, C, F and G:
  - the Einstein–scalar de Sitter-sliced flow programme of Kiritsis, Nitti and collaborators ("Kiritsis" does not occur in any Markdown, text, Python or JSON file under D-Blast 3; zip archives were not searched);
  - catalysed and black-hole-seeded creation;
  - collision and crunch outcomes;
  - holographic reheating;
  - the 2025–2026 detuned-tension inflation paper;
  - swampland context.

## Main findings

1. **The general idea has substantial prior art (context).** Our 4D universe as a wall in 5D, born in a 5D gravitational event, has been published many times. The events proposed include:
   - a smooth 5D instanton (B3, 1998);
   - brane-world creation from nothing, including a thermal instanton with a bulk AdS black hole that "may describe the creation of a hot universe" (B1, 2000);
   - ekpyrotic brane collisions (C1, 2001);
   - colliding 5D bubbles (C8–C10, 2001–2003);
   - a bulk true-vacuum bubble hitting an inflating, tension-detuned brane and reheating it (C11, 2001/2002);
   - a 5D black hole or white hole producing the brane universe (B7, 2013).

   The line is active in 2025–2026:
   - black-hole-catalysed nucleation of a dark bubble, in which the black hole's matter becomes the universe's matter (A11, 2025);
   - dark bubbles realising the dark dimension (A12, June 2026);
   - inflation driven by a detuned brane tension that makes the fifth dimension grow (G6, 2025/2026).

   HDBLAST should therefore not claim novelty for the concept. Its own content is the registered model and its exact, certified and numerical results. A search for HDBLAST or Maldonado found no web-indexed item (query 108).

2. **Where the registered model sits (EXACT consistency check, new to this folder, not new physics).** The model has three identifications with standard literature:
   - **A detuned wall of the DeWolfe–Freedman–Gubser–Karch type (D9).** The first-order flow φ′=W_φ, A′=−W/3 solves the Einstein equations computed from the metric, and φ=−tanh y is the balanced flat wall.
   - **Balanced tension equals the critical tension (D1, D2).** 2W equals the critical Randall–Sundrum/Karch–Randall tension 6/ℓ at both AdS vacua: 2/3 at φ=+1 (ℓ=9) and 10/3 at φ=−1 (ℓ=9/5).
   - **Garriga–Sasaki geometry (B1).** The static shells use dS-sliced AdS closed off by a regular cone, with ρ′(0)=1.

   Two consequences follow:
   - **H² agrees with the thin-brane relation up to a small scalar-profile term.** The standard thin-brane de Sitter relation, H²=(σ/6)²−1/ℓ² (D3; Nihei; Kim and Kim), with σ evaluated at φ=1, reproduces the checkpoint's expansion exactly at O(δ). At O(δ²) it differs by exactly the scalar-profile term −c²δ²/384. Numerically, it equals the checkpoint's `metric_only_H2` to a relative difference of 2.3×10⁻¹⁴.
   - **The 7.87 ppm correction is the scalar profile's effect.** The full scalar-profile branch lies 7.869 ppm below the thin-brane H. The second-order formula predicts 7.849 ppm; the difference is O(δ³). This confirms that the checkpoint's 7.87 ppm correction is the scalar profile's departure from the textbook Randall–Sundrum de Sitter brane.

   These are anchoring checks only. The checkpoint already reported the thin-brane value as `metric_only_H2` and quoted "exact RS 0.6067", and the −c²/384 term is part of its published expansion. What this folder adds is the explicit identification with the published formula and the exact separation. Four deliberate wrong variants all fail, as required: a one-sided junction σ/3, the φ=−1 radius ℓ=9/5 used with the φ=+1 tension, the wrong-sign flow A′=+W/3, and a 1% change in c.

3. **The HDBLAST negative result fits the literature (context, supports the negative result).** Every hot-origin model found adds an ingredient beyond a pure-tension wall:
   - kinetic energy of colliding walls (C11, C12);
   - a bulk black hole (B1 thermal instanton, D5, D6, A11);
   - a holographic hot sector (F6);
   - matter placed on the wall (A3, A7).

   A pure-tension de Sitter brane is empty. In the registered model the bulk is conformally flat, so the projected Weyl term vanishes (Chat 9 moving-shell analysis): there is no dark radiation either. HDBLAST's finding that the vacuum roll-off relaxes to an empty Randall–Sundrum de Sitter brane, with no radiation branch, is therefore the expected outcome for this ingredient set, not an anomaly.

4. **The collapse fate has a concrete literature hypothesis (to be tested).** Several results point the same way:
   - Self-gravitating colliding walls generically form a black brane with the walls trapped behind a horizon (C13).
   - AdS bubble interiors crunch behind an apparent horizon even in potentials unbounded below (B17). HDBLAST's U is unbounded below, with U ~ −(2/27)φ⁶ (checked).
   - Symmetric bubble collisions collapse unless the prior expansion rate is large compared with the momentum transfer in the fifth dimension (C10).

   No trapped-surface or apparent-horizon diagnostic was found in the Chat 13–14 or folder-152 files (text search of `.md`, `.py` and `.json` files for "apparent horizon" and "trapped surface").

5. **Stability: the registered tachyon is typical, and the +1 branch is not covered by any theorem found (constrains).** Tachyonic moduli of de Sitter branes are generic:
   - two-brane radion (E2);
   - "stabilized" inflating branes, confirmed with BraneCode (E3, E12);
   - bulk-scalar power-law branes (C15).

   A single pure-gravity de Sitter brane has no radion (E2), so HDBLAST's tachyon must come from the bulk scalar and its tension coupling, as its effective theory says. The Chat 9 prior-art report already reached this conclusion; this sweep re-confirms it. The fake-supergravity stability theorem covers flat and AdS-sliced walls, but the snippet does not mention de Sitter-sliced walls (E6). Stability of the +1 branch therefore has to be computed; the round's two stability workstreams are doing this. The Kaluza–Klein continuum threshold (3H/2)² = 9H²/4 (B1, E10) agrees with the threshold used in Chat 10. Coupling brane fields to bulk fields can make an infinite ladder of tachyonic modes normalizable (E8). The matter extension therefore needs a coupled stability analysis.

6. **The most directly usable toolset is the Einstein–scalar de Sitter-sliced flow programme (method).** The relevant papers are F1–F5, B10 and B11:
   - existence and classification of regular de Sitter-sliced flows;
   - Coleman–de Luccia bounces read as flows on de Sitter;
   - competition between two branches decided by free energy, with a phase transition as the curvature changes (F4);
   - first-order Hamilton–Jacobi numerics for de Sitter-foliated walls (B11);
   - a no-go theorem for smooth flows from an AdS boundary into a de Sitter interior (F3). This constrains any HDBLAST variant with a positive-energy "blast region" joined smoothly to the AdS bulk.

7. **Challenges for the matter extension.**
   - **Equivalence principle.** In the dark bubble, bulk-mediated couplings give the proton different gravitational and inertial masses, severely violating equivalence-principle tests (A10, 2025). HDBLAST's proposed χ field has a φ-dependent mass, so it couples directly to the bulk scalar. An equivalence-principle and fifth-force estimate is required.
   - **Dark radiation.** Any bulk-Weyl ("dark") radiation is limited by nucleosynthesis and the CMB to −0.41 ≤ ρ_d/ρ_γ ≤ 0.105 (D12). It cannot supply the Standard Model radiation era.

8. **Swampland context (constrains any UV embedding; not a test of the current model).**
   - A late-time empty de Sitter phase would face the de Sitter conjecture and trans-Planckian censorship (G1, G2, G4, G9).
   - The instability of non-supersymmetric AdS that dark-bubble creation needs is contested (G8, B14).
   - HDBLAST has no string embedding, so these are context only.

9. **Observables produced by this class of model.** Each is a reminder that HDBLAST has none yet:
   - micron-scale weakening of gravity and positive spatial curvature (A11, A12);
   - suppressed, blue-tilted large-scale power and oscillatory tensor modes (G6);
   - Kaluza–Klein continuum signatures in primordial correlators (E10);
   - suppressed metric perturbations except on the largest angular scales (B2);
   - a 2.7σ tension with a simple power law (B8).

## Results table (script-backed)

All values come from `theory_checks/rs_thin_brane_crosscheck.json` (SymPy exact algebra; mpmath at 30 digits). The comparison input is the checkpoint's `static_branch/PLUS_BRANCH_RESULTS.json`, read-only, SHA-256 `b03e1c5569afb3567873fa990d16af9321b06e12a80ed89f9088a1bc00cd4fc7`.

| Check | Result | Label |
|---|---|---|
| AdS radius at φ=+1 and φ=−1 (from U=−6/ℓ²) | ℓ = 9 and 9/5 | EXACT |
| Balanced tension 2W equals the critical tension 6/ℓ | 2/3 = 2/3; 10/3 = 10/3 | EXACT |
| dS-sliced AdS5 with ρ=ℓ sinh(y/ℓ) is Einstein with R_AB=−(4/ℓ²)g_AB; regular cone ρ′(0) | True; ρ′(0) = 1 | EXACT |
| First-order flow φ′=W_φ, A′=−W/3 solves R_AB=φ_Aφ_B+(2/3)U g_AB | True | EXACT |
| Flat kink φ=−tanh y solves φ′=W_φ | True | EXACT |
| U for large φ (either sign) | U ~ (-2/27) φ⁶, unbounded below | EXACT |
| Thin-brane H²=(σ(1)/6)²−1/81 minus checkpoint expansion, through O(δ²) | checkpoint − thin = -c**2*delta**2/384 | EXACT |
| Thin-brane H² vs checkpoint `metric_only_H2` at δ=10⁻³ | 0.00005924108026 vs 5.9241080267048024e-05; relative difference 2.3082e-14 | NUMERICAL |
| H/H₀: thin-brane vs full +1 branch | 0.606726506259 vs 0.606721731972 | NUMERICAL |
| Shift of H, full vs thin (ppm); second-order prediction −c²δ²/384 | -7.86893 vs -7.84928 | NUMERICAL |
| KK continuum threshold (units of H²) | 9/4 | EXACT (literature value) |
| Control: one-sided junction σ/3 matches? | False | control, expected False |
| Control: φ=−1 radius with φ=+1 tension matches? | False | control, expected False |
| Control: wrong-sign flow A′=+W/3 solves Einstein eqs? | False | control, expected False |
| Control: c perturbed by 1%, relative difference from checkpoint | 0.0037311 | control, expected ≫10⁻¹² |
| All checks and controls | True |  |

## Opportunities suggested by the literature (ranked by cost and value)

1. **Add a bulk black-hole (Weyl) parameter behind the shell.**
   - What: a static and roll-off version with an AdS–Schwarzschild mass behind the shell (D5, D6, and the B1 thermal instanton). This is the cheapest literature-backed "hot" ingredient.
   - Background: the programme already has the formalism, since Chat 9 notes the μ/a⁴ term of the AdS–Schwarzschild case, and the 22 Sept matter extension carries the Weyl scalar 𝒲 as an unclosed quantity. What the literature adds is that a bulk black hole is the standard route to a hot brane (B1, D6, A11).
   - Tests: whether the +1 branch persists, what a⁻⁴ term it gives, and whether the D12 bounds can be met.
   - Scope: this is a new model branch under project rules.
2. **Compare the Euclidean on-shell actions of the two static branches.** The two branches are the original shell (cone near −1) and the +1 branch (cone near +1). Methods: B1, B9, F4. This would say which configuration a nucleation event prefers.
3. **Search the Chat 14 collapse runs for trapped surfaces and an apparent horizon** (C13, B17). Also evaluate the C10 criterion: expansion rate against momentum transfer in the fifth dimension.
4. **Calibrate the solver externally.** Test the HDBLAST 5D solver on the exactly solvable Koyama–Takahashi model (E9), or reproduce a BraneCode test (E12). Either is a cheap external calibration control.
5. **Coupled brane–bulk stability for the matter extension** (E8), plus an equivalence-principle and fifth-force estimate for the φ-dependent χ mass (A10).
6. **Classify the +1 branch within the de Sitter-sliced flow framework** (F1, B10), with an independent first-order Hamilton–Jacobi existence check (B11, E7).
7. **Consider a kinetic heating channel.** Test a wall-collision or bulk-bubble heating channel (C11, C12) as a separate model variant, calibrated against T_R ≈ 0.88 g N_b^{1/4}.
8. **Communication.** Any Zenodo or GitHub text should cite the prior art in finding 1 and describe HDBLAST as a specific registered model within an established class.

## Annotated bibliography

Each entry gives the title; the authors (if the snippet named them); the year; the venue or identifier; the URL; what the snippet says; relevance to HDBLAST; the implication label; the access level; and the query numbers. Themes:

- A: dark bubbles.
- B: creation, nucleation and catalysis.
- C: collisions, ekpyrotic and cyclic models.
- D: Randall–Sundrum and Karch–Randall foundations, and brane cosmology.
- E: stability of de Sitter walls and branes.
- F: Einstein–scalar de Sitter-sliced flows and holographic reheating.
- G: swampland and string context.

### A. Dark bubbles: de Sitter on a wall from decaying AdS5

**A1. Emergent de Sitter Cosmology from Decaying Anti-de Sitter Space** (2018)
- Authors: S. Banerjee, U. Danielsson, G. Dibitetto, S. Giri, M. Schillo
- Venue or identifier: Phys. Rev. Lett. 121, 261301; arXiv:1807.01570
- URL: <https://arxiv.org/abs/1807.01570>
- Access: search-snippet-only (queries 1). Implication: **context**.
- *Finding (from snippet):* Proposes that positive-energy FLRW cosmology arises on a brane (bubble wall) that mediates the decay of a non-supersymmetric false AdS5 vacuum to a true vacuum; 4D gravity is confined on the brane and a 4D observer sees an effective positive cosmological constant coupled to matter and radiation, without scale separation or a fundamental de Sitter vacuum.
- *Relevance to HDBLAST:* The best-developed published precedent for 'our universe is an expanding codimension-one wall in 5D AdS with de Sitter induced on it'. Unlike HDBLAST it has different AdS5 vacua on the two sides (no Z2 doubling) and the creation event is a quantum nucleation, not a classical roll-off of a registered shell.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A2. de Sitter Cosmology on an expanding bubble** (2019)
- Authors: S. Banerjee, U. Danielsson, G. Dibitetto, S. Giri, M. Schillo
- Venue or identifier: JHEP 10 (2019) 164; arXiv:1907.04268
- URL: <https://arxiv.org/abs/1907.04268>
- Access: search-snippet-only (queries 26). Implication: **context**.
- *Finding (from snippet):* Develops the 4D cosmology on a codimension-one bubble wall separating two AdS5 vacua, relates the scenario to the weak gravity conjecture, and interprets the de Sitter temperature as the Unruh temperature of an accelerated observer riding the bubble.
- *Relevance to HDBLAST:* The Unruh-temperature reading applies equally to the HDBLAST static shells (temperature H/2pi of the dS slicing). It is a vacuum temperature of an empty de Sitter brane, not a radiation era, so it does not remove the HDBLAST negative result.
- *Already in the programme record:* not found in the scanned programme files

**A3. Dark bubbles: decorating the wall** (2020)
- Authors: S. Banerjee, U. Danielsson, S. Giri
- Venue or identifier: JHEP 04 (2020) 085; arXiv:2001.07433
- URL: <https://arxiv.org/abs/2001.07433>
- Access: search-snippet-only (queries 58). Implication: **method**.
- *Finding (from snippet):* Adds matter and radiation to the dark bubble, treating backreaction in the bulk and on the brane through brane bending, and computes the backreacted metric on the bent brane and in the 5D bulk.
- *Relevance to HDBLAST:* A worked template for the sector HDBLAST lacks: brane matter that backreacts consistently on a 5D wall geometry. The 22 Sept HDBLAST matter extension has the same goal but has not been evolved.
- *Already in the programme record:* not found in the scanned programme files

**A4. Curing with hemlock: escaping the swampland using instabilities from string theory** (2021)
- Authors: S. Banerjee and colleagues (full author list not in snippet)
- Venue or identifier: arXiv:2103.17121 (snippet names Int. J. Mod. Phys. D)
- URL: <https://arxiv.org/abs/2103.17121>
- Access: search-snippet-only (queries 90). Implication: **context**.
- *Finding (from snippet):* Argues that a dark bubble with strings attached can carry our universe out of the swampland, using the instabilities of apparently hostile corners of the swampland; dark energy and possibly other dark components are presented as features of higher-dimensional physics.
- *Relevance to HDBLAST:* States the dark-bubble response to swampland objections. HDBLAST has no string embedding, so it can neither use nor be excluded by this argument yet.
- *Already in the programme record:* not found in the scanned programme files

**A5. Gravitational waves in dark bubble cosmology** (2022)
- Authors: U. Danielsson and co-authors (snippet)
- Venue or identifier: Phys. Rev. D 106, 024002; arXiv:2202.00545
- URL: <https://arxiv.org/abs/2202.00545>
- Access: search-snippet-only (queries 59, 61). Implication: **method**.
- *Finding (from snippet):* Constructs the 5D uplift of 4D gravitational waves in de Sitter cosmology on a nucleated bubble in AdS5, extends the link between dark bubbles and Vilenkin's quantum cosmology to gravitational perturbations, and explains apparently negative energy contributions in the 4D Einstein equations that distinguish the dark bubble from Randall-Sundrum.
- *Relevance to HDBLAST:* Method for uplifting tensor perturbations on a nucleated wall. A related snippet (query 61) states that bubble-nucleation boundary conditions select the Vilenkin weighting; that is the kind of weighting any HDBLAST 'creation event' would need.
- *Already in the programme record:* not found in the scanned programme files

**A6. Features of a dark energy model in string theory** (2023)
- Authors: S. Banerjee, U. Danielsson, S. Giri
- Venue or identifier: Phys. Rev. D 108, 126009; arXiv:2212.14004
- URL: <https://arxiv.org/abs/2212.14004>
- Access: search-snippet-only (queries 103). Implication: **context**.
- *Finding (from snippet):* Clears up misconceptions about the dark bubble: points out important differences from Randall-Sundrum and explains why gravity neither is, nor needs to be, localized on the dark bubble.
- *Relevance to HDBLAST:* HDBLAST is Randall-Sundrum-like: a Z2 doubled bulk closed off by a regular cone, so the bulk volume is finite. Its 4D gravity mechanism therefore differs from the dark bubble's, and results cannot be transferred between the two without care.
- *Already in the programme record:* not found in the scanned programme files

**A7. Shedding light on dark bubble cosmology** (2023)
- Authors: I. Basile, U. Danielsson, S. Giri, D. Panizo
- Venue or identifier: JHEP 02 (2024) 112; arXiv:2310.15032
- URL: <https://arxiv.org/abs/2310.15032>
- Access: search-snippet-only (queries 45). Implication: **method**.
- *Finding (from snippet):* Incorporates electromagnetic fields: worldvolume fields backreact on the ambient 5D universe, changing the energy-momentum distribution and the effective gravity induced on the brane, and the resulting 4D cosmology consistently contains electromagnetic waves.
- *Relevance to HDBLAST:* Shows that radiation on a wall must source bulk fields for consistency. In HDBLAST the only bulk fields are the metric and phi, so a radiation era would need an analogous consistent bulk source (for example a bulk black hole or Weyl term).
- *Already in the programme record:* not found in the scanned programme files

**A8. Experimental tests of dark bubble cosmology** (2024)
- Authors: U. Danielsson, D. Panizo (snippet)
- Venue or identifier: Phys. Rev. D 109, 026003; arXiv:2311.14589
- URL: <https://link.aps.org/doi/10.1103/PhysRevD.109.026003>
- Access: search-snippet-only (queries 2, 9). Implication: **context**.
- *Finding (from snippet):* The snippet gives only the title, authors and journal; no result is summarized here.
- *Relevance to HDBLAST:* Listed so that the testability line of the dark bubble programme is not missed; content not assessed.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A9. The dark bubbleography** (2024)
- Authors: S. Banerjee, U. Danielsson, M. Zemsch
- Venue or identifier: JHEP 02 (2024) 102; arXiv:2311.16242
- URL: <https://arxiv.org/abs/2311.16242>
- Access: search-snippet-only (queries 8). Implication: **method**.
- *Finding (from snippet):* Presents the holographic construction of the dark bubble and shows, following holographic renormalization, that non-normalizable modes are essential for a vanishing induced graviton mass in any braneworld model; applies this to the propagator on the wall.
- *Relevance to HDBLAST:* The claim is stated for 'any braneworld model'. HDBLAST's compact regular-cone bulk has a normalizable zero mode, and its tensor sector was found mode-stable (Chat 12). How these two statements fit together should be checked before HDBLAST claims 4D graviton behaviour.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A10. Dark bubble cosmology and the equivalence principle** (2025)
- Authors: I. Basile, A. Borys, J. Masias
- Venue or identifier: arXiv:2507.03748; Phys. Rev. D (2026)
- URL: <https://arxiv.org/abs/2507.03748>
- Access: search-snippet-only (queries 6, 2). Implication: **challenges**.
- *Finding (from snippet):* Couples the electroweak and strong sectors to the induced braneworld gravity by the dark-bubble mechanism. The electroweak sector is unaffected, but the proton's gravitational and inertial masses differ significantly, severely violating equivalence-principle measurements.
- *Relevance to HDBLAST:* A sharp negative test that any braneworld with bulk-mediated matter couplings must pass. HDBLAST's proposed matter field chi has a phi-dependent mass, so the bulk scalar couples directly to brane matter; an equivalence-principle and fifth-force check is required before that extension can be called viable.
- *Caveat:* Snippets disagree on the journal reference: one gives Phys. Rev. D 113 (2026) issue 2, another wrongly attaches Phys. Rev. D 109, 026003 (which belongs to A8). The programme's Chat 9 prior-art file records Phys. Rev. D 113, 026009.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A11. Weak gravity at micron scales from dark bubble cosmology and its cosmological consequences** (2025)
- Authors: U. Danielsson, S. Giri
- Venue or identifier: arXiv:2511.21362; Phys. Rev. D (2026)
- URL: <https://arxiv.org/abs/2511.21362>
- Access: search-snippet-only (queries 7, 9). Implication: **context**.
- *Finding (from snippet):* Predicts that gravity becomes weaker, not stronger, than Newtonian at about 1e-5 m, with explicit table-top predictions. The same effect lowers effective gravity at high energy densities and gives early inflation with nothing beyond radiation. The paper also discusses a quantum origin of the universe in which a 5D black hole of critical size catalyses nucleation of the dark bubble at its horizon, and the black hole's matter becomes the matter on the bubble.
- *Relevance to HDBLAST:* This is the closest recent published proposal to 'a five-dimensional event created our universe'. It also supplies what HDBLAST lacks: a source of matter, the catalysing black hole. HDBLAST cannot claim the general concept as new and should cite this line of work.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A12. Dark bubbles, dark dimensions and fat gravitons** (2026)
- Authors: U. Danielsson, S. Giri
- Venue or identifier: arXiv:2606.20942 (18 June 2026)
- URL: <https://arxiv.org/abs/2606.20942>
- Access: search-snippet-only (queries 5, 1). Implication: **context**.
- *Finding (from snippet):* Uses the instabilities behind the de Sitter swampland conjectures to make accelerated expansion inevitable. The dark bubble is presented as a realization of the dark-dimension proposal and of Sundrum's fat graviton. Predictions: a micron-size dark dimension, gravity fading at micron distances, a string scale of tens of TeV, and a measurable positive spatial curvature.
- *Relevance to HDBLAST:* Shows the kind of falsifiable output this class of model can produce. HDBLAST has no such prediction yet. Positive spatial curvature is natural for a closed de Sitter slicing, which HDBLAST's static shells also use; this is a possible future observable, but it is not derived for HDBLAST.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A13. Self-gravitating electromagnetic waves in the dark bubble model** (2026)
- Authors: U. Danielsson and one co-author (snippet)
- Venue or identifier: arXiv:2606.16547 (15 June 2026)
- URL: <https://arxiv.org/abs/2606.16547>
- Access: search-snippet-only (queries 55). Implication: **method**.
- *Finding (from snippet):* Embeds gravitational and electromagnetic waves using AdS5 pp-wave geometries glued across a three-brane. For localized light beams under mixed AdS5 boundary conditions, the gravitational corrections are consistent with 4D gravity weakening at the 5D AdS scale. On the brane, electromagnetic waves source the Kalb-Ramond B-field in the bulk.
- *Relevance to HDBLAST:* Exact wall-plus-radiation solutions of this type could serve as analytic benchmarks for any HDBLAST radiation-on-shell calculation.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**A14. Dynamical dark energy in 0'B braneworlds** (2025)
- Authors: I. Basile and two co-authors (snippet)
- Venue or identifier: Eur. Phys. J. C (Sept 2025); arXiv:2502.20438
- URL: <https://arxiv.org/abs/2502.20438>
- Access: search-snippet-only (queries 31). Implication: **challenges**.
- *Finding (from snippet):* Builds a dark-bubble realization in non-supersymmetric type 0'B string theory, the unique scale-separated option among the simplest models. The dark energy varies logarithmically. The blue-shifted spectrum excludes inflation, the late-time predictions conflict with Standard Model couplings, and the model appears not to be phenomenologically viable.
- *Relevance to HDBLAST:* A concrete example of a 5D-wall cosmology failing when fully specified, published as a negative result. It sets the standard of candour HDBLAST follows.
- *Already in the programme record:* not found in the scanned programme files

### B. Creation, nucleation, catalysis and bubbles of nothing

**B1. Brane-world creation and black holes** (2000)
- Authors: J. Garriga, M. Sasaki
- Venue or identifier: Phys. Rev. D 62, 043523; arXiv:hep-th/9912118
- URL: <https://arxiv.org/abs/hep-th/9912118>
- Access: search-snippet-only (queries 12, 37). Implication: **method**.
- *Finding (from snippet):* An inflating brane-world can be created from nothing with its AdS bulk. The spatial sections are compact and bounded by the brane, and the de Sitter brane instanton has its interior filled with AdS. The paper discusses Nariai-like pair creation of 'black cigars', and finds that 5D and 4D actions and entropies agree when the instanton is much larger than the AdS radius. Thermal instantons with an AdS black hole in the bulk may describe the creation of a hot universe from nothing.
- *Relevance to HDBLAST:* HDBLAST's static geometry is this instanton's Lorentzian section with a bulk scalar added: a dS-sliced AdS region closed off by a regular cone, Z2-doubled, with one brane. The thermal-instanton variant is a documented route to a hot initial state through a bulk black hole. HDBLAST has not tried it.
- *Already in the programme record:* yes, in 5 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/instanton_orchestrator/SHELL_INSTANTON_ENTROPY_IDENTITY.md`

**B2. Brane New World** (2000)
- Authors: S. W. Hawking, T. Hertog, H. S. Reall
- Venue or identifier: Phys. Rev. D 62, 043501; arXiv:hep-th/0003052
- URL: <https://arxiv.org/abs/hep-th/0003052>
- Access: search-snippet-only (queries 28). Implication: **context**.
- *Finding (from snippet):* A Randall-Sundrum domain wall in AdS with a strongly coupled large-N CFT on it. The conformal anomaly acts as an effective tension that gives the wall de Sitter geometry, the analogue of Starobinsky inflation. The graviton correlator is computed from the no-boundary proposal, and the CFT strongly suppresses metric perturbations on all but the largest angular scales.
- *Relevance to HDBLAST:* An alternative quantum origin of a de Sitter brane in AdS, with a holographic matter sector and a perturbation prediction. It is a benchmark for any HDBLAST perturbation spectrum (open item 5).
- *Already in the programme record:* not found in the scanned programme files

**B3. Smooth 'creation' of an open universe in five dimensions** (1998)
- Authors: J. Garriga
- Venue or identifier: arXiv:hep-th/9804106
- URL: <https://arxiv.org/abs/hep-th/9804106>
- Access: search-snippet-only (queries 51). Implication: **context**.
- *Finding (from snippet):* A non-singular 5D instanton for creating an open universe with a compact extra dimension. Its 4D section is a singular Hawking-Turok instanton, and the 'singularity' is a smooth 5D bubble of nothing. Flat space with a compact extra dimension is gravitationally metastable, but long-lived if the extra dimension is much larger than the Planck length.
- *Relevance to HDBLAST:* Conceptual precedent: a smooth higher-dimensional event can look like a 4D singularity. This supports the logic of looking for a 5D cause of a 4D Big Bang, but it is a different mechanism.
- *Already in the programme record:* not found in the scanned programme files

**B4. Catalytic creation of a bubble universe induced by quintessence in five dimensions** (2020)
- Authors: I. Koga, Y. Ookouchi
- Venue or identifier: Phys. Rev. D 104, 126015 (2021); arXiv:2011.07437
- URL: <https://arxiv.org/abs/2011.07437>
- Access: search-snippet-only (queries 47). Implication: **context**.
- *Finding (from snippet):* Studies 5D bubble nucleation catalysed by quintessence, in the decay of a metastable Minkowski vacuum to AdS, and the dynamics of the bubble that carries a 4D expanding universe. Also discusses the trans-Planckian censorship conjecture and proposes mechanisms for inflation and dark energy.
- *Relevance to HDBLAST:* Precedent for catalysed creation of a 4D universe on a 5D bubble. The related 2019 paper 'Catalytic Creation of Baby Bubble Universe with Small Positive Cosmological Constant' (arXiv:1909.03014) appeared in the same result list.
- *Already in the programme record:* not found in the scanned programme files

**B5. dS4 universe emergent from Kerr-AdS5 spacetime: bubble nucleation catalyzed by a black hole** (2022)
- Authors: not captured in snippet
- Venue or identifier: JHEP 05 (2023) 107; arXiv:2209.05625
- URL: <https://arxiv.org/abs/2209.05625>
- Access: search-snippet-only (queries 9, 10). Implication: **context**.
- *Finding (from snippet):* Studies nucleation of a vacuum bubble in Kerr-AdS5, sufficient conditions for nucleation with a rotating black hole, and how the black hole changes the transition rate; the bubble carries an emergent dS4 universe.
- *Relevance to HDBLAST:* Quantitative black-hole-catalysed creation rates for a dS4 wall in AdS5. It is the closest rate calculation to a '5D event creating our universe'.
- *Already in the programme record:* not found in the scanned programme files

**B6. Dark bubbles and black holes** (2021)
- Authors: not reliably captured in snippet
- Venue or identifier: JHEP 09 (2021) 158
- URL: <https://link.springer.com/article/10.1007/JHEP09(2021)158>
- Access: search-snippet-only (queries 10, 58). Implication: **context**.
- *Finding (from snippet):* Only the title and venue are reliable from the snippets. The accompanying summary text mixed content from several dark-bubble and catalysis papers (black-hole- and string-cloud-catalysed vacuum decay; charged Nariai black holes on the dark bubble), so no specific finding is attributed here.
- *Relevance to HDBLAST:* Part of the black-hole-plus-dark-bubble line (see A11, B5).
- *Caveat:* Snippet summary was not attributable to this paper alone.
- *Already in the programme record:* not found in the scanned programme files

**B7. Out of the White Hole: A Holographic Origin for the Big Bang** (2013)
- Authors: R. Pourhasan, N. Afshordi, R. B. Mann
- Venue or identifier: JCAP 04 (2014) 005; arXiv:1309.1487
- URL: <https://arxiv.org/abs/1309.1487>
- Access: search-snippet-only (queries 11). Implication: **context**.
- *Finding (from snippet):* In a DGP braneworld (4D induced gravity plus 5D bulk gravity), the universe emerges as a spherical 3-brane from the formation of a 5D Schwarzschild black hole. The holographic fluid's pressure singularity lies inside the white-hole horizon and need not be physical. A thermal atmosphere at about 20% of the 5D Planck mass can induce scale-invariant curvature perturbations without inflation.
- *Relevance to HDBLAST:* A prominent published model of a 5D gravitational event producing the Big Bang, with its own perturbation mechanism. It is major prior art for the HDBLAST headline idea.
- *Already in the programme record:* not found in the scanned programme files

**B8. Cosmological Perturbations in the 5D Holographic Big Bang Model** (2017)
- Authors: N. Altamirano, E. Gould, N. Afshordi, R. B. Mann
- Venue or identifier: arXiv:1703.00954
- URL: <https://arxiv.org/abs/1703.00954>
- Access: search-snippet-only (queries 57, 93). Implication: **constrains**.
- *Finding (from snippet):* Computes the exact curvature power spectrum from a thin atmosphere accreting onto the 3-brane. It is scale-invariant on small scales, red on intermediate scales and blue beyond the atmosphere height, and is marginally disfavoured against a simple power law (2.7 sigma). A 2025 snippet (query 93) reports a best-fit nucleation temperature at least three orders of magnitude above the 5D Planck mass, and a 2025 book chapter, 'Cosmic Bubbles and the 5D Holographic Big Bang Model'.
- *Relevance to HDBLAST:* Shows how a 5D-origin model is confronted with CMB data, and that doing so can produce tension. HDBLAST has no spectrum yet (open item 5).
- *Already in the programme record:* not found in the scanned programme files

**B9. Nucleation of de Sitter from the anti de Sitter spacetime in scalar field models** (2024)
- Authors: M. Cadoni, M. Pitzalis, A. P. Sanna
- Venue or identifier: Eur. Phys. J. C (2025); arXiv:2407.10469
- URL: <https://arxiv.org/abs/2407.10469>
- Access: search-snippet-only (queries 76). Implication: **method**.
- *Finding (from snippet):* In Einstein-scalar gravity, gravitational coupling can drive nucleation of de Sitter from AdS through a static, spherically symmetric, metastable scalar lump. Euclidean actions and free energies are compared: the AdS lump is generally less favoured, and the most preferred state is a de Sitter vacuum.
- *Relevance to HDBLAST:* A worked Euclidean-action comparison in Einstein-scalar gravity. The same comparison between HDBLAST's two static branches (cone near phi=-1 versus near phi=+1) has not been done and is cheap.
- *Already in the programme record:* not found in the scanned programme files

**B10. Revisiting Coleman-de Luccia transitions in the AdS regime using holography** (2021)
- Authors: J. K. Ghosh, E. Kiritsis, F. Nitti, L. T. Witkowski
- Venue or identifier: JHEP (Sept 2021); arXiv:2102.11881
- URL: <https://arxiv.org/abs/2102.11881>
- Access: search-snippet-only (queries 106). Implication: **method**.
- *Finding (from snippet):* Coleman-de Luccia AdS-to-AdS decays in Einstein-scalar theories are interpreted as vev-driven holographic RG flows of a QFT on de Sitter space. These flows do not exist for generic potentials; that is the holographic statement that gravity can stabilise false AdS vacua. Existence of the tunnelling solutions is tied to exotic RG flows, and explicit potentials admitting them are constructed.
- *Relevance to HDBLAST:* HDBLAST's static solutions belong to exactly this class: Einstein-scalar solutions with de Sitter slicing and a regular endpoint where the slice shrinks. The existence and classification criteria bear directly on the open existence and uniqueness question for the +1 branch.
- *Already in the programme record:* not found in the scanned programme files

**B11. Searching for Coleman-de Luccia bubbles in AdS compactifications** (2022)
- Authors: G. Dibitetto, N. Petri
- Venue or identifier: Phys. Rev. D 107, 046020 (2023); arXiv:2207.02172
- URL: <https://arxiv.org/abs/2207.02172>
- Access: search-snippet-only (queries 102). Implication: **method**.
- *Finding (from snippet):* Fully backreacted gravitational instantons are obtained by numerically integrating first-order Hamilton-Jacobi equations: smooth domain walls with de Sitter foliations connecting a supersymmetric and a non-supersymmetric AdS vacuum, showing a nonperturbative instability of the latter.
- *Relevance to HDBLAST:* The same mathematics as HDBLAST, a superpotential (Hamilton-Jacobi) description of dS-foliated walls, applied in truncations of string theory. It is a numerical method and a possible cross-check for HDBLAST's static-branch solver.
- *Already in the programme record:* not found in the scanned programme files

**B12. Nothing really matters** (2020)
- Authors: G. Dibitetto, N. Petri, M. Schillo
- Venue or identifier: JHEP 08 (2020) 040
- URL: <https://ui.adsabs.harvard.edu/abs/2020JHEP...08..040D/abstract>
- Access: search-snippet-only (queries 107). Implication: **context**.
- *Finding (from snippet):* Constructs gravitational instantons that generalise Witten's bubble of nothing to AdS with a cosmological constant in any dimension. Their expansion is described by a lower-dimensional de Sitter geometry within a non-compact foliation. Covariantly constant spinors are discussed as a possible topological obstruction to such decays.
- *Relevance to HDBLAST:* An alternative 'higher-dimensional blast': a bubble of nothing in AdS whose wall is de Sitter. The spinor obstruction matters if HDBLAST's W is ever derived from a supersymmetric theory.
- *Already in the programme record:* not found in the scanned programme files

**B13. Instability of the Kaluza-Klein vacuum** (1982)
- Authors: E. Witten
- Venue or identifier: Nucl. Phys. B 195, 481
- URL: <https://ui.adsabs.harvard.edu/abs/1982NuPhB.195..481W/abstract>
- Access: search-snippet-only (queries 71). Implication: **context**.
- *Finding (from snippet):* The Kaluza-Klein ground state is unstable to semiclassical barrier penetration (the bubble of nothing), because the positive-energy conjecture fails for Kaluza-Klein theory; elementary fermions would stabilise it.
- *Relevance to HDBLAST:* The foundational example of a higher-dimensional vacuum event with a de Sitter-like expanding wall.
- *Already in the programme record:* not found in the scanned programme files

**B14. Nonperturbative Instability of AdS5 x S5/Zk** (2007)
- Authors: G. T. Horowitz, J. Orgera, J. Polchinski
- Venue or identifier: Phys. Rev. D 77, 024004 (2008); arXiv:0709.4262
- URL: <https://arxiv.org/abs/0709.4262>
- Access: search-snippet-only (queries 66). Implication: **context**.
- *Finding (from snippet):* For a freely acting, supersymmetry-breaking orbifold there are no tachyons at large 't Hooft coupling, yet a nonperturbative bubble-of-nothing instability analogous to Kaluza-Klein vacuum decay exists.
- *Relevance to HDBLAST:* Non-supersymmetric AdS5 backgrounds are generically suspect in string theory. Any future string embedding of HDBLAST's bulk vacua would face this.
- *Already in the programme record:* not found in the scanned programme files

**B15. Gravitational Effects on and of Vacuum Decay** (1980)
- Authors: S. Coleman, F. De Luccia
- Venue or identifier: Phys. Rev. D (1980); journal details not in snippet
- URL: <https://www.semanticscholar.org/paper/Gravitational-Effects-on-and-of-Vacuum-Decay-Coleman-Luccia/a60b863946ca2d85da9d8b48b72257e8ec05e3ed>
- Access: search-snippet-only (queries 99). Implication: **context**.
- *Finding (from snippet):* As summarised in the snippet: gravity can prevent the decay of a Minkowski or AdS false vacuum; in the thin-wall limit this happens when the wall tension exceeds a bound proportional to the difference of the square roots of the two vacuum energy densities.
- *Relevance to HDBLAST:* HDBLAST's balanced tension sigma=2W is the critical, BPS-like tension at which a flat wall exists (checked exactly in theory_checks: 2W = 6/ell at both vacua). The registered detuning delta>0 is supercritical, which is why the static shells are de Sitter branes.
- *Caveat:* Year from general knowledge; the snippet gave title and authors only.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/God_Plays_Dice_Architect_Research_Handoff_and_100_Explore_Questions_2026-09-08.md`

**B16. Generalized surface tension bounds in vacuum decay** (2017)
- Authors: A. Masoumi, S. Paban, E. J. Weinberg
- Venue or identifier: arXiv:1711.06776
- URL: <https://arxiv.org/abs/1711.06776>
- Access: search-snippet-only (queries 101). Implication: **method**.
- *Finding (from snippet):* Some thin walls have a curvature radius that changes substantially across the wall, so the Coleman-de Luccia tension is ill-defined. The paper proposes a tension definition for this regime, shows it obeys a CDL-like bound, and derives a general bound for all bounces with Minkowski or AdS false vacua.
- *Relevance to HDBLAST:* Gives a principled way to define an effective shell tension when the scalar is nonconstant across the wall, as on HDBLAST's +1 branch.
- *Already in the programme record:* not found in the scanned programme files

**B17. Crunch from AdS bubble collapse in unbounded potentials** (2024)
- Authors: not captured in snippet
- Venue or identifier: arXiv:2411.07692 (Nov 2024)
- URL: <https://arxiv.org/abs/2411.07692>
- Access: search-snippet-only (queries 100). Implication: **method**.
- *Finding (from snippet):* Classical evolution of an AdS bubble nucleated from a Minkowski false vacuum in a potential with an infinitely deep true vacuum along a quartic slope. The interior collapses and forms a spacelike curvature singularity behind an apparent horizon. The scalar's kinetic energy overtakes the negative potential energy, the core density turns positive and trapped surfaces form. This requires no lower bound on the potential.
- *Relevance to HDBLAST:* HDBLAST's potential is also unbounded below (U ~ -(2/27) phi^6 at large |phi|, checked exactly in theory_checks). Its Chat 13 drafts wrongly blamed the reversal on this; the corrected cause is sub-balanced tension. This paper supplies the right diagnostic for the collapse fate: search for trapped surfaces and an apparent horizon.
- *Already in the programme record:* not found in the scanned programme files

**B18. End-of-the-World Branes and Inflationary Predictions for Rocky and Swampy Landscapes** (2024)
- Authors: B. Hassfeld, A. Hebecker, A. Westphal
- Venue or identifier: JHEP 03 (2025) 196; arXiv:2411.11944
- URL: <https://arxiv.org/abs/2411.11944>
- Access: search-snippet-only (queries 89). Implication: **constrains**.
- *Finding (from snippet):* A framework for anthropic predictions assuming de Sitter vacua have finite-dimensional Hilbert spaces. Even with eternal inflation, predictions depend on the rates of creating universes from nothing, and these rates are highly sensitive to the existence of end-of-the-world branes. Distinguishes 'swampy' from 'rocky' landscapes.
- *Relevance to HDBLAST:* Any claim that a 5D event is a typical origin of our universe depends on creation-rate measures. HDBLAST has not addressed measures.
- *Already in the programme record:* not found in the scanned programme files

**B19. Bubbles of cosmology in AdS/CFT** (2023)
- Authors: A. Sahu, P. Simidzija, M. Van Raamsdonk
- Venue or identifier: JHEP (Nov 2023); arXiv:2306.13143
- URL: <https://arxiv.org/abs/2306.13143>
- Access: search-snippet-only (queries 105). Implication: **context**.
- *Finding (from snippet):* The typical big-bang/big-crunch cosmologies of holographic effective theories are not asymptotically AdS, but arbitrarily large spherical bubbles of them can be embedded in asymptotically AdS spacetimes with a Schwarzschild-AdS exterior.
- *Relevance to HDBLAST:* A holographic framing in which a cosmological region is bounded by a shell with a black-hole exterior. It is a possible interpretation of a future HDBLAST shell with a bulk black-hole mass.
- *Already in the programme record:* not found in the scanned programme files

### C. Brane and bubble collisions; ekpyrotic and cyclic models

**C1. The Ekpyrotic Universe: Colliding Branes and the Origin of the Hot Big Bang** (2001)
- Authors: J. Khoury, B. A. Ovrut, P. J. Steinhardt, N. Turok
- Venue or identifier: Phys. Rev. D 64, 123522; arXiv:hep-th/0103239
- URL: <https://arxiv.org/html/hep-th/0103239v2>
- Access: search-snippet-only (queries 4). Implication: **context**.
- *Finding (from snippet):* The hot big bang is produced when a bulk brane collides with a bounding orbifold plane, starting from a cold, empty, static universe. The scenario is claimed to address the horizon, flatness and monopole problems and to give a nearly scale-invariant spectrum without inflation.
- *Relevance to HDBLAST:* Canonical prior art for 'a higher-dimensional event makes the hot Big Bang'. HDBLAST must cite it and state how it differs: a single Z2 shell in a registered Einstein-scalar bulk, with no second brane or orbifold plane.
- *Already in the programme record:* not found in the scanned programme files

**C2. M theory model of a big crunch/big bang transition** (2004)
- Authors: N. Turok, M. Perry, P. J. Steinhardt
- Venue or identifier: Phys. Rev. D 70, 106004
- URL: <https://ui.adsabs.harvard.edu/abs/2004PhRvD..70j6004T/abstract>
- Access: search-snippet-only (queries 82). Implication: **context**.
- *Finding (from snippet):* Treats a big crunch/big bang transition as the collision of two empty orbifold planes and shows that p-brane states winding around the extra dimension propagate smoothly across the collision.
- *Relevance to HDBLAST:* Shows how a 5D collision can be continued through; relevant if HDBLAST's collapse fate is ever continued past the singularity.
- *Already in the programme record:* not found in the scanned programme files

**C3. Cosmological perturbations in a big crunch/big bang space-time** (2004)
- Authors: A. J. Tolley, N. Turok, P. J. Steinhardt
- Venue or identifier: Phys. Rev. D 69, 106005; arXiv:hep-th/0306109
- URL: <https://arxiv.org/pdf/hep-th/0306109>
- Access: search-snippet-only (queries 52). Implication: **method**.
- *Finding (from snippet):* A prescription for matching general-relativistic perturbations across collision singularities of orbifold planes, of the kind met in ekpyrotic and cyclic scenarios.
- *Relevance to HDBLAST:* A matching method for perturbations across a 5D collision event.
- *Already in the programme record:* not found in the scanned programme files

**C4. Solution of a Braneworld Big Crunch/Big Bang Cosmology** (2005)
- Authors: P. McFadden, N. Turok, P. J. Steinhardt
- Venue or identifier: Phys. Rev. D 76, 104038 (2007); arXiv:hep-th/0512123
- URL: <https://arxiv.org/abs/hep-th/0512123>
- Access: search-snippet-only (queries 111). Implication: **constrains**.
- *Finding (from snippet):* Solves for perturbations around two separating or colliding boundary branes as an expansion in the collision speed (V/c). The 4D effective description fails at the first nontrivial order, (V/c)^2, where the growing and decaying 4D modes mix nontrivially.
- *Relevance to HDBLAST:* A warning for HDBLAST: its closed-form 4D effective theory agrees with 5D to about 1e-6 for the static tachyon, but may fail near fast events (turnaround, collapse). This supports keeping full 5D evolution as the reference.
- *Already in the programme record:* not found in the scanned programme files

**C5. Bouncing Negative-Tension Branes** (2007)
- Authors: J.-L. Lehners, N. Turok
- Venue or identifier: Phys. Rev. D 77, 023516 (2008); arXiv:0708.0743
- URL: <https://arxiv.org/pdf/0708.0743>
- Access: search-snippet-only (queries 82). Implication: **context**.
- *Finding (from snippet):* The snippet gives the bibliographic record only.
- *Relevance to HDBLAST:* Listed for completeness of the brane-collision line.
- *Already in the programme record:* not found in the scanned programme files

**C6. A new kind of cyclic universe** (2019)
- Authors: A. Ijjas, P. J. Steinhardt
- Venue or identifier: arXiv:1904.08022
- URL: <https://arxiv.org/abs/1904.08022>
- Access: search-snippet-only (queries 32). Implication: **context**.
- *Finding (from snippet):* Ekpyrotic contraction combined with a non-singular classical bounce: H, energy density and temperature oscillate while the scale factor grows exponentially from cycle to cycle. A companion snippet says the earlier colliding-brane models have given way to bouncing models built from ordinary 4D scalar fields.
- *Relevance to HDBLAST:* The ekpyrotic programme itself moved from 5D brane collisions to 4D bounces. This is context for how strong the 5D evidence would need to be.
- *Already in the programme record:* not found in the scanned programme files

**C7. Robustness of slow contraction to cosmic initial conditions** (2020)
- Authors: A. Ijjas, W. G. Cook, F. Pretorius, P. J. Steinhardt, G. N. Davies
- Venue or identifier: arXiv:2006.04999
- URL: <https://arxiv.org/abs/2006.04999>
- Access: search-snippet-only (queries 62). Implication: **method**.
- *Finding (from snippet):* Numerical-relativity simulations vary all freely specifiable initial data (shear, curvature, scalar field and velocity profiles), including data far outside the perturbative regime, and find slow contraction robust. A follow-up (arXiv:2104.12293) finds robustness enhanced with multiple modes and reduced symmetry.
- *Relevance to HDBLAST:* A methodological template for testing whether HDBLAST's two fates are robust to non-perturbative initial disturbances. The balanced two-bump seeds so far cover a narrow family.
- *Already in the programme record:* not found in the scanned programme files

**C8. A Braneworld Universe From Colliding Bubbles** (2001)
- Authors: M. Bucher
- Venue or identifier: Phys. Lett. B 530, 1 (2002); arXiv:hep-th/0107148
- URL: <https://ar5iv.arxiv.org/html/hep-th/0107148>
- Access: search-snippet-only (queries 95). Implication: **context**.
- *Finding (from snippet):* The Big Bang hypersurface is the collision locus of two bubbles that nucleate in 5D flat space with AdS inside; the resulting domain wall localizes gravity as in Randall-Sundrum.
- *Relevance to HDBLAST:* Prior art for a 5D nucleation-and-collision origin of an RS brane universe.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/HDBLAST_ARCHITECT_MASTER_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md`

**C9. Cosmological Perturbations Generated in the Colliding Bubble Braneworld Universe** (2001)
- Authors: J. J. Blanco-Pillado, M. Bucher
- Venue or identifier: Phys. Rev. D 65, 083517 (2002); arXiv:hep-th/0111089
- URL: <https://arxiv.org/abs/hep-th/0111089>
- Access: search-snippet-only (queries 95, 52). Implication: **method**.
- *Finding (from snippet):* Computes the perturbations in the colliding-bubble braneworld, where quantum fluctuations develop on the AdS-filled bubbles before they collide to form a Randall-Sundrum brane.
- *Relevance to HDBLAST:* A method for primordial perturbations from a 5D origin event.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/HDBLAST_ARCHITECT_MASTER_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md`

**C10. When do colliding bubbles produce an expanding universe?** (2003)
- Authors: J. J. Blanco-Pillado and colleagues (snippet); Gratton and Turok are acknowledged, not authors
- Venue or identifier: Phys. Rev. D 69, 103515 (2004); arXiv:hep-th/0306151
- URL: <https://arxiv.org/abs/hep-th/0306151>
- Access: search-snippet-only (queries 68, 67). Implication: **method**.
- *Finding (from snippet):* A symmetric collision leaves a collapsing universe on the final brane unless the bulk expansion rate just before the collision is large compared with the momentum transfer in the fifth dimension. That prior expansion can come from negative spatial curvature or a positive 5D cosmological constant. Thick-wall numerical collisions confirm the thin-wall result.
- *Relevance to HDBLAST:* A published analogue of HDBLAST's two fates, expansion versus reversal and collapse, with an explicit criterion (expansion rate against fifth-dimension momentum transfer) that could be evaluated on the Chat 14 trajectories.
- *Already in the programme record:* not found in the scanned programme files

**C11. Brane Big-Bang Brought by Bulk Bubble** (2001)
- Authors: not captured in snippet
- Venue or identifier: Phys. Rev. D 66, 023519 (2002); arXiv:hep-th/0110286
- URL: <https://arxiv.org/pdf/hep-th/0110286>
- Access: search-snippet-only (queries 96, 94). Implication: **supports**.
- *Finding (from snippet):* A brane first inflates because the bulk vacuum energy and the brane tension are mismatched, with the bulk in a false vacuum. Bulk false-vacuum decay nucleates a negative-energy true-vacuum bubble that expands, hits the brane and produces a hot Randall-Sundrum big-bang brane. Inflation ends because the bulk vacuum energy changes, and the bubble's kinetic energy heats the universe.
- *Relevance to HDBLAST:* Very close to the HDBLAST narrative: a detuned, inflating brane plus a bulk event that ends inflation and reheats it. It shows that a 5D-event reheating mechanism has been proposed before, and it offers a concrete heating channel (wall kinetic energy deposited on the brane) that HDBLAST's vacuum roll-off lacks.
- *Already in the programme record:* not found in the scanned programme files

**C12. Collision of Domain Walls and Reheating of the Brane Universe** (2004)
- Authors: Y. Takamizu, K. Maeda
- Venue or identifier: Phys. Rev. D 70, 123514; arXiv:hep-th/0406235
- URL: <https://arxiv.org/abs/hep-th/0406235>
- Access: search-snippet-only (queries 20, 94). Implication: **method**.
- *Finding (from snippet):* Particle production when two domain walls collide in 5D Minkowski space, proposed as the reheating mechanism of ekpyrotic or cyclic brane universes. As rendered in the snippet, the created energy density is about 20 g^4 N_b m^4 and the reheating temperature T_R is about 0.88 g N_b^(1/4), for coupling g, number of bounces N_b and wall mass scale m. A follow-up treats collisions in asymptotically AdS5 (arXiv:hep-th/0603076).
- *Relevance to HDBLAST:* A calibrated wall-collision reheating estimate against which HDBLAST's formal production proxy (22 Sept matter extension) could be compared, once physical scales are chosen.
- *Caveat:* The formula's symbols were garbled in the snippet ('20g^4 N_mn^4'); the reading above is a best interpretation.
- *Already in the programme record:* not found in the scanned programme files

**C13. Dynamics of colliding branes and black brane production** (2007)
- Authors: Y. Takamizu, H. Kudoh, K. Maeda
- Venue or identifier: Phys. Rev. D 75, 061304; arXiv:gr-qc/0702138
- URL: <https://arxiv.org/abs/gr-qc/0702138>
- Access: search-snippet-only (queries 53). Implication: **challenges**.
- *Finding (from snippet):* Self-gravitating colliding domain walls, with initial data from a BPS wall in 5D supergravity, generically form a spacelike curvature singularity covered by a horizon: a black brane with trapped walls. Exceptions are non-relativistic weak-field cases, where the walls pass through each other or bounce several times.
- *Relevance to HDBLAST:* Suggests a definite hypothesis for HDBLAST's collapse fate: black-brane formation with the shell trapped. It is testable with an apparent-horizon search on the Chat 14 runs.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS.md`

**C14. Born-Again Braneworld** (2002)
- Authors: S. Kanno, M. Sasaki, J. Soda
- Venue or identifier: Prog. Theor. Phys. 109, 357 (2003); arXiv:hep-th/0210250
- URL: <https://arxiv.org/abs/hep-th/0210250>
- Access: search-snippet-only (queries 34). Implication: **context**.
- *Finding (from snippet):* Two inflating branes collide and re-emerge with tensions of opposite sign. The radion acts as the gravitational scalar; in the Einstein frame the scenario resembles pre-big-bang cosmology; vacuum gravitational waves have a very blue spectrum.
- *Relevance to HDBLAST:* An example where radion (modulus) dynamics of inflating branes, rather than matter, drives the transition, as in HDBLAST's tachyonic roll-off.
- *Already in the programme record:* not found in the scanned programme files

**C15. Regular collision of dilatonic inflating branes** (2005)
- Authors: E. Leeper, K. Koyama, R. Maartens
- Venue or identifier: Phys. Rev. D 73, 043506 (2006); arXiv:hep-th/0508145
- URL: <https://arxiv.org/abs/hep-th/0508145>
- Access: search-snippet-only (queries 42). Implication: **supports**.
- *Finding (from snippet):* A two-brane system with a bulk scalar driving power-law inflation on the branes has a radion instability that can lead to collision; brane quantities such as the scale factor stay regular at the collision; a low-energy expansion reproduces the exact solution.
- *Relevance to HDBLAST:* Bulk-scalar inflating branes with a radion instability are established; HDBLAST's tachyonic shell mode is of this general type (in a single-shell, regular-cone geometry).
- *Already in the programme record:* not found in the scanned programme files

**C16. The big bang as a higher-dimensional shock wave** (2000)
- Authors: P. S. Wesson, H. Liu, S. S. Seahra
- Venue or identifier: Astron. Astrophys. 358, 425; arXiv:gr-qc/0003012
- URL: <https://arxiv.org/pdf/gr-qc/0003012>
- Access: search-snippet-only (queries 110, 108). Implication: **context**.
- *Finding (from snippet):* An exact solution of the 5D field equations describes a shock wave moving in time and the extra coordinate, suggesting the 4D big bang was a 5D shock wave. Different observers can measure different ages. The induced 4D spacetimes match standard FRW matter- and radiation-era models.
- *Relevance to HDBLAST:* The earlier work whose wording is closest to 'higher-dimensional blast'. It uses a different framework, the induced-matter (space-time-matter) Kaluza-Klein approach, not a brane in an Einstein-scalar bulk, and should be cited to avoid a novelty claim on the phrase or the idea.
- *Already in the programme record:* not found in the scanned programme files

**C17. Simulating the universe(s) II: phenomenology of cosmic bubble collisions in full General Relativity** (2014)
- Authors: Johnson, Peiris, Lehner and collaborators (per query context)
- Venue or identifier: arXiv:1407.2950 (also arXiv:1112.4487, arXiv:1508.03641)
- URL: <https://arxiv.org/pdf/1407.2950>
- Access: search-snippet-only (queries 43). Implication: **method**.
- *Finding (from snippet):* Full-GR simulations of bubble collisions, computing cosmological observables directly from the scalar Lagrangian; the first fully relativistic predictions for an ensemble of eternal-inflation models, which differ significantly from non-relativistic approximations.
- *Relevance to HDBLAST:* Methodology for turning a violent, collision-type origin event into CMB templates, if HDBLAST ever reaches that stage.
- *Caveat:* Exact author list per paper not given in snippet.
- *Already in the programme record:* not found in the scanned programme files

### D. Randall–Sundrum / Karch–Randall foundations and brane cosmology

**D1. An Alternative to Compactification** (1999)
- Authors: L. Randall, R. Sundrum
- Venue or identifier: Phys. Rev. Lett. 83, 4690
- URL: <https://ui.adsabs.harvard.edu/abs/1999PhRvL..83.4690R/abstract>
- Access: search-snippet-only (queries 35). Implication: **context**.
- *Finding (from snippet):* A single 3-brane in 5D AdS reproduces 4D Newtonian and general-relativistic gravity even without a Kaluza-Klein mass gap, because the warped bulk has finite volume.
- *Relevance to HDBLAST:* The base model that HDBLAST's late-time static branch reduces to. The +1 branch is an RS de Sitter brane plus a small scalar-profile correction.
- *Already in the programme record:* yes, in 5 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS.md`

**D2. Locally localized gravity** (2001)
- Authors: A. Karch, L. Randall
- Venue or identifier: JHEP 0105, 008
- URL: <https://iopscience.iop.org/article/10.1088/1126-6708/2001/05/008>
- Access: search-snippet-only (queries 3). Implication: **method**.
- *Finding (from snippet):* An AdS4 brane in AdS5 has a bound-state graviton when the AdS4 cosmological constant is small, and the graviton has an ultrasmall mass. The brane can be Minkowski, de Sitter or AdS depending on its tension.
- *Relevance to HDBLAST:* The tension classification behind HDBLAST: tension above the critical value gives a de Sitter brane. The theory_checks script verifies that the registered balanced tension 2W equals the critical value 6/ell at both AdS vacua, so delta(1+c phi)>0 is precisely the supercritical detuning.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/God_Plays_Dice_Architect_Research_Handoff_and_100_Explore_Questions_2026-09-08.md`

**D3. Bent Domain Walls as Braneworlds** (1999)
- Authors: N. Kaloper
- Venue or identifier: Phys. Rev. D 60, 123506; arXiv:hep-th/9905210
- URL: <https://arxiv.org/pdf/hep-th/9905210>
- Access: search-snippet-only (queries 19, 83). Implication: **method**.
- *Finding (from snippet):* Studies 3-brane domain walls in a curved 5D bulk that couple gravitationally to the bulk, including curved (de Sitter/AdS) walls, as braneworld models.
- *Relevance to HDBLAST:* One of the original derivations of the de Sitter brane in AdS5. A snippet for query 83 quotes the thin-brane relation H^2 = (4 pi G5 sigma/3)^2 + Lambda5/6, attributing the related inflation analyses to Nihei (Phys. Lett. B 465 (1999) 81) and Kim and Kim (Phys. Rev. D 61, 064003). No direct URL to those two papers was obtained. In HDBLAST units this relation is H^2 = (sigma/6)^2 - 1/ell^2, which reproduces the checkpoint's H^2 expansion through O(delta), and through O(delta^2) apart from the -c^2 delta^2/384 scalar-profile term (theory_checks).
- *Already in the programme record:* not found in the scanned programme files

**D4. The Einstein Equations on the 3-Brane World** (1999)
- Authors: T. Shiromizu, K. Maeda, M. Sasaki
- Venue or identifier: arXiv:gr-qc/9910076 (journal details not in snippet)
- URL: <https://arxiv.org/abs/gr-qc/9910076>
- Access: search-snippet-only (queries 98). Implication: **method**.
- *Finding (from snippet):* Derives the effective 4D Einstein equations on a Z2-symmetric brane in 5D Einstein gravity. The non-local bulk effect enters through the projected 5D Weyl tensor, a symmetric traceless tensor constrained by the Bianchi identities.
- *Relevance to HDBLAST:* The standard way to see the bulk Weyl ('dark radiation') and bulk-scalar contributions to HDBLAST's shell Friedmann equation.
- *Caveat:* Year inferred from the arXiv number.
- *Already in the programme record:* yes, in 5 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS.md`

**D5. Dynamics of Anti-de Sitter Domain Walls** (1999)
- Authors: P. Kraus
- Venue or identifier: JHEP 12 (1999) 011; arXiv:hep-th/9910149
- URL: <https://arxiv.org/abs/hep-th/9910149v1>
- Access: search-snippet-only (queries 64). Implication: **method**.
- *Finding (from snippet):* Moving domain walls in the Randall-Sundrum universe with the bulk patched together from AdS5 black-hole solutions; the junction equations fix the wall motion, which observers on the wall see as cosmological expansion or contraction.
- *Relevance to HDBLAST:* An exact way to add a bulk black-hole mass behind a shell and obtain an a^-4 ('dark radiation') term. It is the simplest literature-backed extension HDBLAST has not tried.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/moving_shell/MOVING_SHELL_FRIEDMANN_THEOREM.md`

**D6. CFT and Entropy on the Brane** (2001)
- Authors: I. Savonije, E. Verlinde
- Venue or identifier: Phys. Lett. B (2001) (not in snippet)
- URL: <https://www.semanticscholar.org/paper/CFT-and-Entropy-on-the-Brane-Savonije-Verlinde/7e51da58584a6054a76ef7e1d6c126e0c797e242>
- Access: search-snippet-only (queries 65). Implication: **supports**.
- *Finding (from snippet):* A brane with fine-tuned tension moving in an AdS black-hole background carries a radiation-dominated closed FRW universe, the radiation being the thermal CFT dual to the bulk black hole; when the brane crosses the horizon the Friedmann equation coincides with the Cardy(-Verlinde) entropy formula.
- *Relevance to HDBLAST:* Gives the holographic meaning of the a^-4 term from a bulk black hole (item D5). It is a candidate 'hot' ingredient for HDBLAST, with the caveat that this radiation is a hidden CFT sector, not Standard Model radiation.
- *Caveat:* Journal and year not given in snippet; year from general knowledge.
- *Already in the programme record:* not found in the scanned programme files

**D7. Gravity in the Randall-Sundrum Brane World** (1999)
- Authors: J. Garriga, T. Tanaka
- Venue or identifier: arXiv:hep-th/9911055
- URL: <https://arxiv.org/abs/hep-th/9911055>
- Access: search-snippet-only (queries 109). Implication: **method**.
- *Finding (from snippet):* A single positive-tension wall recovers Einstein gravity at leading order, with Kaluza-Klein corrections. With two branes of opposite tension, gravity is linearized Brans-Dicke; on the positive-tension wall the BD parameter exceeds 3000 if the separation is more than 4 AdS radii.
- *Relevance to HDBLAST:* Sets the benchmark for gravity on HDBLAST's late-time shell, and indicates what a scalar (radion-like) admixture would do to 4D tests.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/stability/independent_GN/GN_DERIVATION.md`

**D8. Brane-World Gravity** (2010)
- Authors: R. Maartens, K. Koyama
- Venue or identifier: Living Rev. Relativ. 13, 5; arXiv:1004.3962
- URL: <https://arxiv.org/abs/1004.3962>
- Access: search-snippet-only (queries 39). Implication: **context**.
- *Finding (from snippet):* A review of the geometry, dynamics and perturbations of warped 5D brane-worlds, mainly Randall-Sundrum based.
- *Relevance to HDBLAST:* A standard reference for conventions (junctions, Weyl term, perturbations) when comparing HDBLAST with the literature.
- *Already in the programme record:* yes, in 5 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS.md`

**D9. Modeling the fifth dimension with scalars and gravity** (1999)
- Authors: O. DeWolfe, D. Z. Freedman, S. S. Gubser, A. Karch
- Venue or identifier: Phys. Rev. D 62, 046008 (2000); arXiv:hep-th/9909134
- URL: <https://arxiv.org/pdf/hep-th/9909134>
- Access: search-snippet-only (queries 16). Implication: **method**.
- *Finding (from snippet):* A first-order (superpotential) method for 5D scalar-plus-gravity solutions, inspired by gauged supergravity but not requiring supersymmetry; applied to a full nonlinear treatment of interbrane stabilization; any effectively compactified fifth dimension contains a massless graviton.
- *Relevance to HDBLAST:* HDBLAST's registered model uses exactly this structure: U = (1/2)W_phi^2 - (2/3)W^2, phi' = W_phi, A' = -W/3, balanced tension 2W. The theory_checks script verifies that the first-order flow solves the Einstein equations from the metric, and that the kink phi = -tanh y solves phi' = W_phi. The registered model is therefore a detuned DeWolfe-Freedman-Gubser-Karch wall.
- *Already in the programme record:* yes, in 7 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md`

**D10. Cosmological evolution with brane-bulk energy exchange** (2002)
- Authors: E. Kiritsis, G. Kofinas, N. Tetradis, T. N. Tomaras, V. Zarikas
- Venue or identifier: JHEP 2003; arXiv:hep-th/0207060
- URL: <https://arxiv.org/abs/hep-th/0207060>
- Access: search-snippet-only (queries 40). Implication: **method**.
- *Finding (from snippet):* Energy exchange between brane and bulk produces a rich variety of brane cosmologies depending on the transfer mechanism, the equation of state and the spatial topology; an accelerating era is generic.
- *Relevance to HDBLAST:* Classification template for HDBLAST's proposed exchange law rho_dot + 3H(rho+p) = j phi_dot.
- *Already in the programme record:* not found in the scanned programme files

**D11. Reheating the Universe in Braneworld Cosmological Models with bulk-brane energy transfer** (2008)
- Authors: not captured in snippet
- Venue or identifier: arXiv:0805.1792
- URL: <https://arxiv.org/abs/0805.1792>
- Access: search-snippet-only (queries 24). Implication: **method**.
- *Finding (from snippet):* Analyses how the cosmological composition emerges (the reheating era) after inflation in 5D braneworld models with brane-bulk energy exchange.
- *Relevance to HDBLAST:* A direct precedent for the reheating calculation HDBLAST still needs (open item 4).
- *Already in the programme record:* not found in the scanned programme files

**D12. Observational Constraints on Dark Radiation in Brane Cosmology** (2002)
- Authors: not captured in snippet
- Venue or identifier: Phys. Rev. D 66, 043521; arXiv:astro-ph/0203272
- URL: <https://arxiv.org/abs/astro-ph/0203272>
- Access: search-snippet-only (queries 87). Implication: **constrains**.
- *Finding (from snippet):* Big Bang nucleosynthesis limits brane dark radiation (the bulk Weyl a^-4 term) to -1.23 <= rho_d/rho_gamma <= 0.11 at BBN; adding CMB constraints narrows this to -0.41 <= rho_d/rho_gamma <= 0.105. The same snippet also quotes a range of -12.1% to +6.2% of the total density at 10 MeV, possibly from another paper.
- *Relevance to HDBLAST:* Any HDBLAST 'hot' ingredient built from a bulk black hole or Weyl term must satisfy these bounds. By itself it cannot be the radiation era, which must be mostly Standard Model radiation.
- *Already in the programme record:* not found in the scanned programme files

### E. Stability of de Sitter-sliced walls and branes

**E1. Perturbations on domain walls and strings: A covariant theory** (1991)
- Authors: J. Garriga, A. Vilenkin
- Venue or identifier: Phys. Rev. D 44, 1007
- URL: <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007>
- Access: search-snippet-only (queries 33). Implication: **method**.
- *Finding (from snippet):* A covariant formalism for perturbations of walls and strings, including nucleated true-vacuum bubbles and walls in de Sitter space. The perturbation field often has a tachyonic mass and a nonminimal coupling to the worldsheet curvature. Ignoring gravity, one scalar describes the wall's translations. Local observers see perturbations grow; external observers see them freeze.
- *Relevance to HDBLAST:* Tachyonic masses on de Sitter walls can partly reflect wall translation, which is not a physical instability. HDBLAST's certified tachyon (mu^2 between -7.71788 and -7.71786 in units of H^2) was analysed with the full linearized equations in Chats 11-12, which also clarified the status of the 'mu^2 = -4' harmonics. This paper is the classic reference for separating wall-translation modes from physical instabilities.
- *Already in the programme record:* not found in the scanned programme files

**E2. Radion on the de Sitter brane** (2000)
- Authors: U. Gen, M. Sasaki
- Venue or identifier: arXiv:gr-qc/0011078
- URL: <https://arxiv.org/pdf/gr-qc/0011078>
- Access: search-snippet-only (queries 13). Implication: **supports**.
- *Finding (from snippet):* With two de Sitter branes of opposite tension, the radion has a negative mass squared proportional to the brane curvature. With a single positive-tension de Sitter brane there is no radion, and ordinary Einstein gravity is recovered with massive KK corrections.
- *Relevance to HDBLAST:* Explains why HDBLAST's tachyon cannot be a pure-gravity radion: HDBLAST has one shell, and the regular cone is not a second brane. The instability must come from the bulk scalar and its tension coupling, which agrees with the Chat 9 effective theory.
- *Already in the programme record:* yes, in 4 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**E3. Can Inflating Braneworlds be Stabilized?** (2003)
- Authors: A. Frolov, L. Kofman (from the programme's Chat 9 prior-art file, which read the paper; not in this snippet)
- Venue or identifier: arXiv:hep-th/0309002
- URL: <https://arxiv.org/abs/hep-th/0309002v1>
- Access: search-snippet-only (queries 18). Implication: **supports**.
- *Finding (from snippet):* For de Sitter branes the radion mass squared is typically negative, a strong tachyonic instability, so the parameters of 'stabilized' inflating braneworlds are constrained; BraneCode numerics confirm the instability; the flat-brane limit recovers the usual positive radion mass. A related result (Phys. Rev. D 67, 063515) stabilizes dS and AdS braneworlds with the Casimir force.
- *Relevance to HDBLAST:* In the literature, tachyonic moduli of de Sitter branes are generic. HDBLAST's registered-shell tachyon is typical rather than anomalous, and its stable shell S_8/5 (curved tension, d = 8/5) is an instance of evading the generic instability.
- *Already in the programme record:* yes, in 8 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/00_READ_FIRST.md`

**E4. Wave function of the radion in a brane world** (1999)
- Authors: C. Charmousis, R. Gregory, V. A. Rubakov
- Venue or identifier: Phys. Rev. D 62, 067505 (2000); arXiv:hep-th/9912160
- URL: <https://arxiv.org/abs/hep-th/9912160>
- Access: search-snippet-only (queries 50). Implication: **method**.
- *Finding (from snippet):* Computes the linearized radion in the two-brane Randall-Sundrum model and its couplings to matter on each brane, in agreement with Garriga-Tanaka.
- *Relevance to HDBLAST:* A reference construction of scalar metric modes; useful when interpreting any HDBLAST shell-scalar coupling to brane matter.
- *Already in the programme record:* yes, in 1 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/stability/independent_GN/GN_DERIVATION.md`

**E5. Thick Brane Worlds and Their Stability** (2001)
- Authors: S. Kobayashi, K. Koyama, J. Soda
- Venue or identifier: Phys. Rev. D 65, 064014 (2002); arXiv:hep-th/0107025
- URL: <https://arxiv.org/abs/hep-th/0107025>
- Access: search-snippet-only (queries 14). Implication: **method**.
- *Finding (from snippet):* Non-singular thick Poincare, de Sitter and AdS branes with a dilaton and potential; the effective potentials of the scalar-perturbation master equations are positive definite, so these systems are stable.
- *Relevance to HDBLAST:* A counterpoint: smooth thick de Sitter walls with no thin shell are stable, whereas HDBLAST's thin shell with detuned tension has a tachyon. The master-variable method is a candidate for the +1 branch stability test being run in this round.
- *Already in the programme record:* yes, in 4 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/GOD_PLAYS_DICE_ARCHITECT_RESEARCH_HANDOFF_AND_100_EXPLORE_QUESTIONS_20260908.md`

**E6. Fake supergravity and domain wall stability** (2004)
- Authors: D. Z. Freedman, C. Nunez, M. Schnabl, K. Skenderis
- Venue or identifier: Phys. Rev. D 69, 104027
- URL: <https://ui.adsabs.harvard.edu/abs/2004PhRvD..69j4027F/abstract>
- Access: search-snippet-only (queries 17). Implication: **method**.
- *Finding (from snippet):* A generalized Witten-Nester spinor argument proves stability of flat domain walls with no supersymmetry and in any dimension, and is extended to AdS-sliced walls (Janus is shown stable).
- *Relevance to HDBLAST:* Covers HDBLAST's balanced flat wall (delta=0) but, per the snippet, only flat and AdS-sliced walls: no comparable positive-energy theorem is cited for de Sitter-sliced walls. That fits the registered de Sitter shell having a tachyon, and it means stability of the +1 branch must be computed directly.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/VECTOR STORE/HDBLAST VECTOR STORE/DBlast 4 vector store/God_Plays_Dice_Architect_Research_Handoff_and_100_Explore_Questions_2026-09-08.md`

**E7. Hidden Supersymmetry of Domain Walls and Cosmologies** (2006)
- Authors: K. Skenderis, P. K. Townsend
- Venue or identifier: Phys. Rev. Lett. 96, 191301
- URL: <https://www.osti.gov/etdeweb/biblio/20777232>
- Access: search-snippet-only (queries 48). Implication: **method**.
- *Finding (from snippet):* Every domain-wall solution of gravity plus scalars with Minkowski or AdS worldvolume admits Killing spinors and first-order equations with a solution-determined superpotential; by analytic continuation, flat or closed FLRW cosmologies obey similar equations via pseudo-Killing spinors.
- *Relevance to HDBLAST:* A first-order (fake-superpotential) description of curved flows. HDBLAST's de Sitter-sliced static solutions need the curved or pseudo version; this is a possible route to an existence proof for the +1 branch.
- *Already in the programme record:* not found in the scanned programme files

**E8. Coupled bulk and brane fields about a de Sitter brane** (2006)
- Authors: A. Cardoso, K. Koyama, A. Mennim, S. S. Seahra, D. Wands
- Venue or identifier: Phys. Rev. D 75, 084002 (2007); arXiv:hep-th/0612202
- URL: <https://arxiv.org/abs/hep-th/0612202>
- Access: search-snippet-only (queries 97, 92). Implication: **challenges**.
- *Finding (from snippet):* A bulk scalar in AdS linearly coupled to a scalar on a de Sitter brane has zero, one or two bound states plus continuum resonances, depending on the masses and coupling. A light radion with m = sqrt(2) H satisfies the two-brane boundary conditions but is not normalizable with a single brane. With test-field matter on the brane, this mode, the zero mode and an infinite ladder of discrete tachyonic modes become normalizable.
- *Relevance to HDBLAST:* Directly relevant to HDBLAST's matter extension, which couples the brane field chi to the bulk phi. Adding brane fields can make tachyonic ladders normalizable, so a coupled-sector stability check is required before trusting any matter-driven trajectory.
- *Already in the programme record:* not found in the scanned programme files

**E9. Exactly solvable model for cosmological perturbations in dilatonic brane worlds** (2003)
- Authors: K. Koyama, K. Takahashi
- Venue or identifier: arXiv:hep-th/0307073 (also hep-th/0312087)
- URL: <https://arxiv.org/pdf/hep-th/0307073>
- Access: search-snippet-only (queries 38). Implication: **method**.
- *Finding (from snippet):* A bulk scalar with an exponential potential and an exponential coupling to the brane tension gives power-law inflation on a single brane. The background with scalar backreaction and the perturbations obeying the junction conditions are solved exactly.
- *Relevance to HDBLAST:* An exact benchmark for validating HDBLAST's 5D perturbation and evolution codes on a problem with a known answer, a cheap calibration control.
- *Already in the programme record:* yes, in 2 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/literature/PRIOR_ART_SHELL_STABILITY.md`

**E10. Primordial Correlators from a Kaluza-Klein Graviton Continuum** (2026)
- Authors: not captured in snippet
- Venue or identifier: arXiv:2608.01762 (3 Aug 2026)
- URL: <https://arxiv.org/abs/2608.01762>
- Access: search-snippet-only (queries 41, 37). Implication: **method**.
- *Finding (from snippet):* A spectral representation for inflationary correlators mediated by a continuum of states, realized in an RS2-like braneworld with the inflaton on a de Sitter brane in AdS5. The tensor sector contains a localized zero mode and a KK graviton continuum with a gap of order H; the query-37 snippet puts the threshold at 3H/2, the principal-series edge.
- *Relevance to HDBLAST:* HDBLAST's certified 'no scalar bound state below 9H^2/4' (Chat 10) uses exactly this threshold, (3H/2)^2. This paper shows how such a continuum could become an observable (cosmological-collider-like signals), if HDBLAST ever has an inflating stage with perturbations.
- *Already in the programme record:* not found in the scanned programme files

**E11. Complex scalar field thick branes: stability of linear perturbation and evolution of scalar Kaluza-Klein modes coupled with gravity** (2026)
- Authors: not captured in snippet
- Venue or identifier: arXiv:2609.21421 (18 Sept 2026)
- URL: <https://arxiv.org/abs/2609.21421>
- Access: search-snippet-only (queries 92, 49). Implication: **method**.
- *Finding (from snippet):* Minkowski, de Sitter and AdS thick branes from a complex scalar show no instability in the scalar and vector sectors. Factorizing the tensor equation excludes tachyonic tensor modes. The scalar zero mode is localized on Minkowski and de Sitter branes, with a gapless continuous scalar KK spectrum.
- *Relevance to HDBLAST:* A current example of a sector-by-sector stability analysis of de Sitter branes, comparable in scope to HDBLAST Chat 12.
- *Already in the programme record:* not found in the scanned programme files

**E12. BRANECODE: A Program for Simulations of Braneworld Dynamics** (2004)
- Authors: J. Martin, G. Felder, A. Frolov, M. Peloso, L. Kofman (per query; snippet names 'the authors you mentioned')
- Venue or identifier: arXiv:hep-ph/0404141; companion paper 'Brane world dynamics with the BraneCode', arXiv:hep-th/0309001 (also in the query-88 results)
- URL: <https://arxiv.org/abs/hep-ph/0404141v1>
- Access: search-snippet-only (queries 88, 18). Implication: **method**.
- *Finding (from snippet):* A C++ code for fully nonlinear evolution of 5D braneworlds with bulk scalar fields, described as the first of its type; it was used to confirm the tachyonic instability of 'stabilized' de Sitter branes (item E3).
- *Relevance to HDBLAST:* Prior numerical art for HDBLAST's rolloff5d evolutions. Reproducing a published BraneCode test case would calibrate HDBLAST's solver externally.
- *Already in the programme record:* yes, in 11 scanned programme file(s), e.g. `D-Blast 3/untitled folder 146/00_READ_FIRST.md`

### F. Einstein–scalar de Sitter-sliced flows; holographic reheating

**F1. Holographic RG flows on curved manifolds and quantum phase transitions** (2017)
- Authors: J. K. Ghosh, E. Kiritsis, F. Nitti, L. T. Witkowski
- Venue or identifier: JHEP 05 (2018) 034; arXiv:1711.08462
- URL: <https://arxiv.org/abs/1711.08462>
- Access: search-snippet-only (queries 73). Implication: **method**.
- *Finding (from snippet):* Einstein-dilaton flows with dS_d, AdS_d or S^d slicing and a general potential: UV and IR asymptotics, new flows at finite curvature with no flat counterpart, and bouncing flows that persist at finite curvature.
- *Relevance to HDBLAST:* A classification framework for HDBLAST's static de Sitter-sliced solutions. The +1 branch, whose cone sits near phi=+1, may be one of the finite-curvature flows with no flat counterpart.
- *Already in the programme record:* not found in the scanned programme files

**F2. De Sitter and Anti-de Sitter branes in self-tuning models** (2018)
- Authors: J. K. Ghosh, E. Kiritsis, F. Nitti, L. T. Witkowski
- Venue or identifier: JHEP 11 (2018) 128; arXiv:1807.09794
- URL: <https://arxiv.org/abs/1807.09794>
- Access: search-snippet-only (queries 72). Implication: **method**.
- *Finding (from snippet):* Maximally symmetric curved-brane solutions in dilatonic self-tuning braneworlds: no de Sitter or AdS brane vacua exist unless the near-boundary asymptotics are modified, after which curved solutions with a stabilized brane position generically exist.
- *Relevance to HDBLAST:* The closest published model class: a brane with scalar-dependent tension in an Einstein-dilaton bulk with de Sitter solutions. HDBLAST's bulk is compact (regular cone, no UV boundary), so the existence obstruction differs, but the matching-condition machinery carries over.
- *Already in the programme record:* not found in the scanned programme files

**F3. de Sitter versus Anti de Sitter flows and the (super)gravity landscape (and Part II)** (2019)
- Authors: E. Kiritsis, C. Tsouros (Part I); Part II authors per snippet: Kiritsis, Morales-Tejera, Rosen
- Venue or identifier: arXiv:1901.04546 (JHEP 08 (2023) 126 per listing); Part II arXiv:2510.12373, JHEP 04 (2026) 111
- URL: <https://arxiv.org/abs/1901.04546>
- Access: search-snippet-only (queries 69, 75). Implication: **constrains**.
- *Finding (from snippet):* Part I extends the classification of regular Einstein-scalar solutions to the de Sitter regime, including cosmic clocks that reverse direction without a curvature singularity, and finds no regular solutions interpolating between de Sitter and AdS extrema for generic potentials. Part II finds no regular ('Centaur') solutions interpolating between an AdS boundary and a de Sitter interior for d>2, even with several scalars and a nontrivial field-space metric.
- *Relevance to HDBLAST:* Constrains any HDBLAST variant in which a positive-energy 5D region, the 'blast', is joined smoothly to the AdS bulk. Such a junction would need a shell or a singularity, not a smooth scalar flow.
- *Caveat:* Journal attribution of Part I and author list of Part II come from search listings and were not cross-checked.
- *Already in the programme record:* not found in the scanned programme files

**F4. Holographic confining theories on space-times with constant positive curvature** (2025)
- Authors: J. Kastikainen, E. Kiritsis, F. Nitti
- Venue or identifier: JHEP 09 (2025) 139; arXiv:2502.04036
- URL: <https://arxiv.org/abs/2502.04036>
- Access: search-snippet-only (queries 77). Implication: **method**.
- *Finding (from snippet):* Two branches of solutions compete, with a phase transition as the curvature varies: the low-curvature phase has the flat-space IR geometry, the high-curvature phase a regular interior; the transition is first- or higher-order depending on the potential's leading exponent.
- *Relevance to HDBLAST:* HDBLAST also has two static branches (cone near -1 and near +1). The same free-energy (on-shell action) comparison would say which is preferred at the registered H, a cheap, well-defined next calculation.
- *Already in the programme record:* not found in the scanned programme files

**F5. Holographic self-tuning of the cosmological constant; Brane cosmology and the self-tuning of the cosmological constant** (2017)
- Authors: C. Charmousis, E. Kiritsis, F. Nitti (2017); A. Amariti, C. Charmousis, D. Forcella, E. Kiritsis, F. Nitti (2019)
- Venue or identifier: JHEP 09 (2017) 031, arXiv:1704.05075; JCAP 10 (2019) 007, arXiv:1904.02727; bulk-black-hole version arXiv:2003.05767
- URL: <https://arxiv.org/abs/1704.05075>
- Access: search-snippet-only (queries 74, 78). Implication: **method**.
- *Finding (from snippet):* A Standard Model brane separates an infinite-volume UV region from a finite-volume IR region in an Einstein-scalar bulk; for generic brane vacuum energy, regular flat-brane solutions exist. The cosmology of such branes, including matching conditions in several coordinate systems and bulk black holes, is worked out in the follow-ups.
- *Relevance to HDBLAST:* The time-dependent PDE-plus-junction problem is the same type as HDBLAST's roll-off; its coordinate choices and solution strategies are directly reusable.
- *Already in the programme record:* not found in the scanned programme files

**F6. A dynamical inflaton coupled to strongly interacting matter** (2023)
- Authors: C. Ecker, E. Kiritsis, W. van der Schee
- Venue or identifier: Phys. Rev. Lett. 130, 251001; arXiv:2302.06618
- URL: <https://arxiv.org/abs/2302.06618>
- Access: search-snippet-only (queries 84, 81). Implication: **supports**.
- *Finding (from snippet):* Self-consistently couples the Einstein-inflaton equations to a strongly coupled QFT described holographically, and finds inflation, a reheating phase, and finally a universe dominated by the QFT in thermal equilibrium.
- *Relevance to HDBLAST:* A controlled demonstration of reheating into a hot (holographic) sector, the step HDBLAST has not achieved. It suggests treating part of the 5D bulk as the hot sector's dual.
- *Already in the programme record:* not found in the scanned programme files

**F7. Gravitational reheating at strong coupling** (2023)
- Authors: A. Buchel
- Venue or identifier: JHEP 07 (2023) 159; arXiv:2304.11195
- URL: <https://arxiv.org/abs/2304.11195>
- Access: search-snippet-only (queries 86). Implication: **method**.
- *Finding (from snippet):* Uses gauge/gravity duality to estimate the maximal reheating temperature of strongly coupled theories after a rapid exit from de Sitter; reheating is most efficient when H is much larger than the conformal-breaking scale and the breaking operators are nearly marginal.
- *Relevance to HDBLAST:* HDBLAST's shell moves between de Sitter rates (H/H0 falling from about 0.64 to 0.607). A bound of this type could say how much gravitational particle production such a mild change can give, probably very little.
- *Already in the programme record:* not found in the scanned programme files

**F8. Self-sustained, out-of-equilibrium inflation** (2025)
- Authors: J. Casalderrey-Solana and colleagues (snippet)
- Venue or identifier: arXiv:2512.18079
- URL: <https://arxiv.org/abs/2512.18079>
- Access: search-snippet-only (queries 80, 81). Implication: **context**.
- *Finding (from snippet):* Holographic de Sitter-invariant states of non-conformal strongly coupled QFTs on dS4: out-of-equilibrium effects can sustain exponential inflation with H far below the QFT scale and the species scale; fine-tuning scales only logarithmically; apparent horizons with growing area signal growing comoving entropy; the regime can be the late-time limit of an FRW start.
- *Relevance to HDBLAST:* A de Sitter attractor with a coupled hot sector. It contrasts with HDBLAST's attractor, an empty RS de Sitter brane, and suggests one way a coupled matter sector could change the endpoint.
- *Already in the programme record:* not found in the scanned programme files

**F9. Brane Cosmology from AdS/BCFT** (2025)
- Authors: Fujiki, Kanda, Kohara, Takayanagi (surnames only, per snippet)
- Venue or identifier: JHEP 03 (2025) 135; arXiv:2501.05036
- URL: <https://arxiv.org/abs/2501.05036>
- Access: search-snippet-only (queries 29). Implication: **method**.
- *Finding (from snippet):* An end-of-the-world brane in AdS with a brane-localized scalar: the brane equation becomes a Friedmann-like equation, the model can describe creating a universe via a big bang, the near-hyperplane effective action is Liouville gravity with scalar matter, and a timelike g-theorem is proven from the null energy condition (AdS3/BCFT2).
- *Relevance to HDBLAST:* A recent formal analogue of a shell carrying a scalar in AdS that produces a big-bang-like creation; the null-energy-condition monotonicity is a structural tool HDBLAST could test on its shell trajectories.
- *Already in the programme record:* not found in the scanned programme files

**F10. Effective Dynamics of Inflationary End-of-the-World Branes in AdS3** (2026)
- Authors: K. Fujiki and others (snippet)
- Venue or identifier: arXiv:2609.11643 (Sept 2026)
- URL: <https://arxiv.org/abs/2609.11643>
- Access: search-snippet-only (queries 30). Implication: **method**.
- *Finding (from snippet):* Integrates out the bulk to obtain a Liouville-like effective theory for a 2D cosmological EOW brane with a localized scalar; constructs slow-roll inflating trajectories, regular Euclidean geometries continuing to Lorentzian de Sitter, and the semiclassical on-shell action.
- *Relevance to HDBLAST:* A method for deriving a brane-scalar effective action and a Euclidean continuation in a lower-dimensional toy model, comparable to HDBLAST's closed-form 4D effective theory (Chat 9).
- *Already in the programme record:* not found in the scanned programme files

**F11. Cosmology at the end of the world; Cosmology from the vacuum; Accelerating cosmology from a holographic wormhole** (2020)
- Authors: S. Antonini, B. Swingle (2020); S. Antonini, P. Simidzija, B. Swingle, M. Van Raamsdonk (2022-2024)
- Venue or identifier: Nature Phys. 16 (2020); arXiv:2203.11220 (CQG 41, 2024); Phys. Rev. Lett. 130, 221601 (2023)
- URL: <https://www.nature.com/articles/s41567-020-0909-6>
- Access: search-snippet-only (queries 27). Implication: **context**.
- *Finding (from snippet):* Cosmologies on end-of-the-world branes moving in charged AdS black-hole spacetimes, realized microscopically in AdS/CFT, and follow-ups obtaining cosmology, including accelerating cosmology, from holographic constructions.
- *Relevance to HDBLAST:* Microscopic (holographic) control of brane cosmologies inside AdS black holes: a possible long-term UV framing for a shell-plus-bulk-black-hole version of HDBLAST.
- *Already in the programme record:* not found in the scanned programme files

### G. Swampland and string-theory context

**G1. De Sitter Space and the Swampland** (2018)
- Authors: G. Obied, H. Ooguri, L. Spodyneiko, C. Vafa
- Venue or identifier: arXiv (2018)
- URL: <https://www.semanticscholar.org/paper/De-Sitter-Space-and-the-Swampland-Obied-Ooguri/599c99078a502b7d462d0b8783ff4a2c4436dc19>
- Access: search-snippet-only (queries 21). Implication: **constrains**.
- *Finding (from snippet):* Conjectures that |grad V| >= c V for the scalar potential of any consistent quantum-gravity theory, which forbids de Sitter vacua.
- *Relevance to HDBLAST:* HDBLAST's late-time attractor is an empty de Sitter brane. If the model were embedded in string theory, it would be judged against this family of conjectures; brane-induced de Sitter (dark bubble) is argued to evade them, and HDBLAST's version has not been assessed.
- *Caveat:* Year from general knowledge; not stated in the snippet.
- *Already in the programme record:* not found in the scanned programme files

**G2. de Sitter Bubbles and the Swampland** (2020)
- Authors: A. Bedroya, M. Montero, C. Vafa, I. Valenzuela
- Venue or identifier: arXiv:2008.07555
- URL: <https://arxiv.org/abs/2008.07555>
- Access: search-snippet-only (queries 46). Implication: **constrains**.
- *Finding (from snippet):* Effective theories of rolling scalars dual to a cascade of short-lived de Sitter spaces decaying by bubble nucleation; the trans-Planckian censorship conjecture (TCC) essentially incorporates the weak gravity and distance conjecture constraints on the dual potentials.
- *Relevance to HDBLAST:* TCC-type bounds would limit how long any HDBLAST de Sitter phase (the static shell or the +1 plateau) can last in a UV-complete version.
- *Already in the programme record:* not found in the scanned programme files

**G3. Lectures on the string landscape and the Swampland** (2022)
- Authors: N. B. Agmon, A. Bedroya, M. J. Kang, C. Vafa
- Venue or identifier: arXiv:2212.06187
- URL: <https://arxiv.org/pdf/2212.06187>
- Access: search-snippet-only (queries 54). Implication: **context**.
- *Finding (from snippet):* Lecture notes on the string landscape and on the Swampland programme's constraints for EFTs with a quantum-gravity UV completion.
- *Relevance to HDBLAST:* Background reference for the swampland criteria cited here.
- *Already in the programme record:* not found in the scanned programme files

**G4. On the Origin and Fate of Our Universe** (2025)
- Authors: C. Vafa
- Venue or identifier: Gen. Relativ. Gravit. (2025); arXiv:2501.00966
- URL: <https://arxiv.org/abs/2501.00966>
- Access: search-snippet-only (queries 56). Implication: **constrains**.
- *Finding (from snippet):* A short review of swampland bounds on positive potentials, the de Sitter conjecture, TCC and its relation to the species scale, with implications for inflation and the fate of the universe.
- *Relevance to HDBLAST:* The current mainstream quantum-gravity view on origins: any HDBLAST claim about an inflating or de Sitter origin should be stated against it.
- *Already in the programme record:* not found in the scanned programme files

**G5. The dark dimension scenario (Montero, Vafa, Valenzuela 2022; popular account)** (2022)
- Authors: M. Montero, C. Vafa, I. Valenzuela (per snippet; arXiv:2205.12293 named in snippet, not opened)
- Venue or identifier: Quanta Magazine article (2024) and follow-up EPJC paper arXiv:2309.09330
- URL: <https://www.quantamagazine.org/in-a-dark-dimension-physicists-search-for-missing-matter-20240201/>
- Access: search-snippet-only (queries 22). Implication: **context**.
- *Finding (from snippet):* One micron-size extra dimension with a Kaluza-Klein scale of order meV, motivated by the tiny cosmological constant (about 1e-122 in Planck units) and swampland arguments.
- *Relevance to HDBLAST:* A concrete single-extra-dimension scenario with laboratory tests; A12 ties the dark bubble to it. HDBLAST's extra-dimension scale (AdS radius ell = 9 in model units at phi=+1) has no physical normalization yet (open item 7).
- *Already in the programme record:* not found in the scanned programme files

**G6. Inflation with a Growing Fifth Dimension** (2025)
- Authors: Harvard group (snippet gives affiliation only)
- Venue or identifier: JHEP 05 (2026) 169; arXiv:2512.04177
- URL: <https://arxiv.org/abs/2512.04177>
- Access: search-snippet-only (queries 85). Implication: **supports**.
- *Finding (from snippet):* Inflation with a finite initial time in warped AdS5 with UV and IR branes. The inflaton potential detunes the brane tension, so the fifth dimension grows: a two-field (inflaton plus radion) hyperbolic model with early radion fast-roll and late inflaton slow-roll. The earliest modes are radion-sourced and give a suppressed, blue-tilted scalar spectrum and oscillatory tensors, in principle visible in the CMB. The model links to the dark dimension.
- *Relevance to HDBLAST:* Structurally the closest 2025-2026 match to HDBLAST's core mechanism: a detuned brane tension drives the dynamics of the extra dimension. It shows how such a mechanism yields observables (large-scale suppression), and is a template for HDBLAST open item 5.
- *Already in the programme record:* not found in the scanned programme files

**G7. De Sitter space constraints on brane tensions and couplings** (2024)
- Authors: S. Hassan, G. Obied, J. March-Russell
- Venue or identifier: arXiv:2411.14529
- URL: <https://arxiv.org/abs/2411.14529>
- Access: search-snippet-only (queries 44). Implication: **context**.
- *Finding (from snippet):* Festina-Lente-type arguments using Nariai de Sitter black holes bound p-brane tensions in de Sitter by the Hubble rate and by Chern-Simons-like worldvolume couplings; D-branes satisfy them; axion domain walls evade them.
- *Relevance to HDBLAST:* Tangential: HDBLAST's shell has no gauge couplings, so these bounds would matter only if gauge fields were added on the shell.
- *Already in the programme record:* not found in the scanned programme files

**G8. Stability of non-supersymmetric vacua from calibrations** (2025)
- Authors: not captured in snippet
- Venue or identifier: JHEP 11 (2025) 070; arXiv:2507.02787
- URL: <https://arxiv.org/abs/2507.02787>
- Access: search-snippet-only (queries 70, 23). Implication: **context**.
- *Finding (from snippet):* Calibrations bound D-brane energies and forbid, in the probe approximation, nucleation of D-brane bubbles in several type II AdS4 and AdS5 non-supersymmetric vacua, many of which resisted all decay channels tested. Query 23 snippets summarise the opposite view: the Ooguri-Vafa conjecture that non-SUSY AdS is at best metastable via Brown-Teitelboim brane nucleation.
- *Relevance to HDBLAST:* The instability of non-supersymmetric AdS, which the dark-bubble origin mechanism needs, is contested. HDBLAST's bulk vacua come from a real superpotential and are perturbatively stable; whether a 5D 'blast' event is available depends on that unresolved UV question.
- *Already in the programme record:* not found in the scanned programme files

**G9. Dark energy from string theory: an introductory review** (2026)
- Authors: D. Andriot
- Venue or identifier: arXiv:2603.25797 (Mar 2026, rev. Apr 2026)
- URL: <https://arxiv.org/abs/2603.25797>
- Access: search-snippet-only (queries 104). Implication: **context**.
- *Finding (from snippet):* A review of obtaining dark energy from string theory, either as a cosmological constant (de Sitter solution) or dynamical (quintessence), including historical no-go constraints and attempts to evade them.
- *Relevance to HDBLAST:* The current reference point for whether any late-time de Sitter phase, such as HDBLAST's +1 plateau, can be UV-complete.
- *Already in the programme record:* not found in the scanned programme files

**G10. Our universe: An expanding bubble in an extra dimension (press release on A1)** (2018)
- Authors: Uppsala University press release via ScienceDaily
- Venue or identifier: ScienceDaily, 28 Dec 2018
- URL: <https://www.sciencedaily.com/releases/2018/12/181228164824.htm>
- Access: search-snippet-only (queries 67). Implication: **context**.
- *Finding (from snippet):* Popular account of the 2018 dark bubble proposal: the universe rides an expanding bubble where two five-dimensional spaces meet. A 2026 search found no newer comparable press story.
- *Relevance to HDBLAST:* Shows that 'our universe is a 5D bubble' has already been widely publicised; HDBLAST communication should cite this line to avoid implying the concept is new.
- *Already in the programme record:* not found in the scanned programme files


## Seen but not used

- Query 60: 'What's the (Dark) Matter with Cosmological Bubbles?' is about 4D first-order phase transitions; off-topic.
- Query 79: 'Scalar stars and lumps with (A)dS core' concerns 4D compact objects; only indirectly relevant, covered by B9.
- Query 108: No web-indexed item mentions HDBLAST or Ricardo Maldonado; the query surfaced C16 and standard ekpyrotic items.
- Query 63: Double-holography / Karch-Randall information-theory papers (e.g. arXiv:2504.21856) surfaced; judged not relevant to the HDBLAST dynamics question.
- Query 25: Popular articles on 'universe from a 4D-bulk black hole' and dark-dimension primordial black holes (arXiv:2506.14874); popular framing of B7, not new physics.

## Limitations

- **Snippet-level evidence only.** No abstract or full text was opened, and snippets may misattribute, garble or overstate. Several such cases are flagged in the caveats.
- **Search coverage is incomplete.** It depends on the search engine's index and on the queries chosen. A sweep cannot establish that prior art is absent.
- **Unverified bibliographic details.** Several entries omit authors or journal details because the snippet omitted them; nothing was filled in by guesswork. The few details taken from general knowledge are flagged.
- **Relevance judgements are the sweep's own.** They were not checked by the cited authors or by external review.
- **The cross-check is a consistency check.** The thin-brane comparison uses the checkpoint's floating-point results. It adds no new existence, stability or dynamical claim.

## Reproduction

From `/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/literature`:

```
python3 theory_checks/rs_thin_brane_crosscheck.py   # exact + numeric cross-check, writes rs_thin_brane_crosscheck.json (a few seconds, 1 core)
python3 theory_checks/prior_mention_scan.py         # marks sources already in the programme record, writes prior_mentions.json (~20 s)
python3 theory_checks/build_theory_md.py            # validates the source records and regenerates theory.md and theory_sources.json
python3 theory_checks/builder_negative_control.py   # confirms the builder rejects corrupted records (temporary copy only)
```

Requirements: python3, sympy 1.14 and mpmath (numpy is not needed). The builder fails if any source lacks a URL, year, finding, relevance or valid label; if any source cites a query number that does not exist; if fewer than 15 queries are recorded; or if URLs are duplicated.

## Appendix: queries run (28 Sept 2026)

| # | Query |
|---|---|
| 1 | Emergent de Sitter cosmology from decaying anti-de Sitter space Banerjee Danielsson Dibitetto Giri Schillo |
| 2 | "dark bubble" braneworld cosmology Danielsson 2024 OR 2025 |
| 3 | Karch Randall "locally localized gravity" de Sitter brane AdS |
| 4 | ekpyrotic universe colliding branes big bang Khoury Ovrut Steinhardt Turok |
| 5 | "Dark bubbles, dark dimensions and fat gravitons" |
| 6 | "Dark bubble cosmology and the equivalence principle" 2507.03748 |
| 7 | "Weak gravity at micron scales from dark bubble cosmology" Danielsson Giri |
| 8 | "The dark bubbleography" Banerjee Danielsson Zemsch non-normalizable modes |
| 9 | dark bubble nucleation catalyzed by five-dimensional black hole quantum origin of the universe |
| 10 | "Dark bubbles and black holes" Danielsson Dibitetto Giri AdS5 nucleation |
| 11 | "Out of the white hole" holographic origin for the Big Bang Pourhasan Afshordi Mann 5D black hole |
| 12 | Garriga Sasaki "Brane-world creation and black holes" de Sitter brane AdS instanton |
| 13 | Gen Sasaki "Radion on the de Sitter brane" negative mass squared |
| 14 | thick de Sitter brane stability scalar field domain wall Kobayashi Koyama Soda |
| 15 | stability of de Sitter domain walls bulk scalar superpotential perturbations curved domain wall 2023-2025 |
| 16 | DeWolfe Freedman Gubser Karch "Modeling the fifth dimension with scalars and gravity" |
| 17 | "Fake supergravity and domain wall stability" Freedman Nunez Schnabl Skenderis |
| 18 | "Can inflating braneworlds be stabilized" de Sitter brane radion |
| 19 | Kaloper "Bent domain walls as braneworlds" de Sitter brane AdS5 |
| 20 | Takamizu Maeda collision of domain walls reheating of the brane universe colliding branes numerical |
| 21 | de Sitter swampland conjecture review 2024 2025 status Obied Ooguri Spodyneiko Vafa trans-Planckian censorship |
| 22 | "dark dimension" Montero Vafa Valenzuela micron extra dimension cosmological constant 2022 |
| 23 | non-supersymmetric AdS vacua unstable Ooguri Vafa conjecture bubble nucleation brane Brown-Teitelboim |
| 24 | holographic reheating braneworld particle production brane energy transfer bulk 2023 OR 2024 OR 2025 |
| 25 | universe emerged from higher-dimensional black hole brane 2025 study Big Bang fifth dimension |
| 26 | "de Sitter cosmology on an expanding bubble" Banerjee Danielsson Dibitetto Giri Schillo JHEP 2019 |
| 27 | Antonini Swingle Van Raamsdonk end-of-the-world brane cosmology from AdS |
| 28 | Hawking Hertog Reall "Brane new world" de Sitter brane AdS created from nothing CFT |
| 29 | "Brane cosmology from AdS/BCFT" end-of-the-world brane scalar field big bang |
| 30 | "Effective Dynamics of Inflationary End-of-the-World Branes" |
| 31 | "Dynamical dark energy in 0'B braneworlds" Basile |
| 32 | cyclic universe brane collision 2024 2025 new work Ijjas Steinhardt bounce higher dimensional |
| 33 | Garriga Vilenkin perturbations on domain walls de Sitter wall fluctuation tachyonic mode nucleated bubble |
| 34 | "Born-Again Braneworld" Kanno Sasaki Soda radion de Sitter brane |
| 35 | Randall Sundrum "An alternative to compactification" 1999 single brane AdS5 |
| 36 | Binetruy Deffayet Langlois "Non-conventional cosmology from a brane-universe" rho squared Friedmann dark radiation |
| 37 | de Sitter brane Kaluza-Klein graviton spectrum mass gap (3/2)H Randall-Sundrum inflating brane |
| 38 | Koyama Takahashi bulk inflaton de Sitter brane scalar cosmological perturbations exactly solvable model |
| 39 | Maartens Koyama "Brane-World Gravity" Living Reviews in Relativity |
| 40 | Kiritsis Kofinas Tetradis Tomaras Zarikas "Cosmological evolution with brane-bulk energy exchange" |
| 41 | "Primordial Correlators from a Kaluza-Klein Graviton Continuum" |
| 42 | "Regular collision of dilatonic inflating branes" |
| 43 | numerical relativity bubble collisions observational signatures Johnson Peiris Lehner Wainwright eternal inflation CMB |
| 44 | "De Sitter space constraints on brane tensions and couplings" |
| 45 | "Shedding light on dark bubble cosmology" radiation brane Maxwell |
| 46 | "de Sitter bubbles and the swampland" 2008.07555 |
| 47 | "Catalytic creation of a bubble universe induced by quintessence in five dimensions" |
| 48 | Skenderis Townsend "Hidden supersymmetry of domain walls and cosmologies" curved domain walls superpotential |
| 49 | Giovannini "gauge-invariant fluctuations of scalar branes" thick brane stability |
| 50 | Charmousis Gregory Rubakov "Wave function of the radion in a brane world" brane bending ghost |
| 51 | Garriga "Smooth creation of an open universe in five dimensions" |
| 52 | "Cosmological perturbations in the 5D Big Bang" Tolley Turok |
| 53 | Takamizu Kudoh Maeda "Dynamics of colliding branes and black brane production" AdS bulk |
| 54 | "The String Landscape, the Swampland and the Observed Universe" Agmon Bedroya Kang Vafa |
| 55 | "Self-gravitating electromagnetic waves in the dark bubble model" |
| 56 | "On the Origin and Fate of Our Universe" 2501.00966 |
| 57 | "Cosmological Perturbations in the 5D Holographic Big Bang Model" 1703.00954 |
| 58 | "Dark bubbles: decorating the wall" string cloud matter on bubble |
| 59 | "Gravitational waves in dark bubble cosmology" Danielsson Giri Panizo |
| 60 | "What's the (Dark) Matter with Cosmological Bubbles" |
| 61 | dark bubble Vilenkin tunneling wave function quantum cosmology nucleation of brane AdS5 Danielsson |
| 62 | ekpyrotic slow contraction numerical relativity robustness initial conditions Ijjas Cook Pretorius Steinhardt |
| 63 | double holography de Sitter brane Karch-Randall braneworld 2023 2024 2025 induced gravity dS brane AdS bulk |
| 64 | Kraus "Dynamics of anti-de Sitter domain walls" brane in AdS-Schwarzschild bulk black hole mass dark radiation |
| 65 | Savonije Verlinde "CFT and entropy on the brane" Friedmann equation Cardy-Verlinde brane moving AdS black hole |
| 66 | Horowitz Orgera Polchinski "Nonperturbative instability of AdS5 x S5/Zk" bubble of nothing |
| 67 | our universe bubble expanding in fifth dimension new study 2026 cosmology string theory |
| 68 | "When do colliding bubbles produce an expanding universe" Gratton Turok |
| 69 | "De Sitter versus Anti de Sitter flows and the (super)gravity landscape" |
| 70 | "Stability of non-supersymmetric vacua from calibrations" brane nucleation AdS |
| 71 | Witten "Instability of the Kaluza-Klein vacuum" bubble of nothing 1982 |
| 72 | Ghosh Kiritsis Nitti Witkowski "De Sitter and Anti-de Sitter branes in self-tuning models" |
| 73 | Ghosh Kiritsis Nitti Witkowski "Holographic RG flows on curved manifolds and quantum phase transitions" |
| 74 | Charmousis Kiritsis Nitti "holographic self-tuning of the cosmological constant" brane Einstein-dilaton |
| 75 | Kiritsis Tsouros 2025 2026 Einstein-scalar dS AdS flows "Centaur" no regular solutions interpolating |
| 76 | "Nucleation of de Sitter from the anti de Sitter spacetime in scalar field models" |
| 77 | "Holographic confining theories on space-times with constant positive curvature" Kiritsis |
| 78 | "Brane cosmology and the self-tuning of the cosmological constant" Kiritsis Nitti Amariti Ghosh |
| 79 | "Scalar stars and lumps with (A)dS core" Einstein-scalar |
| 80 | "Self-sustained, out-of-equilibrium inflation" holographic |
| 81 | holographic reheating strongly coupled sector de Sitter thermalization end of inflation AdS/CFT 2024 2025 arXiv |
| 82 | Lehners Turok "Bouncing negative-tension branes" OR "M-theory model of a big crunch/big bang transition" |
| 83 | Nihei OR "Kim and Kim" inflation on the brane de Sitter brane AdS5 tension exceeds critical Hubble rate formula 1999 |
| 84 | "A dynamical inflaton coupled to strongly interacting matter" holography reheating |
| 85 | "Inflation with a growing fifth dimension" dark dimension |
| 86 | "Gravitational reheating at strong coupling" 2304.11195 |
| 87 | brane world dark radiation BBN constraint bulk Weyl C/a^4 2023 2024 N_eff braneworld constraint |
| 88 | BraneCode numerical code brane world dynamics bulk scalar five-dimensional simulations |
| 89 | "End-of-the-World Branes and Inflationary Predictions for Rocky and Swampy Landscapes" |
| 90 | "Curing with hemlock" escaping the swampland using instabilities from string theory dark bubble |
| 91 | "Colliding bubble worlds" Gen Ishihara Sasaki braneworld big bang bubble collision 5D |
| 92 | de Sitter brane bulk scalar field linear stability perturbations 2025 OR 2026 arXiv warped braneworld tachyon radion single brane |
| 93 | Afshordi holographic big bang 5D black hole collapse brane 2024 OR 2025 update |
| 94 | particle production nucleated bubble wall brane universe expanding bubble reheating quantum fields on bubble wall |
| 95 | "The global structure of the colliding bubble braneworld universe" OR "A braneworld universe from colliding bubbles" Bucher |
| 96 | "Brane big bang brought on by a bulk bubble" |
| 97 | "Coupled bulk and brane fields about a de Sitter brane" radion sqrt(2)H normalizable single brane |
| 98 | Shiromizu Maeda Sasaki "The Einstein equations on the 3-brane world" projected Weyl tensor |
| 99 | Coleman De Luccia "Gravitational effects on and of vacuum decay" bubble nucleation AdS true vacuum crunch |
| 100 | "Crunch from AdS bubble collapse in unbounded potentials" |
| 101 | "Generalized surface tension bounds in vacuum decay" Espinosa superpotential tunneling potential |
| 102 | "Searching for Coleman-de Luccia bubbles in AdS compactifications" |
| 103 | "Features of a dark energy model in string theory" dark bubble 2212.14004 |
| 104 | "Dark energy from string theory: an introductory review" 2603.25797 |
| 105 | "Bubbles of cosmology in AdS/CFT" Sahu Simidzija Van Raamsdonk |
| 106 | "Revisiting Coleman-de Luccia transitions in the AdS regime using holography" |
| 107 | "Nothing really matters" JHEP 2020 bubble of nothing dark bubble |
| 108 | "higher-dimensional blast" Big Bang hypothesis five-dimensional brane Maldonado |
| 109 | Garriga Tanaka "Gravity in the Randall-Sundrum brane world" linearized gravity brane bending |
| 110 | "The big bang as a higher-dimensional shock wave" Kaluza-Klein |
| 111 | "Solution of a Braneworld Big Crunch/Big Bang Cosmology" McFadden Turok Steinhardt five-dimensional dynamics |
