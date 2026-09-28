# HDBLAST literature sweep 1: theory

Research round of 27 September 2026, written 28 September 2026. **Work in progress.** This is an annotated bibliography with a short synthesis. It adds no new physics. Nothing here is a claim of the HDBLAST programme until it appears in a dated checkpoint or a Zenodo record.

## Question

What published work (2022–2026, plus foundational papers) bears on the HDBLAST hypothesis that a five-dimensional gravitational event sparked our Big Bang? Specifically: which results support, constrain or challenge the registered model and its current status? That status is: a de Sitter shell in a Z2-doubled Einstein–scalar bulk; a certified tachyonic roll-off with two fates (relaxation to an empty Randall–Sundrum de Sitter brane, or reversal and collapse); no radiation era; and the corrected static "+1 branch" of 22 Sept 2026.

## How the sweep was done, and what "read" means here

- **{{N_QUERIES}} distinct web searches** were run on 28 Sept 2026. All are listed in the appendix. They cover:
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

Totals: **{{N_SOURCES}} annotated sources.** By implication: {{IMPLICATION_COUNTS}}. By period: {{PERIOD_COUNTS}}.

**Relation to the programme's own earlier literature record.**
- The Chat 9 package (16 Sept 2026) contains `literature/PRIOR_ART_SHELL_STABILITY.md`. That report read papers through ar5iv at the time. It already covers: Frolov–Kofman; BraneCode; the bulk-inflaton papers; Koyama–Takahashi; DeWolfe–Freedman–Gubser–Karch; Garriga–Sasaki; Gen–Sasaki; and the dark-bubble papers up to June 2026.
- The Chat 9 moving-shell analysis already uses the Kraus/Ida moving-wall method and the Maeda–Wands projected equations.
- `theory_checks/prior_mention_scan.py` searched the programme files for each source's arXiv number and title: {{PRIOR_FILES}} Markdown/text files in D-Blast 3 folders 146–152, the DBlast 4 vector-store handoffs and the 22 Sept package. Result: **{{N_PRIOR}} of {{N_SOURCES}} sources were already mentioned there, and {{N_NEW}} were not found.** Each entry below says which.
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

All values come from `theory_checks/rs_thin_brane_crosscheck.json` (SymPy exact algebra; mpmath at 30 digits). The comparison input is the checkpoint's `static_branch/PLUS_BRANCH_RESULTS.json`, read-only, SHA-256 `{{INPUT_SHA}}`.

{{RESULTS_TABLE}}

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

{{BIBLIOGRAPHY}}

## Seen but not used

{{EXCLUDED}}

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

{{QUERY_TABLE}}
