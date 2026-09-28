# HDBLAST literature sweeps, 27 September 2026 round

**Work in progress.** Nothing in this folder is a claim until it appears in a dated checkpoint or a Zenodo record. Other sweeps (observations, methods) may add their own files and sections here. This section covers **sweep 1: theory**.

## Sweep 1: theory (`theory.md`)

### Question

What published work, 2022–2026 plus foundational papers, supports, constrains or challenges the HDBLAST hypothesis and the registered model's current status? That status is: a certified tachyonic de Sitter shell; two roll-off fates (an empty Randall–Sundrum de Sitter brane, or collapse); no radiation era; and the corrected static "+1 branch".

### Method

- 111 distinct web searches on 28 Sept 2026.
- Papers could not be opened: arxiv.org, zenodo.org, Springer, DiVA and uu.se were blocked, and fetches returned `EGRESS_BLOCKED`. **Every entry is based on search-result snippets only**, and each is labelled so.
- Records are kept in `theory_checks/theory_sources.py`. The builder validates them (URLs, labels, query numbers, arXiv-id consistency) and renders `theory.md` and `theory_sources.json`.
- A scan of the programme's own files marks which sources were already in the HDBLAST record.
- One exact and numerical cross-check places the registered model inside the standard Randall–Sundrum, Karch–Randall and DeWolfe–Freedman–Gubser–Karch formulas.

### Results

| Item | Result | Label |
|---|---|---|
| Sources annotated | 95 (31 from 2022–2026); 26 already in the programme record, 69 not found there | inventory |
| Implication labels | supports 7, constrains 8, challenges 4, method 41, context 35 | judgement |
| Balanced tension 2W = critical RS/KR tension 6/ℓ at both vacua (ℓ = 9 at φ=+1, 9/5 at φ=−1) | holds | EXACT |
| First-order flow φ′=W_φ, A′=−W/3 solves the Einstein equations computed from the metric; φ=−tanh y is the flat wall | holds | EXACT |
| Thin-brane H²=(σ(1)/6)²−1/81 vs the checkpoint's H² expansion | equal at O(δ); at O(δ²) they differ by exactly −c²δ²/384 | EXACT |
| Thin-brane H² vs checkpoint `metric_only_H2` (δ=10⁻³) | relative difference 2.3×10⁻¹⁴ | NUMERICAL |
| Full +1 branch vs thin brane, shift in H | −7.869 ppm (second-order prediction −7.849 ppm) | NUMERICAL |
| U for large φ (either sign) | ~ −(2/27)φ⁶, unbounded below | EXACT |
| Controls: one-sided junction, wrong vacuum radius, wrong-sign flow, c perturbed by 1% | all fail as they should | control |
| Builder negative control: corrupted arXiv id and query number | build refused | control |

**Main conclusions.** All are labelled in `theory.md`, and none claims novelty or confirmation.

1. **Prior art.** "Our universe as a wall born in a 5D event" has substantial prior art (1998–2026). HDBLAST's distinct content is its registered model and results.
2. **Negative result expected.** No published hot-origin model gets heat from a pure-tension wall, so the HDBLAST no-radiation result is the expected outcome.
3. **Collapse fate testable.** The literature suggests black-brane or crunch formation behind an apparent horizon for the collapse fate.
4. **Tachyon typical; +1 branch open.** The registered tachyon is typical of de Sitter branes with bulk scalars. No stability theorem found covers de Sitter-sliced walls.
5. **Most useful tools.** The Einstein–scalar de Sitter-sliced flow programme (Kiritsis, Nitti and collaborators) is new to the programme and is the most directly useful toolset.
6. **Matter extension must pass two tests.** It must pass equivalence-principle and dark-radiation bounds.

### Limitations

- Snippet-level evidence only; snippets can misattribute or garble. Conflicts are flagged per entry.
- Search coverage is incomplete, and absence of prior art cannot be established.
- The cross-check is a consistency check against published formulas. It adds no existence, stability or dynamical result.

### Reproduction

Run from this folder, with python3, sympy 1.14 and mpmath:

```
python3 theory_checks/rs_thin_brane_crosscheck.py    # writes theory_checks/rs_thin_brane_crosscheck.json (a few seconds, 1 core)
python3 theory_checks/prior_mention_scan.py          # writes theory_checks/prior_mentions.json (~20 s; reads programme files read-only)
python3 theory_checks/build_theory_md.py             # validates records, writes theory.md and theory_checks/theory_sources.json
python3 theory_checks/builder_negative_control.py    # confirms the builder rejects corrupted records (uses a temp copy)
```

### Files

- `theory.md`: annotated bibliography, synthesis, results table and query list (generated; edit the template or records instead).
- `theory_checks/theory_sources.py`: the source records and query list.
- `theory_checks/theory_template.md`: the narrative template.
- `theory_checks/rs_thin_brane_crosscheck.py` and `.json`: the exact and numerical cross-check with controls. `run_log.txt` is its console output.
- `theory_checks/prior_mention_scan.py` and `prior_mentions.json`: which sources were already in the programme record.
- `theory_checks/build_theory_md.py` and `theory_sources.json`: the validator and renderer, and the machine-readable bibliography.
- `theory_checks/builder_negative_control.py` and `.json`: the builder control.

## Sweep 3: methods (`methods.md`)

**The web search did not run.** All four queries attempted for this sweep were refused because the session's web-search budget was exhausted (200 of 200 calls). A page fetch was blocked (`EGRESS_BLOCKED`). The required minimum of 15 queries was **not met**, so the sweep is a procedural negative result. `methods.md` is therefore built only from sources already on record:
- 26 URLs from sweep-1 search records (snippet-level, second-hand);
- 16 URLs written in the programme's own archive files (not re-verified).

It also carries 23 general-knowledge leads without URLs (unverified, not cited) and 20 prepared queries for a later run. The builder machine-checks every URL against its claimed origin.

### Question

Which mathematical and numerical methods would most directly advance the open problems?
- certification of the +1 branch;
- negative modes and stability;
- constraint growth and far-boundary contamination in long 5D evolutions;
- the collapse fate;
- shell and junction numerics;
- particle production with consistent backreaction.

### Main results

| Item | Result | Label |
|---|---|---|
| Searches executed | 0 of ≥15 | procedural negative |
| The programme's existing kv/Krawczyk certification pipeline (M461/M462, original shell) | most direct route to certifying the +1 branch | judgement |
| Regular linear cone germ on the +1 side: η ∝ C₁₄⁽²⁾(cosh u), growth exponents 14 / −18 | holds | exact |
| Germ amplification at δ = 0.001 | 6.288×10¹⁸; a direct φ-formulation would need 39 digits, so the η-formulation is mandatory | interval / numerical |
| Leading-order model vs checkpoint nonlinear branch | O(δ) relative defects (−0.049δ in η_b, +0.0157δ in H²); a usable centre for a proof | numerical |
| The evolution chart's far end is a Rindler-type horizon (z → −∞) | any finite far boundary is artificial; a horizon-penetrating or characteristic chart is the remedy | exact observation / judgement |
| Verified 2020–2026 lattice work with renormalized backreaction | none found (no search) | negative for search goal |
| "Gregory–Ruth" reference | not identified | open |

Full question, method, results table, limitations and file list: `methods_checks/README.md`.

### Reproduction

```
python3 methods_checks/certification_pilot.py
python3 methods_checks/build_methods_md.py
python3 methods_checks/builder_negative_control.py
```
