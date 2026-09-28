# Sweep 3: mathematical and numerical methods (HDBLAST round of 27 September 2026)

**Work in progress.** Nothing here is a claim of the programme until it appears in a dated checkpoint or a Zenodo record. The deliverable is `../methods.md`. This folder holds its records, builder, controls and one pilot calculation.

## The main limitation, stated first

**No web search ran in this sweep.** The four queries attempted were all refused because the session's web-search budget was exhausted (200 of 200 calls, used by earlier work in the same session). A page fetch to journals.aps.org returned `EGRESS_BLOCKED`. The required minimum of 15 distinct queries was therefore **not met**. This is a procedural negative result, stated at the top of `methods.md`.

The bibliography uses only:
- **26 URLs from sweep-1 search records** in this session (`../theory_checks/theory_sources.json`). These are snippet-level and second-hand.
- **16 URLs written verbatim in the programme's own archive files.** These were not re-verified here.

It also lists 23 general-knowledge **leads** with no URL, which are not verified and not cited, and **20 prepared queries** for a later run.

## Question

Which methods would most directly advance the programme's open problems?

- **P1:** certify the static +1 branch.
- **P2:** stability and negative modes beyond the computed sectors.
- **P3:** constraint growth and far-boundary contamination in long 5D evolutions.
- **P4:** the collapse fate.
- **P5:** shell and junction numerics.
- **P6:** particle production with consistent backreaction.

## Method

1. Attempted web search; it was refused (see above).
2. Selected method-relevant items from the sweep-1 records and re-annotated them for methods. Searched the programme's archive (`new-files/D-Blast 3`, folders 13x–15x, read-only) for method references with URLs. Found:
   - the programme's own validated-numerics bibliography (M488G);
   - its computer-assisted BVP proof (M461/M462, which used the `kv` validated integrator and a Krawczyk test);
   - its particle-production benchmarks (Chat 7 and Chat 8).
3. `methods_sources.py` holds every record. `build_methods_md.py` validates it and renders `../methods.md` and `methods_sources.json`. The build fails unless:
   - each URL is byte-identical to its sweep-1 record, or appears verbatim in the named programme file;
   - access labels match provenance;
   - leads carry no URL;
   - the pilot JSON is current.
4. `certification_pilot.py` is a small validated-numerics pilot for P1, with exact, interval and numerical parts and controls.

## Results

| Item | Result | Label |
|---|---|---|
| Web searches executed | 0 of the ≥15 required (budget exhausted); 1 fetch blocked | procedural negative |
| Sources with URLs | 42 (26 sweep-1 search records, 16 programme records); implications: method 31, context 6, constrains 5 | inventory |
| Leads without URLs | 23, including the task-named CCZ4, Z4c, generalized harmonic, hyperboloidal, CDL negative modes, CosmoLattice, kv/CAPD/Arb | unverified |
| "Gregory–Ruth" | not identified in any available record | open |
| Regular linear cone germ on the +1 side: η ∝ C₁₄⁽²⁾(cosh u), with U''(1)/k² = 252 and growth exponents 14 / −18 | holds | exact |
| Evolution chart ends at a Rindler-type horizon: z = ln tanh(y/18) + const → −∞, ρ/y → 1 | holds; so any finite far boundary is artificial | exact |
| Germ amplification G(u_b) at δ = 0.001 | 6.288479×10¹⁸; enclosure relative width about 10⁻⁵⁹ | interval |
| Unwanted-mode suppression outward over [0.25, u_b] | about 10⁴³·³ | conditional (linear regime) |
| Digits for a direct φ-formulation (16 significant digits of η_h) | 39; binary64 gives 1+η_h = 1 | numerical |
| LO truncated model vs checkpoint nonlinear solutions | relative defects −0.0493δ (η_b), −0.0587δ (transport ratio), +0.01572δ (H², matching the predicted 27c²/(384(1+c)) = 0.015717); slopes 1.000–1.002 | numerical |
| LO metric-only H² and LO η_b vs the checkpoint's own linear estimates | ≤ 6×10⁻¹⁴ and ≤ 4×10⁻¹⁶ relative | numerical reproduction |
| Controls: wrong Gegenbauer degree, index or dimension; wrong degree in the amplification; c + 1 %; wrong-sign junction; containment self-test | all detected | control |
| Builder negative control: fabricated programme URL, altered sweep-1 URL, lead with URL, upgraded access label | all refused; unmodified records build | control |

**One control needed correcting.** The wrong-sign-junction control was first written to expect a sign flip. The sign does not flip; η_b becomes about 8× too large instead. The criterion was changed to a magnitude test, and the original outcome is kept in `certification_pilot.json` under `control_design_note`.

## Limitations

- No fresh literature search. All evidence is second-hand snippet evidence or unverified programme records.
- The methods recommendations are judgements and were not implemented.
- The pilot encloses only a closed-form leading-order model. It does **not** prove existence, uniqueness or stability of the nonlinear +1 branch.
- The checkpoint values it compares with are floating-point, with solver tolerances 2×10⁻¹² and 8×10⁻¹⁴.

## Reproduction

Run from `research/HDBLAST_CHECKPOINT_20260927/literature/`, with python3, SymPy 1.14 and mpmath 1.3. One core, a few seconds in total:

```
python3 methods_checks/certification_pilot.py        # writes methods_checks/certification_pilot.json (reads the 22 Sept checkpoint JSON read-only, SHA-256 checked)
python3 methods_checks/build_methods_md.py           # validates records; writes methods.md and methods_checks/methods_sources.json
python3 methods_checks/builder_negative_control.py   # writes methods_checks/builder_negative_control.json; exits 0 only if all corruptions are refused
```

## Files

- `methods_sources.py`: records, attempted and planned queries, problem list, leads.
- `methods_template.md`: narrative template for `../methods.md`.
- `build_methods_md.py` and `methods_sources.json`: validator and renderer, and the machine-readable bibliography.
- `certification_pilot.py`, `certification_pilot.json` and `certification_pilot_run_log.txt`: the pilot.
- `builder_negative_control.py` and `builder_negative_control.json`: the builder control.
- `build_log.txt`: console output of the last build.
