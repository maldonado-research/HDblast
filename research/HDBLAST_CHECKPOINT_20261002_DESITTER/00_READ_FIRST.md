# HDBLAST: massive de Sitter quantum sources and classical shell susceptibility

Ricardo Maldonado · 2 October 2026 · [ORCID 0009-0009-3937-6527](https://orcid.org/0009-0009-3937-6527)

This checkpoint constructs actual one-loop vacuum stress and a paired scalar source for a free, minimally coupled, massive scalar in the Euclidean/Bunch–Davies de Sitter state. A separately resolved classical shell response supplies the derivatives needed for a future stationary quantum correction. The two calculations have not yet been combined into a physical quantum-corrected shell solution.

The source calculation extends the inherited finite-action convention to an explicit convergent sphere action. With `z=H²`, `x>0` the effective mass squared, and a fixed positive reference `r`, both sources follow from that same action density `W`:

```text
Q = 2 W_x,   rho = W − z W_z/2,   p = −rho,   J_phi = x_phi Q/2.
```

The reference stays fixed during every variation. The scalar source retains the unspecified mass-law derivative `x_phi`. Stress is obtained by metric variation, including the sphere volume. The trace identity is checked afterward.

## What was obtained

All four registered mathematical points pass. The table reports dimensionless quantities; `r=1` sets the calculation's units and does not select a physical HDBLAST mass.

| x/r | H²/r | rho/r² | Q/r |
|---:|---:|---:|---:|
| 1 | 0.5 | 2.902429739754e−4 | 4.221715985097e−3 |
| 1 | 1 | 2.119232911866e−3 | 2.071990800425e−2 |
| 2 | 0.5 | −4.118184790862e−5 | −1.351351456310e−4 |
| 2 | 1 | 3.694001486960e−4 | 2.110857992549e−3 |

The negative entries are retained finite-scheme vacuum expectation values. They do not describe particle occupation numbers or heating. See [the numerical report](outputs/summary/NUMERICAL_RESULTS.md) for the action values and error evidence.

A separate exact-mode calculation at `x=r=2H²` directly integrates energy, pressure and variance against the inherited subtraction orders. It obtains

```text
rho = 11 H⁴/(960 pi²),   p = −11 H⁴/(960 pi²),   Q = H²/(12 pi²).
```

Pressure is calculated separately in this check. This establishes a curved-state subtraction bridge at one point, beyond agreement between two determinant evaluations.

The classical calculation resolves the previously unresolved tiny derivative on the inherited corrected `+1` branch at `delta=0.001`: `E1_ell ≈ +9.63825e−14`, where `ell=log10(abs(phi_h−1))`. Its empirical absolute error scale is `3.34e−18`. The archived negative finite-difference estimate was at a numerical floor. Nonsingularity survives; small response coefficients change. For independently supplied infinitesimal sources `s=K rho` and `t=K J_phi`, the resolved response is

```text
Delta phi_b = −7.3550543e−6 s −0.140658540041 t + O(source²),
Delta H_b²  =  0.0371257908752 s +2.1785725e−10 t + O(source²).
```

No physical `K`, mass law, or actual source amplitude is supplied in that response calculation. Its error scales and validity conditions are in [the sensitivity report](sensitivity/STATIC_SENSITIVITY.md).

## Evidence and qualifications

| Component | Recorded outcome |
|---|---|
| Primary spectral calculation | 57 gates pass; 3 wrong-source controls detected |
| Symbolic common-action derivation | 37 exact assertions pass; 6 wrong-formula controls detected |
| Independent proper-time computation | 8 evaluations; 12 refinement comparisons pass |
| Independent audit | 12 cross-method comparisons and 8 trace gates pass; 5 corrupted-record controls rejected |
| Classical sensitivity | 19 integrations; 11 numerical gates pass; 3 wrong-formula controls detected |

The separate derivations, implementations and review were produced within this AI-assisted research workflow. They are internal cross-checks, not external peer review or independent human replication. The largest refined proper-time versus spectral difference is below `5.99e−31`. Analytic bounds control specified series, harmonic and infrared tails; the combined quadrature, ultraviolet remainder, special-function evaluation and arithmetic have no global certified enclosure. Agreement and refinement do not certify every displayed digit.

The original proper-time aggregator saved trace residuals but omitted them from aggregate acceptance. The later audit restores enforcement of all eight originally registered trace gates and passes them. This post-run validation correction and the original evidence are retained. The sensitivity's first symbolic verifier failed; the corrected expansion step passes, with the failed source and log preserved. Primary source evaluation had no failed run or gate repair.

Primary registration was frozen locally before evaluation and publicly committed as `85aea9955b99bba911a0e66e5869c184dd86a660` before the primary run. Independent proper-time and sensitivity registrations were local prospective records; their computations began before that public commit. The exact-mode bridge and later audit are post-registration diagnostics. [The review](independent/INDEPENDENT_SCIENTIFIC_REVIEW.md) and [prospective record](PROSPECTIVE_RECORD.md) distinguish these claims.

This checkpoint provides static vacuum polarization and classical local susceptibility. It establishes no causal in-in response kernel, quantum-corrected shell solution, time evolution, stability theorem, particle production, thermalization, heating, or hot Big Bang mechanism. Physical mass and gravitational matching remain unspecified; the massless endpoint is excluded.

Start with [claims and next steps](CLAIMS_AND_NEXT_STEPS.md), [reproduction](REPRODUCE.md), [the derivation](theory/DESITTER_COMMON_ACTION_DERIVATION.md), and [the literature review](literature/LITERATURE_REVIEW.md). The next priority is a constrained stationary quantum correction with explicit dimensionless physical matching, followed by a separately constructed causal in-in evolution. A static Euclidean vacuum action alone does not supply a time-dependent FRW source.

`CITATION.cff` and `zenodo_metadata.json` prepare citation and deposition metadata. No new DOI or completed Zenodo deposition is claimed. Reports and data use CC BY 4.0; source code uses MIT, as specified in [LICENSE.md](LICENSE.md).
