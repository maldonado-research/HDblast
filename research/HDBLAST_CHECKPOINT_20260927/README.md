# HDBLAST research checkpoint, 27 September 2026

Ricardo Maldonado's HDBLAST program; prepared with AI assistance. Calculations and independent audits were completed on 28 September 2026. **Start with [`00_READ_FIRST.md`](00_READ_FIRST.md).** A one-page overview for visitors is in [`PUBLIC_SUMMARY.md`](PUBLIC_SUMMARY.md).

## Question

The 22 September checkpoint (Zenodo record 22922928) left three questions open:

1. Is the corrected static "+1 branch", the proposed late-time state of the registered 5D Einstein–scalar shell model, stable?
2. Can any mechanism turn the five-dimensional "blast" into a hot, radiation-dominated Big Bang?
3. Can long 5D evolutions be made constraint-accurate?

## Method

Eight workstreams. The six computational ones each have an independent adversarial audit (`VERIFICATION.md` plus a `verify/` folder); the literature sweeps and the program map are self-checked only:

| Folder | Question | Audit |
|---|---|---|
| `stability_gauge_invariant/` | Linear stability of the +1 branch, all dS₄ harmonics, scalar and tensor sectors | yes |
| `stability_time_domain/` | The same by an independent discrete-operator and time-domain method | yes |
| `analytic_structure/` | Exact small-detuning series and solvable structure of the branch | yes |
| `preheating/` | Mode-function particle production along the archived roll-off | yes |
| `evolution_constraints/` | Source of, and remedies for, the constraint stall in 5D evolutions | yes |
| `mechanisms/` | Quantitative screen of routes to a radiation era, including coupled 5D pilots | yes |
| `literature/` | 2022–2026 sweeps: theory, observations, methods (snippet-level) | self-checked only |
| `program_map/` | Chronology and claims ledger of the whole program | self-checked only |
| `synthesis/` | Cross-method reconciliation, quoted-number checks, builders for this report | this report |

## Results

| Result | Status |
|---|---|
| The +1 branch is linearly stable (scalar with brane bending and both junctions, and tensor) at six δ values from 0.0003 to 0.1. Two methods agree where they overlap (the second covers only the homogeneous sector), and both reproduce the known tachyon (growth 1.6571936) | numerical, verified |
| Exact series for φ_b, H² and ρ_b through δ⁸ (δ⁹ from the audit); the degree-14 Gegenbauer polynomial is explained by Δ = 3W″/W = 18 | exact-verified |
| Instant preheating yields ≤ 1.7×10⁻⁴ of the energy a radiation era needs, at the M₅ cutoff (screen-passing couplings, data-end energy, fixed archived trajectory; ≤ 5.9×10⁻⁴ over the stored profile) | conditional negative, verified |
| The 4D energy budget allows ≤ 0.59 e-folds of radiation domination; about 22 are needed | conditional, verified |
| Tuned-vacuum model change (pilots at δ = 0.1): radiation reaches 69% of H², but bulk dark radiation is 1.23× the radiation at the end of the reliable window, above the N_eff limits so far | inconclusive |
| Constraint stall traced to an initial-data defect; discrete projection restores convergence of order ≈ 3 | numerical, verified |
| Three claims refuted by the audits and corrected here (δ⁹ locality, c-family continuation, a stale decay number) | corrected |

The full claim list is in `RESULTS_SUMMARY.json` (66 claims) and `CLAIM_REVIEW.md`. The literature is in `LITERATURE_2022_2026.md` and the next steps are in `NEXT_TESTS.md`.

## Limitations

- **Precision.** All numerics are floating point; nothing is interval-certified.
- **Stability scope.** The vector sector, nonlinear stability, tunnelling and dynamical attraction were not studied. Method B covers only homogeneous perturbations; Method A's real-μ² scan covers [−400, 2.2499], and regularity at the cone is an assumed criterion.
- **Coupling to the bulk.** The preheating and mechanism results use a fixed trajectory or a phenomenological coupling, and the coupled bulk problem is not solved.
- **Literature.** All sources are search snippets; no paper was opened. Several references cited in workstream READMEs have no record in the literature files; they are listed in `CLAIM_REVIEW.md` ("Citations") and are unverified.
- **Scope of claims.** Nothing here is observational support for the hypothesis, and no discovery or novelty is claimed beyond this project.

## Reproduction

See [`REPRODUCE.md`](REPRODUCE.md), with `requirements.txt` for numpy 2.3.5, scipy 1.16.3, sympy 1.14.0 and mpmath 1.3.0. The quick checks of this report take under 30 s:

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927
python3 synthesis/reconcile_stability.py && python3 synthesis/check_quoted_numbers.py
python3 critic_corrections.py --check-only
python3 synthesis/make_manifest.py --check
```

---

**Note for this public copy:** the raw simulation arrays (`.npz`, 113 files, about 195 MB) are not included in the public repository because of their size. `MANIFEST.sha256.json` still lists them, so a manifest check will report them as missing. Every JSON result, script, log and report is included. The arrays are kept in the author's archive and are available on request.
