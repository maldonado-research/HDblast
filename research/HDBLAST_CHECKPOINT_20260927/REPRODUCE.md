# Reproducing the HDBLAST checkpoint of 27 September 2026

## Environment

- Linux, Python 3.11 (tested with 3.11.15).
- numpy 2.3.5, scipy 1.16.3, sympy 1.14.0, mpmath 1.3.0. Install them with `pip install -r requirements.txt`.
- Run single-threaded, so that timings and results match the recorded runs:
  `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`.
- Read-only inputs from outside this folder:
  - the 22 Sept 2026 package at `new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/`, in particular `frozen/registered_solver.py`, `static_branch/solve_plus_branch.py`, `static_branch/PLUS_BRANCH_RESULTS.json` and `source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip`;
  - for the program map and literature scans, files under `new-files/D-Blast 3/`.

  The scripts check the SHA-256 of the frozen solver on import where stated.

**Re-running scripts overwrites the saved outputs.** To compare a fresh run with the published one, copy the workstream folder first (`cp -r stability_gauge_invariant /tmp/sgi && cd /tmp/sgi`), run there, and compare the JSON with the originals or with `MANIFEST.sha256.json`. The copy still needs the read-only inputs above at their absolute paths. `analytic_structure/run_all.sh` ends with `make_manifest.py`, which regenerates that folder's own manifest, so run it only in a copy.

Runtimes are single-core wall times measured on the 4-core machine used for this round.

**Audit corrections come from the audit scripts.** Several numbers quoted in `00_READ_FIRST.md` and `CLAIM_REVIEW.md` are the auditors' corrections (for example ρ_b agreement 8.5×10⁻¹⁴, the B* margin, −1.500 ± 0.001, the preheating profile maximum 5.9×10⁻⁴, the damping factors 25–200, the exact δ⁹ coefficients, the c* family and the d* budget). They are produced by the scripts in each workstream's `verify/` folder; the full audit command block is at the end of each `VERIFICATION.md`. The steps below list the producer scripts plus the audit scripts that generate numbers quoted in the report.

## 1. Stability of the +1 branch, Method A (`stability_gauge_invariant/`, about 14 min)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant
python3 derive_linearized.py      # ~12 s  -> DERIVATION_RESULTS.json (55 checks, 4 wrong-formula controls)
python3 reproduce_backgrounds.py  # ~6 s   -> BACKGROUNDS.json (run first: the next three read it)
python3 spectrum.py               # ~4 min, 2 processes -> SPECTRUM_RESULTS.json
python3 winding.py                # ~4 min, 2 processes -> WINDING_RESULTS.json
python3 controls.py               # ~5.5 min -> CONTROLS_RESULTS.json
```

The audit commands are listed in `VERIFICATION.md`, and the audit scripts are in `verify/`. The B* stability margin quoted in the report comes from `verify/indep_spectrum.py` (~5.5 min) → `verify/indep_spectrum.json`.

## 2. Stability of the +1 branch, Method B (`stability_time_domain/`, about 45 min)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/stability_time_domain
./run_spectra.sh                  # dense spectra, ~15 min
./run_td.sh & ./run_td2.sh; wait  # time-domain batches, ~23 min each on 1 core each
./run_shift.sh                    # shift-invert scans, ~2 min
python3 analyze.py                # -> results/ANALYSIS.json
python3 summarize.py              # -> results/KEY_RESULTS.json; prints 24 pass/fail checks (all pass)
# audit numbers (Richardson order sensitivity, shell excitation ~8.5e-9); needs the audit re-runs of VERIFICATION.md first:
(cd verify && python3 v4_compare_and_extract.py)   # -> verify/v4_compare_and_extract.json
```

## 3. Exact analytic structure (`analytic_structure/`, about 20–25 min)

Run this in a copy of the folder (see the warning above about `make_manifest.py`).

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/analytic_structure
./run_all.sh   # series_expansion.py 8, gegenbauer_structure.py, five run_bvp_scan.py scans,
               # analyze_series_vs_bvp.py, resonance_probe.py, spectrum_leading_order.py, make_manifest.py
```

The individual commands are in `run_all.sh`. The auditor's independent solver and exact δ⁹ coefficients are in `verify/v1_series_independent.py` (28 s) and `verify/v2_bvp_independent.py`; the arguments are listed in `VERIFICATION.md`.

## 4. Instant preheating (`preheating/`, about 10 min)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/preheating
python3 run_preheating.py     > run_preheating.log      # ~7 min -> PREHEATING_RESULTS.json
python3 run_controls.py       > run_controls.log        # ~2 min -> CONTROLS.json
python3 run_correction_fit.py > run_correction_fit.log  # ~1 min -> CORRECTION_FIT.json
python3 build_summary.py      > build_summary.log       # seconds -> SUMMARY.json
# audit aggregates quoted in the report (decay maximum 2.34e-4, profile maximum 5.9e-4):
(cd verify && python3 check_claims.py)                  # -> verify/CHECK_CLAIMS.json
```

## 5. Constraint control (`evolution_constraints/`, about 4–5 h on one core)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/evolution_constraints
./reproduce_all.sh
# or single pieces, for example:
python3 -B run_case.py --hmin 1e-4 --out runs/final/base_h1e-4                    # baseline
python3 -B run_case.py --hmin 1e-4 --project --out runs/final/proj_h1e-4          # remedy A
python3 -B run_case.py --hmin 1e-4 --project --kappa 10 --damp-off 0.02 0.05 --out runs/final/k10p_h1e-4
python3 -B controls.py   # -> CONTROLS.json (17 controls)
python3 -B analyze.py    # -> ANALYSIS.json
# audit recomputation behind the corrected damping factors (25-200) and the field-order caveats:
(cd verify && python3 -B recompute_saved.py)   # -> verify/recompute_saved.json
```

The finest grid (h = 5×10⁻⁵) needs about 18–22 min and several GB of memory per run. Do not run two such jobs at once on a 15 GB machine.

## 6. Mechanism screen (`mechanisms/`, about 15 min plus about 5–10 min per pilot run)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms
python3 m1_weyl_along_chat14.py; python3 m2_energy_budget_and_scales.py; python3 m3_exact_checks.py
python3 m4_dissipative_junction_toy.py; python3 m5_observational_screen.py
python3 m7_initial_shell_vs_c.py; python3 m8_quadratic_tension_tuning.py   # m8 took ~8 min in the audit
# coupled 5D pilots: the full list of 14 commands is in README.md, "Reproduction"; for example
(cd pilot_5d && python3 rolloff5d_matter.py --Y 1 --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 \
    --tag runs/pilot_dstar_Y1_dc1e-2_L17_dzf2e-3)
python3 m6_pilot_analysis.py
```

The corrected results (the signed c-family, the c* pilot, and the d* run with the φ_b cutoff lifted) come from the audit scripts. Run them from `verify/`:

```bash
python3 v1_independent_derivations.py
python3 v3_c_family_through_zero.py
python3 v3b_continue_to_cstar.py
(cd pilot && python3 v5_cstar_pilot_signed.py 12)
python3 v5_analysis.py
# d* pilot re-run and the extended run with the phi_b cutoff lifted (~4-8 min each), then the budget/Weyl analysis:
(cd pilot && python3 rolloff5d_matter.py --Y 1 --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 --tag rerun_dstar_Y1_dzf2e-3 \
  && python3 rolloff5d_matter.py --Y 0 --tag rerun_pilot_Y0 \
  && python3 rolloff5d_matter.py --Y 1 --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 --phimax 1.3 --tag ext_dstar_Y1_dzf2e-3_phimax1.3)
python3 v2_pilot_consistency.py   # -> V2_PILOT_CONSISTENCY.json (69% radiation share, Weyl/radiation 1.756 and 1.23, budget fractions)
python3 v4_arithmetic_checks.py   # -> V4_ARITHMETIC_CHECKS.json
```

`v2_pilot_consistency.py` reads every `*_summary.json` in `../pilot_5d/runs/` and `pilot/`, so run the pilots first.

## 7. Literature sweeps (`literature/`, under 1 min each; no network needed to rebuild)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/literature
python3 theory_checks/rs_thin_brane_crosscheck.py
python3 theory_checks/prior_mention_scan.py
python3 theory_checks/build_theory_md.py
python3 theory_checks/builder_negative_control.py
python3 observations_checks/knee_and_scale_checks.py
python3 observations_checks/project_scan_1yr.py
python3 observations_checks/build_observations_md.py
python3 observations_checks/build_observations_md.py --selftest
python3 methods_checks/certification_pilot.py
python3 methods_checks/build_methods_md.py
python3 methods_checks/builder_negative_control.py
```

The web searches themselves cannot be replayed. Their queries are recorded in the source files.

## 8. Program map (`program_map/`, about 1 s)

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/program_map
python3 -B verify_program_map.py   # -> program_map_checks.json (34 checks)
```

## 9. Synthesis for this report (`synthesis/`, under 30 s in total)

These scripts only read saved outputs. Run them in this order:

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927
python3 synthesis/reconcile_stability.py     # -> synthesis/RECONCILIATION.json (20 checks incl. 4 must-fail controls)
python3 synthesis/check_quoted_numbers.py    # -> synthesis/QUOTED_NUMBER_CHECKS.json (56 numbers + 3 controls)
python3 synthesis/build_results_summary.py   # -> RESULTS_SUMMARY.json (validator refuses bad records)
python3 synthesis/build_claim_review.py      # -> CLAIM_REVIEW.md
python3 synthesis/build_literature.py        # -> LITERATURE_2022_2026.md, synthesis/LITERATURE_MERGED.json
python3 critic_corrections.py                # re-applies the critic corrections to RESULTS_SUMMARY.json, re-renders
                                             # CLAIM_REVIEW.md, runs the critic checks -> CRITIC_CHECKS.json
python3 synthesis/make_manifest.py           # -> MANIFEST.sha256.json (run last)
python3 synthesis/make_manifest.py --check   # verifies every file against the manifest
```

## What reproduction does and does not establish

Re-running the scripts reproduces floating-point results, with the tolerances recorded in each JSON. It does not turn any numerical result into a proof: nothing in this checkpoint is interval-certified. The independent audits (`*/VERIFICATION.md`) re-ran most audited scripts and obtained bit-identical results apart from runtime fields. Exceptions: the long evolutions, the winding contours and the time-domain batches were re-run only in subsets; some `evolution_constraints` re-runs agree to ≤ 10⁻⁹–10⁻⁸ relative rather than bit for bit, with a platform-dependent initial-data difference; the `analytic_structure` BVP scans were replaced by an independent solver rather than re-run; and `mechanisms/m7_initial_shell_vs_c.py` was not re-run (its last row is reproduced by an audit control). The audits also used independent code for every central claim.

**Order matters for the top-level summary files.** `synthesis/build_results_summary.py` and `synthesis/build_claim_review.py` regenerate the pre-critic versions of `RESULTS_SUMMARY.json` and `CLAIM_REVIEW.md`. Always run `python3 critic_corrections.py` after them (it is idempotent); `python3 critic_corrections.py --check-only` verifies the published files without rewriting them.
