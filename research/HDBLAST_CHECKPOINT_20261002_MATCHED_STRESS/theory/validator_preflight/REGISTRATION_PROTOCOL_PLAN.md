# Prospective independent matched-stress validation

These files perform no registered source, response or mode evaluation on import. Before numerical execution, publicly commit and remotely verify the complete registration, including EXPERIMENT.json, both producers, both independently implemented source jets, this validator and raw helper, independent/MANIFEST.json, the inherited common-action W2/W4 protocol, and the stress/derivative-tail library with pulse-norm proof inputs. Pin FULL_REGISTRATION.json by SHA256 and pass the same complete public freeze commit to both producers and the validator. Do not change cutoffs, sources, time grid, tolerances, bounds or control definitions after outcomes.

Place validate_stress.py, raw_stress_audit.py and test_validator_guards.py in the checkpoint code directory. Include them and this plan in the full prospective manifest. Run syntax and the synthetic guard checks normally and with Python -O before freezing. These synthetic checks evaluate arbitrary algebra at an unregistered time, not either registered pulse or a physical response. They never call either numerical producer.

After the public freeze is announced, produce fresh primary and independent results outside the frozen checkpoint. Run the validator twice into distinct fresh output directories:

```
python code/validate_stress.py --primary PRIMARY/results.json --modes INDEPENDENT/results.json --registration-sha256 FULL_REGISTRATION_SHA256 --public-freeze-commit PUBLIC_FREEZE_COMMIT --output CHECK_NORMAL
python -O code/validate_stress.py --primary PRIMARY/results.json --modes INDEPENDENT/results.json --registration-sha256 FULL_REGISTRATION_SHA256 --public-freeze-commit PUBLIC_FREEZE_COMMIT --output CHECK_OPTIMIZED
```

The validator requires exact frozen dependency, source, normalization, grid, resource and manifest evidence, authenticates every independent NPZ hash, checks full archive names/shapes/dtypes, independently rederives all source jets, and reconstructs all eight finite-band quantities at each archived observation from the actual forced complex modes. It compares saved bare stress, full W2/W4 variations and combined integrands to that reconstruction. It derives neither density nor pressure from trace or conservation.

Registered report blocks contain 180 core cross-route comparisons (36 finite-band points times q, q-prime, q-second, density and pressure), 108 separate baseline/anomaly/current comparisons, 72 independent direct-stress-versus-closed comparisons, 36 finite-K trace checks, 72 archive reconstructions (both resolutions), 72 direct-history Simpson Ward endpoints, 36 Ward refinement checks and 180 core refinement comparisons. The continuum block contains 288 comparisons; each uses the actual directed-interval omitted-band bound plus separately recorded numerical estimates. The finite-K trace and Ward source use actual Q0,K and anomaly,K rather than continuum replacements.

The fine-resolution Ward endpoints are acceptance-gated at the registered endpoint tolerance; coarse endpoint residuals are reported. Fine-minus-coarse ledgers are separately acceptance-gated at the registered Ward-refinement tolerance.

The Ward ledger is recomputed using the complete archived direct density and pressure history and independently derived source derivatives, then compared with the direct density endpoint. Observation Wronskians are independently reconstructed. The all-time global Wronskian maxima are authenticated producer declarations: full modes are not archived at every intermediate time. The source-free history and saved source-free observation modes must remain exactly zero.

Exactly twenty controls are registered. Five use the same archived modes to change explicit bare density/pressure mass operators, omit pressure W4, or add the xi=1/6 stress improvement to density/pressure without corresponding dynamics/subtraction changes. The improvement pressure is -[A-second - 3 L A-prime + 3 L-squared A]/(6 epsilon); it is an explicit operator mutation, not a complete alternative theory. Ten algebraic diagnostics omit or flip density/pressure local contacts, omit or flip the current quadratic mass-law contact, replace finite-K Q0 or anomaly with continuum contacts, flip a variance-derivative sign, or change a stress conformal factor. Five synthetic guard mutations change the reference-scale declaration, make a raw value nonfinite, replace extended precision with float64, alter initial data or violate the Wronskian. These last controls are format/state-data controls, not numerical excited-state experiments.

Every algebraic/operator control records its actual maximum absolute residual and whether it exceeds the associated operational cross-route tolerance anywhere. A nonzero algebraic sensitivity witness is separate from operational gate rejection. In particular continuum-Q0 substitutions can lie below the registered cross-Q0 tolerance: the validator does not claim those wrong formulas all fail the physical cross-route gate. The five synthetic guard controls must fail explicit exceptions in normal and optimized execution.

Failure output is retained separately, with traceback and Python optimization mode. Never overwrite a run or relax a gate after failure. Passing this calibration establishes only fixed-geometry linear scalar-to-minimal-stress/current susceptibility under the inherited fixed-positive-r common-action prescription. It does not supply metric-scalar/metric-metric kernels, shifted-root kernels, bulk/boundary perturbation or state/initial data and cannot unlock coupled shell evolution, stability, particle yield or heating claims.

The exact twenty control identifiers are:

1. `missing_density_mass_operator` — same archived modes, operator mutation.
2. `missing_pressure_mass_operator` — same archived modes, operator mutation.
3. `missing_pressure_W4` — same archived modes, subtraction mutation.
4. `improved_stress_density` — same archived modes, xi=1/6 operator addition only.
5. `improved_stress_pressure` — same archived modes, xi=1/6 operator addition only.
6. `omit_density_local_contact` — algebraic contact diagnostic.
7. `omit_pressure_local_contact` — algebraic contact diagnostic.
8. `wrong_density_contact_sign` — algebraic contact diagnostic.
9. `wrong_pressure_contact_sign` — algebraic contact diagnostic.
10. `omit_current_quadratic_mass_law_contact` — algebraic current diagnostic.
11. `wrong_current_contact_sign` — algebraic current diagnostic.
12. `continuum_Q0_at_finite_K` — exact finite-contact sensitivity diagnostic.
13. `continuum_anomaly_at_finite_K` — exact finite-contact sensitivity diagnostic.
14. `wrong_variance_derivative_sign` — algebraic derivative diagnostic.
15. `missing_stress_conformal_factor` — algebraic normalization diagnostic.
16. `moving_reference_scale_declaration` — synthetic reference-parameter guard.
17. `nonfinite_raw_value` — synthetic raw-value guard.
18. `insufficient_raw_precision` — synthetic raw-dtype guard.
19. `altered_initial_data` — synthetic initial-data guard.
20. `altered_raw_Wronskian` — synthetic raw-Wronskian guard.
