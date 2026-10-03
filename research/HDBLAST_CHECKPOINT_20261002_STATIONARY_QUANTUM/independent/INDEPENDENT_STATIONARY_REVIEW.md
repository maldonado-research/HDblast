# Independent review of the registered stationary quantum shell family

The independent calculations support the reported stationary numerical roots for the declared mathematical model. All twelve independent radial solves pass; all eight direct proper-time source evaluations pass 48 explicit gates; a separate audit of the producer's raw output passes 564 gates across 48 source/Hessian rows, six independently solved endpoints, and three deliberately wrong model roots. This is evidence for bounded numerical stationary closure. It does not establish a physical parameter choice, controlled quantum gravity, time evolution, or quantum stability.

## Prospective record and computational separation

The complete source/model registration was publicly frozen at commit `b4f77f5826e0a603400513832aa53b0673df7533` before any new radial root or quantum-source evaluation. `PREFLIGHT.json` records successful verification of all 17 entries in `FULL_REGISTRATION.json` before execution. The independent protocol hash is `60813c5b4ecc881263915c7caed35bfc1e90ca9c0183bf60448c5b8b6da2bb0a`. Its exact source bytes and all numerical gates were frozen before execution.

The model fixes N=1, delta=.001, c=.5975949350280132, eta_ref=-.00008405268308670538, z_ref=.00005924014794328956, and r_ref=.00011848029588657912. The mass is x=r_ref*f^2 with f=1+b*(eta-eta_ref)/2>0, x_phi=b*r_ref*f, x_phiphi=b^2*r_ref/2. This is the inherited generic quadratic interaction, specialized to explicit benchmark parameters. No reference is recomputed during a root search.

The independent subset was selected prospectively: b=±1 and gamma={0,100,1000000}, at two integration resolutions. The four nonzero finer roots supplied the proper-time endpoint subset, each evaluated at two source resolutions. There was no result-dependent endpoint selection.

The root code imports no producer. It integrates the acceleration system (R,R',eta,w), including R''=-R(w^2/4+U/6), rather than imposing the Hamiltonian square root after every step. Its separately evaluated Hamiltonian defect is

    C_H=R'^2-1-R^2(w^2/12-U/6).

Symbolic differentiation gives C_H'=0 under these equations. The numerical monitor therefore independently tests cone initialization and propagation. The exact polynomial in eta retains the tiny inherited cone displacement, which would disappear if phi=1+eta were rounded first.

The root source route uses the analytic resolvent finite part for Q and direct radius variation of the fully normalized spectral action for rho. It does not reconstruct rho from the trace equation. The separate proper-time implementation evaluates the scalar and metric integrals directly from the S4 heat kernel and its first spectral moment. The trace identity is imposed only as a subsequent acceptance check. The proper-time implementation is the frozen independent DESITTER source code, not a producer import.

## Numerical evidence

`ROOTS.json` contains all twelve runs, objective evaluations, diagnostics and control checks. Every optimizer returned success; every verified root passed both normalized junction gates |E_i|/delta<=2e-12 and the Hamiltonian gate. No failed root, integration, or negative-control case was omitted.

| Independent radial diagnostic | Largest observed value | Registered gate |
|---|---:|---:|
| Either normalized junction residual | 1.029e-13 | 2e-12 |
| Relative Hamiltonian defect throughout integration | 3.477e-15 | 2e-10 |
| H^2 refinement difference, relative | 3.615e-13 | 2e-8 |
| Shooting-coordinate refinement difference | 1.642e-12 | 2e-7 |
| Phi refinement difference | Zero after double rounding | 2e-9 |

A rounded zero is not an exact equality or an error certificate. Both independent resolutions use floating-point arithmetic, and the constraint maxima are measured at integrator nodes.

`PROPER_TIME.json` contains eight source evaluations and 48 enforced gates. These include source refinement and cross-method checks, matched trace checks at both settings, explicit spectral-plus-IR tail-bound gates for W/rho/Q, and wrong-W/metric-sign controls.

| Proper-time diagnostic, normalized source units | Largest observed value |
|---|---:|
| Rho/Q discrepancy from independent analytic sources | 3.089e-18 |
| Rho/Q refinement spread | 7.348e-27 |
| Matched trace residual | 2.969e-26 |
| Analytic spectral plus IR tail bound | 1.687e-28 |

Density and action differences use z^2 units; Q differences use z units. The comparison gate is 1e-9 absolute plus 1e-7 relative. The trace gate and tail gates use the corresponding explicitly registered normalized scales. The analytic tail estimates bound those contributions only: the UV heat-series truncation, quadrature error, and finite arithmetic are not certified interval enclosures. The results retain these limitations explicitly.

`PRIMARY_CROSS_AUDIT_V2.json` is a post-run aggregator that compares raw producer rows against the already frozen independent formulas and gates, without trusting the producer validator status. It checks every one of the 48 source rows and all source derivatives rho_z, rho_eta, j_z, j_eta, rho_x, Q_x and Q_z. The largest absolute rho discrepancy is 5.655e-27 and the largest absolute Q discrepancy is 2.118e-22. All derivative checks pass the predeclared 1e-20 absolute plus 1e-8 relative gate.

Across the six selected endpoints, the largest independent-to-producer differences are 1.111e-16 in phi, 1.243e-17 in H^2, and 9.380e-13 in either shooting coordinate. These are far inside the prospective independent comparison gates. The audit also reconstructs the correct residual at each wrong-model root:

| Wrong model | Correct C residual | Correct B residual |
|---|---:|---:|
| Reverse scalar current | -2.78e-22 | 2.98045e-5 |
| Omit metric source | -1.51263e-7 | -3.76e-19 |
| Freeze both sources | -1.15702e-9 | 8.62128e-8 |

All three defects exceed the registered detection threshold. The independent endpoint controls also detect reversing/omitting the scalar current, changing the static pressure sign, and varying the reference with x. At gamma=0 those defects vanish by construction and are marked not applicable. The proper-time wrong-W and curvature-sign controls are detected at all four nonzero endpoint cases.

## Scientific interpretation and limitations

Both actual shell junctions and the regular interior are required. The squared endpoint constraint alone admits solutions without the selected regular scalar profile. The independent acceleration integration, both junction residuals, and the endpoint comparisons supply distinct checks of this issue.

The fixed-reference action must be varied with r held constant. Its exact static stress is rho=W-zW_z/2, p=-rho; the volume variation is essential. The paired scalar current is J=x_phi Q/2, and its finite-feedback derivative includes x_phiphi Q/2. It is generally incorrect to replace J by partial_phi rho. Curvature-dependent and field-dependent loop terms are already included in the exact action; omitting their scalar variations or adding them a second time changes the model.

Exact de Sitter invariance and tracelessness force the projected Weyl tensor to vanish in this warped symmetry sector. An arbitrary dark-radiation constant is incompatible with the ansatz. A regular interior still has to be integrated. Once symmetry is relaxed, this argument supplies no dynamical Weyl closure or perturbation evolution.

The large-gamma roots are formal finite-amplitude solutions of the declared one-loop semiclassical action. Root convergence does not imply a controlled gravitational EFT. In particular the interior curvature scale is of order 1/9 in the inherited length units, despite small shell H. The separately reported Planck hierarchy screen includes that bulk scale; a screen failure is not a numerical-root failure. Conversely a screen pass is a diagnostic, not a proof that omitted quantum-gravity, higher-loop, or additional local operators are negligible.

The matter state is a de Sitter invariant massive Euclidean vacuum. Its static stress is vacuum polarization. It supplies neither a produced-particle population nor radiation transfer, thermalization, reheating, or a trajectory from an unstable initial shell. The selected local positive-mass factor domain excludes the global massless zero. No massless Euclidean vacuum is inferred by continuity.

The numerical evidence establishes neither uniqueness nor a maximal continuation branch. The finite-amplitude roots and static Jacobians do not determine the nonlocal response needed for dynamical stability. Observed refinement agreement is not a rigorous existence theorem or total error enclosure. Physical masses, matching data, and gravitational normalization remain benchmark choices.

## Preserved implementation failures and reproduction

The frozen independent root and proper-time programs completed without failure or post-result edits. Their `roots.exit` and `proper_time.exit` files both contain zero.

The first separately written post-run audit aggregator reached its JSON output step but could not serialize a NumPy boolean. Its exact script, traceback log and exit status are retained as `audit_primary.py`, `primary_cross_audit.log` and `primary_cross_audit.exit`. `audit_primary_v2.py` casts numerical output metadata to ordinary JSON types and provides portable input/output CLI arguments. No source formula, parameter, gate, or root output changed. `DEVELOPMENT_FAILURES.json` records this mechanical correction. Its v2 output passes 564 gates and exits zero.

The producer's originally frozen validator also had a separate metadata keyword collision, producing 48 validation exceptions and aggregate FAIL. That original source/output must remain retained. Its separately documented v2 repair changes metadata naming and is not a change to a scientific equation, root, parameter or gate. This independent raw-row audit did not depend on that validator's status; the complete primary v2 validation outcome is documented in the primary validation evidence.

For a new-location cross-audit:

    python -B audit_primary_v2.py --checkpoint /path/to/checkpoint \
      --primary /path/to/primary/results.json --roots /path/to/independent/ROOTS.json \
      --output /new/path/PRIMARY_CROSS_AUDIT.json

Use the frozen root/proper-time CLIs and frozen registration hash for full independent replay. Never overwrite archived evidence. The original independent execution recorded its actual NumPy/SciPy/mpmath versions in `ROOTS.json`; a separately labeled replay in the repository-pinned environment is a distinct portability check.
