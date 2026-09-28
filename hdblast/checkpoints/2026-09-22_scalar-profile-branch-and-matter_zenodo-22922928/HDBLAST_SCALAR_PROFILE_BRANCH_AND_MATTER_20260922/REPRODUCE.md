# Reproduction and scope of the checks

Use Python 3.12 with NumPy 2.3.5, SciPy 1.16.3 and SymPy 1.14.0. These versions were used through existing local installations; no global environment was changed. Exact scripts also run with compatible earlier Python interpreters when their dependencies are available. requirements.txt lists the versions.

## Verify the delivered package

From a fresh extraction:

```sh
python3 -B verify_package.py --receipt ../REPLAY_RECEIPT.json
```

For hash verification without numerical dependencies, use `--hashes-only`. The normal verification checks the manifest first, creates a temporary copy, and runs the listed check stages there. The original extraction remains unchanged. Historic absolute paths in provenance are descriptive and are not required for these replays.

This replay includes independent radial integrations and reconstruction of initial data. It does not rerun the full PDE trajectories. The matched-time refinement check advances one archived state by ten small RK4 steps to recover a common comparison time; that bounded operation is stated explicitly in its report. The receipt names all ten executed stages and preserves output/errors.

## Regenerate the static scalar-profile branch

Run in a disposable package copy because producers overwrite their own adjacent result files:

```sh
python3 -B static_branch/solve_plus_branch.py
python3 -B static_branch/independent_branch_checks.py
```

The producer solves six detunings at two tolerances with a first-integral formulation. The checker uses the archived shooting parameters in 18 independent second-order Einstein integrations with three cone starts. It does not refit the roots, and does not produce interval enclosures or a stability spectrum. The scalar η=φ−1 is retained directly, including its ~10⁻²³ cone displacement at δ=.001.

## Regenerate the balanced-data PDE controls

The frozen bulk solver uses exactly the previous checkpoint's continuum equations and spatial operators. The new wrapper samples a compensated seed and a zero-bump reference independently at the same shell-anchored conformal z, then replaces all background coefficients consistently. Numerical dissipation is zero in all new runs.

Example:

```sh
python3 -B evolution/evolve_balanced.py --epsilon .01 --hmin .0001 --stretch 2 --L 3 --tf 1 --output recomputed/balanced_epsp01_wide_h1
```

The completed matrix is:

| Output stem | ε | hmin | stretch | L | tf |
|---|---:|---:|---:|---:|---:|
| balanced_epsp01_h2 | .01 | .0002 | .05 | 6 | 2.5 |
| balanced_epsp01_h1 | .01 | .0001 | .05 | 6 | 2.5 |
| balanced_epsp01_wide_h2 | .01 | .0002 | 2 | 3 | 1 |
| balanced_epsp01_wide_h1 | .01 | .0001 | 2 | 3 | 1 |
| balanced_epsp001_wide_h2 | .001 | .0002 | 2 | 3 | 1 |
| balanced_epsm001_wide_h2 | −.001 | .0002 | 2 | 3 | 1 |
| balanced_epsp01_wide_h05 | .01 | .00005 | 2 | 3 | .5 |

Use `--epsilon=-.001` for the negative value. The model detuning remains δ=.001 throughout; ε is an initial-data parameter and does not change the theory. The first two runs were produced before the wrapper gained its optional `--stretch` argument. That code version is preserved as frozen/evolve_balanced_v1.py; the default executable behavior is the same. Each run records the hash of the wrapper actually executed. `init_plus_h2` is an initial-only diagnostic, not an additional time evolution.

The grid is z_j=−s[exp(hmin*j/s)−1] with s=stretch. Thus equal hmin does not mean equal bulk resolution. RK4 uses a step no larger than .4*hmin. A C2 taper on approximately [−.99L,−.85L] modifies the distant part of the seed; requested times stay below .85L, the earliest continuum taper signal at the shell. The core is nevertheless not the whole domain: outgoing errors can be large elsewhere before any far-boundary signal reaches the shell.

Each NPZ preserves initial/final six-field states, reference profiles and sampled full-state/constraint snapshots. Each JSON stores shell observables and diagnostics about every .025 time units. Logs preserve progress. No coupled matter or quantum field is evolved by this wrapper. The source is not a validated late-time solver.

## Numerical and symbolic reviews

- source_audit/recompute_chat14.py reads the unchanged external Chat 14 ZIP and recalculates event times, comparisons and integrity checks. It never imports the external evolution solver.
- evolution_review/check_wrapper_initialization.py rebuilds the zero and ε=.01 seeds, verifies coordinate/background consistency, and compares numerical accelerations against the initial-profile equations.
- evolution_review/verify_sampling_and_difference.py checks exact paired-profile difference identities and causal scales from copied reference profiles. Those difference equations are a proposed improvement; the present wrapper still subtracts separately integrated profiles.
- evolution_review/analyze_evolution_residuals.py reconstructs the full nonlinear constraint terms from saved fields and velocities, retaining both original and actual-term normalizations. It is deliberately able to report large residuals; successful execution is not a claim those residuals are acceptable.
- evolution_review/matched_time_refinement.py advances the saved coarse state by ten original RK4 steps to coordinate t=0.5 and compares three exactly nested grids at that same time. It reports field self-convergence and the remaining constraint plateau separately.
- evolution_review/check_inverse_tolerance.py shares two dense radial solutions across six initialization-only cases to isolate the effect of coordinate-inversion tolerance. It modifies the sampling function only inside the checking process, never the producer files, and performs no PDE evolution.
- matter/verify_matter_extension.py checks action-derived factors, energy exchange, Gauss/Weyl identities, dimensions, the static expansion, formal occupation integrals and conditional radiation thresholds. Wrong-sign/factor controls must remain nonzero. The occupation integral does not become a validated instantaneous renormalized stress during a nonadiabatic event.
- matter/analyze_crossing_eligibility.py extracts proper-time crossing derivatives from three archived grids. Its free-coupling screening conditions are necessary approximation tests, not sufficient conditions for particle production, acceptable backreaction or reheating.
- build_results_summary.py accounts for the completed runs and their exact recorded values. It does not solve new dynamics.

The main report and independent reviews state remaining scientific limitations. A passing hash check only establishes byte identity. A passing algebra check establishes the stated identity under its assumptions. A successful floating root or integration does not supply an interval proof, a dynamical attractor, or observational evidence.
