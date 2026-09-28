# Reproduction and interpretation

Python 3.12 is recommended. Numerical runs used NumPy 2.3.5 and SciPy 1.16.3; exact symbolic work used SymPy 1.14.0. The symbolic checks also ran in an existing Python 3.9 environment. Standard-library programs include the saved-run analysis, proper-clock check, and package verifier. Dependencies are listed in requirements.txt. No user-wide environment was modified for this work.

## Fast replay

From an extracted package, run:

```sh
python3 -B verify_package.py --receipt ../REPLAY_RECEIPT.json
```

To use an already installed separate SymPy interpreter, add `--symbolic-python /absolute/path/to/python`. To check only payload identity, add `--hashes-only`. No private source path is required by this replay. The verifier checks every manifest payload, then copies the package to a temporary directory and reruns eight bounded check stages there. It does not alter the supplied payload or its original receipts.

These stages analyze the saved runs, check exact coordinate identities, inspect archived reference data, check the new solver's stencils and Jacobian against saved modes, verify exact constraint transport, check initial-data algebra, independently verify seed compatibility identities and saved candidates, and differentiate saved initial profiles. Passing them is not the same as regenerating the numerical science below. The delivery receipt states precisely what was replayed.

## Regenerate the spectra

The four main registered controls use δ=.001, L=6, stretch=.05, and shell spacings .0004, .0002, .0001, .00005. For example:

```sh
python3 -B solver/registered_solver.py spectrum --tdet .001 --hmin .0001 --L 6 --stretch .05 --shift 1.6572 --output recomputed/spectrum_reg_h1
```

Change spacing and output name to regenerate each resolution. The larger-domain run uses L=8 at spacing .0001. The larger-detuning control uses δ=.1, spacing .0004, L=6 and target shift near 1.627. Sparse eigensolver rounding can differ across platforms and repeated runs; the physical comparison uses convergence and residuals, not byte-identical eigenvectors. Only the selected near-target mode was validated. No full-spectrum claim is made.

## Regenerate the six short evolutions

Use saved input modes or freshly generated modes with the identical grid:

```sh
python3 -B solver/registered_solver.py evolve --tdet .001 --hmin .0001 --L 6 --stretch .05 --dc 0 --mode runs/spectrum_reg_h1.npz --amplitude 1e-8 --tf 2.5 --ko 0 --output recomputed/reg_mode_h1
```

The other three mode runs use spacing .0002 and `runs/spectrum_reg_h2.npz`: amplitudes +1e−8 and −1e−8 with KO 0, and +1e−8 with KO .02. All end at t=2.5. Use `--amplitude=-1e-8` for the negative case.

The two shifted-tension controls use no `--mode`, `--dc 1e-7`, spacings .0002 and .0001, KO 0, and `--tf 4.5`. They reproduce growth but their constraints worsen under refinement. They are intentionally retained as a failed nonlinear initial-data acceptance test.

All evolutions use RK4 with step at most 0.4 times the minimum mapped derivative spacing. The producer retains all six final fields plus sampled constraint arrays in each NPZ; JSON stores shell observables and monitors approximately every .025 coordinate-time unit. Recorded runs took approximately 45–287 seconds each on the original machine; timings are not portable guarantees.

The normalized residual denominator is `1+6Hc(z)^2+phi_z(z)^2`; diagnostics use `z>−0.8L` and exclude the first six nodes. Raw constraint arrays accompany snapshots. This normalization is not a fractional error bound on small perturbations. The independent solver review also reports amplitude-scaled and one-sided boundary checks, which are less flattering and must be retained.

`analyze_controls.py` analyzes the supplied `runs/` filenames and rewrites its adjacent analysis outputs. Run it in a copy if preserving the manifest. For fresh runs, copy their chosen outputs into the corresponding filenames in a disposable package copy first. Fit windows and all numerical comparisons are explicit in the script and JSON.

## Regenerate the initial-data experiments

Work in a copy of `initial_data/`, since these producers write adjacent output files:

```sh
python3 -B initial_data/constraint_seed.py
python3 -B initial_data/balanced_constraint_seed.py
python3 -B initial_data/check_sampled_constraints.py
```

The first program constructs seven one-bump examples at ε=0, ±10⁻⁸, ±10⁻⁶, ±10⁻⁴ and tolerance comparisons. Nonzero examples demonstrate the second-time junction obstruction. The second constructs compensated two-bump candidates at ε=0, ±10⁻⁶, ±10⁻⁴ using the second-order Hamiltonian and separate mass-transport formulations. Its zero-reference compensation coefficient is conventional because no disturbance exists at ε=0.

The compact bump is an initial-data parameterization, not a source in the physical PDE. Saved profiles do not by themselves supply a validated high-order interpolation onto an evolution grid. No evolved compensated-seed run is included in this checkpoint.

## Archive and source identity

MANIFEST.sha256.json lists every frozen payload other than itself. The ZIP and its checksum and verification receipt are separate delivery artifacts. The copied Chat 13 reference ZIP is preserved byte-for-byte, including its known omission of five intermediate files relative to its original manifest. Our quick check reads a specified member that is present; it does not claim a complete replay of Chat 13.

Original v1/v2 code snapshots are references only. The original folder 151 v2 was changing externally during inspection; source-audit snapshots and the frozen input hash identify which version each finding concerns. The executable new solver is `solver/registered_solver.py`. A final comment-only correction labels its quintic taper C2 instead of C-infinity-like; no executable line changed after the saved numerical runs. The independent solver receipt was refreshed afterward.

The current package contains no publication action, communication with others, global cosmological proof, interval existence proof for the new seeds, late-time registered roll-off, or matter-production calculation.
