# Numerical first-field input and initial-profile checkpoint

This new directory contains no changes to the accepted public capsules.
All mathematical comparisons use exact Fraction decisions and outward rational
intervals. Display-only decimal fields are not used by the proofs.

The new result is a continuous scaled-C3 approximation-error certificate for
the actual retained first master and forced gauge on the complete causal
triangle, an explicit common background approximation with scaled-C4 errors,
and complete-cell evaluation of the actual second initial profiles and full
spacetime one-jet. This is input progress toward the numerical second response.
No numerical second waveform or common-gauge quadratic source has been evaluated.

For each script, use a new output path:

```sh
python -I -S -B check_approximation_errors.py --output /tmp/new-approx.json
python -I -S -B evaluate_background_approximation.py --output /tmp/new-background.json
python -I -S -B evaluate_initial_profiles.py --output /tmp/new-initial.json
```

Repeat with `-O` and compare complete reports. Published local golden filenames
are APPROXIMATION_ERRORS_{NORMAL,OPTIMIZED}.json,
BACKGROUND_APPROXIMATION_{NORMAL,OPTIMIZED}.json and
INITIAL_PROFILES_{NORMAL,OPTIMIZED}.json. Source/control documentation is in
APPROXIMATION_ERROR_PROOF.md and INITIAL_PROFILE_PROOF.md.

`INPUT_PINS.json` binds five small copied inputs by SHA-256 and bytes and
records original absolute paths. The pinned earlier certificates remain
conditional premises; this run does not re-prove their PDE existence results.
No earlier code is imported or executed. The background script imports only
this directory's new initial interval-arithmetic module.

`run_numeric_controls.py --output NEWPATH` verifies complete normal/optimized
byte identity and runs 14 wrong-formula/input processes. Intended rejection
reasons are part of NUMERIC_CONTROLS.json. All rejection tests require a failed
process, matching reason and absence of an output report. The independent
review lives in the adjacent numeric-review directory.

For plotting, INITIAL_PROFILES_NORMAL.json contains 498 complete `profile_cells`
covering r in [1,2], including two analytic endpoint tails. `scalar_mean` is
an interval for the actual probability-normalized sphere mean on the entire
cell. The degree2/degree4 entries use the retained unnormalized scalar angular
basis; they are not orthonormal modal amplitudes. `component_caps` are <=
upper bounds, not claims of strict mathematical inequalities.

The affine background approximation retains inherited shell uncertainty.
The master/gauge comparisons have complete derivative-error budgets and
nonzero stabilized initial jets but have not been improved by numerical
Picard iteration. Actual common-gauge quadratic source values/errors,
orthonormal angular projections, a numerical second-response reconstruction
and its independent continuous residual/error certificate remain required.
