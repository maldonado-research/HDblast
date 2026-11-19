# Independent energy controls

These scripts use the prior public capsule's `native/FORK.json` and
`native/WALLS.zip` directly. They import no preceding scalar response or spectral
solver. Dependencies are Python 3.10+, NumPy and SciPy. Paths can be anywhere;
`NATIVE` below is the directory holding those two files.

From this capsule's root, after generating `solver/BUDGETS.json`, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python check/verify.py --native "$NATIVE" --full
```

The `--full` flag recomputes all independent outputs and the matching single-mode
primary run. Omitting it verifies the supplied numerical outputs. The default
solver paths are relative to this script, not the working directory; overrides
are `--solver /path/to/budget.py` and `--budgets /path/to/BUDGETS.json`. The final
report is `check/CHECKS.json`. No network or repository checkout is needed.

`native_spectrum.py` solves the uncoupled complex parent modes and conformal
clock with proper-time DOP853. The amplitudes are integrated independently on a
4,096-interval proper-time Simpson grid. The frequency rule uses 96 Gauss nodes
on each `[q,Λ]`. Since `dη a u = dτ u`, each amplitude is directly
`∫dτ f u exp(±iωη)`. The positive spectral density is
`(ω²−q²)/(8π²)` and energy is `ω/a_final` times the count. The native radial
weight is `q²dq/(2π²)`. Initial occupations and correlations are reconstructed
from the archived canonical `X,C,P`; no reported prior pair energy is input.
`SPECTRUM.json` covers the three preparations, both grids, both windows and all
matching cutoffs. `SPECTRUM_REFINED.json` doubles the time grid and increases the
frequency rule to 128 nodes at the largest cutoff.

`finite_bath.py` evolves the full canonical transfer matrix of the exact finite
quadratic oscillator model at `y=0,0.05,0.025`. It also evolves the coefficients
through order `y²` by polynomial multiplication of each canonical map. Every
midpoint kick-drift-kick step is symplectic. It uses the native fine-grid
`q=3.7437724285495775`, 24 bath nodes and `Λ=60`, and includes bare and static
counterterm controls for each preparation. Two proper-time grids, 4,096 and
8,192 steps, establish second-order convergence. Richardson extrapolation in
`y²` checks the exact finite-coupling calculation against the series. Separate
Richardson extrapolation in the time step compares the series endpoint and
explicitly integrated work with `SINGLETON_SOLVER.json` from the primary causal
hierarchy. Bath energy is also compared with the independent amplitudes at the
same finite frequency nodes.

The 30 full native spectral energy comparisons agree with the causal solver to
`3.3×10⁻⁹` relative or better. The finite-bath endpoint comparison agrees to
`6×10⁻⁹` absolute; all four grouped work terms agree to `1.5×10⁻⁹`. Canonical
symplectic errors stay below `3×10⁻¹⁴`. Exact finite-`y` extrapolation agrees with
the same discretized series to `9×10⁻⁸` absolute. Numeric limits and per-case
differences are in `CHECKS.json`, with a failing exit status on disagreement.

The exact finite-`y` model is an auxiliary Gaussian bath realization of the
specified order-`y²` Yukawa noise and response kernel. Its higher orders do not
represent a nonperturbative fermionic Yukawa theory. The finite cutoff, native
momentum band and supplied wall remain part of the experiment's definition.

All `.py`, `.md` and `.json` files in this directory are reproducibility assets.
Exclude only `__pycache__/` and `*.pyc`; there are no other scratch files.
