# Evaluated actual first jets and geometric quadratic sources

Portable standard-library-only numerical source checkpoint. The actual first
master/gauge/background PDE constructions remain authenticated conditional
premises. The new result is evaluated enclosure of their common-gauge metric
and scalar two-jets, full geometric quadratic source rows and physical shell
source rows, with continuous cell coverage and useful signs.

The public interface is:

```python
from fractions import Fraction as F
from actual_first_jets import I, actual_jet
m, s, b = actual_jet(I(F(1,4), F(1,4)+F(1,4096)),
                     I(F(5,4), F(5,4)+F(1,4096)))
```

`m` contains frame-sector amplitudes A,B,C,D,E,L,M and their ordinary physical
y,t,yy,yt,tt derivatives. `s` contains s,sy,st,syy,syt,stt. `b` contains the
actual background rho,c,theta,a,ay,p,py,Vp,Vpp,Vppp,ell. The geometric adapter
adds frame connections. Mode `comparator` is a separate diagnostic and is
never used as a replacement for the actual solution.

Every source run retains 81 boxes covering the complete causal triangle,
64 whole cells covering an evolved radial profile, two narrow evolved bulk
witness cells, and six cells covering the full shell time interval. The
boxes meeting the diagonal are intersected with x+r<=8 explicitly. All
enclosures are complete-cell results in the registered exact angular basis.

The bulk witness has a strictly negative yy Einstein mean-sector coefficient.
Every cell of the evolved radial profile has strictly negative yy mean-sector
coefficients for both QE and QEred. The reflected witness has a strictly
negative Qphi mean-sector coefficient. A shell witness has strictly positive
one-copy scalar affine Hessian mean coefficient. These are coefficient signs;
they are not angular pointwise sign claims. The source report also derives
the positive-bundle bulk forcing C0 bound F_H<=17387513156573, compared with
the retained 171610068644000 bound. No numerical second waveform is claimed.

Run each command with a new output path:

```sh
python -I -S -B actual_first_jets.py --output /tmp/fresh-jets.json
python -I -S -B check_actual_jets.py --output /tmp/fresh-jet-check.json
python -I -S -B exterior_error_caps.py --output /tmp/fresh-exterior.json
python -I -S -B evaluate_actual_sources.py --output /tmp/fresh-sources.json
python -I -S -B run_numeric_controls.py --output /tmp/fresh-controls.json
```

Repeat the first four with `-O`; their complete bytes must match the retained
normal and optimized goldens. The default source-map input is the sibling
`projection/GEOMETRIC_COEFFICIENT_MAPS.json`. Its final bytes/SHA256 are in the
source report. A coarse replay can use `--profile-cells 8 --cover-divisions 1`;
it remains a complete causal cover and an evolved radial profile, but its
report differs from the tighter default golden.

`run_numeric_controls.py` first verifies complete golden identity, then runs
ten direct checkpoint formula/input rejection processes and four separately
identified inherited rejection processes. Every process must fail with its
specific intended reason and must not write an output report. Direct CLI
cases and expected reasons are in `NUMERIC_CONTROLS.json` and the replay
roster. These controls supplement the independent numeric and geometric
reviews; they are not proofs of the retained PDE existence premises.

`ACTUAL_FIRST_JET_AND_SOURCE_PROOF.md` gives the reconstruction, derivative
conversion, analytic flat tails, exterior causal restriction, exact collected
source substitution, shell normalization, whole-domain coverage, and limits.
`ORACLE_FREEZE.json` records the stable oracle used by the independent shell
invariant calculation. The final producer freeze and replay roster bind the
complete public numerical artifact set.

No new first-field Picard iteration, numerical second-response quadrature,
second waveform reconstruction, or full continuous second-response residual
has been certified. The bulk C0 improvement alone does not change the retained
operator, shell Green map or the complete response error budget.
