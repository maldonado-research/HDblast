# A past–future scalar comparison for a timelike-shell packet

4 November 2026, America/Los_Angeles. Conditional linear Einstein–scalar research note.

For one registered first packet, a graph-energy estimate bounds the difference
between the actual past and future shell scalar coefficients. Applying that
bound to retained future interval enclosures gives two past lobe signs and a
narrower past zero bracket. Earlier papers and capsules remain frozen.

## Results

Use the registered width `ell=10^-6`, initial master displacement `Q(z/ell)`
and nonzero velocity `exp(-z)Q(z/ell)`. The physical first scalar is
`Phi_1(t,Omega)=S(t)(n1^2-n2^2)`. Under the retained branch, reconstruction,
operator and linear-development premises,

```
sup_{0<=t<=3ell} |S(-t)-S(t)| < B,
B = 264000000001/13738809600 < 20.
```

The paired proper shell times are `-rho_b*t` and `+rho_b*t` at the same anchor.
This is an absolute difference of **first linear coefficients in model units**.
It does not bound either individual amplitude by 20 or give a relative error.

The earlier actual future waveform certificate uses the same `ell=10^-6`
packet, coupled field and clock. Its `ell_S` export is `ell*S`; this conversion
must precede interval widening by `B`. Its lobe witnesses imply

```
S(-5ell/4) < -37376,     S(-7ell/4) > 37379.
```

Retained **whole-cell** enclosures yield stronger endpoint signs:

```
S(-6145ell/4096) > 63,
S(-6143ell/4096) < -53.
```

Continuity therefore establishes at least one past amplitude zero in

```
(-6145ell/4096, -6143ell/4096).
```

Its width is `ell/2048`, exactly **1024 times narrower than the earlier coarse
half-ell past bracket**. This comparison concerns those two past brackets.
It does not improve or transfer the old finest future zero brackets, establish
uniqueness or optimality, or preserve the small future tail's sign.

## Why the difference estimate is smaller

The actual master equation and its complete dynamic boundary law are even in
time. Reflection identifies the current past solution with a future development
having the reversed initial velocity. The difference of these two developments
has zero displacement and velocity `-2ell*exp(-ell*r)Q(r)` in scaled coordinates.

A complete energy retains the bulk, endpoint mass and endpoint spring. Its
higher graph energy gives `||A0*Delta Z||<33ell`, while the local boundary
trace satisfies `|Delta q|<ell`. The physical reconstruction then gives the
exact bound `B` above. The nonlocal graph field is not claimed to have compact
support. The original master stays inside `r<=5`; the auxiliary operator agrees
with the actual coefficients through `r<=6`, inside the validated collar `r<=8`.

Time reflection, positive product-space operators, higher graph energies and
widening an interval by an absolute error are established methods. The result
here is their quantitative application to this registered packet and the
explicit transfer of retained actual-field bounds. It is not a new general
mathematical principle or a first amplitude calculation for this packet.

## Materials and verification

Read the [primary proof](primary/PROOF.md),
[independent derivation](independent/PROOF.md),
[continuous-cell refinement](refinement/PROOF.md), and
[closest-prior comparison](REFERENCES.md).
The working note is available as
[PDF](manuscript/PAST_FUTURE_SCALAR_COMPARISON.pdf) and
[editable TeX](manuscript/PAST_FUTURE_SCALAR_COMPARISON.tex).

Selected runtime inputs are included. With Python 3.12.14, from this capsule,
choose fresh external output directories:

```sh
python -I -S -B checks/run_checkpoint.py --output-dir /tmp/hdblast-past-future-fresh
python -I -S -B checks/run_controls.py --output-dir /tmp/hdblast-past-future-controls-fresh
```

The first command authenticates the closed inventory before and after three
mathematical normal/optimized pairs and one manuscript audit pair, requiring
identical canonical bytes. The second exercises six actual source/CLI rejection
cases and two authentication-only positives; these are not mathematical proof
mutants. The new arithmetic uses standard-library exact fractions and selected
source checks. No old producer campaign, numerical PDE, background trajectory,
or field/source/bump oracle is executed.

The old actual waveform and continuous-cell theorems remain **retained premises**.
Selected-file authentication and new interval postprocessing do not reproduce
their full numerical or analytic campaigns. Internal mathematical review and
exact arithmetic are not external peer review or proof-assistant certification.

## Remaining questions

This concerns a conditional first linear scalar coefficient with zero angular
mean. It supplies no physical injected-energy law, radiation or heating,
finite nonlinear amplitude window, incoming source mechanism, SI calibration
or evidence for a higher-dimensional Big Bang.

An earlier project result already constructs a conditional formal complete
second Taylor coefficient through future contact. The outstanding quantitative
particular past preparation, nonlinear remainder and finite-amplitude observable
must be distinguished from a blanket absence of second-order existence. See
[next tasks](NEXT_TASKS.md). Broader novelty and significance still require
independent specialist assessment. No new Zenodo edition or journal submission
accompanies this GitHub working note.
