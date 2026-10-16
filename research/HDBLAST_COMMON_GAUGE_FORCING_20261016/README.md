# Constructive second-order continuation for one aligned quadrupole

16 October 2026. Conditional research checkpoint for the retained smooth
regular Einstein–scalar background, fixed bulk potential, one fixed shell
tension, and the selected registered aligned quadrupole.

A convergent causal integral series constructs the full second metric–scalar
coefficient for this state and its retained compatible second initial data.
The construction uses the actual transformed first metric and scalar in one
relative-harmonic gauge. It includes the full quadratic bulk forcing, every
physical shell row, the normal gauge row, and the initial and subsidiary
constraint joins. Taking identical copies supplies the copy-even matching
conditions. This is a second perturbative coefficient of the equations,
with the parameter convention

\[
g_\epsilon=g_0+\epsilon q+\tfrac12\epsilon^2k,\qquad
\phi_\epsilon=\phi_0+\epsilon\eta+\tfrac12\epsilon^2w_2.
\]

The controlled causal domain is

\[
\ell=10^{-6},\qquad x=t/\ell,\quad r=z/\ell,\qquad
0\le x\le3,\quad r\ge0,\quad x+r\le8,
\]

times the entire unit three-sphere, through first shell contact. The first
fields remain the actual solutions of the retained master and gauge integral
equations. Interval boxes bound those fields; selections inside the boxes
do not define replacement solutions.

## The convergent second coefficient

The complete generated angular tensor-section space has dimension **287**.
It retains all fifteen symmetric metric fibre components and the scalar,
including the generated transverse-traceless tensor sectors. The number 287
counts the finite section space, not 287 distinct harmonic degrees. No
second-order scalar-only metric ansatz is imposed.

On a declared norm of angular Hilbert values and first scaled derivatives,
the exact physical half-line integral equation is

\[
U_2=\mathcal B+\mathcal T U_2,\qquad U_2=(k,w_2),\qquad
\|\mathcal T\|\le q_*=\frac{2109}{10^6}=0.002109<1.
\]

Here \(\mathcal B\) contains the retained second initial fields and the
actual bulk and boundary forcing. The operator includes the separate
Dirichlet and Robin rows; its reflected kernels and boundary lift preserve
the physical boundary and flat initial corner. Thus

\[
U_2=\sum_{j=0}^{\infty}\mathcal T^j\mathcal B,\qquad
\left\|U_2-\sum_{j=0}^{N-1}\mathcal T^j\mathcal B\right\|_X
 \le\frac{q_*^N}{1-q_*}\|\mathcal B\|_X.
\]

This defines a unique actual solution of the forced second-coefficient
problem, rather than an arbitrarily chosen element of source intervals.
The scalar equation, physical shell laws and quadratic relative-harmonic
constraint are retained together. The subsidiary system has zero compatible
initial and boundary data, so its constraint field vanishes. The natural
second Einstein–scalar equations and copy-even junctions then follow.

The proof also establishes qualitative smoothness from the smooth packet,
background, characteristic equations and boundary trace ODEs, including the
master’s dynamic boundary law. This smoothness is separate from the orders
whose numerical derivatives have been bounded.

The explicit coarse comparison-norm cap is

\[
\|U_2\|_X<1.02\times10^{36},\qquad
\|U\|_X=\max\{\sup|U|_{\mathcal H}/3,
                         \sup|U_x|_{\mathcal H},\sup|U_r|_{\mathcal H}\}.
\]

The angular Hilbert norm uses a positive Frobenius comparison norm and
probability-normalized sphere measure; it is not gravitational energy.
The bound is dominated by the inherited second initial scalar one-jet.
It is an upper bound, not an evaluated response amplitude or a physical
signal. Finite-dimensional evaluation also gives every adapted-frame
component below \(6\times10^{37}\) and each scaled frame-component
derivative \(\partial_x,\partial_r\) below \(2\times10^{37}\).
The numerical response bounds cover values and these first derivatives.
No numerical second-response waveform or candidate truncation has been
evaluated, and no numerical higher-derivative response bound is asserted.

## Full forcing in the same first gauge

The first gauge vector transforms both the metric and scalar. The calculation
retains the mixed and anisotropic angular metric components, shell clock and
angular transport, and perturbed unit-normal normalization. The continuous
strict upper bounds include:

| Affine quadratic source | Each adapted background-frame component |
|---|---:|
| Natural Einstein-equation residual, all fifteen independent rows | \(9.444\times10^{20}\) |
| Registered reduced Einstein-equation residual, all fifteen rows | \(9.444\times10^{20}\) |
| Scalar bulk row | \(3.855\times10^{15}\) |
| Quadratic relative-harmonic vector | \(5971\) |
| First covariant derivative of that vector | \(2.684\times10^9\) |
| Physical metric shell residual, all ten tangential rows | \(1.4512\times10^9\) |
| Physical scalar shell residual | \(4.3539\times10^9\) |

These are absolute individual component bounds in the stated frames and
model units. They are not energy norms or measured amplitudes. The shell
numbers use the outward one-copy residuals
\(K_{ab}-\sigma\gamma_{ab}/6\) and
\(n\phi+\sigma'/2\). The reports give the explicit factors for the two-copy
normal sums and coordinate boundary equations. The first physical shell
rows eliminate normal derivatives from the quadratic physical source values
while retaining the full anisotropic tangential tensor. First intrinsic
tangential covariant derivative bounds are also supplied for those sources.

Every quadratic source is a second parameter derivative along the affine
first-pair path. The coefficient of \(\epsilon^2\) is one half the second
jet, and the second equation has forcing minus the residual Hessian. The
direct source calculation uses the certified first master and gauge jets
through order three, without separately bounding a covariance expression
that would expose fourth gauge derivatives.

## Proof and verification map

- `evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`: complete finite-angular
  operator, causal series, bounds, qualitative regularity and subsidiary joins.
- `bulk-producer/COMMON_GAUGE_BULK_BOUND_PROOF.md`: transformed first pair,
  full natural/reduced source recipes and continuous exact bounds.
- `shell-producer/COMMON_HARMONIC_PHYSICAL_SHELL_FORCING.md`: registered
  normal/tension signs, all ten physical tensor rows, scalar row and adapters.
- `bulk-review/`: separate field-jet, tensor-source and bound checks.
- `boundary-review/BOUNDARY_COMPATIBILITY_AUDIT.md`: physical/gauge rows,
  normal gauge forcing and compatibility identities.
- [Result status](RESULT_STATUS.md): how this local construction advances
  the earlier task description and differs from the inspected literature.
- `literature/COMMON_GAUGE_FORCING_PRIMARY_COMPARISON.md`: bounded
  primary-reference comparison, preserving its original task-stage context.

The source programs use Python’s standard library. Exact finite checks and
deliberately wrong controls supplement the all-point derivations. Separate
internal AI-assisted review is not external specialist peer review.

## Reproduce

`REPLAY_ROSTER.json` declares the exact checker commands, modes and complete
expected reports. The frozen `reproduce.py` authenticates the source inventory
and writes every replay result outside the source tree. Obtain the manifest
SHA256 from the separate
`hdblast/continuation/2026.10.16-common-gauge-forcing/PUBLICATION_PINS.json`
or its dedicated CI workflow, then run

```sh
python -I -S -B reproduce.py --manifest-sha256 THE_EXTERNALLY_AUTHENTICATED_HASH --output-dir /tmp/hdblast-common-gauge-replay
```

Choose a new output directory. The separate publication verifier authenticates
the ZIP and all member bytes before extraction and replay. These checks treat
the authenticated earlier background, master, gauge and second-initial-data
certificates as declared premises; they do not silently rerun all earlier
research.

## Scope

The result applies to this selected registered aligned state, its retained
second initial correction, and the stated causal collar. It is conditional
on the retained smooth regular branch and actual first master/gauge
solutions. It constructs the formal second coefficient of the equations;
it does not establish an exact finite-amplitude nonlinear family, a
nonlinear Taylor remainder, an estimate uniform over unrestricted angular
frequencies, an invariant force, or an observable cosmological response.
Internal checks establish neither external novelty nor journal readiness.
No higher-dimensional Big Bang origin is established. Earlier frozen
capsules and the manuscript are preserved; this checkpoint does not include
a journal submission or a new Zenodo version.
