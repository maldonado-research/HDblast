# A pointwise scalar-trace bound through past shell contact

4 November 2026, America/Los_Angeles. Conditional linear Einstein–scalar
research continuation. Earlier source identifiers retain their UTC dates.

This checkpoint continues the registered first packet through its first
possible past shell contact and bounds the actual first physical scalar
coefficient on the shell. It retains the specified background, initial packet,
dynamic boundary law and coupled reconstruction. It does not specify a remote
physical source or complete nonlinear preparation.

## The quantitative result

Use `ell=10^-6`, `x=t/ell`, `r=z/ell`, and the same normalized first master
`U(x,r)=w(ell*x,ell*r)`. The shell is at `r=0`. Its actual reduced system is

```
U_xx - U_rr + V(x,r) U = 0,
-U_r + Lambda*q_xx + b(x)*q = 0,  q(x)=U(x,0).
```

The actual shell enclosures give positive boundary mass and spring,
`1000<Lambda<2500` and `180ell<b<300ell`. The static collar extends through
`r=8`, and the directly derived potential satisfies
`5500ell^2<V<103000ell^2`. On the chosen past interval `-3<=x<=0`, the
constructed first solution stays within `r<=5`, so an auxiliary coefficient
extension outside the collar does not affect it.

The estimating energy retains both boundary terms:

```
E = (1/2) integral [U_x^2 + U_r^2 + V*U^2] dr
                  + (Lambda/2)*q_x^2 + (b/2)*q^2.
```

The common time dependence of the bulk and boundary gives

```
(1 - 63ell^2/4400)*E(0) <= E(x) <= E(0),  -3<=x<=0.
```

This is a positive mathematical norm with time-dependent coefficients. It is
not identified with physical injected or Noether energy.

**The main refinement controls the physical scalar trace.** The basic energy
alone does not bound its required pointwise normal derivative. A higher graph
energy of the static operator instead gives `||A0 Z||_H<417` and `|q|<2`.
The exact boundary/conjugation cancellation then yields

```
Phi_1(t,Omega) = S(t)*Y(Omega),
|Phi_1(t,Omega)| <= |S(t)| < 250000000,
                          -3ell<=t<=0.
```

Here `Y=n1^2-n2^2`, with `|Y|<=1`, is the registered quadrupole. The bound is
coarse and applies to the actual **first linear coefficient**, in the fixed
model units. It is not a measured amplitude, an SI calibration or a bound for
a finite nonlinear perturbation. The quadrupole has zero angular mean; this
is not a monopole expansion, radiation or heating result.

The higher graph field `A0^(1/2)Z` is generally nonlocal. Its estimate uses the
abstract graph energy; no compact-support claim is made for that field. Finite
propagation is used for the original first master.

## Contact, earlier work and limits

For backward elapsed time `s=-x`, the original packet has support in
`[1-s,2+s]` before `s=1`, and in `[0,2+s]` through `s=3`. The actual dynamic
shell law is retained after contact. The physical first scalar trace is zero
before the front and is nonzero somewhere in each interval immediately beyond
it. That qualitative activation argument adapts an earlier
registered result; its auxiliary contradiction bound is not an actual numerical
lower bound for the reflected response.

The endpoint-mass Hilbert space, positive-form construction, boundary law and
qualitative activation method already occur in earlier HDBLAST work and the
primary literature. They are credited rather than claimed as new principles.
The contribution is the explicit past-time application and the graph-based
pointwise upper bound for this fixed packet. Broader originality and significance
require specialist assessment.

The analytic branch, validity of its shell enclosure, coupled reconstruction and
relevant gauge/uniqueness premises remain explicit. Finite arithmetic and source
authentication do not reprove them. The earlier near-half principal-channel
shares are not extended through reflection by this energy estimate.

A compatible second-order state, nonlinear past development, physical incoming
source law, radiation production, effective nonlinear parameter interval and
cosmological-origin evidence remain unresolved. Earlier papers and capsules are
preserved. No new Zenodo edition or journal submission accompanies this GitHub
working checkpoint.

## Materials and reproduction

Read the [continuous proof](primary/PROOF.md), the
[independent derivation](independent/PROOF.md), the
[closest-prior comparison](REFERENCES.md), and the [next tasks](NEXT_TASKS.md).
The working paper is available as [PDF](manuscript/PAST_SHELL_CONTACT.pdf) and
[editable TeX](manuscript/PAST_SHELL_CONTACT.tex) for author and specialist review.

All selected runtime input copies are included. With Python 3.12.14, from this
capsule directory, choose fresh output directories outside it:

```sh
python -I -S -B checks/run_checkpoint.py --output-dir /tmp/hdblast-past-contact-fresh
python -I -S -B checks/run_controls.py --output-dir /tmp/hdblast-past-contact-controls-fresh
```

The first command authenticates the complete inventory before and after two
mathematical verifier pairs and one manuscript display/diagram audit pair,
requiring identical normal/optimized report bytes. The second checks six
deliberately invalid source/CLI cases and two authentication-only positives.
Those rejection controls are not mathematical proof mutants.

New computation is standard-library rational arithmetic and source checks.
No numerical PDE solution, background trajectory, bump sample or source-oracle
evaluation is used. The continuous proof and retained analytic premises remain
explicit; authentication is not an independent proof of earlier theorems.
