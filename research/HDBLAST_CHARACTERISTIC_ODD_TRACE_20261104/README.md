# Characteristic bounds for the past scalar response of a timelike shell

4 November 2026, America/Los_Angeles. Conditional Einstein–scalar research
continuation. Previous papers and capsules remain frozen.

For the same registered `ell=10^-6` packet, the absolute difference of the
actual first scalar coefficients now satisfies

```
|S(-t)-S(t)| < 1/40 = 0.025,               0<=t<=3ell,
|S(-t)-S(t)| < 1/40000 = 0.000025,   49ell/25<=t<=3ell.
```

These are bounds in model units at paired background proper times
`-rho_b*t` and `+rho_b*t`. They concern a difference of **first linear
coefficients**, not either individual amplitude or a finite-amplitude
nonlinear waveform. The earlier graph estimate was below 20; the new
characteristic estimate uses the difference's zero displacement and small
integrated initial velocity mass.

## Two past zeros and a signed coefficient tail

With the retained actual future waveform intervals, the new bounds establish
at least one temporal amplitude zero in each disjoint open past interval:

```
(-1.49999ell, -1.49998ell),
(-1.97ell,    -1.96ell).
```

The result gives **at least two distinct zeros**. It does not prove uniqueness,
optimal localization, a total zero count or monotonicity. On the complete
paired post-pulse interval the same physical coefficient has the signed tail

```
-1/1500 < S(-t) < -1/4000,                 2ell<=t<=3ell.
```

The angular perturbation is `Phi_1=S(t)(n1^2-n2^2)`. The angular factor
changes sign and has zero mean; a negative coefficient tail is not a globally
negative scalar field, isotropic signal, force or radiation result.

## A past first-field quadratic norm

Define `H_-(x)=ell*(w_z(-ell*x,0)-gamma*w(-ell*x,0))`. The incoming trace
estimate gives `||H_--H_+||_L2(0,3)<3ell`. Retaining the cross term yields

```
|integral H_-^2 - integral H_+^2| < 18ell + 9ell^2.
```

With the same positive background prefactor and proper-clock/cosh cancellation,
the retained future certificate therefore gives

```
6.06e11 < C_minus < 7.44e11,
C_minus = (K_background/ell)*integral_0^3 H_-(x)^2 dx.
```

This is a quadratic norm of the first scalar field, not physical injected
energy, a computation of the second scalar coefficient, or a finite-amplitude
nonlinear response law.

## What changed mathematically

Time reflection identifies the registered past solution with a forward
development having reversed initial velocity. Their difference has zero
displacement and unscaled velocity `-2exp(-z)Q(z/ell)`. Its integrated velocity
mass is small even though its pointwise velocity is not small.

The primary proof uses a velocity-only free comparator, its boundary inertia,
and the complete actual-error history. A separate proof controls the actual
incoming characteristic and the damped boundary equation directly through
two positive majorants. Both retain the endpoint mass, spring and physical
reconstruction. Neither infers a pointwise normal derivative from an ordinary
bulk energy norm or substitutes a free-comparator value for the actual field.

Time reflection, characteristic transport, positive Volterra estimates and
interval widening are established methods. The quantitative fixed-packet
application and transferred past conclusions are the focused contribution.
The old actual future waveform/tail/norm theorem remains a retained premise;
its full producer campaign is not rerun here.

## Read and reproduce

- [Working paper PDF](manuscript/CHARACTERISTIC_ODD_TRACE.pdf) and
  [editable TeX](manuscript/CHARACTERISTIC_ODD_TRACE.tex).
- [Primary proof](primary/PROOF.md) and
  [independent direct-characteristic proof](independent/PROOF.md).
- [Primary-source comparison](REFERENCES.md).
- [Separate conditional nonlinear bridge](nonlinear/NONLINEAR_BRIDGE.md) and
  [nonlinear gap ledger](nonlinear/GAP_LEDGER.md).

With Python 3.12.14, choose fresh output directories outside this capsule:

```sh
python -I -S -B checks/run_checkpoint.py --output-dir /tmp/hdblast-odd-trace-fresh
python -I -S -B checks/run_controls.py --output-dir /tmp/hdblast-odd-trace-controls-fresh
```

The first command authenticates the complete closed inventory before and after
two mathematical normal/optimized pairs and one manuscript audit pair, with
identical canonical report bytes required. The second exercises six source/CLI
rejections and two authentication-only positives. These are transport controls,
not mathematical proof mutants. The nonlinear supplement is an independently
reviewed analytic note, not a third numerical producer pair.

Selected public inputs and source pins are included. New arithmetic uses
standard-library fractions and source checks, with no old campaign, numerical
PDE, background trajectory, or field/source/bump oracle. Finite arithmetic and
hashes do not reprove the retained analytic theorems. Internal review is not
external peer review or proof-assistant certification.

## Nonlinear status and remaining work

Earlier project results already supply conditional exact nonlinear families
and qualitative family-wise signs. The separate supplement specifies an
actual-area quadrupole moment and a sufficient central-zero persistence guard,
`|epsilon|*K(a)<106`. The second-derivative bound `K(a)` and effective amplitude
radius are **not computed**. Its proper-time slab covers the central bracket;
it does not establish a nonlinear late tail or second late zero.

For selected equivariant completions with fixed affine amplitude normalization,
the supplement also gives an odd moment and vanishing quadratic coefficient
for that one moment. Other second-order sectors can remain nonzero. This is
not a new nonlinear PDE theorem, a statement about arbitrary completions or
parameters, or a cubic remainder from only `C^2` regularity.

The most concrete next target is an effective observable/remainder bound for
the chosen past completion and physical observers. Incoming preparation,
matter/radiation production, scale and observational discrimination remain
separate. These results do not establish a higher-dimensional Big Bang origin.
No new Zenodo edition or journal submission accompanies this working checkpoint.
