# Exact-packet initial data and conditional response-error bounds

31 October 2026. A separately defined approximation, `W1_exact_packet`,
removes the old artificial displacement cutoff and displacement convolution.
The actual formal second-response target U, original bulk source F and F_x,
shell source G, Robin matrix M, angular Gram norm and physical-side domain
are unchanged. `W1_quadratic`/`W1_dagger` and the distinct source-smoothed
`Z_epsilon` remain frozen.

## Result

At the initial point (x,r)=(0,65/64), the new radial second-derivative
error obeys

```
||U0_rr - W1_exact_packet_rr||_H / ||U0_rr||_H < 4.369e-13.
```

The actual target denominator is separately bounded away from zero, near
0.6942781. The new absolute error is below 3.033e-13. The old candidate's
error divided by the same target norm exceeds 0.999. Comparing the new
upper error with the old certified lower error proves an actual improvement
at this point. This is an initial coordinate derivative in the stated
positive all23 H norm, not a physical or observational relative error.

Complete interval coverage also gives the following conditional upper
certificates. Values below are conservative readable roundings; the JSON
reports retain exact rational inequalities.

| Quantity | Retained old upper certificate | New upper certificate |
| --- | ---: | ---: |
| Initial radial second-derivative error, H supremum | 6.495e26 | 1.388e23 |
| Complete C1 error, X norm | 3.541e25 | 1.603e23 |
| Evolved xx and xr errors, each H supremum | 6.643e26 | 2.135e25 |
| Evolved rr error, H supremum | 6.643e26 | 2.144e25 |

All three new second-derivative upper certificates are more than 30 times
tighter than their retained exact counterparts. This comparison concerns
upper bounds for different approximations. It does not establish a
factor-30 reduction of the unknown global true error. The unchanged
velocity-error charges remain in every complete estimate.

The comparison domain and norm are

```
D+ = {0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x},
X(V) = max(sup ||V||_H/3, sup ||V_x||_H, sup ||V_r||_H),
ell = 10^-6.
```

The same retained global target scale gives a second-seminorm ratio below
6.597e-9. That normalization uses a large outer-endpoint target derivative;
it cannot replace the separate local comparison. Regionwise absolute
initial errors are reported, with relative ratios left undefined or
unresolved wherever a positive target denominator is not proved.

## Candidate and proof structure

The [pre-execution freeze](reference/FREEZE.json) fixes the candidate,
target, gates and bounded refinement budgets. On 1<r<=7/4 use

```
q(r)=exp(4-1/((r-1)(2-r))), E(r)=exp(-ell*r),
```

and set q=0 for r<=1. Insert their exact jets and the retained nominal
background polynomial in the completed registered all23 displacement
formulas, including the five mean overrides. The background difference and
validated derivative remainder channels remain. The shared exact q and E
differences vanish by functional identity before interval evaluation.
Exponential enclosure midpoints never define the candidate.

At r=7/4, continue each fully reconstructed displacement component with
its own exact value, first derivative and second derivative. This gives C2
data; it does not assert global C3 regularity. Continuing only the raw packet
quadratically would fail a required reconstructed scalar trace. The complete
stored smoothed P1 velocity is unchanged and has the required regularity.

The initial and local gates pass at the first 128-bit campaign, with 49
complete producer leaves and no adaptive splits. The second, 192-bit
campaign is therefore unnecessary. A separate evaluator uses its own
arithmetic, direct all23 formulas and 48 active leaves, plus the exact-zero
collar. Its wider bounds independently pass the initial gates. Exact tail
polynomial inequalities cover the flat edge; no tail interval is omitted.

After those gates, define the new finite response by the same retained
causal operators:

```
W0 = P[Bs + K(L0 Bs) + LG],
W1_exact_packet = P[Bs + K(La W0 + F) + LG].
```

The new response ledger recomputes absolute candidate jets and paired
target-error jets separately. It re-establishes the support wedge, Robin
inverse, W0-Bs cancellation, initial acceleration charges and exact r1(0)=0
from the common first Cauchy jet. The original unfiltered F_x bound is used.
A fresh whole-domain bound replaces old candidate-dependent shell histories;
old shell point or quadrature values are not reused.

## Reproduce and assess

Use a complete checkout and Python 3.12.14. Distributed computations use
the standard library and explicit runtime validation. From the repository
root, choose a fresh output directory:

```sh
python -I -S -B research/HDBLAST_EXACT_PACKET_INITIAL_20261031/checks/run_checkpoint.py \
  --repo-root . --output-dir /tmp/hdblast-exact-packet-fresh
```

The runner authenticates the new closed capsule and its selected inherited
source closure, reproduces canonical reports in normal and optimized
Python, and executes the documented deliberate-corruption controls.
Eight complete canonical reports reproduce in both modes. All 98 deliberate
corruption subprocesses are rejected: 92 formula or proof-ledger cases and
six source-authentication cases. These counts describe executed checks,
not independent theorem counts.
Input authentication is distinct from replaying every inherited proof.
Independent evaluators and reviews are internal AI-assisted work, not an
external specialist assessment or an originality determination.

The [bounded literature comparison](literature/METHODS_COMPARISON.md)
explains the established dependency-preserving arithmetic principles and
why unrelated wave-equation results are not imported as shell theorems.

## Limits and next target

These are conditional coordinate error bounds for an auxiliary finite
approximation to the retained formal second perturbative coefficient.
Absolute global bounds remain large. The new auxiliary initial data are
not presumed to satisfy the physical constraints exactly. No new
natural-equation or constraint residual, shell-point value, invariant
observable or cosmological fit is produced here.

Existing analytic target/reconstruction, regularity, physical shell and
subsidiary premises remain explicit. Earlier conditional exact initial-data
families and protected precontact Cauchy evolution retain their own
hypotheses. Full finite-amplitude nonlinear evolution through shell contact
and a higher-dimensional origin of the Big Bang remain unestablished.

The next quantitative task is to isolate the remaining velocity and
coefficient/correction contributions, then define any further candidate as
a new comparison with its own error ledger. No changed velocity or extra
bulk iteration is included in this frozen experiment. Existing manuscripts,
source capsules, websites and other projects are preserved. No Zenodo
version or journal submission accompanies this checkpoint.
