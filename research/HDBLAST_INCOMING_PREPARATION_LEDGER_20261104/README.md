# Incoming preparation constraints and quadratic source budgets

4 November 2026. This checkpoint makes the next physical preparation question
more concrete. It audits the registered first packet, bounds the actual saved
quadratic forcing on the target future cone, and derives constraint and inverse
relations that a specified past preparation would have to satisfy. It does not
supply a physical blast source or solve the past preparation problem.

## What the calculation establishes

**The registered packet is not a purely one-way packet on the present slice.**
Its normalized master has `w(0,z)=q_ell(z)` and `w_t(0,z)=exp(-z)q_ell(z)`.
Both principal characteristics are nonzero. Their conventional positive quadratic
weights are almost balanced: their relative difference is at most `ell^2/4`,
where `ell=10^-6`. These are normalized-master principal weights, not physical
Noether energies or radiation fluxes. Potential scattering or earlier shell
interaction could create mixed present data from an incoming past state; that
history has not been reconstructed.

**The actual retained quadratic source has a continuous integrated budget.**
Exact clipping of the existing continuous source cells to
`D_T={0<=x<=5/4,0<=r<=5/4-x}` and its active part `x+r>=1` gives:

| Model-normalized squared forcing norm | Safe upper bound |
|---|---:|
| Bulk integral `integral_D_T |S0|_H^2 dx dr` | `8526418129275428520152526` |
| Shell integral `integral_0^T |G0|_H^2 dx` | `3708320490000` |
| Combined source budget `B_src` | `4263209064642351515776089` |

Five clipped source cells cover active area `9/32`; the shell source is active
only on the last time interval of length `1/4`. The mean norm retains all five
channels with `H=diag(1,2,1,3,1)`. This postprocessing uses saved continuous
interval bounds, not new source-oracle samples.

For the compatible reduced mean problem, a tilted energy with a boundary-value
state obeys

```
E' <= (301/100)E + (1/2)||S0||_L2_r^2
                    + (5002001/4000000)|G0|_H^2,
E(x) <= 44(E(0)+B_src), 0<=x<=5/4.
```

It retains the actual nonsymmetric Robin rows. The energy is an estimating norm,
not physical injected energy. The auxiliary zero-Cauchy forced component has a
scalar-trace upper bound `27392057158547`. It is not a physical solution by itself:
zero absolute second Cauchy data generally violate the nonzero quadratic
constraints. This trace ceiling is looser than the retained Duhamel bound; no
improved observable prediction or new absolute proper-volume sign is claimed.

**A partial past-null constraint can be solved explicitly.** With a common affine
null generator and declared free second scalar trace and corner data, the
Raychaudhuri constraint determines a second expansion correction by a weighted
integral of first expansion, shear and scalar-gradient squares. A nonzero
integral of those first-order squares forces nonzero second geometry even when
the free second scalar datum vanishes.
It solves one trace constraint, not the full characteristic/shell problem.
Null expansion and logarithmic area are not the present shell expansion or
actual second area; the latter also contains the square of the first variation.

**The current preparation is conditionally identifiable from a full scalar history.**
For a homogeneous difference between admitted families sharing their first
perturbation, within the retained compact metric-mean fiber, the normalized
physical shell scalar-mean history differs from its free initial scalar profile
by a C1 operator
of norm less than `0.000376`. The resulting causal inverse is stable on the
admissible image and determines only the initial segment reached by that history,
not exterior preparations. This is not surjectivity onto arbitrary histories, a past-state
existence theorem, or an inference from measured cosmological data.

## Physical meaning and reproduction

The current shell action contains pure scalar-dependent tension, with surface
stress `S_ab=-sigma*gamma_ab`, density `sigma` and pressure `-sigma`. Its flux
balance changes vacuum tension; it does not by itself describe a thermal radiation
bath. A matter interaction, state, production history and backreaction would need
additional specified physics. The standard conservation identities and literature
comparison keep those ingredients separate from the present model.

Read the [derivation and limitations](DERIVATION.md),
[primary-source comparison](REFERENCES.md), and [next tasks](NEXT_STEP.md).
The working paper is available as [PDF](manuscript/INCOMING_PREPARATION_LEDGER.pdf)
and [editable TeX](manuscript/INCOMING_PREPARATION_LEDGER.tex) for author and
specialist review.

Use Python 3.12.14 and a fresh output directory outside the checkout:

```sh
python -I -S -B research/HDBLAST_INCOMING_PREPARATION_LEDGER_20261104/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-incoming-preparation-ledger-fresh
```

The replay authenticates the new inventory, all 30 inherited closed capsules and
selected legacy inputs before scientific imports and afterward. New arithmetic
uses the standard library; normal/optimized reports, actual rejected subprocesses
and successful controls are counted separately. No new background, evolved-field,
numerical bump or source-oracle queries or Einstein–scalar PDE solves are used.

Existing continuous source enclosures, operator caps, finite propagation,
constraint inverse, linear uniqueness and exact-family development assumptions
remain retained inputs. Authentication does not reprove them. The null data choice
is explicit and does not extend future operator bounds to the past. No physical
preparation law, full past solution, uniform nonlinear amplitude window, physical
scale, radiation production, observation or higher-dimensional Big Bang origin is
established. Earlier capsules and papers remain unchanged; no new Zenodo edition
or journal submission accompanies this GitHub checkpoint.
