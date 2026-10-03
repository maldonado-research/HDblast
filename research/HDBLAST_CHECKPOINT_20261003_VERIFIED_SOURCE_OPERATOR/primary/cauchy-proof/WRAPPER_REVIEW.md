# Primary wrapper mathematical audit

Verdict: **GO for the reviewed wrapper mathematics and fabricated entry.**

Reviewed wrapper:
`/workspace/hdblast-research-work/verified-integration-20261003/primary-design/primary_route.py`

Wrapper SHA256:
`b32e1a515cfa3f62cd23f520fd201c502a856a82dc81f405b0c0daf6a292e6a6`

The original mathematics audit reviewed SHA256
`334a1cdee4992493231c4646ea2dcf242196ca67542560054a3ddd45424b456f`.
The subsequent authorized change touches only `run_fabricated` and adds a
manufactured operation rehearsal helper. Production module nodes are unchanged
by AST comparison. `MANUFACTURED_REHEARSAL_DIFF_REVIEW.md` records that scope,
zero physical callbacks, context restoration, and normal/-O byte equality.

Engine SHA256:
`de995d2edeca3944400fd90b892521f05c941b6d8a1f75c2a820740b82d995c0`

Scope: source-model and phase-model budget multipliers; additive global
composition; local drift; source-work margins; exact-real projection at k=0;
meaning of the arithmetic baseline; outward interval-width acceptance. Root
owns bootstrap, public freeze, source authorization, exact event counting,
protocol, and the final physical-run decision. This audit invokes only the
wrapper's `run_fabricated()` entry; no physical source callback was executed.

## Budget multipliers

For a normalized panel polynomial p(x)=sum a_j*x^j with half width H and a
proved uniform source tail r, the local model disks are correct:

\[
 R_{M0}=2Hr,\quad R_{Mexp}=2Hr,\quad R_{Mu}=2H^2r.
\]

The exponential phase has modulus one on the real path. The drift bound uses
|Phi(d)|<=d, so integral_0^(2H) d dd=2H^2. The independently audited engine
implements these same factors inside the complete enclosures.

Let A=sum_j |a_j|. The wrapper's `polynomial_l1` uses exact rational outward
absolute upper bounds of the coefficient balls, hence its computed A is a
valid upper bound on all exact target coefficients. With the integrated E/Q
kernel tails e_tail and q_tail, the panel phase disks H*A*e_tail and
H^2*A*q_tail are therefore correct.

For a whole interval of length one, the source disks r, r, r/2 are correct:
the drift bound integrates b-s over the interval, giving 1/2. For the phase
disk in Mu, each panel's local Q tail receives H^2*A*q_tail. The additional
small-phase stable-drift tail is d*3/98! times an M0 polynomial bound 2H*A.
It is included exactly for the branch |2kd|<=1. The large-phase branch uses a
validated full exponential divided by a nonzero exact represented frequency,
so it contributes arithmetic uncertainty but no finite-series model tail.

The rotations multiplying local E/Q moments have exact-target modulus one;
their arithmetic ball radii are retained in the complete output. There is no
missing model-tail multiplier from those rotations.

## Source-work margin

The engine constructs Lg by multiplying the normalized forcing series by the
normalized L series. Its source-work majorant 32 follows from
|L*g| <= (8/27)*64 = 512/27 < 32 on the proved complex disks. Its Taylor
remainder is half the g remainder. The local source-work integral budget
2H*r_work and complete-domain budget r_work are correct. Source work uses M0
and has no exponential phase model error. Its coefficient and arithmetic
enclosures are retained as real Arb balls throughout summation and export.

## Global composition and k=0

The engine's global source moments use the identities

\[
 e^{\lambda(b-s)}=e^{\lambda(b-r)}e^{\lambda(r-s)},
 \quad\Phi(b-s)=\Phi(b-r)+e^{\lambda(b-r)}\Phi(r-s),
 \qquad\lambda=2ik.
\]

These prove the additive whole-domain E and drift sums. The wrapper reads
these source moments with zero fabricated incoming modes and epsilon=1; it
does not infer an integral from an endpoint discrepancy.

At k=0 every source coefficient and integration kernel is real, exp=1 and
Phi(d)=d. Thus Mexp and Mu have exactly zero imaginary components. The
wrapper's projection sets only those source-moment imaginary intervals to
exact zero. M0 and source work are real already. It does not narrow incoming
complex modes or project nonzero-momentum moments. This is a proved target
identity, not a numerical acceptance adjustment. At k=0, reporting zero
phase-model error is correct even though conservative tail padding remains
inside the complete output intervals.

## Baseline and acceptance semantics

The arithmetic baseline suppresses the source Taylor tails and local E/Q
model-tail inflation, and retains coefficient-ball uncertainty, finite
polynomial arithmetic, rotations, and the tiny stable-drift series enclosure.
It is a separately defined finite-polynomial computation. Its enclosure
width is not an upper bound on every rounding effect after model tails are
added, and does not form an additive partition of the final interval width.
The budget explicitly states both limitations. No subtraction, cancellation,
or negative attribution is used.

The acceptance quantities are the complete output radii reconstructed from
outward exact rational endpoints. The sum of real and imaginary half widths
is a conservative L1 radius about the rectangle midpoint, hence also bounds
Euclidean distance from that midpoint. Every whole source moment and whole
source-work moment is included in the 1e-26 width gate. Integers and zero are
exported with explicit reduced n/d syntax, including 0/1.

## Fabricated checks

Command:

```
PYTHONDONTWRITEBYTECODE=1 /workspace/hdblast-cloud-setup/venv-arb/bin/python \
  /workspace/hdblast-research-work/verified-integration-20261003/primary-design/cauchy-proof/audit_wrapper_fabricated.py
```

Result: `FABRICATED_WRAPPER_EXACT_MOMENT_AND_BUDGET_CHECKS_PASS`.

The checks cover all 1170 moment rows and 130 work rows. For the fabricated
piecewise real polynomials, exact rational integration independently confirms
every M0, every k=0 E/drift moment, and every source-work integral. Every local
source and phase multiplier and every whole source and phase multiplier agrees
with an independently assembled exact rational expression. The fabricated
complete-width gate passes. The fabricated entry refuses any callback name
other than `fabricated_source`; no production callback is invoked.

The earlier engine audit separately tested global E/drift orientation against
independent high-precision closed polynomial integrals at zero, tiny, moderate,
and maximal fabricated momentum. This wrapper audit introduces no outstanding
mathematical finding. It is not a physical result, a bootstrap approval, or a
claim that a registered scientific run has passed.
