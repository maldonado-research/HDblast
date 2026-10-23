# An evaluated conditional exterior forced-response mean

This is a source-interval and causal-integral certificate for one component of
the full auxiliary forced increment. It is not a computed global waveform, an
exact nonlinear family, a physical force, or a claim of mathematical novelty.
All assumptions on the background, actual first solution, source projection,
relative-harmonic reduction, and inverse are the retained conditional premises.

Write `W=(I-T)^(-1)(KF+LG)` in the registered23-amplitude selected space. The
first solution determining F and G retains its dynamic material-shell boundary.
The complete second operator includes every bulk coupling and its Robin lift.
The increment W has zero initial data; its sum with the auxiliary homogeneous
comparison gives the registered total second response. Neither summand is
claimed separately to satisfy the physical quadratic subsidiary constraints.

## Closed exterior domain of dependence

Let `t0=1/4`, `r0=5/4`, and

    C={0<=x<=t0, |r-r0|<=t0-x}.

The lower radial boundary is `r=1+x`, so `r-x>=1` throughout C. The upper
characteristic bound is `x+r<=3/2<8`. Thus C lies strictly away from the shell
and inside the controlled domain. If `(s,v)` lies in the backward cone of
`(x,r)` in C, then

    |v-r0|<=|v-r|+|r-r0|<=(x-s)+(t0-x)=t0-s.

Therefore `(s,v)` also lies in C. No reflected source or shell trace can enter
this cone. Both the affine `LG` and the Robin feedback `L(M PN W|b)` vanish
there as exact terms of the original retarded equation. This changes no shell
data and is not a replacement of the full operator by different boundary rows.
The restriction of the actual global W satisfies

    W|C=KF|C+K(A W_x+B W_r+D0 W)|C.

## Full-state coupling bound

Use the same fixed positive probability-sphere/Frobenius norm H and

    ||V||X,C=max(supC |V|H/3, supC |V_x|H, supC |V_r|H).

Endpoint differentiation of the unreflected K gives each derivative bound
`t0 sup|f|H`, and its value bound is `t0^2 sup|f|H/2`. Since `t0^2/6<=t0`,
`||K||C0(H)->X,C<=t0`. The retained individual full-state coefficient bounds
are `||A||<=1e-4`, `||B||<=1e-4`, `||D0||<=1e-6`. Hence

    sigma=203/1000000,
    qC=t0 sigma=203/4000000<1,
    BC=t0 FH,C/(1-qC) >= ||W||X,C.

This is a bound for all23 amplitudes together, including both TT orientations.
No modal decoupling is assumed. At the apex the full feedback value, rather
than its X norm, has the stronger bound

    |(W-KF)(t0,r0)|H <= (t0^2/2) sigma BC =: E.

In particular E bounds the error in any H-dual-norm-one mean functional.

## Actual source enclosure and exact-area integration

The evaluator requires the previously frozen169-member source capsule with
manifest SHA256
`3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2`.
It authenticates every member before and after execution. Its input is the
same actual-source ordinary2-jet oracle, retaining every authenticated scaled
C3 first-field error. It does not substitute the comparator into the source.
Only the exact registered reduced-metric and scalar forcing maps are evaluated.

Set `h=t0/n`, with n=4 in the final calculation. In row i=0,...,n-1 take
rectangles

    [ih,(i+1)h] x [r0-t0+kh,r0-t0+(k+1)h],
    k=i,...,2n-i-1.

Their intersections with C have disjoint interiors and cover C including its
boundary and apex. The two end cells in each row have triangular intersection
area `h^2/2`; every intermediate cell has area `h^2`. There are n(n+1)=20
positive-area cells and total area `t0^2=1/16`. The rectangular source boxes
slightly overbound their clipped intersections, but remain inside the original
controlled exterior source domain. This only widens the source enclosures.

Let lambda take the covariant yy component's probability S3 mean in the
registered frame/gauge. The scalar angular basis is `(1,Z,Fang)` with means
`(1,0,0)`, checked exactly. Consequently `lambda(F)` is exactly the constant
metric-forcing yy coefficient, with normalization `-2 ell^2 rho^2 QEred`.
The fixed functional has dual H norm one. For each cell j the actual source
enclosure supplies `lj<=lambda(F)<=uj` throughout the cell. Therefore

    (1/2) sum_j area_j lj
       <= lambda(KF)(t0,r0)
       <= (1/2) sum_j area_j uj.

There is no floating-point quadrature estimate in this step. The ranges include
both variation inside cells and the actual first-field uncertainty. Every area
and decision is rational. The factor one half is the registered retarded-wave
normalization, not a perturbation-factorial convention.

The full-state FH,C bound uses all forcing coefficient magnitudes from all
cells, exact probability-sphere Gram matrices, and multiplicity two for the
mixed metric blocks. It may be reduced by the accepted full-domain FH bound
`17387513156573`. Summing this norm contribution with the analytic feedback
bound includes every retained bulk amplitude and every feedback iteration.

The evaluated N4 result encloses `lambda(KF)` approximately between
`5.4235689e7` and `6.3966679e10`. The point-value feedback allowance is
approximately `2.7567907e7`. Thus the full auxiliary W yy mean is strictly
positive, approximately enclosed by `[2.6667782e7,6.3994248e10]`.
The exact rational endpoints in the machine receipt govern; these decimals are
explanatory and are not used for the sign decision.

## Continuous rectangle corollary

Set `epsilon=1/1000000` and let R be the closed rectangle

    |x-1/4|<=epsilon, |r-5/4|<=epsilon.

It is contained in the enlarged backward cone with apex
`(tplus,r0)=(1/4+2epsilon,5/4)`: for every point of R,
`|r-r0|<=epsilon<=tplus-x`. This larger cone is still strictly exterior,
with `r-x>=1-2epsilon>0`, and remains inside the controlled domain. Thus the
same full-state causal estimate, now using the accepted global source cap,
gives

    Bplus=tplus*17387513156573/(1-tplus*sigma).

Each of the two mean derivatives of W has magnitude at most Bplus there.
The straight segment from the evaluated point to every point of R is inside
this convex enlarged cone, with L1 displacement at most2epsilon. Hence

    |lambda(W)(x,r)-lambda(W)(1/4,5/4)|<=2epsilon Bplus.

Adding this exact rational Lipschitz allowance to the point interval preserves
strict positivity on the entire stated rectangle. No new source derivative,
sampled-residual assertion, or evaluated waveform is used. The derivative norm
weights are one; its value-only weight three does not enter this corollary.

## Limits

This certifies one coordinate/frame-dependent mean of the full auxiliary W,
at the evaluated point and on its explicitly bounded rectangle,
before shell causal contact. It does not compute the shell-driven response,
certify a complete full-domain numerical residual, construct the physical total
second waveform, identify an invariant observable, or establish a cosmological
origin mechanism. The large interval width is reported rather than hidden.
The result is a standard retarded-integral and norm consequence of the accepted
model-specific source certificates; no general new wave theorem is asserted.
