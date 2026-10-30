# Charging a curvature-matched extension to a new finite candidate

The earlier finite reconstruction W1 is preserved. Define a new auxiliary
candidate `W1_dagger` by changing only the nominal displacement's extension
past the right endpoint R=7/4, before applying its same fixed even smoothing
kernel of width `delta=2^-30`. Write h=r-R. The old outer extension is
affine. The new one is

```
u_ext_dagger(R+h)=u(R)+u'(R)h+(1/2)kappa*h^2,  h>=0,
kappa=u_nom''(R-).
```

Kappa denotes the exact fixed one-sided curvature of the specified nominal
profile. It is enclosed by the authenticated component intervals; it is not
chosen as their midpoint. The interior nominal profile, P1 velocity, kernel,
actual coefficient functions, actual first field and actual affine sources
F,G are unchanged. The physical target U is unchanged. This is a new finite
candidate, not a change of the archived W1 or a new nonlinear solution.

## Exact half-kernel charges

The kernel is `k_delta(s)=35/(32delta)*(1-(s/delta)^2)^3` for |s|<=delta,
and zero elsewhere. For `r=R-d` with `0<=d<=delta`, the displacement change
inside the retained initial footprint is

```
Delta u_s(r)=(1/2)kappa integral_d^delta k_delta(s)(s-d)^2 ds,
(Delta u_s)'(r)=kappa integral_d^delta k_delta(s)(s-d) ds.
```

The scalar factors are nonnegative and increase as d decreases to zero.
Their maxima follow from the exact positive-half moments

```
integral_0^delta k_delta=1/2,
integral_0^delta s*k_delta=35delta/256,
integral_0^delta s^2*k_delta=delta^2/18.
```

Thus, for `kappa_H>=||kappa||_H`,

```
sup ||Delta u_s||_H <= kappa_H*delta^2/36,
sup ||(Delta u_s)'||_H <= kappa_H*35delta/256.
```

These are half-kernel charges. Using the full second moment would miss the
factor two. The parity extension at the negative endpoint is obtained by
the same registered H-isometric reflection. On D+ one has x+r<=R, so both
free-wave initial arguments lie in [-R,R]; no unbounded exterior quadratic
tail is used in this norm estimate. The unchanged velocity and the one-half
weights in the free wave give

```
C_B=Delta Bs_X <= kappa_H*max(delta^2/108,35delta/256)
               =kappa_H*35delta/256.
```

## Propagation through the finite constructions

Let `Bdry=L M PN`, `P=(I-Bdry)^-1`, `k=||KLa||_X`,
`k0=||KL0||_X`, `b=||Bdry||_X<1`. The same actual M,F,G and frozen L0
are used in both constructions. The affine terms cancel exactly:

```
Delta W0=P(I+KL0)Delta Bs,
Delta W1=P[Delta Bs+KLa Delta W0].
```

These identities preserve operator order; no bulk or boundary operator is
commuted through another. Consequently

```
Delta W0_X <= (1+k0)/(1-b)*C_B,
Delta W1_X <= [1+k(1+k0)/(1-b)]/(1-b)*C_B.
```

The frozen comparison constant is recomputed as
`k0=(3/2)(nA0+nB0+(3/4)nD0)`, after independently checking all twenty
frozen bulk/boundary block Frobenius inequalities against the same full23 H Gram. Actual
k and b are authenticated from the earlier full-bulk update certificate.
The unchanged conditional target then obeys the triangle bound

```
X(U-W1_dagger) <= E_old_exact+Delta W1_X < 3.588e25.
```

The numerical ledger has `kappa_H` approximately `3.68937076612e33`,
`C_B` approximately `4.69764843052e23`, and the general-inverse
`Delta W1_X` cap approximately `4.69948131057e23`. The resulting error
upper certificate is approximately `3.58732997845e25`. This compares upper
certificates for two different auxiliary candidates; it does not measure a
change in their unknown true errors. The higher-curvature motivation and
the endpoint error reduction are established separately by the extension
calculation, not by this C1 triangle bound.

## Sharper direct reproof of the same global C1 upper certificate

The triangle comparison above is a valid independent fallback. It is not
the best certificate for the new candidate. Every component of kappa is
bounded by the old global nominal second-derivative cap M2_i. The quadratic
continuation is C1 at R and has constant second derivative kappa outside it;
with the same reflected continuation it remains globally C1,1 with those
same component M2_i bounds. The checker also verifies
`sum H_i M2_i^2 <= M2_H^2`. Therefore the original pointwise smoothing bounds

```
||u_s_dagger-u_nom||_H <= M2_H delta^2/18,
||(u_s_dagger)'-u_nom'||_H <= 35 M2_H delta/128
```

hold on the retained footprint for the new extension as well. No global
exterior value or first-derivative bound is required: the global free-wave
majorant bounds the unchanged interior u_nom,u_nom' and adds these smoothing
errors. This distinction matters because the quadratic extension itself is
not globally bounded on an unbounded exterior interval.

The checker reconstructs each input of the earlier global C1 majorant:

* The actual-minus-interior-nominal first-jet initial errors are unchanged.
  The M2_H smoothing charges are identical. P1 velocity, its approximation
  error, its slope jumps and its exact smoothing charge are unchanged. Hence
  the initial free-wave error and smooth free-wave size majorants are equal
  to the old numerical inputs.
* The complete admitted shell history ends at x=5/4, before the changed
  collar can arrive. Its previous Bs, first-feedback and V traces are
  unchanged exactly. The same M,G and those trace bounds give the same
  global boundary-correction majorant.
* The global first-feedback bound depends only on the same frozen coupled
  operator cap and the same free-wave size cap. It is therefore also valid
  for the new construction. The coefficient-error and feedback residual
  majorants are reconstructed from these caps; F,G and all actual/frozen
  coefficient norms are unchanged. The same bulk residual response upper
  bound R remains valid.
* The actual boundary inverse, retained contraction q and boundary-inverted
  bulk contraction beta have not changed. The formula
  `(initial_error+beta*R)/(1-q)` is re-evaluated for the new candidate and is
  exactly equal as a rational number to the old numerical majorant.

Thus the sharper direct argument proves

```
X(U-W1_dagger) < 3.541e25.
```

This is the **same numerical upper certificate for a different candidate**,
not equality of the two unknown actual errors. The report retains both the
direct proof and the weaker `<3.588e25` triangle comparison. The same direct
residual ledger also re-proves the old numerical reduced-defect majorant
for W1_dagger, rather than asserting the residual field itself is unchanged.

In the separate comparison route, the old exact reduced forcing-defect cap
needs an explicit charge. Since both constructions have the same actual F,

```
Delta(W1-W0)=P K[La Delta W0-L0 Delta Bs],
X(Delta(W1-W0)) <= [k Delta W0_X+k0 C_B]/(1-b),
sup ||r1_dagger-r1||_H <= (k/(3/2))*X(Delta(W1-W0)).
```

The additional forcing-defect cap is approximately `6.477e16`. Adding it to
the old exact cap gives approximately `3.57395021547e22`, still strictly
below the old readable threshold `3.574e22`. The old exact cap is not
silently reused. These are reduced forcing defects, not full natural
Einstein or constraint residuals.

## Causal preservation of previous point results

Inside the retained initial footprint the changed data are supported only
on `[R-delta,R]` and their reflected negative copy. The free-wave change
therefore vanishes when `x+r<R-delta`. Local coefficient multiplication,
the physical reflected causal propagator K and the causal boundary lift L
preserve this zero region. Each finite term and the convergent Volterra
inverse preserve it as well. This is a finite-propagation argument, not an
inference from the small C1 change.

The selected shell event has x+r=9/8<R-delta with an open gap. The old and
new W0,W1 fields agree on its neighborhood, so the PR98 candidate first jets
and the PR99 direct harmonic-constraint values are unchanged exactly. The
first field and QC are unchanged too. Existing selected-point target boxes
need not be widened merely because the new global triangle bound is larger.

In fact all admitted shell times satisfy x<=5/4<R-delta. The perturbation
cannot reach that boundary segment, so P acts as the identity on these
specific differences throughout D+. The report separately records the
slightly smaller support-refined caps
`Delta W0_X<=(1+k0)C_B`,
`Delta W1_X<=[1+k(1+k0)]C_B`. The general-inverse caps above remain valid
and are used for the main readable certificate.

## Verification and scope

The numerical comparison authenticates the initial second-jet report, old
bulk bound, raw full23 response and profiles, and uses exact rational norm
and transfer arithmetic. Normal and optimized outputs match byte for byte.
An independent manufactured checker proves whole-interval half-kernel
charges using Bernstein positivity, exercises noncommuting finite matrices
and common-source cancellation, and checks a causal filtration. Its proof
and limitations are recorded in `CURVATURE_TRANSFER_REVIEW.md`.

No complete bulk-source derivative, proof of smaller true second-derivative
error throughout the domain, quantitative C2 theorem, full Einstein/Hamiltonian/momentum residual,
self-consistent second-order solution, nonlinear family or observational
support for a higher-dimensional Big Bang follows from this comparison.
