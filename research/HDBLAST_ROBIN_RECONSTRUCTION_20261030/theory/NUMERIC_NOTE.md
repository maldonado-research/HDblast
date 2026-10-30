# Validated evaluation of the Robin-resummed boundary subproblem

This calculation replaces the first truncated Robin feedback of the prior
smooth candidate by an integral-defined reconstruction that obeys the actual
affine boundary law. It retains the **same smoothed nominal initial data** and
the frozen principal bulk source `f0=A0 B0_x+Bmat0 B0_r+D0 B0`, where `Bmat0`
denotes the frozen radial-derivative matrix and `B0` the free wave. It does not solve
the complete coupled bulk problem or insert its omitted affine bulk source.

The actual Robin matrix `M` is fixed and constant in shell time, with its
coefficient uncertainty enclosed by the inherited response report. The actual
affine shell source `G` is the fixed inherited continuous function; it is not
replaced by zero or chosen independently from an interval at each operation.
Only its magnitude bounds, not its pointwise values or signs, are needed here.

## Exact reconstruction and finite numerical enclosure

Let `J` denote normalized integration from zero and let

```
h = PN(B0+K f0)|shell,
Z + M JZ = h - JG,
W = B0+K f0+L(MZ+G).
```

The causal lift is `Lg(x,r)=-Jg((x-r)_+)`. The exact integral-defined trace is

```
Z = sum_(j>=0) (-M)^j J^j h
      - integral_0^x exp[-M(x-s)] G(s) ds.
```

Thus `PN W_r=M PN W+G` and `PD W=0` exactly. The zero residual reported by
the checker is this analytic identity with the same fixed `M,G` throughout.
It is **not** an interval subtraction of independently selected boxes or a
claim that a finite series alone satisfies the boundary law exactly. The
numerical partial sum is used only to enclose the integral-defined candidate,
with a proved remainder. The principal equation remains
`W_xx-W_rr=f0`, whereas the full target's source depends on its full field.

At `(x,r)=(9/8,0)` the active shell history has length `d=x-1=1/8` because
the inherited data, source and fields vanish through the zero collar. The
checker keeps terms `j=0,1,2,3` of the homogeneous trace series.

## Endpoint trace moments

Write `Uj=J^j u_s` and `Vj=J^j v_s` for the smoothed initial displacement and
velocity, with `U0=u_s,V0=v_s`. The derivative of `u_s` is used consistently.
For the physical coupled source ridges

```
F(z)=[(A0+Bmat0)(u_s'(z)+v_s(z))+D0(u_s(z)+V1(z))]/2,
G_ridge(z)=[(Bmat0-A0)(u_s'(z)-v_s(z))+D0(u_s(z)-V1(z))]/2,
```

let `Fk,Gk` denote normalized primitives. The endpoint moments of the reflected
bulk trace are

```
Hj=x F_(j+1)(x) -(j+1) F_(j+2)(x)
      +[G_(j+2)(x)+(-1)^j G_(j+2)(-x)]/2,
hj=J^j h=PN[Uj(x)+V_(j+1)(x)+Hj].
```

The output projection is applied after the full physical matrix action.
Negative initial primitives use their proper component parity inside each
action; no parity projector is moved through a coupled matrix. Independent
polynomial controls in the companion derivation test these moment formulas.

All `U1,...,U5` are enclosed together from 64 interval panels on each of the
two analytic initial branches `[65/64,66/64]` and `[66/64,9/8]`. On a panel
`[l,r]`, the pinned nominal interval producer encloses every value of the
expression-defined displacement. The exact nonnegative integration mass for
`Uj` is `[(x-l)^j-(x-r)^j]/j!`. There is no omitted quadrature tolerance.
The published per-component convolution error `alpha_i=M2_i delta_u^2/18`
adds `alpha_i (x-c)^j/j!`, with `c=65/64-delta_u`, covering the entire
smoothed support, including its tiny extension below `65/64`.

The velocity moments through `V6` are exact rational compact-kernel hinge
convolutions of the frozen P1 velocity. Widths remain `delta_u=2^-30` and
`delta_v=1/2048`. Nominal endpoint displacement and radial derivative use the
inherited interval evaluator with their proved convolution remainders.

## Rational remainder bounds

The independently published local cone report gives a uniform cap on the
whole trace history, `||h||_H<=Hcap<1.046e31`. The new checker validates the
actual-minus-frozen Robin coefficient entry bounds and weighted block operator
bounds, obtaining `||M||_H<=m=||M0||_cap+deltaM_cap`. Put `a=md`; then
`a<3.425e-5`. The homogeneous series remainder after degree `N=3` satisfies

```
||R_h||_H <= Hcap a^(N+1)/(N+1)! / [1-a/(N+2)] < 6e11.
```

This follows from `||J^j h||<=d^j Hcap/j!` and the ratio bound on the
remaining factorial series. No matrix exponential is approximated without a
remainder, and no inverse of `M` is used.

For the exact source convolution, the leading `-JG` term is enclosed rowwise
by `[-d Gi_cap,d Gi_cap]`. Its remaining higher terms have Gram-norm cap

```
m d^2 Gcap/2 / [1-a/3] < 1.043.
```

The whole actual source has the inherited uniform bound `Gcap=487122` in the
comparison norm; the scalar-mean leading source contribution has magnitude
below `11786`. These are conservative source enclosures, not source signs.
Zero rows of actual `M` have exactly zero higher matrix tails. In particular,
the TT trace is the isolated affine-source lift, rather than a falsely zero
candidate tensor trace. The five Dirichlet source rows vanish by the pinned
source factorization.

## Candidate value, derivatives and correction from the old approximation

The reconstruction's shell derivatives follow from

```
q=MZ+G,
W_x=h_x-q,
W_r=PD(B0_r+(K f0)_r)+q.
```

For the registered scalar mean, the detailed intervals give the strict
outward-rounded consequences

```
8.29037e30 < W_w:0 < 8.29041e30,
1.97462e32 < partial_x W_w:0 < 1.97463e32,
2.15406e27 < partial_r W_w:0 < 2.15408e27.
```

The change from the preceding first-feedback candidate is enclosed without
subtracting two large trace boxes. With `Z0=PN(B0|shell)` and
`V=B0+K f0+L(M0 Z0)`, the identity is

```
W|shell - V|shell
 = -(M-M0)JZ0 - M J(PN Kf0|shell)
     +sum_(j=2)^3(-M)^j hj + source contribution + homogeneous remainder.
```

Its scalar-mean interval is strictly positive:

```
7.45e20 < W_w:0 - V_w:0 < 1.526e21.
```

The direct trace-series enclosure and the old trace plus this independent
algebraic difference are intersected; they overlap in every row. The report
gives exact rational bounds for all 23 modes and each retained trace moment.
The sign concerns a correction to the chosen approximate formal coefficient,
not a measured force, radiation, energy scale or cosmological signal.

## Authentication, reproduction and limits

The checker authenticates the prior shell, coarse and source reference reports,
the initial producer, initial interval report, response profiles and response
report. The old producer authenticates its upstream roots before scientific
imports and again after evaluation. Input hashes authenticate recorded results;
they do not independently re-prove the inherited field certificates.

Run `check_robin_resummed.py --repo-root PATH`. The prior capsule defaults to
`PATH/research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030`; `--prior-capsule` permits
a separate local copy. Normal and optimized outputs must match byte for byte.
The default reproducible report uses 64 panels per analytic branch.

The boundary identity needs continuous inherited `G` and a compatible zero
corner; a physical-side C2 assertion also requires the inherited C1 regularity
of `G`. The source magnitude enclosure alone does not supply that regularity.
No numerical C2 error, complete coupled bulk solution, pointwise Einstein or
subsidiary-constraint residual, exact nonlinear family, invariant observable or
higher-dimensional Big Bang evidence is established here. The new candidate's
error against the complete formal response requires the separate residual/error
ledger, rather than being inferred from the exact boundary identity.
