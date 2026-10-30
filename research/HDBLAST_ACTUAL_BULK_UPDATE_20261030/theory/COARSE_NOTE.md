# Independent retained-majorant audit for the first actual-bulk update

This calculation independently checks the conditional norm estimate for a
boundary-preserving actual-bulk update. It imports neither the new update
checker nor a bulk/source evaluator. It uses authenticated JSON and exact
rational arithmetic, recomputes the Gram operator inequalities, smoothing and
initial error, and reconstructs the earlier residual majorant. The same actual
coefficient functions and sources are retained throughout; interval entries do
not describe independently changing coefficients at each use.

## Operators and different initial data

Write `La V=A V_x+B V_r+D V` on the physical half-line. The complete coupled
source is reflected only after matrix action. Let

```
B V = L(M PN V),
P = (I-B)^(-1),
Psi_s(Y) = P[Bs+K(La Y+F)+LG].
```

`Bs` is the same free wave from smoothed nominal data used by the preceding
reconstruction. `M` is the same fixed actual Robin matrix, and `G,F` are the
same actual affine sources, rather than midpoints or zero substitutions.
The first update is `W1=Psi_s(W0)`. Its bulk forcing is evaluated formally at
`W0`, not at `W1`, so a finite update is not yet bulk self-consistent.

The retained formal target `U` uses the actual initial data instead of the
smoothed nominal initial data. Thus it is not the fixed point of `Psi_s`.
One cannot claim `X(U-W1)<=beta X(U-W0)` without an additional initial-data
term. That omission would remove most of the uncertainty after this update.

## Independent norm reconstruction

The compatible zero-wedge comparison space has

```
X(V)=max(sup ||V||_H/3, sup ||V_x||_H, sup ||V_r||_H),
sup ||V||_H <= (3/4)X(V),
sup_shell ||V||_H <= (1/4)X(V).
```

Consequently, using verified actual coefficient H-norm caps,

```
||K La||_X <= k=(3/2)(nA+nB+(3/4)nD),
||B||_X <= b=nM/4,
||P||_X <= 1/(1-b),
Lip(Psi_s) <= k/(1-b).
```

The script checks the positive Gram weights and every weighted block
Frobenius bound for the frozen matrices and coefficient errors. It verifies
all actual entry intervals against their frozen centers and error caps,
including the Robin rows/columns and off-block zeros. No commutation of
`M`, the bulk matrices or the parity projector is used.

Numerically, the exact bounds give approximately

```
k = 0.000321529951748768,
b = 0.0000684876118606234,
k/(1-b) = 0.000321551974075561.
```

For an independent, deliberately conservative route, the checker uses the
older retained majorant `q=3409/8000000`. After explicitly verifying
`k+b<=q<1`, it may assign the extra slack to the bulk term and use

```
beta_bar=(q-b)/(1-b) = approximately 0.000357661883547634.
```

Subtracting `b` from an unrelated total bound would not suffice. Here the raw
coefficient inequalities prove the required relation `k<=q-b` directly.

## Correction and error bounds

Let `r0=La W0+F-L0 Bs`, where `L0 Bs` is the earlier frozen bulk forcing.
The preceding direct ledger bounds `R>=X(K r0)`; this is a response norm,
not a pointwise forcing bound. The script independently reconstructs it from
the frozen free-wave size, the global boundary correction, actual-minus-frozen
coefficient errors, the full-cone feedback bound and the affine response
estimate. In particular `SOURCE` bounds `X(KF+LG)`, so
`X(KF)<=SOURCE+sup||G||_H` is used. No raw pointwise F cap is inferred.

The correction satisfies

```
C1=W1-W0=P K r0,
X(C1)<=D=R/(1-b).
```

Its boundary data are homogeneous: `PN C1_r=M PN C1`, `PD C1=0`.
The affine source cancels in this correction because both reconstructions use
the same actual M and G. The exact values imply

```
D = approximately 1.66731453122847e26.
```

Let `a` be the authenticated free-wave mismatch between actual and smoothed
nominal initial data. The old conservative certificate is
`E0=(a+R)/(1-q)`. Inverting the Robin term before applying the bulk update gives
the independent error recurrence

```
X(U-W1) <= a/(1-b)+beta_bar E0
         = (a+beta_bar R)/(1-q) =: E1.
```

Equivalently, comparison with the smoothed-data fixed point gives

```
E1 = a/(1-q)+beta_bar D/(1-beta_bar).
```

The script checks both forms by exact equality. Its values are

```
a = approximately 3.53286359987874e25,
E1 = approximately 3.54033516534262e25 < 3.541e25,
E0/E1 > 5.70.
```

The number `a/(1-q)`, approximately `3.53436968315997e25`, is a limiting
**upper-certificate** contribution as bulk iterations continue. It is not a
positive lower bound on the true error. The recurrence accounts for the
different actual and nominal initial data at every iteration.

## What the point enclosure does and does not resolve

At the shell sample `(x,r)=(9/8,0)`, the inherited zero-wedge hypotheses give
`||C1||_H<=D/8`. The generic derivative bound is `D`, and for Robin rows

```
|(C1_r)_i| <= (D/8) sqrt(sum_j max|M_ij|^2/H_j).
```

The script widens the existing W0 boxes by these radii to enclose all 23
components of W1. It then charges E1 to enclose U. This is a norm-based
enclosure of the unresolved actual bulk convolution, not a quadrature
evaluation. The report names it accordingly. For zero-M tensor rows the
radial correction is exactly zero, so the old actual-source radial interval
is preserved. Their value and time-derivative corrections need not vanish.

The exact identity

```
D+E1=E0
```

is particularly important. Widening W0 by an unresolved correction and then by
the new error reproduces the old target value/time norm radius. The shared
Robin radial radii obey the same identity. Thus the smaller error certificate
for the new integral-defined candidate does not by itself yield narrower
best pointwise target enclosures. A generic radial box may be narrower than
the old generic radial box because the correction itself satisfies a Robin
law, but it remains weaker than the existing shared-Robin radial enclosure.
A sharper enclosure of the actual causal correction is needed to gain that
numerical information. The script verifies the value/time and shared-Robin
width identities exactly in the complete point table. This is the global
norm-only baseline. The separate [local estimate](LOCAL_NOTE.md) supplies
additional causal-domain information and does narrow the best target boxes,
without evaluating the convolution or determining its sign.

## Reproduction and scope

Run `check_bulk_coarse.py --repo-root PATH` with ordinary and optimized Python;
the exact JSON outputs agree byte for byte. Optional `--initial-capsule` and
`--robin-capsule` select authenticated copies outside the supplied old checkout.
The default paths use the published prior research capsules. No project
scientific code is imported by this independent program.

This result remains conditional on inherited background, first-field,
coefficient/source, zero-wedge and formal integral-equation premises. It does
not evaluate the actual bulk source pointwise, prove a numerical C2 error,
complete the finite update's bulk self-consistency, establish full Einstein
constraint propagation or a nonlinear family, or supply a physical observable.
This is internal mathematical review, not external specialist assessment or
a claim of a new iteration method.
