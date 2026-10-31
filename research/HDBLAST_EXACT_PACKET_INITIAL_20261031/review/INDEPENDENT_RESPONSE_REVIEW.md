# Independent exact-packet finite-response review

31 October 2026. **Pass, conditional on the retained formal-target analytic
premises and the independently checked new initial interval expressions.**
The independent exact-rational ledger reproduces every reported response
node and each retained old C2 specialization. The local true-error gate,
all initial gates, the complete C1 no-increase gate and all three exact
factor-two C2 certificate gates pass. The global factor-two conclusion is a
comparison of complete upper certificates. The separate local comparison
uses a positive lower bound for the actual target norm at the same point.

## Inputs and independence

`independent/check_response_independent.py` imports only Python's standard
library and executes no scientific helper. It builds the ledger afresh as
sums of named nonnegative rational terms. The interval-expression evaluator
and its separate independent evaluation supply the actual new completed
displacement expressions; this response checker additionally recomputes
every leaf's positive-Gram aggregation, exact coverage, channel maxima and
the complete all23 local norm inequalities. It does not infer expression
correctness merely from agreement of two reported decimal bounds.

The selected input/source authentication covers 44 files, including the
checker, initial producer and frozen definition map, response producer,
retained operator/derivative/source reports, registered initial source and
the initial producer's selected closure. Every selected file receives an
exact SHA256 check before and after evaluation. Static retained reports and
maps are read as data; no retained scientific checker is imported.

The relevant immutable SHA256 values are:

| Artifact | SHA256 |
|---|---|
| New `INITIAL.json` | `692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3` |
| Initial producer | `1d7f535ecd5a380cca98408531841ab207b336fd56a249262e7671d06d591896` |
| Frozen initial definition | `aa9a7eae96aae988d3d1a25d83dcd6b436893a9464f7e3813436472a8e137ec1` |
| New `RESPONSE.json` | `991ee9e04d17a879e4e65f3f68ad8f9e76929b4c10d4698d2af06dd6588f3273` |
| Response producer | `c74b48b392d9f3741edc1c80e19aaa1eed3578775d93520c35c7e7fa1126f211` |
| Independent response checker | `48b112a26788f5e715ff439ad10b2166ddebae9c1e65a25c3626d1d3b1c6ad4c` |
| `reference/INDEPENDENT_RESPONSE.json` | `a1d2a554474585ede21f74b64ed22e6036cb0d58041edc2c11ed2f44fb742cbe` |
| Retained `COMBINED.json` | `3e4fb0174df7e88c566346a1b68aedbe4e731fbd062de11eafa28a2cae0a4989` |
| Retained `C2.json` | `85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923` |
| Retained `NORMALIZATION.json` | `e63711829cad2e3dea6ac7e133f754d0e25d04fc40c460c681c7ed1ef3a67a5a` |

Normal and optimized Python produce identical independent canonical bytes.
The same bytes are reproduced when the repository root changes and the
initial/source/report artifacts move inside the published capsule.
Canonical authentication keys use fixed roles for `checker_self`,
`initial_source`, `initial_freeze`, `response_source`, `initial_report` and
`producer_response`; retained sources use repository-relative paths.
Runtime filesystem locations do not enter the canonical result.

## Equations and cancellations checked

Use the unchanged full23 positive H norm on
`D+={0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x}` and
`X(V)=max(sup||V||/3,sup||V_x||,sup||V_r||)`. Write
`T=3/2`, `c=a+b+3d/4`, `c0=a0+b0+3d0/4`,
`dc=delta_a+delta_b+3delta_d/4`,
`comm=cap(A_x)+cap(B_x)+3cap(D_x)/4`,
`h=m/4`, `Pcap=1/(1-h)`, `q=3409/8000000` and
`beta=(q-h)/(1-h)`.

The completed exact displacement and unchanged smoothed P1 velocity vanish
on `r<=1`. The retained causal operators preserve `x+r<=1`, supplying the
same value factor `3/4` and shell history factor `1/4`. Actual coefficients,
actual first field, original unsmoothed `F/F_x`, physical `G` and actual `M`
remain fixed. In particular the derivative input is
`Fx=905405545954284`, never the filtered-source derivative cap.

Let `Dfield=W0-Bs` and `C=W1-W0`. The fresh global C1 route uses

```
Dfield=P[K L0 Bs+L(M PN Bs+G)],
BsX=max((U+T V)/3,U1+V),
I=max((eu+T ev)/3,eur+ev),
D0X=Pcap*((T c0+h)*BsX+Gcap),
R=T*(dc*BsX+c*D0X)+KFcap,
CX=Pcap*R,
E=(I+beta*R)/(1-q),  F0=c*CX.
```

Here the candidate absolute caps are `U,U1,U2,V,V1`, while the separate
actual-minus-candidate caps are `eu,eur,eurr,ev,evr`. The inverse identity
for `Dfield` proves the new bound directly. It requires no inheritance of
old shell samples. The retained `KFcap` is conservative even against the
direct bound `T*FH` for the unchanged original source.

The differentiated calculation retains the cancellation between W0 and
Bs and every initial-acceleration and boundary-source term:

```
Qfree=U2+V1,
f00=a0 V+b0 U1+d0 U,
f0x=(a0+b0)*Qfree+d0*BsX,
D0xX=Pcap*(f00+T*f0x+m*BsX+Gx),
r00=delta_a V+delta_b U1+delta_d U+FH,
r0x=(a+b)*D0xX+(d+comm)*D0X
    +(delta_a+delta_b)*Qfree+(delta_d+comm)*BsX+Fx,
CxX=Pcap*(r00+T*r0x),
F1=(a+b)*CxX+(d+comm)*CX,
w=eurr+a*ev+b*eur+d*eu,
J=max((ev+T*w)/3,evr+w),
Z=(J+T*(F1+comm*E))/(1-q),
Bxx=Bxr=Z,  Brr=Z+c*E+F0.
```

Both corrections have zero initial value and first Cauchy jet, but their
accelerations need not vanish: `Dfield_xx(0)=f0(0)` and
`C_xx(0)=r0(0)`. The positive `f00` and `r00` charges are retained.
Since La is first order, the common first Cauchy jet gives
`r1(0)=La C(0)=0` exactly. The same actual affine `G/M` cancels in the
equations for C and the final target error. Differentiating Dfield still
requires `m*BsX+Gx`. Actual `M_x=B_x=0` is independently checked.

The checker also reproduces the new absolute W1 value/first/second jets
through the prescribed forcing derivative route. Absolute candidate jets
and paired target-error jets remain separate throughout. Substitution of
the authenticated retained old inputs reproduces each old xx/xr/rr cap
exactly, before any new acceptance comparison.

## Results and interpretation

Exact rational decisions precede all displays below.

| Complete quantity | New upper bound, approximate |
|---|---:|
| C1 error X | `1.602148829e23` |
| xx error H | `2.134073127e25` |
| xr error H | `2.134073127e25` |
| rr error H | `2.143224748e25` |

The C1 upper certificate is about 220.97 times smaller than the retained
one. The xx/xr certificates improve by about 31.12 and rr by about 30.99.
Every required factor-two comparison uses its own retained exact rational
cap. The surviving differentiated correction bound is about `7.681e24`,
and the unchanged velocity-gradient error is about `9.672e24`; both remain
charged.

At `P*=(0,65/64)`, the independently aggregated same-point positive target
norm lower is approximately `0.6942780997`; the new paired error upper is
approximately `3.032832001e-13`. Their ratio is at most approximately
`4.368324454e-13`, and the new upper is strictly below half the retained
old true-error lower, approximately `0.6942779566`. This establishes the
predeclared local true-error improvement. It is separate from the global
upper-certificate comparisons.

Only after absolute bounds, the retained outer-endpoint lower scale
`3249000000000000000000000000000000` yields a global rr/max-second-seminorm
upper quotient of approximately `6.596567399e-9`. That endpoint scale is
not used as the local denominator, as a separate actual xx/xr denominator,
or as an invariant or observational relative error.

## Effective controls

`independent/run_effective_controls.py` creates real copied inputs or
checker scripts. Each case executes normal and optimized Python and must
fail before a successful output artifact exists. Fifteen cases, comprising
30 child executions, pass this rejection requirement. Each receipt stores
an explicit nonzero rational witness, failed coverage/role identity or
changed pinned-source hash.

The five initial-report cases omit the unchanged velocity-gradient error,
interchange absolute and error curvature, zero a background channel, skip
the near-edge tail cell, and use the global denominator locally. The four
response-report cases substitute filtered Fx, insert copied old shell
point values, alter a numeric response node, and change the affine shell
error law. The four copied-checker cases omit the frozen initial
acceleration, the bulk initial acceleration, original Fx, or the velocity
gradient in J. The two selected-source cases modify the registered initial
formula file and the actual original-Fx report. Rehashed report fixtures
are rejected semantically, rather than only by their file hashes.

## Limits and retained premises

The original physical-side C2 theorem requires displacement C2, velocity
C1, regular actual coefficients/source and the compatible flat corner.
The componentwise completed-displacement quadratic extension and unchanged
smoothed P1 velocity meet the required data regularity. No global C3
displacement, C2 extension across the artificial reflected shell, or new
actual sixth derivative is claimed.

This ledger-only report produces no response point table or quadrature
values and rejects copied old point tables. Affected response/constraint
tables require fresh evaluation. It supplies no new natural-equation,
frame/connection, clock/history or physical-observable error certificate.

The existing conditional exact initial-data and precontact result remains
available under its stated premises. The auxiliary finite candidate is
not proved to satisfy the physical constraints exactly. A nonlinear
shell-contact family, finite-amplitude remainder and corresponding
existence/compatibility theorem remain unproved. No cosmological or
observational conclusion follows from these bounds.
