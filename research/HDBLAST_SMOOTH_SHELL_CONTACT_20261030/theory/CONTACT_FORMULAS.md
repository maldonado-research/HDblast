# Physical reflection, active-shell formulas and the remaining boundary defect

30 October 2026. This concerns the derivative-consistent smooth nominal
candidate `V=B0+B1` of the frozen constant-matrix 23-mode system. All signs
below use the registered inward scaled derivative `partial_r`. The actual
boundary law is `PN U_r=M_actual PN U+G`, while the candidate uses the frozen
constant matrix `M0`. They are not outward-normal derivatives.
The affine shell source is denoted `G_shell` below to distinguish it from the
one-variable ridge function `G(z)`. No physical force or full nonlinear solution
is inferred.

## Full source first, reflection second

Let `P` be the diagonal mixed reflection signs, `PN=(I+P)/2` the eighteen
Robin rows, and `PD=(I-P)/2` the five Dirichlet rows. The full matrices `A,B,D`
need not commute with `P`. For initial displacement `u`, its own derivative
`u'`, velocity `v`, and normalized primitives `Vj` of `v`, put

```
F(z)=[(A+B)(u'(z)+v(z))+D(u(z)+V1(z))]/2,
G(z)=[(B-A)(u'(z)-v(z))+D(u(z)-V1(z))]/2.
```

On the physical side the coupled source is
`f(s,y)=F(y+s)+G(y-s)` for `y>=0`. Its extension for the half-line Duhamel map
is `Ef(s,y)=P f(s,-y)` when `y<0`. This reflects the **complete coupled
source**. Applying the matrices after reflecting their input fields would
require matrix identities that are not assumed and generally do not hold.
The zero-collar hypotheses make every primitive normalization at zero
unambiguous and preserve the initial corner compatibility.

For `x,r>=0` the exact map is

```
Kf(x,r)=(1/2) integral_0^x ds integral_(r-x+s)^(r+x-s) Ef(s,y) dy.
```

Let `p` and `m` be the integrals of `Ef` along the upper and lower characteristic
edges. Direct differentiation of these integration limits gives
`Kf_x=(p+m)/2` and `Kf_r=(p-m)/2`. For `r<x`, set `t=x-r`, `a=x+r`. The
following split uses only physical `f` and holds on a continuous post-contact
region, rather than only at one sample:

```
2 Kf = integral_0^t ds [integral_0^(a-s) f(s,y)dy
                        +P integral_0^(t-s) f(s,y)dy]
       +integral_t^x ds integral_(s-t)^(a-s) f(s,y)dy,
p = integral_0^x f(s,a-s)ds,
m = P integral_0^t f(s,t-s)ds +integral_t^x f(s,s-t)ds.
```

At `r=x` the reflected triangle has zero area and its edge has zero length.
The formulas join the unreflected map continuously. Under the inherited smooth
profiles and zero collar they have the previously stated physical-side C2
regularity. This is not a numerical C2 error certificate. An odd output source
extension can fail to be smooth across the shell; the physical-side convention
is essential.

## Exact shell reduction

At `r=0`, the lower edge is the reflected upper edge, so `m=Pp`. Using
normalized first and second ridge primitives `F1,F2,G1,G2` gives

```
H(x)=x F1(x)-F2(x)+[G2(x)+G2(-x)]/2,
p(x)=x F(x)+[G1(x)-G1(-x)]/2,
Kf(x,0)=PN H(x),
Kf_x(x,0)=PN p(x),
Kf_r(x,0)=PD p(x).
```

Here the output projectors occur **after** the matrix actions. A direct proof
also avoids primitive cancellation: the physical source triangle is

```
H(x)=integral_0^x z F(z)dz
       +(1/2) integral_-x^x (x-|z|)G(z)dz,
p(x)=integral_0^x f(s,x-s)ds.
```

For the first term, fixing `z=y+s` in `0<=s<=x,0<=y<=x-s` leaves a length
`z` interval in `s`. For the second, fixing `z=y-s` leaves length
`(x-|z|)/2`. This proves all weights and signs. Differentiating the direct
weighted integral gives the displayed `p`.

Only `u(x),u'(x),U1(x)=integral_0^x u`,
`U2(x)=integral_0^x(x-z)u(z)dz`, and exact smoothed P1 velocity primitives
through `V3` are needed for interval evaluation. Negative-argument primitives
obey `u(-z)=P u(z)`, `u'(-z)=-P u'(z)`,
`Uj(-z)=(-1)^j P Uj(z)`, and the corresponding `Vj` rule. These identities
act on the input profile primitives; they do not move `P` through a coupling.

An active physical-shell point is `(x,r)=(9/8,0)`: it belongs to the retained
cone, and its characteristics meet the nonzero initial data beyond the collar.
The inherited zero-wedge value error factor is `x+r-1=1/8` there. This factor
sharpens a value error inherited in the `X` norm; it does not sharpen its first
derivative errors or supply second derivative control.

## Robin lift and the exact defect of the truncated candidate

For `g(x)=M0 PN B0(x,0)`, the half-line lift is

```
Lg(x,r)=-integral_0^max(x-r,0) g(s)ds.
```

Thus its inward `r` derivative is **positive** `g(x-r)` while its `x`
derivative is **negative** `g(x-r)`. With normalized `U1,V1,V2`, define

```
Z0=PN[u(x)+V1(x)],
J0=PN[U1(x)+V2(x)].
B1(x,0)=PN H(x)-M0 J0,
B1_x(x,0)=PN p(x)-M0 Z0,
B1_r(x,0)=PD p(x)+M0 Z0.
```

The registered `M0` has zero Dirichlet output rows. Since `PN B0_r=0`,
`PD B0=0`, and `PD B1=0`, the complete candidate `V=B0+B1` has exact
homogeneous Dirichlet values and obeys `PN V_r=M0 PN B0`. Consequently its
residual against the **frozen** affine Robin operator is the identity

```
PN V_r-M0 PN V-G_shell = -M0 PN B1-G_shell.
```

There is no omitted plus sign or factor of two. The registered actual
coefficient is enclosed by the pinned response report and need not equal `M0`.
The residual against the **actual** affine boundary law is therefore

```
PN V_r-M_actual PN V-G_shell
    =-M0 PN B1-(M_actual-M0) PN V-G_shell.
```

The coefficient-correction term must be included when enclosing the actual-law
defect; the frozen-law defect alone cannot be labeled the actual-law defect.
For constant `M_actual`, its `x` derivative is
`-M0 PN B1_x-(M_actual-M0) PN V_x-(G_shell)_x`. If `M_actual` depends on `x`,
the additional term `-(M_actual)_x PN V` is required. The exact formal target
satisfies `PN U_r=M_actual PN U+G_shell`; a separate boundary-law evaluation
for that target must multiply its trace enclosure by the actual coefficient
box, not the frozen candidate matrix.

A nonzero interval enclosure of the candidate residual concerns a defect of the
specified truncated candidate. It neither contradicts the inherited error
ledger, which already allows coefficient variation, omitted forcing and
feedback, nor proves that the unknown exact solution violates its boundary
equation.

A future alternative can enforce the frozen homogeneous Robin row exactly.
For `h=PN(B0+Kf)|shell`, solve the trace Volterra equation
`Z=h-M0 integral_0^x Z`. For constant `M0` and zero corner trace this is

```
Z(x)=h(x)-M0 integral_0^x exp[-(x-s)M0] h(s)ds,
W(x,r)=B0(x,r)+Kf(x,r)-integral_0^max(x-r,0) M0 Z(s)ds.
```

No inverse of `M0` is required. Then `PN W_r=M0 PN W` and `PD W=0` exactly.
This does not make the frozen bulk source self-consistent and does not insert
`G_shell`. It is a proposed next reconstruction, requiring validated matrix
exponential/trace quadrature and a fresh candidate-transfer bound before any
numerical error claim. The present checkpoint evaluates the existing first
feedback candidate instead.

## Independent algebra audit

`checks/check_contact_identities.py` uses only Python's standard library and
exact fractions. It imports no upstream feedback implementation. It builds a
three-component polynomial free wave, forms its physical coupled source with
deliberately noncommuting matrices, then integrates that source over split
physical triangles and characteristic line segments. It compares value and
both first derivatives against primitive formulas at six points, including
unreflected, contact, post-contact and two active shell points. Further tests
check the shell projectors, direct polynomial lift derivative signs, Dirichlet data, frozen affine Robin defect and its
`x` derivative, and actual-law coefficient-correction identities including a
variable-coefficient derivative term. A plausible incorrect input-reflection implementation is
rejected. All 34 named groups pass identically in normal and optimized Python.
These are finite algebra controls supporting the general derivation above;
they are not numerical HDBLAST source bounds, external specialist review, or a
proof of a nonlinear spacetime.
