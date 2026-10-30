# Whole-cone C1 transfer to the exact Robin subproblem reconstruction

This calculation bounds the change from the earlier smooth candidate
`V=B0+B1` to an integral-defined reconstruction `W` that satisfies the actual
affine Robin law. It does not replace the frozen bulk source by the fully
self-consistent coupled source. The resulting C1 error bound is **slightly
larger** than the inherited bound for `V`; the benefit is an exact boundary
identity, while retaining the earlier `1e29` error target.

All norms here are the inherited positive modal Gram norm, and
`X(v)=max(sup||v||/3,sup||v_x||,sup||v_r||)` on the retained cone
`0<=x<=3/2`, `r>=max(0,x-5/4)`, `x+r<=7/4`. The actual formal response `U`
and all inherited background, first-field, source and formal-response
certificates remain conditional inputs. This is not a new proof of an exact
nonlinear spacetime or a numerical second-derivative error certificate.

## Exact correction equation and sign

Let `M0` be the frozen rational Robin matrix used by `V`; let `M` be the
actual constant shell matrix enclosed by the authenticated response report.
The affine shell source is `G`. Both matrices have zero Dirichlet rows and
columns. The actual-law boundary defect of `V` is

```
R = PN V_r-M PN V-G
  = -M0 PN B1-(M-M0)PN V-G.
```

The reconstruction retains the same `B0` and the same reflected principal
bulk forcing `f0=A0 B0_x+B0_matrix B0_r+D0 B0`, but solves the boundary trace
Volterra equation exactly. Equivalently, write

```
C=W-V=Lq,
(Lq)(x,r)=-integral_0^max(x-r,0) q(s)ds.
```

At the shell, `C_r=q`, `C_x=-q`, `C=-Jq`, where `Jq(x)=integral_0^x q`.
Therefore the condition `PN W_r-M PN W-G=0` is equivalent to

```
q+M Jq=-R.                                             (1)
```

This fixes the sign of the correction. If `Z=PN W|shell` and
`Vtr=PN V|shell`, the equivalent trace identity is
`Z-Vtr=J R-M J(Z-Vtr)`, and `q=M(Z-Vtr)-R`.
No matrix inverse or commutation with the parity projector is used.

## Zero collar and the short history

The inherited actual formal response is C1 and zero on `x+r<1`.
Its exact Robin law then gives `G=U_r-M U=0` on the shell for `x<1`, with
the endpoint included by continuity. Alternatively this follows from the
inherited flat first-field collar and quadratic source structure. The uniform
source magnitude cap alone would **not** prove this support property.

The nominal displacement and velocity profiles are zero through `65/64`;
the two compact smoothing widths still preserve zero data for `|r|<=1`.
Consequently `V`, `B1`, `R`, and the solution `q` of (1) vanish for `x<=1`.
On the retained cone, a lift only sees `t=x-r<=5/4`. Thus the longest active
boundary history has length `L=1/4`, not the full `3/2` bulk time span.

Suppose `m>=||M||`, `Rcap>=sup||R||` on this history. Taking the sup norm
of (1) gives the rational contraction majorant

```
Q = Rcap/(1-m L),     provided m L<1,
sup||q|| <= Q,
X(C) <= Q,           sup||C|| <= L Q.                  (2)
```

Existence and uniqueness also follow by the contraction mapping for (1) on
the continuous Robin-valued functions over the short history. This concerns
the boundary correction of a frozen-source subproblem, not the nonlinear
Einstein system. For a shorter history ending at `b`, replace `L` by `b-1`.

## Independent cap on the candidate boundary residual

For every initial mesh cell meeting `[0,b+max(delta_u,delta_v)]`, the checker
forms component caps from maximum absolute P1 nodal values plus the inherited
continuous interpolation error for the nominal displacement and its radial
derivative. The P1 velocity cap uses convexity directly. Positivity of the
compact convolution transfers these to bounds `U0,U1,Vcap` on the smoothed
profiles. Taking outward exact square roots after the positive Gram sum gives
their H-norm caps. No new point evaluator or scientific project code is
imported by this checker.

Write `a0,b0,d0,m0` for the frozen matrix norm bounds and
`dm>=||M-M0||`. On the physical source cone,

```
Fcap=(a0+b0)(U1+Vcap)+d0(U0+b Vcap),
T0=U0+(b-1)Vcap,
B1cap=((b^2-1)/2) Fcap+(b-1)m0 T0,
Vcap_trace=T0+B1cap,
Rcap=m0 B1cap+dm Vcap_trace+Gcap,
m=m0+dm.                                               (3)
```

The bulk trace factor uses the active **physical** source triangle, whose
area is `(b^2-1)/2`; source reflection is applied after the complete coupled
matrix action. Its mixed projector is an H contraction. The active area is
not `(b-1)^2/2`. The second term in `B1cap` bounds the frozen Robin lift.
`Gcap=487122` comes from the separately authenticated uniform actual source
report and covers the whole retained shell history. It is a magnitude upper
bound, not a new source point value or sign.

The actual coefficient uncertainty `dm` is retained explicitly. The checker
verifies each actual coefficient interval relative to its frozen entry and
each weighted block Frobenius cap. The block decomposition is orthogonal in
the registered Gram norm; the matrices and error boxes vanish off block.

## Exact arithmetic output

For the longest retained history, `b=5/4`, the approximate display values are

```
B1cap                       = 2.315393725863152e28,
Rcap                        = 6.343031472057630e24,
m L                         = 6.848761186062340e-5,
Q                           = 6.343465920889470e24,
L Q                         = 1.585866480222368e24.
```

Exact rational endpoints in the machine report support the strict readable
bound `X(W-V)<6.344e24`. Independently recomputing the old nominal-convolution
transfer gives

```
Eprior = 2.183648025252848...e28,
X(U-W) <= Eprior+Q = 2.184282371844937...e28
                         < 2.185e28 < 1e29.            (4)
```

The triangle inequality in (4) is a new candidate transfer, not a reduction
of the total certified error and not a measured approximation error.
For the shorter history ending at the sample `b=9/8`, the correction cap is
approximately `2.616504782636379e24` in X and its value cap is approximately
`3.270630978295474e23`. This shorter-history bound does not replace the global
bound when discussing the entire retained cone.

## Change of the reduced coupled bulk residual

The principal wave operator annihilates `C=Lq` on the physical side. For the
reduced formal residual convention
`mathcal_R(U)=U_xx-U_rr-A U_x-B U_r-D U-F`, the change is

```
mathcal_R(W)-mathcal_R(V)=(A-B)q+D Jq,
||mathcal_R(W)-mathcal_R(V)|| <= (a+b+L d)Q.             (5)
```

Here `a,b,d` enclose the **actual** lower-order coefficient operators. The
checker authenticates their interval boxes and verifies their coefficient
error block norms as well. The resulting global change bound is approximately
`1.359576007844843e21`. Equation (5) bounds only the change due to the boundary
correction. It is not a bound on the complete coupled residual, a pointwise
Einstein tensor, the subsidiary equations or the full nonlinear residual.
The principal source remains frozen at the free wave, so those remaining
equations require further work.

## Reproduction and inheritance

Run `check_robin_transfer.py --repo-root PATH`; the default source report is
`PATH/research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030/reference/AFFINE_SOURCE.json`.
An optional `--source-report` overrides its location without changing the
required hash. Normal Python and `-O` output must agree byte for byte. The
checker uses only the Python standard library and explicit failures.

The profile, initial and response reports, the prior source report and both
reports underlying that source bound are authenticated. The code does not
replay their upstream field evolution or source-factorization proofs. The
zero wedge, regularity, exact formal Robin equation and prior C1 certificate
remain inherited mathematical premises. There is no claim of external
specialist assessment, new physical detection or higher-dimensional origin
of the Big Bang.
