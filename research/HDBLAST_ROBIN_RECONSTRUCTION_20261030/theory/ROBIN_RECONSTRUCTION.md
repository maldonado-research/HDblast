# A boundary-resummed reconstruction with the affine Robin source retained

30 October 2026. This note constructs an integral-defined approximation that
satisfies the registered actual affine Robin law exactly, while retaining the
bulk source frozen at the free smooth wave. It therefore solves a specified
principal-wave boundary subproblem, not the complete coupled second-order
Einstein–scalar problem. The construction is the classical Volterra resolvent;
no originality or nonlinear-existence claim is made.

## Fixed inputs and conditional regularity

Let `PN` project onto the eighteen Robin rows and `PD=I-PN` onto the five
Dirichlet rows. The free reflected smooth wave is denoted `B0`. Define

```
f0 = A0 partial_x B0 + Bmat0 partial_r B0 + D0 B0,
S = B0 + K f0,
h(x) = PN S(x,0).
```

Here `A0,Bmat0,D0` are the frozen rational bulk matrices. The physical source
is coupled first and its completed output is then reflected, as derived in
the preceding shell-contact checkpoint. No commutation of bulk matrices with
the mixed parity projector is assumed. The auxiliary reflected propagator has
`PD S(x,0)=0` and `PN partial_r S(x,0)=0`.

The actual Robin matrix `M` is a single constant matrix with `M=PN M PN`.
Its entries lie within the authenticated actual-coefficient boxes. It need
not equal the previous candidate's frozen matrix `M0`. Throughout the
construction the same fixed actual matrix is used at every time; independent
interval choices at different times do not define the intended reconstruction.

Write the actual affine shell source as `G_shell(x)`, with
`PD G_shell=0`. The inherited zero-collar hypotheses give
`h=G_shell=0` through `x=1`. For the C1 construction, continuity of the source
and the stated regularity of the base suffice. Physical-side C2 regularity of
the lift additionally needs `G_shell` to be C1 and the relevant corner jets to
match. This qualitative regularity is an inherited model assumption; a uniform
source magnitude bound alone does not prove it. The earlier full23 shell and
subsidiary construction supplies the zero-wedge/corner hypotheses, rather
than this note deriving a new first-field regularity theorem.

The retained shell history is `0<=x<=5/4`, so its active length is at most
`L=1/4`. The symbol `L` here is a time interval length, not the physical scale
`ell=1e-6` and not a modal field coefficient.

## Exact trace and boundary reconstruction

For a time function, let `J a(x)=integral_0^x a(s)ds`. The trace is defined by

```
Z = h - J(M Z + G_shell).                                      (1)
```

Because `M` is constant, the unique continuous solution is

```
Z = sum_(k>=0) (-M)^k J^k h
    - sum_(k>=0) (-M)^k J^(k+1) G_shell,                       (2)
```

or, equivalently,

```
Z(x) = h(x) - M integral_0^x exp[-(x-s)M] h(s)ds
             - integral_0^x exp[-(x-s)M] G_shell(s)ds.          (3)
```

No inverse of `M` is needed. Formula (2) converges uniformly on every finite
time interval: if `m>=||M||` and the inputs vanish before `x=1`, then
`||J^k h(x)|| <= d^k H/k!` and
`||J^(k+1)G_shell(x)|| <= d^(k+1) Gcap/(k+1)!`, where
`d=max(x-1,0)`, `H>=sup||h||` and `Gcap>=sup||G_shell||` on that history.
These factorial bounds require no smallness of `m*d`. Termwise integration
in the uniformly convergent series proves (1). If two traces solve (1), their
difference satisfies a homogeneous Volterra equation; iteration gives
`||difference|| <= C (m*d)^n/n!` for every n, hence uniqueness.

With the inward radial derivative convention, set

```
q = M Z + G_shell,
W(x,r) = S(x,r) - integral_0^max(x-r,0) q(s)ds.                (4)
```

The same zero collar makes the lift and its available derivatives vanish
across the initial characteristic corner. At the shell, (1) gives
`PN W=h-Jq=Z`. The inward derivative of the lift is `+q`, while its time
derivative is `-q`. Therefore

```
PD W(x,0) = 0,
PN W_r(x,0) = q(x) = M PN W(x,0) + G_shell(x),
PN W_x(x,0) = h'(x) - M Z(x) - G_shell(x).                    (5)
```

These are exact analytic identities for the integral-defined function, not
small numerical residuals obtained by subtracting independently widened
intervals. A finite series without its tail does not generally satisfy (1)
or (5) exactly. The numerical intervals enclose the function in (2)–(4).

The lift is a function of `x-r` and solves the homogeneous principal wave
equation on the physical side, including the corner under the stated jets.
Consequently

```
W_xx - W_rr = f0.                                            (6)
```

The initial data remain those of `B0` since the lift and `Kf0` vanish initially.
Equation (6) is precisely the remaining limitation: it does not replace the
frozen source by `A W_x+B W_r+D W+F_affine` from the actual coupled system.
The bulk residual against that registered equation is

```
f0 - A W_x - B W_r - D W - F_affine.                          (7)
```

No claim that (7) vanishes, or that all Einstein/constraint residuals vanish,
follows from enforcing the boundary law. Unresolved bulk, subsidiary and
second-derivative errors require their own bounds.

## Endpoint temporal moments without nested integration

Let `U_j=J^j u`, `V_j=J^j v` for the smoothed initial displacement and
velocity, with `U_0=u`, `V_0=v`. Their zero collar normalizes every primitive
at zero. The physical ridge functions in the previous construction are

```
F(z) = [(A0+Bmat0)(u'(z)+v(z)) + D0(u(z)+V_1(z))]/2,
G_ridge(z) = [(Bmat0-A0)(u'(z)-v(z)) + D0(u(z)-V_1(z))]/2.
```

Let `F_j,G_j` be their normalized j-fold primitives. Define

```
H_j(x) = x F_(j+1)(x) - (j+1) F_(j+2)(x)
       + [G_(j+2)(x) + (-1)^j G_(j+2)(-x)]/2,
h_j(x) = PN[U_j(x) + V_(j+1)(x) + H_j(x)].                   (8)
```

Then `H_j=J^j H_0`, where `PN H_0` is the reflected bulk trace, and
`h_j=J^j h`. The identity follows by repeated integration of `x F_1-F_2`:
`J^j(x F_1-F_2)=x F_(j+1)-(j+1)F_(j+2)`. Each integration of a negative
argument contributes one minus sign. All primitive constants vanish at zero.

If `P=PN-PD`, the negative ridge primitives for k>=1 can be evaluated as

```
G_k(-x) = (-1)^(k-1) / 2 * {
    (Bmat0-A0) P [U_(k-1)(x)+V_k(x)]
    - D0 P [U_k(x)+V_(k+1)(x)] } .                          (9)
```

The parity in (9) stays inside each matrix action. It is not moved through a
coupling. Equations (8)–(9) reduce a finite trace series to endpoint values,
positive weighted displacement integrals and exact smoothed velocity
primitives. The derivative uses (5), not differentiation of an interval
remainder. The odd/Dirichlet radial derivative remains `PD S_r`.

## Rational tails and the unresolved affine source

For a homogeneous trace series truncated at k=N and `a=m*d`, the tail satisfies

```
T_h <= H a^(N+1)/(N+1)! / [1-a/(N+2)]                       (10)
```

whenever `a<N+2`. This is a rational majorant for the factorial series: the
ratio of successive omitted terms is at most `a/(N+2)`. This inequality is a
convenience for the enclosure, not an existence restriction on (1).

The full affine source term in (3) is retained in the analytic candidate.
It need not have been evaluated pointwise to be enclosed. Its leading term
is `-J G_shell`, whose ith component has magnitude at most `d G_i` when
`G_i>=sup|(G_shell)_i|`. The remaining source terms obey

```
T_G,higher <= m d^2 Gcap/2 / [1-a/3]                         (11)
```

for `a<3`, because they sum from k=1 with denominator `(k+1)!`. Thus the
numerical method may add the symmetric row interval for `-JG_shell` and a
Gram-norm radius from (11). This uses a uniform inherited source enclosure;
it does not newly determine the source's point value, sign or correlations.
If an actual Robin matrix row vanishes identically, every term containing a
positive matrix power has that row zero, so both relevant tails are exactly
zero there. The affine source itself can remain nonzero in such a row.

## Correction of the earlier boundary mismatch and C1 transfer

Let `Z0=PN B0(x,0)` and let the earlier first-feedback approximation be
`V=S+L(M0 Z0)`, with `Lq=-Jq` along the shell. Its actual-law defect is

```
R = PN V_r - M PN V - G_shell
  = -M0 PN B1 - (M-M0)PN V - G_shell.                        (12)
```

Writing `C=W-V` and `Q=MZ+G_shell-M0 Z0` gives

```
C = LQ,
C_trace = J R - M J C_trace,
Q + M J Q = -R.                                             (13)
```

In particular, omitting the actual-minus-frozen matrix term in (12) would
misidentify the forcing being corrected. If a uniform trace defect bound is
`Rcap`, (13) gives `sup||Q|| <= exp(m*L) Rcap`. The simpler rational bound
`Rcap/(1-m*L)` is also valid when `m*L<1`. On the retained cone the active
lift length is at most `L=1/4`, so for the earlier norm
`X(C)=max(sup||C||/3,sup||C_x||,sup||C_r||)` one has

```
X(W-V) <= sup||Q||,
X(U-W) <= X(U-V) + sup||Q||.                                (14)
```

For example, with `m0>=||M0||`, `deltaM>=||M-M0||`, and existing uniform
trace caps,
`Rcap=m0*B1tracecap+deltaM*Vtracecap+Gcap` suffices. The companion transfer
checker supplies the model-specific rational arithmetic. Equations (13)–(14)
explain why correcting the boundary mismatch preserves a conditional global
C1 approximation estimate; they do not create a numerical C2 error bound.

## Independent algebra controls

`checks/check_robin_resolvent.py` imports no numerical HDBLAST evaluator.
It checks (8)–(9) through j=5 by direct exact polynomial primitives with
noncommuting bulk matrices, then manufactures an exact polynomial solution Z,
a nonzero affine G and h from (1). A bivariate characteristic lift verifies
the exact affine Robin identity and unchanged principal bulk source.
Controls reject use of M0 in place of M, omission of G and the wrong lift sign.
Finite-series residuals telescope exactly and are nonzero in the fixture,
demonstrating the distinction between the exact reconstruction and its finite
enclosure. Rational factorial/source bounds are checked against the independent
known polynomial solution for N=0 through5. All44 named groups pass identically
under normal and optimized Python. These are finite algebra controls supporting
the displayed proof, not external specialist or proof-assistant verification.

All conclusions retain the registered background, first-field, shell source,
coefficient and formal-response assumptions. No physical observable, exact
nonlinear spacetime or higher-dimensional origin of the Big Bang is established.
