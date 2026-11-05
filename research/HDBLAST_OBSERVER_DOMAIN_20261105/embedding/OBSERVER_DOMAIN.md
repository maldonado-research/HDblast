# Global embedding and spacelike cuts for shell observers

5 November 2026, UTC. Conditional finite-dimensional geometry.

## 1. Physical identification and premises

Fix `rho_b>0`, `0<a<a_adm`. On a fixed reference shell product region
`(-T_coord,T_coord) x S^3`, with an open margin around the image, supply
the complete physical fields `h_e=iota_e^*g_e`, `varphi_e=iota_e^*phi_e`.
The physical shell map `iota_e` is an embedding on this region, with all
parameter dependence included. The identification with the reference shell
is fixed before parameter differentiation. Represent the specified anchor
section in this gauge by `sigma(n)=(0,n)` with the specified labels; its
future unit normal is the actual `h_e` unit normal `u_e(0,n)`. A fixed
coordinate representation of the section does not fix its normal or its
physical image in `e`. If the original section varies in this gauge, its
derivatives belong in the initial data and a change to this fixed-section
gauge needs justification. It cannot silently be reset to the background.

Assume a sufficiently regular intrinsic normal-geodesic observer solution
exists jointly on an open neighborhood of

```
|e|<=a,  q=-tau/rho_b in [0,T_q],  n in S^3,
X_e(q,n)=(t_e(q,n),omega_e(q,n)),   |t_e|<=T_coord.
```

This full-cylinder flow-existence and tube-retention premise does **not**
assume noncrossing or spacelike cuts. One can solve the parametrized initial
value equations before knowing whether different clocks or labels collide.
All derivative bounds below hold on this whole cylinder and parameter
interval, including each segment from zero to `e`. Require continuous first
domain derivatives and continuous parameter derivatives of the velocity and
label Jacobi operator, so the displayed commutations and integrations hold.
A separate existence/bootstrap proof must supply these inputs without
assuming the conclusion. Applying the retained observer budget only on an
already noncrossing tube would not remove its domain assumption.

The exact background identification is

```
h_0=rho_b^2[-dt^2+cosh^2(t) kappa_round],
X_0(q,n)=(-q,n),   u_0=partial_tau X_0=(1/rho_b,0).
```

Indeed the background future unit normal is `rho_b^(-1)partial_t`; its
normal geodesics keep the angular labels fixed and have `t=tau/rho_b=-q`.
This is not a definition of coordinate time at nonzero parameter. Parameter
derivatives hold actual proper time and transported labels fixed.
`T_q=T_obs/rho_b` and `T_coord` are distinct; the background containment
premise entails that the coordinate tube contains `[-T_q,0]`.

Use the standard global embedding `S^3 subset R^4` for angular positions.
Derivative differences are ambient differences in this fixed embedding.
The auxiliary target metric is `g_+=dt^2+kappa_round`; the domain metric
is `dq^2+kappa_round`. Finite local charts may obtain the norms, charging
every conversion constant. No global chart or orthonormal sphere frame
is asserted. The full shell fields include the physical embedding already;
the intrinsic flow does not reconstruct missing embedding derivatives.

## 2. A global chord criterion

Set `M=[0,T_q] x S^3`, `Psi_e=(-t_e,omega_e)`,
`Psi_0(q,n)=(q,n)`, and `R_e=Psi_e-Psi_0`, all in `R x R^4`.
If, uniformly,

```
||D_(q,n)R_e||op <= kappa < 2/pi,                    (1)
```

then for every pair `x=(q,n),y=(r,m)` in `M`,

```
|Psi_e(x)-Psi_e(y)| >= (1-pi*kappa/2)|x-y|.           (2)
```

Proof: let `theta in [0,pi]` be the round angle between `n,m`. A minimizing
great-circle arc (any semicircle for antipodes) combined with linear clock
interpolation is a minimizing product geodesic lying inside `M`, of length
`d_M=sqrt((q-r)^2+theta^2)`. The background chord has squared length
`c^2=(q-r)^2+4sin^2(theta/2)`. Concavity of sine on `[0,pi/2]` gives
`sin z>=2z/pi`, hence `theta<=pi sin(theta/2)` and `d_M<=(pi/2)c`
(also using `pi/2>1`). Integrate `DR_e` along this geodesic on the cylinder:
`|R_e(x)-R_e(y)|<=kappa d_M`. The reverse triangle inequality proves (2).
The ambient straight chord need not lie on `M` and is never used as an
integration path. All pairs with different clocks are included.

The coefficient is positive, proving global injectivity and a quantitative
inverse chord bound. Since the background inclusion is isometric on tangent
vectors, `|D Psi_e v|>=(1-kappa)|v|>0` for `v!=0` (`2/pi<1`). Thus
`Psi_e` is an immersion, including the boundary in the manifold-with-boundary
sense. A continuous injection of compact `M` into a Hausdorff space is a
homeomorphism onto its image, proving a `C^1` embedding. Time reflection is
a diffeomorphism, so `X_e` is embedded as well. Composing with the assumed
physical shell embedding `iota_e` gives the physical embedded observer
domain. Distinct clocks or labels cannot cross in this domain. This does
not claim surjectivity onto a chosen coordinate tube. Higher regularity is
inherited from the flow. If `T_q=0` the argument embeds the initial sphere;
a four-dimensional domain then has no duration.

## 3. Projected augmented first-variation inputs

In this fixed ambient representation let `u_e=partial_tau X_e` be the full
proper velocity and `P_e=D_nX_e:T_nS^3 -> R x R^4` the label Jacobi map.
Suppose uniform first parameter bounds give

```
||partial_e u_e|| <= U1,    ||partial_e P_e||op <= P1. (3)
```

Here `P1` is a projected Jacobi-variation budget, not the old observer
note's unrelated `P1` for `||F_e||`. Certified projections of its augmented
first variation may supply these bounds only after charging scales,
embedding/trivialization and chart constants, and all initial normal and
section data. Bounds at parameter zero alone are insufficient. `rho_b U1`
and `P1` have the same inverse-parameter units.

For the isometry `J=diag(-1,I_4)`, `Psi_e=J X_e` and
`partial_qX_e=-rho_b u_e`. Parameter integration at fixed `(q,n)` yields

```
|partial_q R_e| <= |e| rho_b U1,
||D_nR_e||op <= |e| P1,
b := sup ||D_nX_e-D_nX_0||op <= a P1.               (4)
```

For `(alpha,v)`, triangle inequality and Cauchy-Schwarz give
`|DR_e(alpha,v)|<=|e|[rho_b U1|alpha|+P1|v|]`
`<=|e|sqrt((rho_b U1)^2+P1^2)sqrt(alpha^2+|v|^2)`. Therefore one may take

```
kappa=a sqrt((rho_b U1)^2+P1^2).                    (5)
```

A rational sufficient guard for nonnegative rational inputs is

```
C=a(rho_b U1+P1)<1/2,                              (6)
```

since `kappa<=C` and `pi<4` gives `1/2<2/pi`. Then the lower coefficient
in (2) may also be replaced by the positive rational `1-2C`, with a
non-strict inequality valid even for `C=0`. These are sufficient guards;
neither threshold is claimed sharp. The `pi/2` geodesic-to-chord comparison
is attained for equal clocks and antipodal labels, which by itself proves
no sharpness of an injectivity threshold.

## 4. Spacelike cuts from a full metric gap

Suppose throughout the image tube

```
||h_e-h_0||_(g_+) <= rho_b^2 eta,   eta>=0.          (7)
```

This is the symmetric-bilinear operator norm relative to the positive
`g_+`; its quadratic error at `w` is at most `rho_b^2 eta |w|_+^2`.
A bulk coordinate metric restricted to an unperturbed shell is insufficient
to supply this gap for the complete physical shell metric.

For a unit round label tangent `v`, write `D_nX_e(v)=(s,A)`. When `b<=1`,
(4) gives `|s|<=b`, `|A-v|<=b`, and `1-b<=|A|<=1+b`. The vectors `A,v`
may be tangent at different angular points; the ambient comparison and
their round lengths are well defined. For `kcut_e=X_e(q,.)^*h_e`,

```
h_0(D_nX_e v,D_nX_e v)/rho_b^2
 = -s^2+cosh^2(t_e)|A|^2 >= (1-b)^2-b^2=1-2b,
|D_nX_e v|_+^2 <= b^2+(1+b)^2=1+2b+2b^2.
```

Using `|t_e|<=T_coord` and discarding `-s^2` for an upper bound proves,
as quadratic forms relative to the round metric,

```
Lcut I <= kcut_e/rho_b^2 <= Ucut I,
Lcut=1-2b-eta(1+2b+2b^2),
Ucut=cosh^2(T_coord)(1+b)^2+eta[(1+b)^2+b^2].        (8)
```

Thus `Lcut>0` proves spacelikeness without presuming it. This local
positivity does not itself prove noncrossing; section 2 supplies that
separate global conclusion. Also, `eta<1` preserves the shell metric's
Lorentzian signature: relative to `g_+`, `h_0/rho_b^2` has eigenvalues
`-1,cosh^2t,cosh^2t,cosh^2t`, whose distances from zero are at least one.
A symmetric perturbation of norm less than one cannot change the inertia.

## 5. Exact compatible fixed-background spectrum

The background cut is `k_0(q)=rho_b^2 cosh^2(q) kappa_round`. Define
`G_e=k_0^(-1/2) kcut_e k_0^(-1/2)` in background-orthonormal frames, or as
the corresponding endomorphism; no global frame is needed. Suppose

```
0<=b,eta<=1/8,    0<=T_q,T_coord<=1/4.              (9)
```

These are symbolic sufficient thresholds, not evaluated physical budgets.
`Lcut` decreases in nonnegative `b,eta`, so `Lcut>=151/256>0`.
For `0<=z<=1/4`, the continuous series argument

```
cosh z = sum_(j>=0) z^(2j)/(2j)!
       <= sum_(j>=0) (z^2/2)^j
       = 1/(1-z^2/2) <= 32/31
```

uses `(2j)!>=2^j`; hence `cosh^2z<=1024/961`. Lower normalization divides
by at most `cosh^2(T_q)` and upper normalization by at least one. Therefore

```
G_e >= (151/256)(961/1024) I
     = (145111/262144)I > (1/2)I,
G_e <= [(1024/961)(9/8)^2+(1/8)((9/8)^2+(1/8)^2)]I
     = [1296/961+41/256]I
     = (371177/246016)I < 2I.                      (10)
```

The strict margins persist even at the endpoints in (9). The retained
area-density interface can thus use `lambda=1/2,Lambda=2` on this same
admitted domain. No `G'_e,G''_e`, scalar, or remainder inputs are supplied
by these spectral bounds; they require their own full-pullback estimates.

## 6. Scope and relation to retained notes

Given independently justified full-cylinder flow existence, fixed physical
identification, uniform projected velocity/Jacobi bounds, coordinate
containment, and the metric gap, the sufficient guards prove global
noncrossing and spacelike cuts. These two conditions can then be removed
as independent assumptions. Existence and retention are not proved here.
This note supplies no nonlinear field theorem, evaluated constant,
effective numerical amplitude, or moment-zero result by itself.

The source README states the observable and full-pullback scope. The source
observer note sections 1--5 provide the augmented flow interface and its
norm and domain caveats. The area note sections 1 and 3 require the full
cut metric and all pullback terms. The public source-dependency note
separates retained analytic variation premises from missing numerical and
observer-domain inputs. Copies and byte pins identify precisely those sources.
The new argument uses standard
product geometry, elementary calculus and compactness; it is an explicit
conditional interface, not a claim of a new general embedding theorem.
