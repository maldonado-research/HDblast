# Conditional observer-flow input budget for the actual-area moment

5 November 2026, UTC. Working note. This is a finite-dimensional
parameter-ODE estimate on an assumed common compact proper-clock tube.
It does not prove a new Einstein–scalar or shell PDE theorem, certify that
tube, evaluate any nonlinear constant, or produce a numerical amplitude radius.

## 1. Source interface and physical identification

The physical observer and moment definition is the immutable
[nonlinear bridge](https://github.com/maldonado-research/HDblast/blob/fe1d84d420eefcad60156bdfeff1c3cbde9f66a5/research/HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/nonlinear/NONLINEAR_BRIDGE.md).
The current full-pullback metric, area and normalized quotient interface is
the companion [area-density bound](../geometry/AREA_DENSITY_BOUND.md).
The normalized-moment estimates used below are ordinary quotient-rule
bookkeeping.

They fix the normal-geodesic observers, transported round-sphere labels
`Omega`, actual proper time `tau`, template `Y=n1^2-n2^2`, and actual-area
moment. In particular, the quotient budget needs full pulled-back scalar
bounds `E1,E2`, not scalar derivatives at fixed coordinate time. The retained
PDE result is conditional and its parameter interval is not evaluated.

Use a fixed reference shell manifold `S` of dimension four. Its complete
physical fields must already be supplied as

```
h_e = iota_e^* g_e,                 varphi_e = iota_e^* phi_e,
```

where `iota_e` is the actual physical shell embedding, with all parameter
dependence included. These are fields on the reference shell after a fixed
identification. A coordinate restriction to an unperturbed shell does not
supply them. The ODE below adds intrinsic observer transport to these full
fields; it does not reconstruct missing embedding derivatives.

Let `sigma_e(Omega)` be the specified initial section and `u_e(0,Omega)`
its future unit normal in `h_e`, with signature `(-,+,+,+)`. The initial
normal, its label derivatives, and its parameter derivatives are part of
the initial inputs. Setting them to their background values is not generally
valid, even if the section is represented by a fixed coordinate surface.

Fix `0<a<a_adm`. Assume the initial-value problem and the actual cuts are
admitted on an open neighborhood of

```
|e| <= a,   -T_obs <= tau <= 0,   Omega in S^3.
```

Assume its states stay in a compact tube inside the regular shell domain,
the cuts are spacelike and noncrossing, and the fixed label maps remain
valid. The entire interval from the anchor to the cuts is needed, not just
the short interval on which the moment zero is sought. These are premises,
not consequences of the derivative estimates below. A sufficient regularity
input is a `C^2` parameter curve of full shell fields in a uniform shell `C^4`
topology, a uniformly bounded inverse metric on this tube, and correspondingly
regular initial section and normal data. More minimally, require the
augmented ODE field described below to be jointly `C^2` in `(e,X)`, with
continuous clock dependence and the displayed derivative bounds.

`T_obs` is a proper duration. For the central interval in the bridge it may
be chosen as `rho_b*6145 ell/4096` only if the complete connecting tube up
to that proper duration has been admitted. It is distinct from the bridge's
coordinate-slab number `T_coord=13 ell/8`. At nonzero parameter one cannot
set the coordinate component `x^0(-s)=-s/rho_b` by definition.

## 2. Geodesic and label-Jacobi equations, with past orientation

In a fixed shell chart let `x=x_e(tau,Omega)` and retain the future-directed
proper velocity `u=partial_tau x`. Write the Christoffel symbols of `h_e`
as `Gamma_e`. For a fixed local label direction `partial_A`, set

```
P_A^mu = partial_A x^mu,            Q_A^mu = partial_A u^mu.
```

The proper-time system is

```
partial_tau x^mu = u^mu,
partial_tau u^mu = -Gamma^mu_ab u^a u^b,
partial_tau P_A^mu = Q_A^mu,
partial_tau Q_A^mu = -(partial_nu Gamma^mu_ab) P_A^nu u^a u^b
                     -Gamma^mu_ab(Q_A^a u^b+u^a Q_A^b).
```

The last line is the differentiated geodesic equation. Its coordinate
form is sufficient here; no curvature or extra dynamical equation is added.
For Levi-Civita symbols its last term equals `-2 Gamma^mu_ab u^a Q_A^b`.
The components `x,u,P,Q` total `4+4+12+12=32` in such a chart.

To integrate into the past, put `s=-tau >= 0` and keep the definition
`u=partial_tau x`, `Q=partial_A u`. The system `X_s=F(e,s,X)` is then

```
x_s = -u,
u_s = +Gamma[u,u],
P_s = -Q,
Q_s = +(partial_x Gamma)[P,u,u] +2 Gamma[u,Q].       (A)
```

Thus every component of the proper-time vector field changes sign.
The alternative convention `v=partial_s x=-u` gives the usual spray signs
but has past-directed initial velocity and reverses the corresponding Q;
the two conventions must not be mixed. The parameter derivatives below
hold `s` fixed, equivalently hold the actual `tau=-s` fixed. No inverse
coordinate-time derivative enters. Along admitted geodesics the metric
compatibility identity preserves `h_e(u,u)=-1`, so this is the proper clock.

Because the augmented field contains `partial_x Gamma`, its `F_XX` may
contain three spatial derivatives of `Gamma`, hence four of `h_e`.
Mixed bounds include, for example, `partial_e^2 partial_x Gamma`.
Merely assuming a twice spatially differentiable metric would not justify
the augmented estimates. The sufficient `C^2(e;C^4)` hypothesis above
accounts for these derivatives, without claiming bounds for them.

## 3. A precise way to use finitely many charts

All norms below are fixed, parameter-independent, scaled norms; operator
and bilinear norms are induced by them. `x,u,P,Q` have different physical
units. Fix positive reference scales before combining their components.
An unscaled sum of length, velocity and label derivatives is not a physical
norm. With this convention `L*T_obs` is dimensionless, and `V` and `W`
carry respectively one and two inverse parameter units. The parameter need
not have been declared dimensionless.

There are two equivalent bookkeeping choices:

1. Use finite shell/label charts and a finite subdivision in `s`, chosen
   on local parameter/label neighborhoods so every whole segment stays in
   its chart. The switch times are fixed, independent of `e`. Compactness
   allows a finite collection of such chains. Apply the bounds in section 4
   on each segment, and charge every transition as stated below.
2. Choose once an auxiliary smooth embedding of the compact shell tube
   into a finite-dimensional Euclidean space, together with the induced
   representation of its velocities and label derivatives. Use fixed
   smooth label frames (or a finite trivialization). A tubular extension of
   the lifted vector field gives a `C^2` ODE in a neighborhood of the compact
   embedded state set. Its derivative bounds include the embedding,
   trivialization and extension constants. Ordinary Euclidean estimates
   then apply over the entire duration. This auxiliary embedding is not the
   physical `iota_e`, is parameter independent, and asserts no global shell
   coordinate chart. Scalar fields can be extended using the same fixed
   retraction. Equivalent choices change the bounds, not the observable.

For choice 1, let `X_+=Psi(X_-)` be a parameter-independent augmented chart
transition at a fixed switch time. `Psi` includes the transformation of
velocities and of P,Q, as well as any fixed label-frame changes. Then

```
V_+ = D Psi V_-,
W_+ = D Psi W_- +D^2 Psi[V_-,V_-].
```

If `||D Psi||<=alpha`, `||D^2 Psi||<=beta`, the next initial budgets are

```
v_next <= alpha v_end,
w_next <= alpha w_end +beta v_end^2.                (B)
```

These transition charges are necessary. Taking the maximum chartwise L
and ignoring transition Hessians is not a valid global second-derivative
bound. A parameter-dependent exit time would add further terms; it is not
used here. For choice 2, all these effects are included in the single
ambient field F and its fixed norm. The simple duration-wide formulas below
refer to that coherent representation or to an individual chart segment.

## 4. Parameter variational system and its uniform budget

Set `V=partial_e X`, `W=partial_e^2 X`. For fixed `s,Omega`, ordinary
differentiation and symmetry of mixed derivatives give exactly

```
V_s = F_X V+F_e,
W_s = F_X W+F_XX[V,V]+2F_Xe V+F_ee.                (C)
```

Here `F_e` means explicit parameter dependence at fixed state, while V,W
include the motion of that state. The term `2F_Xe V` has coefficient two.
Assume throughout the compact state/parameter tube

```
||F_X|| <= L,       ||F_e|| <= P1,
||F_XX|| <= H,      ||F_Xe|| <= C,      ||F_ee|| <= P2,
||V(0)|| <= v0,     ||W(0)|| <= w0,                (D)
```

where every bound is nonnegative and uniform in the label. `v0,w0` refer to
the full augmented initial state, including normal velocity and all P,Q
data. They are not only initial position bounds.

For `t>=0` define

```
phi_L(t) = (exp(Lt)-1)/L       if L>0,
phi_0(t) = t.                                   (E)
```

The variation-of-constants/Gronwall inequality from the first line of (C)
gives, on `[0,T]`,

```
||V(s)|| <= exp(Ls)v0 +phi_L(s)P1,
V_* = exp(LT)v0 +phi_L(T)P1.                     (F)
```

All coefficients are nonnegative, so the last expression bounds every
earlier s. The second line of (C) then gives

```
||W(s)|| <= exp(Ls)w0
            +phi_L(s)[H V_*^2+2C V_*+P2],
W_* = exp(LT)w0 +phi_L(T)[H V_*^2+2C V_*+P2].    (G)
```

For a finite chart chain use `T` equal to each segment duration, propagate
its endpoint budgets by (B), and take maxima over all segment bounds and
all finitely many chains. Bounds on an earlier segment may exceed the final
endpoint after a change of norm, so include all segments in that maximum.
For one ambient representation use `T=T_obs` directly. A clock-dependent
F is allowed; no `F_s` bound is needed because parameter differentiation is
at fixed clock.

For completeness, (F) follows from
`||V(s)||<=v0+integral_0^s(L||V(r)||+P1)dr`.
The forcing bound in (G) follows from (C), (D), and (F).
The identity `phi_L(s)=integral_0^s exp(L(s-r))dr` proves both cases of (E)
and avoids division by zero. When `T=0`, `(V_*,W_*)=(v0,w0)`; when `L=0`,
`V_*=v0+TP1` and `W_*=w0+T[H V_*^2+2C V_*+P2]`.
If `v0=P1=0`, V vanishes; its nonzero second variation can still be driven
by `w0` or `P2`. Neither an exponential numerical value nor an exponential
interval certificate is needed until actual L,T inputs are supplied.

## 5. Output to the full scalar and cut-metric budgets

Write the position components of the parameter derivatives as

```
Z = partial_e x,      W_x = partial_e^2 x.
```

Commutation with the fixed label differentiation gives

```
partial_A Z   = partial_e P_A,
partial_A W_x = partial_e^2 P_A.                  (H)
```

Thus the x and P components of V,W control all four quantities needed by
the induced cut-metric calculation. If fixed projection constants satisfy

```
||Z|| <= k_x V_*,      ||partial_Omega Z|| <= k_P V_*,
||W_x|| <= k_x W_*,    ||partial_Omega W_x|| <= k_P W_*,
```

use these values as inputs to that calculation. The constants include
reference scales and chart/frame conversion. They equal one only for a
deliberately compatible scaled product norm. In an ambient realization
these are ambient derivatives; converting to a shell chart uses its first
and second derivative bounds, exactly as in (B). Do not silently identify
ambient second derivatives with chart second derivatives. Likewise W_x
alone would not supply `partial_Omega W_x`; that is why P,Q were included.

In a shell chart, the complete observer scalar is

```
A_e(tau,Omega) = varphi_e(x_e(tau,Omega))-phi_b.
```

For the constant background `phi_b`, its parameter derivatives are exactly

```
partial_e A = varphi_e,par +D_x varphi Z,
partial_e^2 A = varphi_e,par,par +2D_x varphi_e,par Z
                 +D_x^2 varphi[Z,Z]+D_x varphi W_x.       (I)
```

The notation `par` denotes explicit parameter differentiation at fixed
shell point, to distinguish it from the subscript indexing the family.
These are derivatives of the complete shell scalar, not of an unembedded
bulk-coordinate restriction. In an ambient realization use its fixed
scalar extension and ambient Z,W_x; the same identity holds and gives the
same physical A. Individual bounds depend on that representation.

Suppose common compact-tube bounds in the chosen representation are

```
|varphi_par| <= b_e,        ||D_x varphi|| <= b_x,
|varphi_par,par| <= b_ee,   ||D_x varphi_par|| <= b_ex,
||D_x^2 varphi|| <= b_xx,
||Z|| <= z,                ||W_x|| <= w.
```

Then the quotient note may use

```
E1 <= b_e+b_x z,
E2 <= b_ee+2b_ex z+b_xx z^2+b_x w.               (J)
```

Choose z,w from (F)–(H), or their finite-chart propagation. Dimensions match:
each summand of E1 has scalar/parameter units and each summand of E2 has
scalar/parameter-squared units. The term `D_x varphi W_x` and the mixed
term in (I) cannot be omitted. The full position includes its coordinate
time component, so clock motion is charged as well as angular motion.

For clarity the corresponding induced cut metric is
`q_AB=h_e(x_e)(P_A,P_B)`. Its differentiated bounds require h and its first
two parameter/space derivatives, P, Z,W_x, and the label derivatives in
(H). The separate area calculation converts these to J1,J2 and positivity;
the ODE estimates alone do not certify a determinant lower bound.

Combining complete scalar inputs (J), complete area inputs J1,J2, the
background anchors, and the same common proper-clock domain with the
quotient note yields, whenever `a J1<=1/2`,

```
K(a) <= H_mom(a),
H_mom(a)=6E2+32E1 J1+16a E1 J2+64a E1 J1^2.
```

If also `a H_mom(a)<106` and the retained linear endpoint signs apply on
that same proper interval, there is an actual-area moment zero for each
admitted `0<|e|<=a`. This is only a conditional input chain. No value of
L,P1,H,C,P2, v0,w0, scalar/metric bounds, J1,J2, a_adm or an effective a is
evaluated here. Oddness does not eliminate the terms in (I), (C), or the
area bounds away from zero. Joint diffeomorphism transport of shell fields,
initial section and labels preserves the observable; atlas, norm and
auxiliary-embedding choices only change these sufficient constants.

## 6. Reproducibility and limits of the checker

`check_budget.py` uses only Python's standard library and exact Fractions.
It checks (C) coefficient by coefficient for a basis of parameter/state
monomials through degree four, the scalar chain identity, past-orientation
signs, chart-transition second-derivative terms, and zero-L/zero-time edge
cases. It also checks the algebraic positivity and dimensionless rewrite
of the phi function by its rational power-series coefficients. This is a
finite algebra check, not a proof of the PDE premises, existence of the
tube, analytic Gronwall lemma, or any numerical physical estimate.

From this directory, choose a fresh report directory outside the closed
capsule, then run

```
report_dir="$(mktemp -d)"
python3 check_budget.py --output "$report_dir/normal.json"
python3 -O check_budget.py --output "$report_dir/optimized.json"
cmp "$report_dir/normal.json" "$report_dir/optimized.json"
```

The required `--output` argument must name a fresh absolute non-symlink
file whose parent directory exists. Exclusive creation forbids overwrites.
Reports may be outside the capsule so replay preserves its manifest.
Both modes emit the same canonical report bytes. `CHECK_REPORT.md` records the result;
`MANIFEST.json` hashes the deliverables and immutable bridge source.
