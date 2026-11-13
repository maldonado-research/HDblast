# Physical shell traces and the complete proper second scalar at zero

13 November 2026 UTC. Conditional derivation with internal mathematical review. All amplitude jets are
derivatives at zero, with expansion `F_e=F_0+e F_1+e^2 F_2/2+o(e^2)`.
This note is conditional on the retained smooth background, selected exact
constraint curve, relative-harmonic shell development and uniqueness
premises. Applying its quantitative bounds on the **past** observer also
requires the separate complete same-state reverse-time join. The older
future-only evolution statements are not altered here.

## 1. Sources, fixed norms and the physical identification

Source paths in this note are relative to `research/` of the unmodified
base tree `7fb54e385e5250deae408910da55bf000a53f5a2`:

* `HDBLAST_SECOND_WAVE_EVOLUTION_20261112/MIXED_SOURCE.md`, section 2,
  supplies the full mixed first-wave norm, rather than only N14.
* `HDBLAST_SECOND_WAVE_EVOLUTION_20261112/SECOND_WAVE.md`, sections 1 and
  8, supplies the actual mixed second-wave X13 norm and its identification
  by lower-energy uniqueness.
* `HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md`, section 4,
  proves the boundary-inclusive cylinder H3 embedding with constant 10^7.
* `HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md`, “Exact local
  reconstruction,” fixes the literal reconstructed metric and scalar.
* `HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/REGISTERED_EXACT_FAMILY_BINDING.md`,
  sections 3–4 and 6, fixes the actual common physical shell chart, selected
  canonical curve, initial open zero collar, harmonic completion and lower
  C2 parameter topology.
  `HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/inputs/NONLINEAR_BRIDGE.md`, sections 2–5,
  fixes the original anchor, proper clock, transported labels and moment.
* `HDBLAST_INITIAL_TRANSPORT_20261106/INITIAL_TRANSPORT.md` specifies the
  fixed positive bundle norm and the exact initial boundary-chart transport.
* `HDBLAST_FIRST_PROPER_GEOMETRY_20261106/FIRST_CUT_BOUND.md`, section 2,
  records the difference between harmonic and longitudinal observer
  components; its longitudinal coordinate bounds are not used below.

Set

```
ell=10^-6,  T=13ell/8=26ell/16,  Q=25ell/16,
B1=10^4000000 exp(10^152017),
B2=10^10200000 exp(10^152019).
```

The bulk spatial norm uses `dz^2+kappa`, full measure `dz dmu_kappa`,
product connection and full derivative square sums. The fixed bundle is
`R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D`; the symmetric
block has its full Frobenius norm, including both off-diagonal entries.
All tensor derivative slots are included. The actual first and second
amplitude fields are `u=V1` and `v=V2`; their mixed norm bounds are

```
sum_(j=0)^14 ||partial_t^j u||H^(14-j) < B1,
sum_(j=0)^13 ||partial_t^j v||H^(13-j) < B2.               (1)
```

These follow on the already proved future slab. Whenever (1) is used at
`t=-q`, the separate reverse-time theorem, including the tensor time-slot
signs and the actual nonzero initial harmonic mixed velocity, is an
explicit prerequisite. Relabelling the actual future solution is not that
theorem.

For actual parameter-dependent observers on the whole interval below,
retain the qualitative exact-family prescribed-slab application on
`-T<=t<=0` as well. The original binding's rectangle ends at `3ell/2`;
it is not silently enlarged. NONLINEAR_BRIDGE section 2 gives the longer
prescribed-slab application with strict causal margins, conditional on
the same analytic state-flow, inverse, support, regularity and existence
premises. The reverse-time join and that conditional application are
both part of the hypothesis here; zero-amplitude wave estimates alone
would not construct a nonlinear observer neighborhood beyond `3ell/2`.

On the reference shell let `p=dt^2+kappa` be the positive product metric
and `D` its product connection. Use p-Hilbert norms for every tensor and
derivative slot. Normalize the actual Lorentz shell metric by

```
g_e=h_e/rho_b^2,
g_0=-dt^2+cosh(t)^2 kappa,
p_a=partial_e^a g_e|0,   a=1,2.
```

The symbol `p_a` is a metric jet; the undecorated `p` is the positive norm
metric. The background radius is fixed, `129<=rho_b<=131`.

The retained physical shell embedding is `z=0` in the common chart.
The literal reconstruction writes `F_z(V)=(v,Htilde,w,Dcross)` and, in
the coframe `(dz,dt,cosh(t)theta^A)`,

```
G_ab=rho(z)^2 J4_ab+rho_b^2 Htilde_ab,
phi=phi0(z)+w.
```

Here `F_0=identity` for all states. Restricting the actual physical metric
and scalar to this shell, and then differentiating in amplitude, gives
the exact identities

```
(h_a)_ab=rho_b^2 (V_a)_H,ab|_(z=0)
                    in the displayed de Sitter coframe,
varphi_a=(V_a)_w|_(z=0),                    a=1,2.         (2)
```

This is the full physical pullback in the specified fixed-shell gauge.
It uses both the physical identification and the exact boundary flow; it
is not a replacement of an arbitrary moving shell by a coordinate surface.
The tangential reconstruction and scalar sector in (2) are linear, so
there are no missing quadratic reconstruction terms. The radial root and
cross block do not contribute to an intrinsic shell restriction. The
background **shell** scalar is `varphi_0=phi_b`, constant; the bulk scalar
`phi0(z)` is not constant.

## 2. Evaluated product-frame C4 trace budgets

Apply the proven H3 embedding to
`partial_t^j D_sp^k V_a`, with the same tensor-slot norms as in (1).
It gives, including both cylinder ends,

```
|partial_t^j D_sp^k u| <= M1:=10^7 B1,   j+k<=11,
|partial_t^j D_sp^k v| <= M2:=10^7 B2,   j+k<=10.         (3)
```

For example the first line needs spatial norm of order `k+3`, which is
available when `j+k+3<=14`. No extension across the shell and no additional
trace constant are used. The original tensor rank plus the field-derivative
slots used here remains within the source embedding's rank-sixteen budget.

Since the product time direction is parallel, for each coefficient block
`U_a` occurring in (2) the full spacetime product derivative satisfies

```
|D^m U_a|^2
 =sum_(j=0)^m binom(m,j)|partial_t^j D_S^(m-j) U_a|^2
 <=2^m M_a^2,                   m<=4.                    (4)
```

This expansion counts every placement of time versus angular derivative
indices; it does not commute two angular derivatives.

Put `A(t)=diag(1,cosh(t) Id_3)`. Converting the two covariant tangential
metric slots from the displayed coframe to a fixed p-orthonormal frame
is `p_a=P(t)(V_a)_H`, where `P=A tensor A` on the symmetric block.
The invariant orthogonal subblocks of P have eigenvalues `1,cosh(t),
cosh(t)^2`. On `|t|<=T<1/100`, for `0<=i<=4`,

```
||partial_t^i P||op <=2^(i+1).                            (5)
```

Indeed `cosh(t)^2=(1+cosh(2t))/2`, and every positive derivative of
its second term is at most `2^(i-1) cosh(2|t|)`. The elementary series
bound `exp(x)<=1/(1-x)` for `0<=x<1` gives `cosh(2|t|)<2` here;
the same bound covers derivatives of cosh(t) and i=0. Angular derivatives
of this multiplier vanish because the time/spatial splitting and the
identity on T*S3 are product-parallel.

The full covariant Leibniz rule, (4) and (5) therefore give

```
|D^m p_a|
 <=sum_(i=0)^m binom(m,i) 2^(i+1) 2^((m-i)/2) M_a
 =2(2+sqrt(2))^m M_a
 <512 M_a,                       0<=m<=4.                 (6)
```

The last inequality follows from `2+sqrt(2)<4` and `m<=4`. Define the
explicit common metric budgets

```
H_a=512 M_a=512*10^7 B_a.
max_(m<=4)|D^m p_a| <=H_a,
max_(m<=4)|D^m h_a| <=rho_b^2 H_a
 <=87864320000000 B_a <10^14 B_a.                         (7)
```

The physical scale in the last line uses `131^2*512*10^7` exactly.
The scalar requires no coframe multiplier:

```
max_(m<=4)|D^m varphi_a| <=4 M_a,
|varphi_a|<=M_a,      |d varphi_a|<=sqrt(2) M_a.           (8)
```

All suprema in (3)–(8) are on the common shell cylinder; their use on its
past side remains conditional on the same-state join. These are separate
bounds on identified derivatives at zero. The retained nonlinear curve
is still only C2 into X10; it is not asserted to be C2 into X13 or X14.

## 3. First harmonic observer, actual initial data and clock

Let `X_e(q,n)` be the original intrinsic shell geodesic, with its original
physical initial anchor and transported round label n. The parameter
`q=-tau/rho_b` is held fixed in amplitude differentiation. Its velocity
`V_e=partial_q X_e` satisfies `g_e(V_e,V_e)=-1` and is past directed.
On the background,

```
X0(q,n)=(-q,n),      V0=-partial_t,      0<=q<=Q.
```

The exact selected initial canonical data agree with the background on
the entire open shell collar `z<ell/2`. The prescribed lapse, zero shift,
and its relative-harmonic completion preserve this initial collar. Thus
the full initial shell metric, anchor and its past unit normal agree for
every parameter, not just for the first seed. Consequently

```
xi1=partial_e X_e|0,
xi1(0,n)=0,             D_q xi1(0,n)=0.                   (9)
```

Let `nabla=LC(g0)`, a fixed connection used for covariant comparison,
and set

```
eta1=nabla_e V_e|0=nabla_q xi1.
```

The second equality uses its torsion-free property. This is a covariant
velocity variation. It is not merely the derivative of coordinate velocity.

Write `A0=nabla-D`. In the p-product splitting, for angular U,W,

```
A0(partial_t,U)=tanh(t) U,
A0(U,W)=cosh(t)sinh(t) kappa(U,W) partial_t,
A0(partial_t,partial_t)=0.
```

Hence along `t=-q`, using `cosh(t)|sinh(t)|<=2|t|`,

```
|A0|_p^2=6tanh(t)^2+3cosh(t)^2 sinh(t)^2
        <=18q^2 <25q^2.
```

For clarity, `cosh(t)<=100/99` and `|sinh(t)|<=100|t|/99` follow
from the same elementary exponential series bound on `|t|<1/100`, so
their product is less than `2|t|`. The first connection variation is
the tensor

```
C1^a_bc=1/2 g0^ad [nabla_b p1_cd+nabla_c p1_bd
                                      -nabla_d p1_bc].    (10)
```

The inverse g0 has p-operator norm at most one. The two covariant metric
slots give

```
|nabla p1|<=|D p1|+2|A0||p1|<=(1+10q)H1,
|C1|<=(3/2)(1+10Q)H1<2H1,                               (11)
```

where `30Q<1` proves the last strict inequality. This includes the
background connection terms omitted by taking raw metric partials.

Use a product-parallel angular frame at the fixed label n and split
`xi1=a partial_t+b`, b tangent to S3. Direct differentiation of the full
harmonic-coordinate shell geodesic equation yields

```
a''(q)=-C1^t_tt(-q,n),
b''(q)+2tanh(q)b'(q)=-C1^S_tt(-q,n),
a(0)=a'(0)=b(0)=b'(0)=0.                                (12)
```

These are also the components of the forced Jacobi equation
`nabla_q eta1+R0(xi1,V0)V0=-C1(V0,V0)`, with the convention in which
the unit de Sitter background has curvature +1. Formula (12) can be
checked without a curvature sign convention: the background angular
connection term is `2 tanh(-q)(-1)b'=2tanh(q)b'`; the time equation
has no linear background spatial-velocity term.

Let `F(q)=C1(V0,V0)(-q,n)`. Solving (12) gives a single product-frame
integral formula

```
xi1'(q)=-int_0^q K(q,s)F(s) ds,
K(q,s)=diag(1,cosh(s)^2/cosh(q)^2 Id_3).
```

Since `0<=s<=q`, the p-operator norm of K is at most one. Therefore,
uniformly in n,

```
|xi1'(q)|<=2qH1,
|xi1(q)|<=q^2 H1,
eta1=xi1'+tanh(q)b,
|eta1(q)|<=2q(1+q^2/2)H1<3qH1.                           (13)
```

The last inequality uses `Q^2<1`. In particular the endpoint budgets
are evaluated exactly as

```
Cxi=Q^2 H1=B1/80,
Ceta=3Q H1=24000 B1,
sup|xi1'|<=2QH1=16000 B1.                                (14)
```

Only xi1 and eta1 in the **same harmonic representative as V2** occur
here. The sharper old longitudinal displacement does not replace (13).
Under the retained first gauge transformation the observer components
shift by `t1 -> t1-b_gauge Y` and `n1 -> n1-d_gauge grad Y`; transforming
the second scalar at the same time would be mandatory if that gauge were
used instead.

The background coordinate margin is

```
T-Q=ell/16=1/16000000>0.                                 (15)
```

The retained C2 lower-topology family and ordinary geodesic dependence,
on the compact background path, permit a qualitative parameter shrinkage
so these actual observers remain in the coordinate cylinder for all
`q<=Q`. Equation (15) is the margin for that application. No numerical
lower bound for the parameter radius follows from (13) or from compactness.

## 4. Complete proper scalar second derivative and evaluated bound

Define the actual observed scalar by
`Phi_e^obs(q,n)=varphi_e(X_e(q,n))`. The derivatives hold actual proper
time and the original transported labels fixed. The full chain rule gives

```
Phi1^obs=varphi1(X0),
Phi2^obs=varphi2(X0)+2 dvarphi1(X0)[xi1].                 (16)
```

The two potential background terms, `D^2 varphi0[xi1,xi1]` and
`dvarphi0[xi2]`, vanish because the full shell scalar `varphi0=phi_b`
is constant. No second observer variation is needed for (16). The same
formula would omit physical embedding terms if applied to a bare bulk
scalar before the identification in (2).

The scalar and observer budgets yield

```
sup_(q<=Q,n)|Phi2^obs|
 <=M2+2sqrt(2) M1 Q^2 H1
 <M2+3 M1 Q^2 H1
 =10^7 B2+375000 B1^2.                                  (17)
```

The last coefficient is an exact evaluation:
`3*(25/(16*10^6))^2*512*(10^7)^2=375000`.
Since `B1^2<B2`, a single explicit numerical envelope is

```
sup|Phi2^obs| <10^8 B2
              =10^10200008 exp(10^152019).               (18)
```

Indeed `8000000<10200000`, `2*10^152017<10^152019`, and
`10^7+375000<10^8`. Formula (17) retains the useful exact dependence
on the two separate bulk inputs; (18) is only its conservative rounding.

## 5. Selected normalized moment and scope

For the retained actual-area moment put
`A_e=Phi_e^obs-phi_b`, `dA_e=j_e dA0`,
`N_e=int Y A_e dA_e`, and `D_e=int Y^2 dA_e`.
Here `A0=0` and the retained exact first scalar is `A1=S(-q)Y`.
At zero,

```
N1=S(-q)D0,
N2=int Y Phi2^obs dA0+2S(-q)D1,
M''0=N2/D0-2N1D1/D0^2
     =int Y Phi2^obs dA0/D0.                             (19)
```

Thus neither j1 nor j2 is needed for this selected second moment.
The background denominator is
`D0=rho_b^3 cosh(q)^3 pi^2/3>0`, and the fixed sphere integrals give
`|M''0|<=(6/pi) sup|Phi2^obs|`. Subject to the actual selected-family
quarter-turn equivariance and uniqueness premises retained in
NONLINEAR_BRIDGE section 5, its stronger exact identity is `M''0=0`.
That conditional oddness is already retained; it is not a new consequence
of (17). For an arbitrary different first mode, the first area/denominator
cross terms in (19) would not cancel this way.

The new evaluated content here is the physical fixed-shell trace ledger,
the first harmonic proper-observer bounds, and the complete second scalar
field bound at zero. It does not provide a numerical nonzero amplitude
interval, a uniform nonlinear second derivative, a finite-amplitude
remainder, a third amplitude derivative, an observational confirmation or
a new physical-origin claim. The quantitative past conclusion remains
conditional on the separate same-state reverse-time theorem and all
retained analytic existence and identification premises.
