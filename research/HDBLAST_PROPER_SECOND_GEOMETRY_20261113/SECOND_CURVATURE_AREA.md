# Intrinsic second curvature and proper area at amplitude zero

13 November 2026 UTC. Conditional mathematical continuation. Every subscript 1 or 2
below is the corresponding actual derivative in the registered affine
amplitude at zero, not its Taylor coefficient. In particular a second
derivative is twice the quadratic Taylor coefficient. This note closes a
new curvature/area interface at that single expansion point; it supplies
no numerical nonzero amplitude interval.

## 1. Physical identification and conditional inputs

Retain the exact selected Einstein–scalar family, original physical
anchor, transported round labels, harmonic fixed-shell identification,
smooth background, compatible data, development/uniqueness, and common
observation-domain premises of the retained sources. In particular the
physical shell restriction is taken before any observer composition.
The [same-state negative-time continuation](REVERSE_TIME.md) is a separate
necessary join: the future-only wave inequalities cannot themselves be
applied at a negative coordinate time. The companion proof retains the
actual nonzero registered velocity and all harmonic tensor time slots;
its reflected positive-time solution is pulled back to the original
negative timeline.

Write

```
ell=10^-6, T=13ell/8, Q=25ell/16, q_D=ell/4,
q=-tau/rho_b, X0(q,n)=(-q,n), V0=partial_q X0=-partial_t,
g_e=h_e/rho_b^2, g0=-dt^2+cosh(t)^2 kappa,
p_a=partial_e^a g_e|0, a=1,2.
```

Here h_e is the full physical shell metric in the actual fixed-shell
gauge, and rho_b is a fixed background constant. Constant rescaling does
not change the connection or Ricci as a covariant two-tensor. Therefore

```
E_e=Ric(h_e)-3h_e/rho_b^2=Ric(g_e)-3g_e, E0 identically 0.
```

Use the positive product metric b=dt^2+kappa for every pointwise norm,
the full Frobenius tensor norm, and the product connection D. This
positive norm is not the Lorentzian metric g0. Let nabla denote the
Levi-Civita connection of g0, also used on derivative indices. All suprema
below include the full S3 and the shell time interval -Q<=t<=0, or the
background observer cylinder 0<=q<=Q as appropriate.

The retained mixed-wave ceilings, after the required negative-time join,
are

```
B1=10^4000000 exp(10^152017),
B2=10^10200000 exp(10^152019).
```

The fixed-frame physical trace and harmonic observer calculation in
[PROPER_SCALAR.md](PROPER_SCALAR.md) gives

```
H_a=512*10^7 B_a,
max_(0<=k<=4) sup |D^k p_a|_b <=H_a,
|xi1(q)|_b <=q^2 H1,
|eta1(q)|_b <=3q H1,
xi1=partial_e X_e|0,
eta1=nabla_e V_e|0=nabla_q xi1, V_e=partial_q X_e.
```

These are bounds in the same harmonic physical shell chart. The
coframe-to-product multiplier and its time derivatives have already
been charged in H_a. The endpoints simplify to

```
Cxi=Q^2 H1=B1/80, Ceta=3Q H1=24000 B1.
```

The covariant velocity variation eta1 differs from the product derivative
D_q xi1. In the product-parallel time/angular splitting it is

```
eta1=D_q xi1+tanh(q) (xi1)_S.
```

The zero observer initial variations follow from exact initial collar
equality and the fixed original anchor; they are not assumptions about
arbitrary initial hypersurfaces. The longitudinal first observer bounds
from the earlier first-cut note have not been combined with harmonic
second fields.

The strict coordinate margin T-Q=ell/16 and retained qualitative C2
dependence permit an unspecified common small parameter neighborhood
around this compact background path. No numerical lower bound for that
neighborhood is inferred here. Only at-zero spatial jets are estimated;
the actual nonlinear family need not be C2 into the upgraded mixed X13
topology for this calculation.

## 2. Explicit conversion to background covariant jets

Let A0=nabla-D. In a b-orthonormal product splitting its only nonzero
components are

```
(A0)^t_AB=alpha(t) delta_AB, alpha=cosh(t)sinh(t),
(A0)^A_tB=(A0)^A_Bt=beta(t)delta^A_B, beta=tanh(t).
```

Because the time/angular splitting and kappa are D-parallel, all D
derivatives of these tensors are their ordinary time derivatives. Thus

```
|A0|_b^2=3alpha^2+6beta^2,
|D A0|_b^2=3cosh(2t)^2+6sech(t)^4,
|D^2 A0|_b^2=12sinh(2t)^2+24sech(t)^4 tanh(t)^2.
```

On |t|<=T<1/100 the elementary hyperbolic bounds give

```
|A0|_b<1, |D A0|_b<4, |D^2 A0|_b<1.
```

For example cosh(2T)<(1-2T^2)^-1, |sinh(2t)|<=2|t|cosh(2T),
and |tanh(t)|<=|t| prove these comparisons without time sampling.
For a rank-r tensor, the connection difference acts on each tensor
slot, so |(nabla-D)Z|<=r|A0||Z|. Applying this also to the added
derivative slots gives, when max_(k<=3)|D^k p|<=H,

```
|nabla p| <= (1+2)H=3H,
|nabla^2 p| <= [1+5+2*4+6]H=20H,
|nabla^3 p| <= [1+9+7*4+26+2+20*4+24]H=170H.
```

For clarity, the unrounded coefficient ledger for the third line is

```
|D^3p|+9|A0||D^2p|
 +(7|D A0|+26|A0|^2)|Dp|
 +(2|D^2 A0|+20|A0||D A0|+24|A0|^3)|p|.
```

It comes from differentiating the second-order ledger
`|D^2p|+5|A0||Dp|+(2|D A0|+6|A0|^2)|p|`
termwise, then adding the four-slot connection action on nabla^2 p.
Every differentiated tensor slot is counted. Convenient evaluated
ceilings are therefore

```
P1=200H1 >=max_(k<=3) sup |nabla^k p1|_b,
P2=20H2  >=max_(k<=2) sup |nabla^k p2|_b.                 (1)
```

## 3. Inverse and connection jets, with the full Ricci quadratic term

Set A=g0^-1. In the b-orthonormal splitting it is
diag(-1,sech(t)^2,sech(t)^2,sech(t)^2), so its operator norm is at
most one, and nabla A=0. The actual inverse-metric amplitude derivatives
at zero are

```
(g^-1)_1=-A p1 A,
(g^-1)_2=2A p1 A p1 A-A p2 A.                            (2)
```

No inverse bound on a numerical nonzero parameter interval is used.
The connection difference C_e=nabla^(g_e)-nabla is the tensor

```
C_e^a_bc=(1/2)(g_e^-1)^ad
            [nabla_b(g_e)_dc+nabla_c(g_e)_db-nabla_d(g_e)_bc].
```

Define S(p)_dbc=nabla_b p_dc+nabla_c p_db-nabla_d p_bc.
Since nabla g0=0, differentiating exactly gives

```
C1=(1/2)A S(p1),
C2=(1/2)A S(p2)-2A p1 C1.                               (3)
```

The second inverse derivative in (2) multiplies S(g0)=0 and disappears
from C2. Its absence in (3) is a cancellation at the background, not an
omitted inverse-metric contribution. In particular the coefficient two
in the second term of C2 is required by the derivative convention.

In the curvature convention Ric(g0)=3g0, the complete Ricci difference is

```
Ric(g_e)_ij-Ric(g0)_ij=L(C_e)_ij+Q(C_e,C_e)_ij,
L(C)_ij=nabla_a C^a_ij-nabla_j C^a_ia,
Q(C,C)_ij=C^a_ab C^b_ij-C^a_jb C^b_ia.
```

Consequently the Eulerian defect jets, defined at fixed shell points,
are exactly

```
E1=L(C1)-3p1,
E2=L(C2)+2Q(C1,C1)-3p2.                                 (4)
```

There is no use of the bulk Einstein equation to replace intrinsic
Ricci, nor an identification of intrinsic Ricci with a projected bulk
Ricci component. The field equations enter through the actual shell
jets p1,p2 already supplied by the wave and physical reconstruction
interfaces.

Here is the full pointwise norm ledger. A permutation of tensor slots
preserves the Frobenius norm, and raising an index by A costs at most
one. Every trace in dimension four costs at most sqrt(4)=2. Hence

```
|C1| <=(3/2)|nabla p1|,
|nabla C1| <=(3/2)|nabla^2 p1|,
|nabla^2 C1| <=(3/2)|nabla^3 p1|,

|nabla C2| <=(3/2)|nabla^2 p2|
             +3|nabla p1|^2+3|p1||nabla^2 p1|,
|L(C)| <=4|nabla C|,
|Q(C,C)| <=3|C|^2.
```

In the last line the first quadratic term uses the trace bound 2,
followed by Cauchy–Schwarz; the second is a contraction between two
arrays and costs 1. This also proves the bound for nonspecial C without
using a gauge cancellation. Applying (4) and differentiating its first
line once yields

```
|E1| <=6|nabla^2 p1|+3|p1| <=9P1,
|nabla E1| <=6|nabla^3 p1|+3|nabla p1| <=9P1,
|E2| <=6|nabla^2 p2|+3|p2|
        +(51/2)|nabla p1|^2+12|p1||nabla^2 p1|
      <=9P2+(75/2)P1^2.                                 (5)
```

Here Ck means the maximum of the supremum norms of derivative orders
zero through k, not their sum. These are the required evaluated E1 C1
and E2 C0 bounds. They require
only first-shell C3 and second-shell C2, well below the available traces.

## 4. Covariantly consistent observer composition

The actual intrinsic geodesics use dimensionless past proper q. Their
normalization is exact, h_e(V_e,V_e)=-rho_b^2, equivalently
g_e(V_e,V_e)=-1. Thus the full intrinsic curvature input is

```
R_e(q,n)=Ric(h_e)_(X_e)(V_e,V_e)
       =-3+(E_e)_(X_e)(V_e,V_e).                         (6)
```

Use the fixed background connection nabla to differentiate tensor fields
along the two-parameter map X(e,q,n). At a fixed q,n,
`nabla_e [E_e(X_e)]=(partial_e E_e)(X_e)+nabla_(partial_e X_e) E_e`.
Differentiate (6) twice. Every term containing E0, nabla E0, or
nabla^2 E0 vanishes because E0 is the zero tensor throughout the shell.
The only surviving terms are

```
R2=E2(V0,V0)
    +2(nabla_xi1 E1)(V0,V0)
    +4E1(eta1,V0).                                      (7)
```

The two copies of the spatial derivative arise from differentiating
the parameter derivative and the spatial pullback separately. The
factor four in the last term is the two amplitude allocations times
the two velocity slots. No second observer position or velocity remains.
The identity is intrinsic: eta1=nabla_e V_e|0, not a coordinate velocity
derivative. In a fixed coordinate chart eta1^a=(partial_e V^a)|0+
Gamma0^a_bc xi1^b V0^c; replacing eta1 by the first term alone would
discard the connection contribution. Torsion freeness gives
eta1=nabla_q xi1 because the parameter coordinate fields commute.

Since |V0|_b=1, (1), (5), and (7) give the completed pointwise bound

```
sup |R2| <= Cr,
Cr=9P2+(75/2)P1^2+18P1 Cxi+36P1 Ceta
  =180H2+[1500000+3600Q^2+21600Q]H1^2.                  (8)
```

All constants are explicit. In exact rational arithmetic,

```
1500000+3600Q^2+21600Q
 =1536000034560009/1024000000 <1500001,
180*(512*10^7)=921600000000,
180*(512*10^7)+1500001*(512*10^7)^2 <4*10^25.
```

The retained B1 and B2 obey B1^2<B2: their base-ten prefactors satisfy
8000000<10200000 and their natural-exponential exponents satisfy
2*10^152017<10^152019. Therefore

```
Cr <4*10^25 B2
   <10^10200026 exp(10^152019).                          (9)
```

The weighted average needed below is immediately controlled, since
the weight Y^2 is nonnegative: sup_q |<R2(q)>_w|<=sup_(q,n)|R2|<=Cr.
Equation (8) preserves the contributions of the true second field,
quadratic first field, first position, and first velocity separately;
(9) is a conservative single ceiling, not a smallness statement.

## 5. Completed denominator and area interface at zero

Use the exact definitions and first-field calculations in the retained
[QUADRATIC_AREA_BOUND.md](../HDBLAST_QUADRATIC_PROPER_AREA_20261106/QUADRATIC_AREA_BOUND.md),
sections 4–5. Thus Y=n1^2-n2^2,
`<f>=(2pi^2)^-1 integral_S3 f dOmega`,
`<f>_w=<Y^2 f>/<Y^2>`, and `<Y^2>=1/6`. The relative actual cut
area density j_e is defined against the fixed background proper cut;
`j0=1` and `D_e/D0=<j_e>_w`. The weight is independent of amplitude.
The background denominator is positive,

```
D0=rho_b^3 cosh(q)^3*pi^2/3, D1=0.
```

Writing the retained first-field functions as F(q), I(q), B(q), and
K(n)=tr(HessY+(8/3)Y Id)^2, the exact Raychaudhuri calculation gives

```
j2=L2_exp+L2_shear+j1^2+L2_R,
L2_R(q,n)=-integral_0^q cosh(s)^2 R2(s,n)
                         [tanh(q)-tanh(s)] ds,
D2/D0=K_D+<L2_R>_w,
-3.65*10^-16<K_D<5.357*10^-16.                           (10)
```

There is no addition of the alternative determinant decomposition;
doing so would count known first-field terms twice. The full-family
quiet interval is only 0<=q<=q_D=ell/4. Exact quiet equality makes
R2 and all area derivatives zero there. Although the first tangent
has the larger quiet interval q<=ell, its gap cannot replace q_D in
the curvature kernel for the actual second field.

Set Cosh=(1-Q^2/2)^-2. Since cosh(s)^2<=Cosh and
0<=tanh(q)-tanh(s)<=q-s, the retained kernel constant is

```
C_R=(1/2)Cosh*(Q-q_D)^2
   =578027520000/671088639998361600000001
   <8.614*10^-13.
```

The now-evaluated curvature input (8) therefore proves

```
|<L2_R>_w| <=C_R Cr,
-3.65*10^-16-C_R Cr < D2/D0 <5.357*10^-16+C_R Cr,
|D2/D0| <5.357*10^-16+8.614*10^-13 Cr
        <10^14 B2
        =10^10200014 exp(10^152019).                    (11)
```

This is a complete explicit bound on the actual second denominator
derivative at zero for every 0<=q<=Q, conditional on the stated analytic
and geometric premises. It is not an assumed curvature projection or
a bound on only its known first-field part.

The same pointwise kernel also closes the full area-density derivative.
For completeness let d=Q-ell=9ell/16, U=23/(2*10^9), Fmax=11/125,
and Bmax=U[1+Qd+(4/3)d^2], as in the retained note. Its first-field
identities and |Y|<=1, 0<=K<=8 imply

```
0<=-L2_exp <=(1/3)Cosh Fmax^2 d^2,
0<=-L2_shear <=(4/3)U^2 d^4,
0<=j1^2 <=9Bmax^2.
```

Exact arithmetic gives
`(1/3)Cosh Fmax^2 d^2+(4/3)U^2 d^4<8.17*10^-16`
and `9Bmax^2<1.191*10^-15`. Consequently

```
-8.17*10^-16-C_R Cr < j2(q,n) <1.191*10^-15+C_R Cr,
sup |j2| <1.191*10^-15+8.614*10^-13 Cr <10^14 B2.         (12)
```

For the total area A_e=integral dA_e, A2/A0=<j2>. Using the exact
moments <Y^2>=1/6 and <K>=40/9 gives the slightly sharper known-part
interval

```
-1.362*10^-16-C_R Cr < A2/A0 <1.985*10^-16+C_R Cr.         (13)
```

Here the lower known budget is
`(1/18)Cosh Fmax^2 d^2+(20/27)U^2 d^4`, and the upper is
`(3/2)Bmax^2`. The curvature average for total area is the unweighted
round average, also bounded by the same pointwise Cr. These area
statements use the same retained hypersurface-orthogonal intrinsic
observer congruence and do not require a second observer/label equation.

## 6. Exact scope

The new evaluated step is the full intrinsic R2 bound and its transfer
to j2, D2/D0, and A2/A0 at zero. It uses actual first and second harmonic
shell jets and first harmonic observer motion, with the inverse,
connection, pullback, and proper-clock terms charged. It does not use
the parity cancellation of the selected scalar quadrupole to set any
of these even second quantities to zero.

The retained development, uniqueness, gauge/physical-shell, and common
observation-domain premises remain conditional. No explicit parameter
radius, uniform bound for E2(e) or R2(e), observer continuation tube,
positive denominator margin at nonzero amplitude, finite-amplitude
Taylor remainder, sign of area change, physical prediction or external
peer review is claimed. These enormous analytic
ceilings are model-normalization bounds rather than physical energy
or evidence for a cosmological origin.
