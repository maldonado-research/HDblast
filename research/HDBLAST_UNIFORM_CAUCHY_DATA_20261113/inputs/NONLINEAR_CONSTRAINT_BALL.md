# A quantitative positive ball for the selected nonlinear initial constraint curve

13 November 2026 UTC. Mathematical derivation with an independent
review recorded in AUDIT.md. This is an initial-constraint statement, conditional on the retained
actual smooth background and selected weighted inverse. It gives no
prescribed-time evolution interval, proper-observer tube, observable
remainder, physical calibration, or cosmological-origin conclusion.

## 1. Statement, exact data and norms

Use the same annulus and explicitly selected weights as the accepted
all-mode inverse checkpoint:

```
ell=10^-6, A={ell/2<z<5ell/2} x S3,
r0=(z-ell/2)(5ell/2-z)/ell^2,
kappa0=10^-8, q=kappa0^2 r0^2, w=exp(-1/r0),
gamma0=dy^2+rho(y)^2 kappa_S3, dy=-rho dz,
rho0=sqrt(det gamma0), Gamma=(Z2)^4.
```

All spatial tensor norms, derivatives and volume measures below use the
fixed positive metric gamma0. The density rho0 is fixed, not varied with
the unknown metric. Let

```
X=(dot H16_(q,w) scalar + dot H14_(q,w) one-form)^Gamma,
Y=(dot H12_(q,w) scalar + dot H12_(q,w) one-form)^Gamma.
```

Use Hilbert direct sums for X,Y and the literal derivative-square weighted
Sobolev norms of the accepted inverse. For u=(N,Z), the correction maps are

```
A_H N=w^2 q^4 P_g^*N,
A_M Z=w^2 q^2 P_m^*Z.
```

The second output is the ordinary symmetric contravariant tensor tau;
its canonical momentum-density correction is rho0*tau. Retain the actual
registered seed in all four canonical sectors. In ordinary fixed-density
variables write the affine physical data map as

```
gamma=gamma0+h,       h=e*h_seed+A_H N,
phi=phi0+psi,         psi=e*psi_seed,
pi=rho0*tau,          tau=e*tau_seed+A_M Z,
varpi=rho0*omega,     omega=e*omega_seed.                 (1)
```

The scalar and scalar-momentum density have no unknown correction.
Neither seed velocity sector is discarded.

For a tensor T and d=0,1,2 define the literal-sum norm

```
E_d(T)=sum_(j=0)^(14-d)
           ||w^-1 q^(j+d-2) nabla0^j T||_L2(dmu0).       (2)
```

Let `Zphys=E0(h)+E0(psi)+E1(tau)+E1(omega)` and equip the parameter and
potential product space with `|(e,u)|=|e|+||u||_X`. The accepted seed
and correction-map certificates give

```
Zphys <= B |(e,u)|,             B=10^200,               (3)
Zphys(direction i) <= B |direction i|.
```

Indeed the seed's four-sector E sum is at most `10^200`, while the
correction sums are bounded by `10^30||N||+10^27||Z||<=10^31||u||_X`.

Define the actual normalized nonlinear map and its preconditioner by

```
Chat=(H,M_density/rho0),
F(e,u)=w^-2 Chat(x0+e*v_seed+A*u),
G(e,u)=L^-1 F(e,u),        L=F_u(0,0).                  (4)
```

The accepted all-mode inverse gives `||L^-1||_(Y->X)<=10^5000`.
This inverse and the family use the same weights, density normalization,
Gamma-invariant domain and fixed correction range.

The quantitative conclusions established below are

```
Q=10^100,
U={|(e,u)|<1/(1024 Q B)},
sup_U ||D^2F||_(R x X)^2->Y <=10^1100,
sup_U ||D^2G||_(R x X)^2->X <=10^6100 <= K=10^6400.      (5)
```

Both maps are C2 on the open set U, which consists of positive initial
metrics. Put

```
a=10^-6402,                 r=10^-6404=K*a^2.           (6)
```

For every `|e|<=a`, the retained selected initial curve has a unique
representative `u(e)` in the closed X ball of radius r, satisfying
`F(e,u(e))=0`. This branch is C2 on an open parameter neighborhood of
`[-a,a]`, agrees with the previously selected local branch near zero,
and obeys

```
u(0)=u'(0)=0,
||u(e)||_X <=K e^2,
sup_(|e|<=a)||u'(e)||_X <=.0202,
sup_(|e|<=a)||u''(e)||_X <=8K.                           (7)
```

These are model-normalization bounds, not estimates of a measured effect.

## 2. Accepted product interfaces and fixed background multipliers

The companion `WEIGHTED_EMBEDDING_PRODUCTS.md` proves, in the norms (2),

```
E_(d+1)(nabla0 T)<=E_d(T),
E_(d1+d2)(T1 tensor T2)<=10^17 E_d1(T1) E_d2(T2),
||w^-2 S||_H12(q,w)<=10^18 E2(S).                       (8)
```

The product holds for `d1+d2<=2` with the specified tensor ranks. Its
low-factor pointwise bound is

```
|nabla0^j T| <=10^12 w q^(-j-d) E_d(T)
                    when j+3<=14-d.                    (9)
```

For the products in this proof the original factors have rank at most
three; intermediate uncontracted factors have rank at most six. Thus the
companion's product input-rank ceiling 25 and pointwise rank ceiling 32
are respected, including the derivative slots. A contraction with fixed
parallel tensors is charged by its Hilbert operator norm. Each Euclidean
trace costs at most two in dimension four; even a ceiling `2^16` is
absorbed in Q. Matrix composition has its usual Frobenius product bound.

Also `E_(d+1)(T)<=E_d(T)` whenever both sides are defined: q<=1 and the
higher grade uses fewer derivatives. Background tensors and the identity
are *not* asserted to belong to E0. They act as separately bounded
multipliers.

Here is an evaluated background ledger. Write `b=1/rho`,
`a_warp=rho_y/rho`, `p=phi0_y`, `eta0=phi0-1`. The only nonparallel
basic multipliers needed below are eta0, `p0=p dy`, and

```
Ric0=(3b^2+p^2-3a_warp^2)dy^2
       +(3b^2+p^2/3-3a_warp^2)gamma_tan.                (10)
```

In the radially parallel product frame their total polynomial coefficient
norm is below 20 and their degree is at most two. The retained actual
polynomial background vector field `D_s=q D_y` has coefficient norm
at most `10^-8` and raises degree by at most eight. Also
`d=q_y` and `beta=q a_warp` have coefficient norm below `10^-8` and
degree at most five.

For `B_j=q^j nabla0^j B`, the exact radial recurrence is
`(D_s-jd)B_j`; the angular recurrence is the fixed product derivative
plus `beta K_(rank+j) B_j`. All generated slot-insertion tensors are
parallel for the round product connection, so the intrinsic angular
derivative of these background expressions is zero. The collected
insertion norm is at most `4(rank+j)`. Polynomial degree is at most
`2+8j`; all original and derivative slots are counted. The loose
stage factor 1000 through j=14 therefore gives

```
sup_(j<=14)|q^j nabla0^j B| <=20*1000^14 <10^50.         (11)
```

The weighted Leibniz rule then bounds background multiplication on E_d
by `(2^(15-d)-1)10^50<10^55`. The same bound applies after the
necessary finite parallel contractions, on enlarging it to Q.
The full background scalar is `phi0=1+eta0`, with its unit part kept
separate. The fixed volume density is parallel for nabla0 and is already
removed in the variables tau,omega. No varying-density multiplier is
inserted.

Consequently Q bounds every graded product, needed fixed contraction,
background multiplier, pointwise constant and source-normalization map
used below. The useful algebra is the algebra of physical E0 images;
the raw weighted potential space is not declared an algebra.

## 3. Exact canonical equations and background subtraction

Identify h with its gamma0-raised endomorphism H when multiplying
matrices; this preserves all norms because gamma0 is parallel. Let

```
Ainv=(I+H)^-1=I+R,
Lambda=det(gamma0)/det(gamma)=det(I+H)^-1=1+s,
p_phi=p0+nabla0 psi.
```

In gamma0-orthonormal tensor notation the *full* constraints are

```
H_constraint=R(gamma)/2
  +2 Lambda[(tr_gamma tau)^2/3-tr(gamma tau gamma tau)]
  -Lambda omega^2/2
  -(gamma^-1)^(ij)(p_phi)_i(p_phi)_j/2
  -Uhat(eta0+psi),

(M_density/rho0)_i
  =-2 gamma_ik nabla0_j tau^jk
   -2 tau^jk nabla0_j h_ki+tau^jk nabla0_i h_jk
   +omega[(p0)_i+nabla0_i psi].                         (12)
```

Here `tr_gamma tau=gamma_ij tau^ij`, and the trace of the four-factor
matrix product is the complete squared momentum contraction. The exact
retained polynomial is

```
Uhat(eta)=-2/27+(14/9)eta^2+(50/27)eta^3
             -(1/6)eta^4-(4/9)eta^5-(2/27)eta^6.        (13)
```

These formulas preserve all metric components and both momentum sectors.
In particular (12) is not the scalar-conformal restriction of the metric
constraint. The momentum formula follows from the canonical density
divergence; rewriting it with nabla0 cancels the derivatives of rho0.

For curvature use the difference tensor of the two Levi-Civita connections:

```
C^k_ij=(gamma^-1)^(kl)
        (nabla0_i h_jl+nabla0_j h_il-nabla0_l h_ij)/2,
Ric(gamma)_ij=Ric0_ij+nabla0_k C^k_ij-nabla0_j C^k_ik
                 +C^k_kl C^l_ij-C^k_jl C^l_ik,
R(gamma)=(gamma^-1)^(ij) Ric(gamma)_ij.                  (14)
```

The retained background satisfies Chat(x0)=0. To define F in Y, subtract
this identity *before* estimating. For example, with `Bric=Ric(gamma)-Ric0`,

```
R(gamma)-R(gamma0)=R:Ric0+(I+R):Bric,
Ainv:(p0+nabla psi)^2-I:p0^2
 =R:p0^2+2(I+R):(p0 tensor nabla psi)
                 +(I+R):(nabla psi tensor nabla psi).
```

The potential difference is its finite polynomial expansion in psi with
coefficients depending on eta0. Every displayed difference belongs to E2
by (8), the grade inclusions and the background multiplier bounds. The
momentum terms already contain tau or omega. Thus the physical constraint
residual belongs to E2, and (8) puts F in Y. This reasoning never assigns
an E2 norm to the nondecaying background Ricci or potential separately.

## 4. Uniform inverse and determinant compositions

On the open set U from (5), put `delta=Zphys`. Equations (3),(5) imply
`Q delta<1/1024`. In particular (9) gives `|H|_op<1/1024`; the initial
metric is positive, with relative eigenvalues in `(1023/1024,1025/1024)`.

For matrix fields use the unital algebra consisting of a parallel identity
plus an E0 field. Multiplication by `I+R` on any admissible E_d tensor
has bound `1+Q E0(R)` when it acts on a single matrix/vector index.
The Neumann series, convergent in this unital algebra, gives

```
E0(R)<=2 E0(h),
E0(D Ainv[h1])<=4 E0(h1),
E0(D^2 Ainv[h1,h2])<=16Q E0(h1)E0(h2).                 (15)
```

For example the exact derivatives are
`D Ainv[h1]=-Ainv h1 Ainv` and
`D^2 Ainv[h1,h2]=Ainv h1 Ainv h2 Ainv+Ainv h2 Ainv h1 Ainv`.
They use the full matrix inverse.

A convenient determinant construction in the same algebra is

```
l=log det(I+H)=sum_(k>=1)(-1)^(k+1) tr(H^k)/k,
Lambda=exp(-l)=1+s.
```

The trace map has norm at most two. The product-series estimates and
`Q E0(h)<1/1024` give

```
E0(l)<=4 E0(h),
E0(Dl[h1])<=4 E0(h1),
E0(D^2l[h1,h2])<=8Q E0(h1)E0(h2),
E0(s)<=8 E0(h),
E0(DLambda[h1])<=8 E0(h1),
E0(D^2Lambda[h1,h2])<=48Q E0(h1)E0(h2).                (16)
```

For the last two lines use
`D Lambda=-Lambda Dl` and
`D^2 Lambda=Lambda[(Dl)(Dl)-D^2l]`. The scalar multiplier
`1+s` has module bound below two. The convergent series agree pointwise
with the actual inverse determinant, since I+H is positive. Their
derivative series converge uniformly on each closed subball of U.

Thus the maps are C2 on an *explicit open domain*, rather than on an
unspecified neighborhood supplied by a qualitative IFT. The graded algebra
also controls spatial derivatives falling on arbitrary state factors in
these series. It does not always place the highest derivative on one of
the two parameter directions.

## 5. Full second-derivative ledger in physical graded norms

For physical directions `z_i=(h_i,psi_i,tau_i,omega_i)`, write

```
D_i=E0(h_i)+E0(psi_i)+E1(tau_i)+E1(omega_i).
```

All estimates in this section hold at every point of U. We use
`Q delta<=1/1024` only to absorb harmless extra state factors, not to
discard any equation term.

### 5.1 Curvature

The three-permutation first-derivative map in (14) has norm at most three.
Equations (8),(15) consequently give

```
E1(C)<=3 delta,
E1(DC[z_i])<=4 D_i,
E1(D^2C[z_1,z_2])<=16Q D_1 D_2.                        (17)
```

In the Ricci difference the two linear traces together cost at most four.
The two quadratic contractions together cost at most eight. Hence

```
E2(Bric)<=13 delta,
E2(DBric[z_i])<=17 D_i,
E2(D^2Bric[z_1,z_2])<=400Q D_1D_2.                     (18)
```

For example the last line is bounded by
`[64+768 Q delta+256]Q D_1D_2<=400Q D_1D_2`.
The scalar trace by the identity costs two; its perturbation has module
cost at most `Q E0(R)`. Differentiate both factors of
`R:Ric0+(I+R):Bric`. The term with fixed Ric0 costs at most
`16Q^2 D_1D_2`. The other terms cost at most
`[208 Q delta+136+1200]Q D_1D_2`. Therefore the curvature contribution
to the second derivative of H_constraint is bounded by

```
E2(D^2[R(gamma)/2]) <=10^4 Q^2 D_1D_2.                 (19)
```

The background multiplier in this estimate is essential.

### 5.2 Scalar gradient and the actual polynomial potential

Multiplication by `p_phi=p0+nabla psi`, increasing the grade by one,
costs at most `2Q`: the background p0 has its multiplier bound and
the perturbation has E1 norm at most delta. The exact second derivative
of `Ainv:p_phi^2` consists of

```
(D^2 Ainv):p_phi^2,
(D Ainv[h_1]):(nabla psi_2 tensor p_phi+p_phi tensor nabla psi_2),
(D Ainv[h_2]):(nabla psi_1 tensor p_phi+p_phi tensor nabla psi_1),
Ainv:(nabla psi_1 tensor nabla psi_2
                        +nabla psi_2 tensor nabla psi_1).
```

Equations (15) and (8) bound these by respectively
`64Q^3`, `16Q^2`, `16Q^2`, and `6Q`, times `D_1D_2`.
The factor one-half and all signs can be bounded by the safe common value

```
E2(D^2[-|dphi|_gamma^2/2]) <=10^4 Q^3 D_1D_2.          (20)
```

The potential's second derivative is exactly
`Uhat''(eta0+psi) psi_1 psi_2`. Its degree is four; the sum of the
absolute coefficients of Uhat is below ten, so that of Uhat'' is below
300. Multiplication by eta0+psi has norm at most `2Q`; its background
part is a multiplier, not an E0 element. The product of the two scalar
directions costs Q in E0 and E0 embeds into E2. Thus

```
E2(D^2[Uhat(eta0+psi)]) <=10^6 Q^5 D_1D_2.             (21)
```

This includes the nonquadratic matter potential at arbitrary states of U.

### 5.3 Both canonical momentum-energy sectors

Put
`S=(tr_gamma tau)^2/3-tr(gamma tau gamma tau)`.
Represent gamma as I+H acting on matrix indices. Its module norm is below
two; traces and fixed contractions of these finite expressions are
covered by 100. Counting its two metric and two momentum factors gives

```
E2(S)<=10^3 Q delta^2,
E2(DS[z_i])<=10^4 Q delta D_i,
E2(D^2S[z_1,z_2])<=10^4 Q D_1D_2.                     (22)
```

To see the last allocation explicitly: there are at most sixteen ordered
two-factor differentiations. When both derivatives hit the momentum
factors, the graded product costs Q and the two metric modules cost at
most four. When one derivative hits a metric factor, the remaining
momentum state contributes an additional `Q delta<=1`. When both hit
metric factors, the two remaining momentum states contribute
`(Q delta)^2<=1`. Thus each term costs at most `4Q D_1D_2` before the
fixed contraction ceiling 100. The same counting proves the first two
lines. No momentum component or trace term is suppressed.

Use (16) and differentiate `2 Lambda S`. The four product-rule terms
are bounded before the factor two by
`2*10^4 Q`, `2*(8*10^4 Q^2 delta)`, and
`48*10^3 Q^3 delta^2`, times `D_1D_2`. Consequently

```
E2(D^2[2 Lambda S]) <=10^6 Q^2 D_1D_2.                 (23)
```

For the distinct scalar momentum-density term, use
`E2(omega^2)<=Q delta^2`,
`E2(D omega^2)<=2Q delta D_i`, and
`E2(D^2 omega^2)<=2Q D_1D_2`, followed by the same determinant derivatives.
This gives

```
E2(D^2[-Lambda omega^2/2])<=10^3 Q^2 D_1D_2.            (24)
```

The normal scalar momentum is not being substituted for the prescribed
canonical scalar density omega.

### 5.4 The full momentum constraint

In the second line of (12), the affine fields have zero second derivative.
The second derivative therefore contains both orders of the pairs
`h_i nabla tau_j`, `tau_i nabla h_j`, and
`omega_i nabla psi_j`, with the exact coefficients and index contractions
in (12). All have grade two. The fixed background `gamma0 nabla tau`
and `omega p0` terms are linear and have zero second derivative.
The contraction costs and numerical coefficients are covered by

```
E2(D^2[M_density/rho0])<=100Q D_1D_2.                  (25)
```

This uses the actual density formula, not a variable-density replacement.

### 5.5 Combined physical and weighted-source bounds

Summing (19)--(25), with the scalar and one-form outputs measured by the
sum of their E2 norms, gives

```
E2(D^2 Chat_H)+E2(D^2 Chat_M)
 <=10^7 Q^5 D_1D_2.                                    (26)
```

For clarity, the sum of the six displayed coefficient ceilings is below
`10^7 Q^5` for every Q>=1. A Hilbert output norm is no larger than this
sum. The physical data map (1) is affine, so there is no additional
second derivative of that map. Equations (3),(8),(26) imply

```
||D^2F(e,u)[xi1,xi2]||_Y
 <=10^7 Q^6 B^2 |xi1||xi2|
 =10^1007 |xi1||xi2| <=10^1100 |xi1||xi2|.              (27)
```

The use of Q for the source normalization deliberately exceeds its
proved `10^18` constant. Applying the accepted inverse gives (5).

## 6. Domain, symmetry and finite-regularity justification

The products, relative-background formulas, matrix Neumann series,
log-determinant series and finite matter polynomial construct F as a C2
map on the explicit open set U. All series estimates are in physical
graded Banach spaces, with separately bounded identity/background
multipliers. For potentials first defined by compact smooth tests,
the estimates extend continuously to X. Approximation can be kept in U
because the point under consideration has a positive distance from its
boundary. The differentiated product and composition estimates extend
in the same way. In particular this does not rely on assigning a known
radius to the older qualitative IFT neighborhood.

The weights, background, exact canonical operators and correction range
commute with Gamma. The registered seed is Gamma-invariant. Products,
contractions, inverse metrics and determinants preserve that equivariance;
therefore the nonlinear source lies in the same invariant Y on which
the selected inverse acts. The fixed density is a density, so reflections
do not introduce an orientation sign into its normalization.

The accepted physical correction maps give general bounded maps into
`H14_0(A)` and `H13_0(A)` and bounded zero extensions. The seed has fixed
compact support inside A. Thus the constructed data equal the background
on the same open shell and vertex collars and have the retained physical
finite regularity. The constraint equations hold in the actual canonical
variables, with no distributional annulus-boundary layer. No claim of one
spatially smooth family of all orders is made.

## 7. Explicit contraction and identification of the selected branch

The exact background and the full linearized seed give

```
G(0,0)=0, G_u(0,0)=I, G_e(0,0)=0.                     (28)
```

For all `(e,u)` in the closed product ball of radii (6), the segment
from the origin lies in U. Integral Taylor and (5) yield

```
||I-G_u(e,u)||<=K(|e|+||u||),
||G(e,0)||<=K e^2/2.                                  (29)
```

The exact finite comparisons are

```
K(a+r)=1/100+1/10000=101/10000<1/2,
K a^2/2=r/2,
Q B(a+r)=1.01*10^-6102<1/2048.                         (30)
```

Hence the closed product ball is strictly inside U. The map
`T_e(u)=u-G(e,u)` has contraction constant at most `.0101` and maps the
closed radius-r ball into a ball of radius at most `.5101 r`.
Banach's fixed-point theorem gives a unique fixed point there, hence
`F(e,u(e))=0` for every `|e|<=a`.

The same radius-r ball works for the larger open parameter interval
`|e|<(101/100)a`: its contraction constant is at most `.0102` and its
self-map bound is at most
`[(1.01)^2/2+.0102]r=.52025 r<r`. It also remains inside U.
The invertible derivative `G_u` and the C2 map on U give C2 parameter
dependence throughout this interval, including an open neighborhood of
the endpoints of the stated closed interval.

The fixed-point estimate sharpens to
`||u(e)||<=K e^2/[2(1-.0101)]<=K e^2`. Differentiating G=0 on this branch,
and using inverse norm at most two, gives

```
||u'||<=2K(a+r)=.0202,
||u''||<=2K(1+||u'||)^2<=8K,
u''(0)=-G_ee(0,0).                                    (31)
```

Equations (28) also give `u(0)=u'(0)=0`. The formula for `u''(0)` is
exactly the previously selected weighted second initial correction.

Finally, this construction uses the same nonlinear equation and correction
range as the retained local implicit family. For sufficiently small e,
that prior branch lies inside the present radius-r ball by its continuity
at zero; fixed-point uniqueness identifies the two. The present argument
therefore supplies an explicit interval for the same selected germ, with
no parameter rescaling or new arbitrary completion. Its extension over
the stated interval is selected by the same ball uniqueness.

## 8. Premises, companion certificates and scope

The retained analytic premises are the actual smooth background on the
annulus, its canonical equations and polynomial potential, the registered
full seed satisfying the linearized constraints, and the selected
Gamma-invariant weighted isomorphism. Its quantitative norm `10^5000`
is the accepted all-mode inverse result. Nothing in this argument replaces
those premises with numerical trajectory samples.

Companion inputs (private reviewed originals; public editorial changes are recorded in AUDIT.md):

* `PHYSICAL_CORRECTION_MAPS.md`, accepted SHA256
  `2422e6fa073cc5878b0949cc596d8087b95e0e62c288a2ffd6df9cfe26a92e7a`:
  literal E0/E1 correction bounds, physical H14/H13 maps, fixed-density
  conversion and zero extension.
* `WEIGHTED_EMBEDDING_PRODUCTS.md`, accepted SHA256
  `9b77705494c5402bb03708cf304d9046e0fcd21401abf31f2d9579d0f4c336e0`:
  tensor embedding, graded products and the literal w^-2 source conversion.
* `SEED_ENVELOPE.md`, accepted SHA256
  `c02abb7c509f5fd7e36e573887592d2629479ffd3a9f14d8b05e2f7202e85b5b`:
  all four actual seed sectors and their common E sum `10^200`.

Retained paths below are relative to
the repository `research/` directory:

| Source | SHA256 |
|---|---|
| `HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/CANONICAL_CONSTRAINTS.md` | `4e25e993f9653b530e30e2a2d3fefcb0a006f04050534190b99a98126ba04f63` |
| `HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md` | `95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b` |
| `HDBLAST_WEIGHTED_REGULARITY_20261106/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json` | `fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0` |

The conclusion is a uniform *initial constraint-family* interval with
fixed finite regularity and positive initial metric. The separate initial
harmonic completion, nonlinear shell evolution, proper-clock domain,
noncrossing and observable remainders still require their own uniform
estimates. The earlier zero-amplitude geometric bounds are not those
uniform estimates. No such later conclusion is inferred from (6).
