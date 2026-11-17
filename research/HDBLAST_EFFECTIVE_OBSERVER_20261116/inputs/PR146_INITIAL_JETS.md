# Differentiated compatible initial jets of the actual parameter curve

16 November 2026. New proof at immutable repository base
`fb2ba28453887056166e9f9b87a5a78b9bc22360`. The source identities are pinned
in `SOURCE_PINS.json`. All retained analytic
background, actual point-coefficient, canonical-family and smooth
frozen-linear hypotheses remain in force. This note proves the initial
parameter interface; evolved derivative existence is proved separately.

## 1. Statement and exact spaces

Put `ell=10^-6`, `Sigma=[0,6ell] x S3`, `T=13ell/8`, and retain the fixed
positive sixteen-field product bundle, full Frobenius symmetric-block
norm, connection `D`, measure, and natural `Gamma=(Z2)^4` action. All norms
below use the specified one-copy or identical-two-copy convention already
chosen for the actual family. Let

```
I=[-a_iter,a_iter], a_iter=10^(-10^500209),
Y14=direct_sum_(j=0)^14 H^(14-j)(Sigma;E),
||(K_j)||Y14=sum_(j=0)^14 ||K_j||H^(14-j).
```

The actual elapsed solution is `v_epsilon(s)=J_E V_epsilon(-s)`.
The time-slot isometry `J_E` is parallel, parameter independent, and
preserves the four/twelve shell split. Its initial pair is

```
p(epsilon)=J_E V_epsilon(0), q(epsilon)=-J_E V_epsilon,t(0).
```

`UNIFORM_CAUCHY_DATA`, Section 6, already proves that this pair is C2 on
an open parameter neighborhood of its larger original closed interval,
with values in `H14 x H13`. In particular this is an ordinary two-sided
initial-pair statement near every point of I. We use that result, including
the actual lapse, current half-action canonical conversion and all four
mixed metric velocities. We do not assert a new initial-pair theorem.

Let `J_j(epsilon)` be the nonlinear initial temporal recursion, with
`J0=p,J1=q`, through `j=14`. Then

```
epsilon -> J(epsilon)=(J0,...,J14) is C2 into Y14
on an open parameter neighborhood of I.                     (1)
```

The entries are the actual initial temporal jets of the retained
`C X14` solution on I. Write

```
W_j=partial_epsilon J_j, Z_j=partial_epsilon^2 J_j.
```

Both tables are continuous into Y14, have finite uniform norms on I,
are Gamma invariant, and are supported in the same closed annulus
`A0=[ell/2,5ell/2] x S3`. Their shell Dirichlet and outer Dirichlet
compatibility holds through `j=13`, and shell Neumann compatibility
through `j=12`, in the ordinary trace spaces at those orders. No ordinary
normal trace of an H1 jet `J13,W13,Z13`, and no ordinary boundary trace
of an L2 jet `J14,W14,Z14`, is asserted.

## 2. Literal principal normalization and full first jets

The exact elapsed equation, after the fixed-connection first-jet
conversion, is

```
v_ss=A^i(x,v) D_i v_s+B^ij(x,v) Hsp_ij v+F(x,v,Q),
Q=Dv=(v_s,D_1v,...,D_4v),
F=C(x,v)+E(x,v)Q+G(x,v)(Q,Q).                         (2)
```

`Hsp` is the complete spatial product Hessian, including its derivative
index connection. In the underlying principal tensor put

```
a=-q^ss>0, b^i=q^si, c^ij=q^ij,
A=2a^-1 b, B=a^-1 c,
C=-a^-1 cbar, E=-a^-1 bbar, G=-a^-1 dbar.             (3)
```

The actual all-state lower bound is `a>=1/18000`. Thus (3) uses the
state-dependent inverse and the minus-complete-lower source; neither
the inverse nor the lower source is frozen at the background. The
coordinate-to-product conversion is the literal one in
`MIXED_INITIAL_JET_PRODUCTS`, (2)--(4), including `P=Q-Omega v`, the
ordered connection products, and all five first-jet directions.
In particular `C(x,0)=0` identically in the base variables.

For completeness, for either principal numerator `k=2b` or `k=c`,
the map `K=a^-1 k` has vertical derivatives

```
K_V[h]=a^-1 k_V[h]-a^-2 a_V[h] k,
K_VV[h,l]=a^-1 k_VV[h,l]
 -a^-2 {a_V[h]k_V[l]+a_V[l]k_V[h]+a_VV[h,l]k}
 +2a^-3 a_V[h]a_V[l]k.                              (4)
```

Consequently both `A_VV` and `B_VV` include the inverse Hessian and
the numerator/inverse cross terms. They must not be assumed zero or
omitted.
The original parameter changes initial data only. Background, potential,
shell generator, reference connection, and coefficient chart are fixed;
there is no explicit external epsilon derivative in (2).

All coefficient derivatives in the sequel are the horizontal/vertical
bundle derivatives of the pinned coefficient proof. A horizontal
derivative includes the induced actions on every output, base, and
state-input slot and the term `-K_V[Omega_i v]`. Spatial and time
derivatives commute for this fixed product connection; no spatial
covariant derivatives are commuted past one another.

## 3. Finite nonlinear Cauchy operator

For `0<=r<=12`, insert a formal time polynomial with jets
`J0,...,J_(r+1)` into the right side of (2), and define

```
R_r(J0,...,J_(r+1))
 =partial_s^r {A(x,v)D_i v_s+B(x,v)Hsp_ij v
                      +F(x,v,Dv)}|s=0,
J_(r+2)=R_r(J0,...,J_(r+1)).                          (5)
```

Here explicit base-time dependence is differentiated as well. Formula
(5) means a finite Taylor coefficient calculation and requires no
temporal solution. In particular its full first-jet substitutions are

```
partial_s^j Q|0=(J_(j+1),D_1J_j,...,D_4J_j).
```

The highest jet on the right is `J_(r+1)`, from a time first-jet or
mixed-principal factor. The coefficient and spatial-Hessian factors
use at most `J_r`. There is no `J_(r+2)` on the right, precisely because
the actual inverse in (3) has already normalized the acceleration.

The retained nonlinear initial-jet proof gives

```
I14(J)<=10^6503000 d, d=||p||H14+||q||H13,
```

under `d<=10^-6503020`; the canonical interval lies strictly inside
that guard. Its Section 6 also states the finite jet map is C2 on the
relative open pair ball `d<2*10^-6503020`. We give the needed graded
justification explicitly, so that differentiating (5) is not a formal
interchange assumption.

For a row with `k` ordered spatial derivatives, `r+k<=12`, the finite
chain/Leibniz expansion of (5) consists of coefficients evaluated at
`(0,x,J0)` times ordered products of factors `D^b J_j`. Replace `C`
by the exact radial mean-value identity

```
C(x,U)=integral_0^1 C_V(x,theta U)[U] dtheta.          (6)
```

Every leaf then has at least one field factor. Give `D^b J_j` intrinsic
order `b+j`. The sum of intrinsic orders in each product is at most
`r+k+2<=14`. Replacing a field factor by either of two Frechet
directions preserves its intrinsic order; a coefficient variation adds
an order-zero state-direction factor. Thus the same total-order bound
holds for both differentials of the row.

Choose a maximal-order field or direction factor for L2. Every other
factor has order at most seven: two factors of order at least eight
would exceed fourteen. For a low factor of intrinsic order `b+j<=7`,
the three further spatial derivatives required by the retained tensor
`H3 -> L-infinity` embedding are available in `H^(14-j)`. The original
tensor slots, its derivatives, and this embedding reserve remain within
the sixteen-slot certificate. Coefficient arrays, including all created
direction slots, use their certified full-array operator norms directly;
they are not embedded as high-rank fields. The maximal factor has all
the required `b` derivatives in L2. This proves continuous multilinear
product estimates for each leaf into `H^(12-r)`.

There are finitely many leaves. On the above relative open pair ball
the triangular retained induction gives `I14<2*10^-20<10^-12` and
strict state-patch membership. At most `r+k<=12` horizontal/vertical
coefficient differentiations arise before parameter variation. Two
Frechet variations require at most fourteen; extraction (6) requires
at most fifteen; one further coefficient derivative controls continuity
of the second differential. All are within the retained mixed-order
sixteen reserve. These counts include explicit base derivatives and
state-input derivatives together; they are not separate allowances of
sixteen for each kind of derivative.

Here is also a direct justification of the differentiability. Apply the
pointwise finite-dimensional first and second Taylor formulas to each
coefficient in a leaf. Its direction/product terms are precisely those
obtained by the product rule above. On any sufficiently small common
neighborhood of a given hierarchy, the state segments stay in the
coefficient patch. The preceding high/low estimates bound the Taylor
remainders in the row's target norm, uniformly for directions in a
bounded ball. For continuity of the Hessian apply the same estimates
to a difference of second-differential leaves, using the additional
coefficient derivative just counted and multilinearity for the other
factors. These estimates imply the Banach-space Frechet remainders,
not only pointwise Gateaux derivatives.

Equivalently one can first perform these finite operations for smooth
collar-supported sections and then pass by strong Sobolev density. The
same product estimates make the operations and their first two
differentials converge in the indicated targets; strict state margin
is preserved for sufficiently close approximants. This proves

```
R_r: product_(j=0)^(r+1) H^(14-j) -> H^(12-r)
is C2 on the relevant open hierarchy neighborhood.              (7)
```

The directions in (7) need not be small: boundedness of a differential
is homogeneous in its one or two direction norms. Finite triangular
composition of the thirteen rows proves that the Cauchy operator
`Jcal:(p,q)->(J0,...,J14)` is C2 into Y14 on the stated relative open
pair ball. Composition with the already proved actual initial curve
establishes (1). This argument uses two parameter derivatives only;
the extra coefficient derivatives certify spatial composition, not a
third derivative of the actual parameter curve.

## 4. Exact first and second parameter recursions

In the following identities `v,w,z` mean formal time jets with entries
`J_j,W_j,Z_j`, respectively; evaluation after `r` time derivatives uses
only already known entries through `r+1`. Let `p_w=Dw` include its
time component. The first differentiated row is

```
W_(r+2)=partial_s^r {
 A^i D_i w_s+B^ij Hsp_ij w
 +F_V[w]+F_Q[p_w]
 +A_V^i[w]D_i v_s+B_V^ij[w]Hsp_ij v}|0.             (8)
```

The second differentiated row is

```
Z_(r+2)=partial_s^r {
 A^i D_i z_s+B^ij Hsp_ij z
 +F_V[z]+F_Q[Dz]
 +A_V^i[z]D_i v_s+B_V^ij[z]Hsp_ij v
 +F_VV[w,w]+F_VQ[w,p_w]+F_QV[p_w,w]+F_QQ[p_w,p_w]
 +A_VV^i[w,w]D_i v_s+B_VV^ij[w,w]Hsp_ij v
 +2A_V^i[w]D_i w_s+2B_V^ij[w]Hsp_ij w}|0,          (9)
```

for every `0<=r<=12`. In particular (9) includes the state Hessians
of both principal maps and both principal cross terms. All principal
terms on the right have positive signs. The `A_V[z]`, `B_V[z]`,
`F_V[z]`, and `F_Q[Dz]` terms are the linear part in z, and are
included once.

For clarity, the exact source derivatives in (8)--(9) are

```
F_V[w]=C_V[w]+E_V[w]Q+G_V[w](Q,Q),
F_Q[p_w]=E p_w+G(p_w,Q)+G(Q,p_w),
F_VV[w,w]=C_VV[w,w]+E_VV[w,w]Q+G_VV[w,w](Q,Q),
F_VQ[w,p_w]=E_V[w]p_w+G_V[w](p_w,Q)+G_V[w](Q,p_w),
F_QV[p_w,w]=F_VQ[w,p_w],
F_QQ[p_w,p_w]=G(p_w,p_w)+G(p_w,p_w).                 (10)
```

The actual G is symmetric in the combined direction/fiber inputs;
(10) displays both allocations without depending on an implicit factor
of two. There are five direction slots in Q and in `p_w` throughout.

For a second explicit way to check temporal multiplicities, let
`A_k=partial_s^k A(x,v)|0`,
`A'_k=partial_s^k(A_V(x,v)[w])|0`, and
`A''_k=partial_s^k(A_V[z]+A_VV[w,w])|0`, and likewise for B.
Then the A contribution to the three recursions is, respectively,

```
sum_(k=0)^r binom(r,k) A_k D_i J_(r-k+1),
sum_(k=0)^r binom(r,k) [A_k D_i W_(r-k+1)
                                      +A'_k D_i J_(r-k+1)],
sum_(k=0)^r binom(r,k) [A_k D_i Z_(r-k+1)
            +2A'_k D_i W_(r-k+1)+A''_k D_i J_(r-k+1)].          (11)
```

The B contributions replace `D_i` and index `r-k+1` by `Hsp_ij`
and index `r-k`. The F contribution is the total time derivative of
the appropriate line of (8) or (9). Formula (11) includes every explicit
time derivative of the maps and every Leibniz multiplicity.

The legitimacy of (8)--(9) now follows from (7) and the ordinary
Banach-space C2 chain rule. For smooth jets their finite coefficient
extraction commutes with parameter differentiation as an algebraic
identity; continuity of the proved row differentials passes that
identity to these actual finite Sobolev jets. No evolution regularity
in epsilon is presumed in this interchange.

## 5. Exact first Cauchy pairs for the variations

The initial pairs for the first and second variation equations are

```
w(0)=W0=p', w_s(0)=W1=q',
z(0)=Z0=p'', z_s(0)=Z1=q''.                          (12)
```

These have no independent second physical velocity. To spell out their
actual initial construction, use the notation of `UNIFORM_CAUCHY_DATA`:
let `P=L I X`, `Qvel=L I Xdot`, and let `Tcal(P)=F_-z(r(P))` be
its root/straightening chart. Primes here are the original epsilon
derivatives of the completed raw pair. The exact elapsed data are

```
W0=J_E Tcal_1[P'],
W1=-J_E {Tcal_1[Qvel']+Tcal_2[P',Qvel]},
Z0=J_E {Tcal_1[P'']+Tcal_2[P',P']},
Z1=-J_E {Tcal_1[Qvel'']+Tcal_2[P'',Qvel]
                 +2Tcal_2[P',Qvel']+Tcal_3[P',P',Qvel]}.        (13)
```

All chart derivatives are evaluated at the same nonzero P. There is no
explicit first time derivative of the normalization at the initial
section, as proved in the retained source; consequently these are the
complete formulas. In particular the factor two and third state
derivative in `Z1` cannot be removed.

The raw `X,Xdot` in (13) are exactly the values and all sixteen
velocities in retained equations (5)--(7). Their inputs have

```
h'=h_seed+A_H N', h''=A_H N'',
psi'=psi_seed, psi''=0,
tau'=tau_seed+A_M Z', tau''=A_M Z'',
omega'=omega_seed, omega''=0,
ell_lapse'=f, ell_lapse''=0, rho0'=rho0''=0.
```

The retained inverse is
`K_ij=2lambda(gamma_ik gamma_jl tau^kl-T gamma_ij/3)` and
`Pi=lambda omega`, with `lambda=(det gamma/det gamma0)^-1/2` and
the prescribed ADM lapse `alpha=rho(1+epsilon f)`. Thus affine
`omega` or scalar canonical field does not imply affine Pi or wave
velocity. Differentiating the retained full raw formulas, rather than
only the diagonal fields, defines `P',P'',Qvel',Qvel''` in (13).
The mixed covector velocity and the metric-dependent reference
contraction remain present. Neither the older doubled-action inverse
nor background lapse is substituted.

## 6. Uniform finite norms, secants and Taylor remainders

Define the finite quantities

```
K_l=max_(epsilon in I) ||partial_epsilon^l J(epsilon)||Y14,
l=0,1,2,
omega_2(delta)=sup_{a,b in I, |a-b|<=delta}
                     ||J''(b)-J''(a)||Y14.                     (14)
```

Compactness and (1) prove `K_l<infinity` and
`omega_2(delta)->0`. This is a uniform parameter modulus, not an
unevaluated numerical coefficient claimed as a new certificate.
For example if `L_1,L_2` denote the finite suprema along the compact
pair curve of the norms of `D Jcal,D^2 Jcal`, the retained initial
pair bounds give

```
K_1<=L_1*10^4200,
K_2<=L_1*10^10400+L_2*(10^4200)^2.                  (15)
```

These merely show explicitly how already proved initial bounds join
the continuous finite jet map. No new numerical value of `L_1,L_2`
is asserted or needed.

For `a,b=a+h in I`, with `h!=0` in the quotient expression, the
Banach-valued fundamental theorem of calculus gives

```
||J(b)-J(a)||Y14<=K_1 |h|,
||J'(b)-J'(a)||Y14<=K_2 |h|,
||J(b)-J(a)-hJ'(a)||Y14<=(K_2/2)|h|^2,
||(J'(b)-J'(a))/h-J''(a)||Y14<=omega_2(|h|),
||J(b)-J(a)-hJ'(a)-(h^2/2)J''(a)||Y14
                           <=(|h|^2/2)omega_2(|h|).             (16)
```

The last identity uses
`h^2 integral_0^1(1-theta)[J''(a+theta h)-J''(a)] dtheta`.
All bounds are uniform over allowed base points and signs of h.
They include all temporal initial jets, hence also every lower initial
mixed norm required by the X13 first-variation and X12 second-variation
estimates. No Lipschitz modulus for `J''`, and no third epsilon
derivative, is inferred. These are the exact initial secant and Peano
inputs for identifying evolved parameter derivatives. The same statements,
with finite maxima and modulus taken on the compact closure of the
slightly larger J constructed in Section 8, hold for secants on that
enclosure. Thus they also give two-sided secants at the endpoints of I.

## 7. Support, symmetry and boundary compatibility

For each epsilon the initial pair and its first two parameter
derivatives vanish outside A0 and on fixed open shell/outer collars,
by the retained local completion and zero-extension proof. These
support conditions define closed linear subspaces of each of the
relevant Sobolev spaces; differentiating a curve valued there preserves
them. Every row (5), (8), (9) is local, `C(x,0)=0`, and its
non-background field/direction factors vanish on these collars.
Induction preserves support in A0 for J,W,Z. The finite products do
not create annulus-edge distributions: they are Sobolev products of
the zero-extended sections, whose validity is part of the retained
initial interface and of the density argument in Section 3.

All maps are natural under Gamma with its action on every tensor slot,
so the same induction preserves Gamma invariance. Nonlinear angular
products remain allowed; no finite harmonic or isotypic restriction
is imposed. Fixed shell projections `Pi_D,Pi_N`, the normal `D_z`,
and outer projection commute with epsilon differentiation wherever
their ordinary trace maps are bounded. For `m=14-j`,

```
Pi_D J_j=Pi_D W_j=Pi_D Z_j=0 at z=0, j<=13,
J_j=W_j=Z_j=0 at z=6ell, j<=13,
D_z Pi_N J_j=D_z Pi_N W_j=D_z Pi_N Z_j=0 at z=0, j<=12. (17)
```

The respective targets are `H^(m-1/2)` and `H^(m-3/2)` on the
boundary, valid for `m>1/2` and `m>3/2`. They are zero because the
sections vanish on a collar. In fact every further spatially
differentiated boundary expression permitted by `m` is zero there.
This does not define an ordinary boundary trace below the trace
threshold. At the top commuted energy level any needed Neumann row
is the separate weak variational/conormal identity from the energy
proof, not an ordinary normal trace of J13. Differentiating the
field-dependent physical shell equations is a separate chart-chain
calculation; (17) concerns the fixed normalized primal four/twelve
rows only.

## 8. Parameter endpoints and admission

The source admission is strict for the actual pair, for a reason
that can be checked without enlarging a threshold. The one-copy
calculation `UNIFORM_CAUCHY_DATA`, (18)--(19), gives the first-pair
derivative at most `10^4128`. Its specified identical two-copy passage
allows at most two for that norm. Therefore in either retained case

```
||pair'(epsilon)||H14xH13 <=2*10^4128<10^4200,
||pair(epsilon)||H14xH13 <=2*10^4128 |epsilon|<Dbar
 for epsilon in I, Dbar=10^4200 a_iter.                         (18)
```

The zero background pair is used in (18). Compactness and the already
open C2 initial-pair curve give an open parameter neighborhood O of I
on which `d(epsilon)<Dbar`. Since `a_iter<a_init`, choose a sufficiently
small unspecified positive eta such that

```
J=(-a_iter-eta,a_iter+eta),
closure(J) subset O intersect (-a_init,a_init).                  (19)
```

The initial curve, all its support/symmetry properties, and its retained
derivative bounds therefore hold on a neighborhood of the compact
closure of J. We may and do take O also inside the finite-jet C2 domain.
In particular the maxima and modulus of Section 6 can be taken on
`closure(J)` and are finite/uniform. The retained generic
reduced finite-data construction in `NONLINEAR_LIMIT_AND_UNIQUENESS`,
Section 1 and Sections 3--7, applies to supported Gamma-invariant
pairs with `d<=Dbar`, and its top-continuity strengthening uses the
same solution class. Thus it supplies the actual normalized solutions
on all of J, with the same common majorant
`M=2 C0 Dbar exp(2HT)<R/4` and the same full slab. Its uniqueness
identifies the restriction to I with the already retained solution
family. The parameter-derivative argument can consequently be carried
out on this open domain. Ordinary two-sided derivatives at each
endpoint of I then mean the derivatives of this extension; they are
also the continuous limits of derivatives from the interior of I.
All uniform estimates claimed only on I retain their original scope.
This argument uses the actual retained sufficient admission `d<=Dbar`;
it does not apply that theorem to data beyond its threshold, does
not declare a factor-two parameter interval, and does not create a
new numerical amplitude ceiling. The `delta<=2Dbar` margin in the
smooth-approximation construction remains its own source-justified
approximation allowance. No physical conclusion for a new numerical
parameter interval is asserted: (19) is a qualitative normalized
evolution neighborhood used to give the precise endpoint meaning.
One could equivalently restrict all conclusions to the interior of I
with continuous one-sided endpoint extensions, using (16) on I alone.

The previous full physical closure remains completed under its
retained hypotheses. This initial-jet proof adds no X15 estimate,
new time slab, new amplitude, physical normalization, observer map,
or third parameter derivative.
