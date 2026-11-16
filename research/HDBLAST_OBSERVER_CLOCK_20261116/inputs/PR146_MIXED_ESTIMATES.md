# Full mixed difference and linearized stability estimates

16 November 2026. New proof at the unchanged merged base
`fb2ba28453887056166e9f9b87a5a78b9bc22360`. All earlier roots and
capsules are read-only. The results below are qualitative finite estimates
on the already admitted full slab. They add no numerical constant,
amplitude restriction or time restriction.

## 1. Actual equation, norms and admitted inputs

Let `Sigma=[0,6 ell] x S3`, `ell=10^-6`, `T=13 ell/8`, and use the fixed
positive sixteen-field product bundle, full Frobenius symmetric-block
norm, fixed product connection, full spatial measure and all angular
modes. The same proof works for the retained identical two-copy Hilbert
union. Put

```
X_r(h;s) = sum_(j=0)^r ||partial_s^j h(s)||H^(r-j)(Sigma),
I_r(h)   = sum_(j=0)^r ||partial_s^j h(0)||H^(r-j)(Sigma),
S_n(f;s) = sum_(j=0)^n ||partial_s^j f(s)||H^(n-j)(Sigma).
```

`S_n` is a full mixed source norm. It includes time derivatives and all
ordered spatial derivative slots, including every normal placement.
It is not a spatial-only norm. For nonsmooth fields these are
distributional derivatives with the indicated representatives.

The actual normalized equation and linearized source are

```
P_v h = h_ss-A^i(x,v)D_i h_s-B^ij(x,v)Hsp_ij h,
P_v v = F(x,v,Q),                       Q=Dv, all five slots,
F = C+E Q+G(Q,Q),
T_v k = F_V[k]+F_Q[Dk]
        +A_V^i[k]D_i v_s+B_V^ij[k]Hsp_ij v,
mathscrL_v = P_v-T_v.                                      (1)
```

All coefficient derivatives are the actual horizontal/vertical ones,
including the state-connection correction and every direction/input slot.
No coordinate first jet is substituted for `Dv`. The background and
other external parameters are fixed; epsilon changes only the admitted
canonical initial curve.

Inputs are the admitted actual family `v_epsilon in C X14`, its common
bound `sup_epsilon,s X14(v_epsilon;s)<=M<R/4`, `R=10^-12`, and the
unchanged interval `|epsilon|<=10^(-10^500209)`. The full initial jets
and their parameter bounds are supplied by the companion initial-jet
proof. No evolved parameter derivative is assumed in proving the
nonlinear difference estimate.

The companion `DIFFERENTIATED_INITIAL_JETS.md`, Section 8, also joins
the retained finite-data theorem to an unspecified open interval `J`
containing this closed interval: the actual initial pair stays strictly
inside the admitted data ball on `closure(J)`, and the same construction
gives the same common bound `M<R/4` and the same slab there. Whenever
ordinary derivatives at the retained endpoints are used below, use
that actual extension. All estimates here hold uniformly on
`closure(J)` by the same proof, with finite constants on that compact
set. The reported canonical interval and its numerical ceiling remain
unchanged.

The retained coefficient certificate supplies ordered mixed
horizontal/vertical derivatives through sixteen. The retained metric
and boundary facts are

```
(1/2)I <= B <= 2I,   19/20 <= B^zz <= 20/19,
|A| <= 1/19,
A^z = B^zA = 0 at z=0.                                    (2)
```

The fixed chart has four Dirichlet and twelve Neumann primal shell
components, and all Dirichlet components at `z=6 ell`. These rows are
linear, independent of epsilon, and preserved by time differentiation
and the six natural rotations. They are the chart rows, not an assertion
that differentiated physical junction laws have no extra terms.

The read-only proofs used below are the exact retained public files
listed in the source ledger in Section 9. Their bytes were checked
against the immutable repository paths in Section 9 and the snapshots
pinned by `SOURCE_PINS.json` (where `MIXED_INITIAL_PRODUCTS.md` names
the same temporal-jet product file).
Their retained numerical constants are not altered or silently reused
as new parameter constants. We repeat the analytic steps needed at the
present orders and finite regularity.

## 2. Exact nonlinear difference and complete five-jet linearization

Write `v_b-v_a=h`, `v_theta=v_a+theta h`, and `Q_theta=Dv_theta`.
The pointwise state segment remains in the certified convex state patch
by the common admitted bound. With the principal operator taken from
the FIRST solution `v_b`, exact subtraction gives

```
P_(v_b) h = F(v_b,Dv_b)-F(v_a,Dv_a)
 +(A(v_b)-A(v_a))^i D_i(v_a)_s
 +(B(v_b)-B(v_a))^ij Hsp_ij v_a
             = T_(b,a) h,                                  (3)

T_(b,a) k = integral_0^1 {F_V(v_theta,Q_theta)[k]
                         +F_Q(v_theta,Q_theta)[Dk]} dtheta
 + integral_0^1 A_V(v_theta)^i[k] dtheta D_i(v_a)_s
 + integral_0^1 B_V(v_theta)^ij[k] dtheta Hsp_ij v_a.
```

The second-derivative factors in the last two terms are those of the
SECOND solution `v_a`. Replacing them by derivatives of `h` or assigning
both principal operators to `v_a` without changing the source is not
the subtraction identity.

The source derivatives in (3) are exactly

```
F_V[k] = C_V[k]+E_V[k]Q+G_V[k](Q,Q),
F_Q[p] = E p+G(p,Q)+G(Q,p).                                (4)
```

Every `Q,p` has all five directions, including elapsed time. Both
quadratic slots are retained. No symmetry shorthand is needed.

For the stability of first-variation candidates, let `w_a,w_b` solve
`mathscrL_(v_a) w_a=mathscrL_(v_b) w_b=0`, with their differentiated
initial jets and homogeneous fixed rows. Put `k=w_b-w_a`, `d=v_b-v_a`.
Then

```
mathscrL_(v_b) k = R_(b,a)[w_a],                            (5)
R_(b,a)[u] = (A_b-A_a)^i D_i u_s+(B_b-B_a)^ij Hsp_ij u
 +(F_V,b-F_V,a)[u]+(F_Q,b-F_Q,a)[Du]
 +(A_V,b-A_V,a)^i[u]D_i(v_a)_s+A_V,b^i[u]D_i d_s
 +(B_V,b-B_V,a)^ij[u]Hsp_ij v_a+B_V,b^ij[u]Hsp_ij d.
```

The high compared-field second derivatives in the first line are those
of the SECOND field `u=w_a`. The terms with `D_i d_s,Hsp d` are separate
and necessary. For complete source accounting the two `F` differences
in (5) have the exact joint-segment formula

```
(F_V,b-F_V,a)[u]+(F_Q,b-F_Q,a)[Du]
 = integral_0^1 {F_VV,theta[d,u]+F_VQ,theta[u,Dd]
                 +F_QV,theta[Du,d]+F_QQ,theta[Dd,Du]} dtheta.
                                                               (6)
F_VV[d,u] = C_VV[d,u]+E_VV[d,u]Q+G_VV[d,u](Q,Q),
F_VQ[u,p] = E_V[u]p+G_V[u](p,Q)+G_V[u](Q,p),
F_QV[p,d] = E_V[d]p+G_V[d](p,Q)+G_V[d](Q,p),
F_QQ[p,q] = G(p,q)+G(q,p).
```

Likewise `A_V,b-A_V,a=integral A_VV,theta[d,.]` and the same for `B`.
Thus (5) contains all state/first-jet cross terms and principal
coefficient Hessians. It is an identity for actual candidate solutions,
not a preliminary assumption that they are derivatives.

## 3. Product estimates with the actual derivative reserve

All constants denoted by `C` below are finite, uniform over the admitted
parameter interval, and may depend on the common `X14` bound, the fixed
coefficients, geometry, finite word sets and the full time slab. Their
values are not asserted to be numerical certificates.

The following estimates are the interfaces used in the proof:

```
S12(T_(b,a) h) <= C X13(h),
S11(T_(b,a) h) <= C X12(h),                                (7)

S11(T_v k) <= C X12(k),
S10(T_v k) <= C X11(k),                                   (8)

S11(R_(b,a)[u]) <= C X13(d) X13(u),
S10(R_(b,a)[u]) <= C X12(d) X12(u).                        (9)
```

The last line holds almost everywhere when `u` has only a bounded top
norm. In particular the first line of (9) requires the compared field's
`X13`, not its `X12`. No `X15` of any nonlinear solution occurs.

Here is a direct proof, including time derivatives. Expand any ordered
time/spatial row of each displayed source using the covariant chain
rule. Horizontal derivatives act on all coefficient slots and state
inputs, and vertical derivatives introduce the differentiated state
factor. Every expansion has finitely many ordered Leibniz leaves;
repeated leaves retain multiplicity. Coefficient maps before state
contraction use their full-array/operator certificate and are never
Sobolev embedded.

For (7), a leaf is linear in one distinguished `h` or `Dh` factor, with
other field factors from `v_theta,v_a`. At degree `n<=12`, each field
factor has order at most fourteen, the distinguished factor has order
at most thirteen, and the sum of intrinsic field derivative orders is
at most `n+2<=14`. The number of factors is finite. At most one factor
has order at least eight. Put a maximal-order factor in `L2` and all
others in `L-infinity`, using the retained `H3` tensor embedding. If
the high factor is from the base family, the distinguished `h` is a
low factor, so its infinity norm is bounded by `C X13(h)` and the
base high factor by the common `X14` bound. If `h` is high, use its
`L2` norm. Every other low factor has order at most seven, and the
three extra spatial derivatives fit within the available mixed norms.
The largest embedded rank is at most `2+7+3=12<16`. The same reasoning
at `n<=11` uses only `X12(h)`, because the distinguished first-jet
factor has order at most twelve. This proves (7) after the finite
spatial-order and temporal sums. It also proves (8), using respectively
degrees eleven and ten.

For (9), use (6) and the corresponding principal coefficient segment
identities first. Each resulting leaf is bilinear in distinguished
`d,u` factors and may contain base-family factors. At degree `n<=11`
its total intrinsic field order is at most thirteen; the only factors
that can reach order thirteen are `D^2 d,D^2 u` with all row derivatives
on them. At most one factor has order at least seven. Put that factor
in `L2`, and all others in infinity norm with the three-spatial-order
reserve. If a base factor is high, both distinguished factors are low;
if either distinguished factor is high, the other is low. This gives
exactly `C X13(d) X13(u)`. Degree ten gives the lower `X12` product.
In particular the principal term `delta A D^2 u` is estimated with
`u` in the high second-derivative position. This calculation does not
exchange that position with an unknown having less regularity.

The largest required original coefficient order here is thirteen:
one vertical derivative plus twelve row derivatives in (7), or two
vertical derivatives plus eleven row derivatives in (6). The order
sixteen input therefore suffices. No high-rank coefficient tensor is
sent into the rank-sixteen field embedding. The fixed full connection,
all normal placements, time derivatives and both quadratic jet slots
remain in these finite expansions.

## 4. Energy at orders thirteen and twelve

For `r=13` or `r=12`, use all ordered words
`W=partial_s^p R_I` of degree at most `r-1`, with the actual natural
rotation action on every tensor and derivative slot. Put

```
E_v[g] = (1/2) integral (|g_s|^2+B^ij<D_i g,D_j g>+|g|^2),
Y_r(h;s)^2 = sum_(p+|I|<=r-1) E_v[partial_s^p R_I h](s).
```

The retained angular bridge gives `Y_r<=C X_r` and
`Y_r(h;0)<=C I_r(h)`. The energy gradient is the complete spatial
gradient. In particular it includes the first normal derivative of
each commuted word.

Consider

```
P_v h = T h+f,                                           (10)
```

where either `(r,T)=(13,T_(b,a))` with `v=v_b`, or
`(r,T)=(12,T_v)`. For this paragraph define the actual word source norm

```
E_(r-1)(f)^2 = sum_(p+|I|<=r-1)||partial_s^p R_I f||2^2.
```

The full mixed source norm controls it: `E_(r-1)(f)<=C S_(r-1)(f)`.
Naturality gives the exact ordered-subset commutator

```
[P_v,W]h = sum_(nonempty ordered subsets J of W)
 {(W_J A)^i D_i partial_s(W_(Jc)h)
  +(W_J B)^ij Hsp_ij(W_(Jc)h)}.                           (11)
```

No rotation generators or spatial covariant derivatives are reordered.
If `q>=1` derivatives hit the coefficient in a degree `n` word, the
unknown order is at most `n-q+2`. Expanding the composed coefficients
covariantly and using the same high/low allocation as Section 3 gives

```
energy aggregate of (11) <= C X_r(h),
normal-recovery rows of degree <=r-2 <= C X_(r-1)(h).     (12)
```

For the energy rows the total intrinsic order is at most `r+1<=14`
and the unknown order at most `r`; for normal rows these are at most
`r` and `r-1`. This is why there is no unestimated top unknown term
in the coefficient-hit normal rows.

The exact retained energy identity has boundary flux

```
integral_boundary [(A.n)|g_s|^2/2+B^ij n_i<D_j g,g_s>].
```

At the shell (2) reduces it to `-B^zz<g_z,g_s>`. The orthogonal fixed
split makes this zero: the Dirichlet component has `g_D,s=0`, and
the Neumann component has `g_N,z=0`. At the outer boundary all
components have `g_s=0`. Only time/rotation words are used here; no
normal derivative is assigned any additional homogeneous boundary row.
The coefficient production terms in the bulk are bounded by `C E_v[g]`
because `D A,D B,B_s` are uniformly bounded. Thus, using (7), (8), (11)
and (12), and regularizing the square root at zero energy,

```
Y_r(h;s) <= C I_r(h)
          +C integral_0^s [X_r(h;t)+S_(r-1)(f;t)] dt.     (13)
```

This is an integrated energy inequality, not a requirement that a
growth coefficient times `T` be small.

## 5. Full normal recovery and explicit absorption

To show exactly what is recovered, let

```
U_j(h)^2 = sum_(p+|I|<=r-j)
                 ||D_z^j partial_s^p R_I h||2^2,
0<=j<=r.
```

The energy and the inverse natural-rotation bridge give
`U_0+U_1<=C Y_r`. This includes the pure rotation words of length `r`
and the time derivatives of length `r`, obtained respectively from
the gradient/zeroth-order energy and from the time energy of words of
length `r-1`. No energy is taken of a normally differentiated unknown.

For `j>=2`, apply `D_z^(j-2) partial_s^p R_I`, of total degree at most
`r-2`, to the UNISOLATED equation (10). Expand its products first.
Only then divide by the zeroth coefficient `B^zz`, whose reciprocal is
at most `20/19`. The unique desired term is
`B^zz D_z^j partial_s^p R_I h`. All zero-hit other principal terms
have smaller normal count: the `h_ss` and angular Hessian terms have
count `j-2`; the `A^z D_z h_s` and `2 B^zA Hsp_zA h` terms have
count `j-1`; the angular mixed time term has count `j-2`.
Interior `A^z,B^zA` are retained. Only their shell traces vanish.

Coefficient-hit terms are precisely the lower rows in (12), and
the differentiated `T h` terms are bounded by the lower lines of
(7) or (8). Put

```
N_(r-2)(f)^2 = sum_(a+p+|I|<=r-2)
                       ||D_z^a partial_s^p R_I f||2^2,
N_(r-2)(f) <= C S_(r-2)(f).
```

After the finite row square sums there is a finite `c>=1` such that

```
U_j <= c sum_(k<j) U_k
          +c [X_(r-1)(h)+N_(r-2)(f)],        2<=j<=r.    (14)
```

This is a triangular full normal-count system. For an explicit
absorption, choose `tau=1/(4c+2)` and set

```
V = U_0+U_1+sum_(j=2)^r tau^(j-1) U_j.
```

On multiplying (14) by `tau^(j-1)` and summing, all terms with
`k>=2` on the right are bounded by

```
[c tau/(1-tau)] sum_(k=2)^r tau^(k-1) U_k < V/4.
```

Terms with `k=0,1` are bounded by `C(U_0+U_1)`. Absorb `V/4` into
the left side. The unweighted normal array is at most a finite
multiple of `tau^(1-r) V`. Finally the product splitting

```
||D_sp^m q||2^2
 = sum_(a=0)^m binom(m,a)||D_z^a D_S^(m-a) q||2^2
```

and the inverse rotation bridge recover every spatial derivative
placement and the literal temporal sum defining `X_r`. The result is

```
X_r(h;s) <= C Y_r(h;s)+C X_(r-1)(h;s)
                             +C S_(r-2)(f;s).           (15)
```

The absorbed terms were top-order terms of lower normal count, with
weights chosen once from finite coefficient bounds. This step uses no
smallness of `h`, `f`, epsilon or `T`. The remaining lower unknown
norm in (15) is instantaneous; it has not been replaced by a spurious
small-time factor. No derivative of `1/B^zz` is used.

## 6. Finite regularity, compatibility and Green's formula

For the nonlinear comparison, actual `v_a,v_b in C X14` provide a
spare order for `r=13`. A time/rotation word of degree at most twelve
therefore satisfies

```
Wh in C H2,       (Wh)_s in C H1,       (Wh)_ss in C L2.
```

Thus the energy identity and the integration by parts are ordinary
finite-Sobolev identities; no smooth nonlinear flow is assumed.

For `r=12` it is enough that a candidate has

```
w in C X12 intersect L-infinity X13.                     (16)
```

For every word of degree at most eleven this implies

```
Ww in C H1 intersect L-infinity H2,
(Ww)_s in C L2 intersect L-infinity H1,
(Ww)_ss in L-infinity L2.                                (17)
```

Distributional differentiation and the bounds in (17) imply
`Ww in W^(1,infinity) H1` and `(Ww)_s in W^(1,infinity) L2`.
Their continuous representatives are those already fixed by (16).
Consequently each commuted energy is absolutely continuous. The
spatial Green formula is legitimate at almost every time, pairing an
`H2` field with an `H1` time derivative. Equations (7)--(12) put the
commuted forcing in `L1 L2` (indeed bounded here), so its integrated
energy inequality is (13). This argument neither assumes nor proves
strong `C X13` for `w`.

The trace rows used in that formula also follow at this regularity.
The original fixed rows may be differentiated in time and by natural
rotations as spacetime distributions. At the highest word in (17),
the spare `X13` gives `Ww in L-infinity H2`, so its Neumann trace
exists in `L-infinity H^(1/2)`; it equals the distributional derivative
of the original zero trace. The same argument gives the Dirichlet
trace of `(Ww)_s` from its `L-infinity H1` bound. Thus the two zero
factors needed for the flux hold almost everywhere. Normal recovery
is an interior distributional identity with all rows in `L2` and
adds no boundary condition. Initial compatibility is the actual full
temporal-jet compatibility supplied by the initial-jet proof; arbitrary
uncompleted Cauchy-pair data have not been substituted.

Every normal source in (15) is continuous at the lower level in the
applications here. For example `S10(R[u])` uses at most twelve
derivatives of `u`, supplied continuously by (16); `S11(R[u])` may
use its thirteenth derivatives and is required only almost everywhere
in the energy integral. Thus the resulting `X12` estimates hold for
the already continuous representatives at every time, not just almost
everywhere. One can alternatively obtain the same conclusion by
taking regular times approaching each endpoint and using lower-level
continuity.

## 7. Closed full mixed estimate; no small-time assumption

The literal mixed fundamental theorem gives, for these classes,

```
X_(r-1)(h;s) <= I_(r-1)(h)+integral_0^s X_r(h;t) dt.      (18)
```

Apply it separately to each actual temporal jet in its own spatial
space. It is not interpolation between unrelated time jets. Combining
(13), (15), (18), and `I_(r-1)<=I_r`, gives

```
X_r(h;s) <= C [I_r(h)+S_(r-2)(f;s)
                  +integral_0^s S_(r-1)(f;t) dt]
                 +C integral_0^s X_r(h;t) dt.            (19)
```

Hence the closed estimate is

```
sup_(0<=s<=t) X_r(h;s)
 <= C exp(C t) [I_r(h)
       +||S_(r-2)(f)||L-infinity(0,t)
       +||S_(r-1)(f)||L1(0,t)],              0<=t<=T.    (20)
```

The infinity norm may include endpoints through the continuous lower
source representative. All constants are finite on the already fixed
slab. Ordinary Gronwall, including its exponential factor, is used;
no condition `C T<1`, shorter slab or reduced amplitude is introduced.

## 8. Consequences for the actual family and first candidates

For the nonlinear difference (3), use `r=13`, `T=T_(b,a)` and `f=0`.
Then (20) proves

```
||v_b-v_a||C X13 <= K I13(v_b-v_a).                      (21)
```

The completed initial-jet curve is `C2` and its first derivative has a
uniform finite `I13` bound on the compact admitted interval. Therefore

```
||v_b-v_a||C X13 <= K_initial K |b-a|.                   (22)
```

This is actual full mixed Lipschitz stability of the evolved canonical
family. In particular all actual first difference quotients are uniformly
bounded in `C X13`. It uses the existing individual `C X14`, not an
unproved top-data continuity statement or an `X15` comparison solution.

For candidates (16) with a common finite bound
`sup_epsilon ||w_epsilon||L-infinity X13<=B`, apply (20) at `r=12`
to (5), with `f=R_(b,a)[w_a]`. Equations (9), (22) give

```
||w_b-w_a||C X12
 <= K [I12(w_b-w_a)
          +(1+T) B ||v_b-v_a||C X13].                   (23)
```

The full differentiated initial-jet curve has a uniform finite second
derivative, so `I12(w_b-w_a)<=K_J |b-a|`. Consequently

```
||w_b-w_a||C X12 <= K_* |b-a|.                           (24)
```

In particular the actual candidate differences divided by `b-a` have
a uniform `C X12` bound. For each fixed nonzero increment their spare
`L-infinity X13` is finite and justifies (17); no bound for that spare
norm divided by the increment enters `K_*`. The constants depend on
the uniform candidate bound `B`, never on its high difference quotient.

At a fixed base, the same estimate with zero residual proves uniqueness
of first-variation candidates in (16) with the same full initial jets.
Existence and identification of candidates with actual parameter
derivatives are the separate quotient-construction task. This note
supplies its actual mixed estimates and finite-regularity justification;
it does not infer differentiability merely from the formal equation.

The entire argument preserves the exact elapsed same-member tensor
reflection, fixed density/lapse binding, fixed four/twelve shell split,
outer Dirichlet row, common canonical amplitude and original time slab.
It makes no observer, response-normalization or cosmological claim.

## 9. Exact retained source ledger

All repository paths below are relative to the unchanged merged base
`fb2ba28453887056166e9f9b87a5a78b9bc22360`. The corresponding files
were read from the immutable repository base and their bytes
matched the private source copies used in the derivation. These pins
identify the retained inputs, including their stated conditional scope;
the present argument supplies the new comparison and finite-regularity
steps rather than treating the earlier smooth-only frozen proof as an
already established parameter theorem.

| Repository path | Bytes | SHA-256 | Git blob |
| --- | ---: | --- | --- |
| `research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/FROZEN_ENERGY_COMMUTATORS.md` | 16215 | `32bae0deb33acda059e4acb55cfa7d511b8699cb5f59affd383d8e463aab746f` | `e10ef82a83cbc0fa9ed8a932fd9b6a87a765ebf4` |
| `research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/INTERIOR_NORMAL_RECOVERY.md` | 16710 | `7be299e660f98cc3c765175f24c598a2367232828cd294482e0beba3875c4267` | `907414e0b854c7f2265f42d4573e2b047b6f989e` |
| `research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/ROTATION_MIXED_PRODUCTS.md` | 16213 | `0299e88a147a5cb88f1b8a56e2b7091a26c9074e1ec7b060bf0618d202f34863` | `66f2b46c22faec46f02b74e29ac77934edeafd19` |
| `research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/NORMALIZED_METRIC_MARGINS.md` | 3917 | `07f637fc7e32cc4ff3ebcef398687df43752e670825a5c4da85941509dc537a2` | `e9f7c1f8dacf840ba96c52c034029846a0dd6ce9` |
| `research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md` | 18417 | `37ea4bbe555eba7aaae1140d778d7e2d6e86ee20f7a489e5e35825750411fcc9` | `85ac189c097c734f40997367b046ec100a851208` |
