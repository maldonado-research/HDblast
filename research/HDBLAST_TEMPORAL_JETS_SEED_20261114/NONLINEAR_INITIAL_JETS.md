# Exact nonlinear initial temporal jets on the retained past slab

14 November 2026. Conditional mathematical checkpoint, retaining the same smooth-background, actual coefficient-chart, and selected canonical initial-family premises as the authenticated source capsules. The result is an evaluated initial-jet map with `I14<=10^6503000 d(a)` on the explicit smaller interval below. This concerns initial temporal jets, not existence of an evolved nonlinear solution. All previous source bytes remain unchanged.

## 1. Norms and retained initial state

Let `ell=10^-6`, `Sigma=[0,6ell] x S3`, and `T=13ell/8`. Use the fixed positive product bundle

```
E=R_v + Sym^2(R dt + T*S3)_H + R_w + (R dt + T*S3)_D
```

with full Frobenius norm in its symmetric block, the fixed product connection, and measure `dz dmu_kappa`. Spatial Sobolev norms are the square root of the sum of squared covariant derivative norms. The symbols `J_j` denote the formal jets forced by the PDE at the initial slice; if a sufficiently regular solution exists, they are its temporal derivatives there. No evolution is constructed in this note. The time hierarchy is the literal sum

```
I14(J)=sum_(j=0)^14 ||J_j||H^(14-j).
```

The construction acts on the full Gamma-invariant space, `Gamma=(Z2)^4`; it imposes no finite angular cutoff or single-isotypic restriction. The identical two-copy norm is charged as in the accepted initial interface.

The actual selected initial curve, with its original lapse `alpha_e=rho(1+e u_profile Y)`, zero shift, and fixed canonical density `rho0`, has its full harmonic completion in `H14 x H13`. This includes the mixed gauge velocities; no physical momentum is rescaled or discarded here. The source normalization is the retained half-action Einstein–scalar system. On `|e|<=a_init=10^-6402`, the completed pair has norm at most `10^-2200` and first parameter derivative at most `10^4200`. At `e=0` it vanishes. Hence for `a<=a_init`,

```
d(a):=sup_|e|<=a (||V_e(0)||H14+||V_e,t(0)||H13)
     <=min(10^-2200,10^4200 a).                                   (1)
```

To represent the same original past state, set `U(s)=J_E V(-s)`. The parallel isometry `J_E` reverses the four odd time slots `H_tA,D_t` and commutes with the spatial product connection and the boundary split. The initial jets are

```
J_0=J_E V(0),             J_1=-J_E V_t(0).                        (2)
```

Their sum norm is unchanged. The coefficients below are the exactly reflected coefficients. Their absolute derivative certificates are unchanged. Equation (2) is not a replacement by an unrelated forward state. The support of both initial fields is contained in the retained annulus, with fixed open zero collars at both ends of Sigma.

## 2. Exact all-state equation

The pinned coordinate equation is `q^{mu nu} partial_mu partial_nu V=f(x,V,P)`, with `f=c+bP+d(P,P)` and `P=partial V`. Write `D_mu=partial_mu+Omega_mu`, and let `gamma` denote the fixed product base connection. In bundle first jets `Q=DV`, the exact source is

```
Fbar(V,Q)=f(V,Q-Omega V)
 +q^{mu nu}[Omega_mu Q_nu+Omega_nu Q_mu-gamma^k_mu_nu Q_k
       +(partial_mu Omega_nu-Omega_nu Omega_mu)V].                (3)
```

The source in (3) contains the true bulk reference Christoffel/Ricci expressions. It includes `-b Omega V` and the quadratic substitution terms. Its constant, linear and quadratic first-jet coefficients are denoted `cbar,bbar,dbar`. The equality `cbar(x,0)=0` holds identically in the base variables because the actual background solves the background-subtracted reduced equation. The product Hessian includes its derivative-index connection. No connection term is counted twice.

After reflection and writing `q^{ss}=-a`, `q^{si}=b^i`, `q^{ij}=c^{ij}` for spatial indices, the exact equation is

```
a(x,U) U_ss=2b^i(x,U)D_iU_s+c^{ij}(x,U)H_ijU-Fbar(x,U,DU),       (4)
```

where `H_ij` is the spatial product Hessian. In particular the sign of the lower source on the right is negative. The reflected geometric source itself has no extra overall minus; the minus in (4) comes from isolating the negative time-time principal coefficient.

Normalizing once by positive `a`, write

```
U_ss=A^i D_iU_s+B^ij H_ijU+C+E^mu D_muU
                         +G^mu nu(D_muU,D_nuU),                 (5)
A=2b/a, B=c/a, C=-cbar/a, E=-bbar/a, G=-dbar/a.
```

Here every coefficient is evaluated at `(s,z,omega,U)`, and `C(x,0)=0`. This normalization is pointwise in state; it does not replace `a` by its background value. Mixed spatial-temporal principal terms and all first-jet time components remain present. The equation is scalar-principal on all sixteen fiber components.

## 3. Triangular nonlinear recursion

For an independent unnormalized check, let `a_k,b_k,c_k` be the total `s` derivatives of the coefficient functions at zero, evaluated on the already constructed jets, and let `F_r=partial_s^r[-Fbar(x,U,DU)]|0`. Differentiating (4) gives

```
a_0 J_(r+2)=2 sum_(k=0)^r binom(r,k)b_k^i D_iJ_(r-k+1)
 +sum_(k=0)^r binom(r,k)c_k^ij H_ijJ_(r-k)
 +F_r-sum_(k=1)^r binom(r,k)a_kJ_(r-k+2),
0<=r<=12.                                                       (6)
```

Every total coefficient derivative includes explicit base-time dependence. For `k>=1`, the acceleration term uses at most `J_(r+1)`. The first-jet source also uses at most `J_(r+1)`. Thus (6) is triangular. Equivalently define `J_(r+2)` by taking `r` formal time derivatives of the right side of (5). This normalized recursion has exactly the same jets and avoids estimating the already isolated acceleration terms twice.

The spatial product connection is time independent, so temporal differentiation commutes with `D_i` and `H_ij`. If moving angular frames are used to evaluate coefficients, their exact time derivatives must first be included in the fixed-product coefficient conversion. Setting those time derivatives to zero is not part of this recursion.

The required target for the `r`th row is `H^(12-r)`. A term with a spatial Hessian has intrinsic derivative order at most `r+(12-r)+2=14`; a linear first-jet term has order at most thirteen, and the quadratic first-jet term at most fourteen. Coefficient chain-rule factors use derivatives of earlier jets and no hidden future jet. The exact zero-state cancellation removes any source term independent of the perturbation.

## 4. Numerical coefficient and finite product ledger

The independently accepted companion [MIXED_INITIAL_JET_PRODUCTS.md](MIXED_INITIAL_JET_PRODUCTS.md), accepted original-proof SHA256 `777ece939f4e34b15f7537585135e1b832abf0b682a8252804aa23b9990b6d9c`, proves that the five normalized maps in (5) have every ordered horizontal/vertical coefficient derivative through total order sixteen bounded as follows, uniformly on the retained state patch. The norm is full Frobenius on every output and base/direction derivative slot and multilinear Hilbert operator norm on state inputs:

```
Cstar=10^500000.                                                 (7)
```

The horizontal derivative is the induced coefficient-bundle derivative, including `-K_V(Omega_i V)` as well as the connection actions on output, direction, and state-input slots. It is not just the background `V=0` slot calculation. The companion's finite numerical ledger is:

* The principal inverse-metric Taylor aggregate through sixteen is `<10^17000`; ordinary derivatives are `<10^18000`.
* Connection derivatives through seventeen have full entry sum `<10^60`; the exact first-jet substitution gives `cbar,bbar,dbar` ordinary derivatives `<10^150200`.
* The scalar reciprocal recurrence uses `a^-1<=18000` at every state and gives ordinary reciprocal derivatives `<10^300000`, rather than only a zeroth reciprocal bound.
* Normalizing the five groups gives ordinary coefficient derivatives `<10^450225`.
* The twenty-one-variable Taylor aggregate starts at `10^450237`; each of sixteen horizontal/vertical conversion steps costs at most `10^80`, including `-K_V(Omega_i V)` and at most forty elementary tensor slots. The final full-array frame conversion costs at most `10^30`, yielding `<10^451547<10^452000<Cstar`.

These entry-sum estimates dominate the complete output/base Frobenius arrays. They leave no direction-array factor to be supplied by a Sobolev embedding of a high-rank coefficient tensor.

The product estimate and numerical induction now follow explicitly from (7). Let

```
S_n=sum_(j=0)^n ||J_j||H^(14-j),       S_n<=10^-12.                (8)
```

The retained global product-cylinder embedding is `||T||infty<=10^7 ||T||H3` for tensor sections with at most sixteen slots. It needs no boundary condition. In particular the full Frobenius-to-star comparison gives

```
||J_0||star,infty <=2*10^7 S_n <=2*10^-5 <10^-3.                 (9)
```

Thus every segment `theta J_0`, `0<=theta<=1`, stays in the actual coefficient domain. The exact cancellation is used in the form

```
C(x,U)=integral_0^1 C_U(x,theta U)[U] dtheta.                     (10)
```

This is an identity of bundle maps before differentiation. Hence its mixed base/time derivatives also contain at least one field factor; it does not leave a hidden nonzero background term. The coefficient in (10) requires at most thirteen coefficient derivatives when the row is differentiated at most twelve times, within (7).

To estimate the `r`th row in `H^(12-r)`, apply `r` time derivatives and `h` ordered spatial covariant derivatives, `0<=h<=12-r`, to (5), using (10). The five groups start with at most two field factors. After `k` differentiations there are at most `k+2` field factors. Differentiating the coefficient gives at most two alternatives: its horizontal derivative, or its vertical derivative applied to a field derivative. Differentiating a field factor gives at most `k+2` alternatives. Thus the number of Leibniz/chain-rule terms is bounded by

```
Nterm <=5 product_(k=0)^11(k+4)=5*15!/3! <10^13.                 (11)
```

This counts multiplicities rather than suppressing binomial factors. Fewer derivatives only reduce the bound. Index rearrangements are tensor permutations of norm one; no covariant derivatives are commuted. Direction/fiber contractions are included in the full coefficient operator norm (7). The horizontal coefficient convention is what makes this chain rule valid for nonzero states. Each `theta` factor in (10) has absolute value at most one.

Every term has at least one field factor, possibly of order zero, with total intrinsic field derivative order at most fourteen. Choose a factor of maximal intrinsic order and put it in `L2`; its norm is at most `S_(r+1)`. Every other factor has intrinsic order at most seven, since two factors of order at least eight would exceed fourteen. Such a factor has the form `D^b J_j`, `b+j<=7`, and has at least three spatial derivatives available within `H^(14-j)`. Its tensor rank together with the three embedding derivatives is at most `2+7+3=12`, within the sixteen-slot certificate. Consequently its `L-infinity` norm is at most `10^7 S_(r+1)<=10^-5`.

There are at most fourteen field factors, but every extra factor costs at most one under (8). Equations (7)–(11) give, for every indicated `h`,

```
||D^h J_(r+2)||L2 <=Nterm Cstar S_(r+1).
```

The square-sum defining the target Sobolev norm has at most thirteen terms. Therefore the deliberately larger evaluated row constant

```
||J_(r+2)||H^(12-r) <=M S_(r+1),
M=10^500200,                     0<=r<=12,                       (12)
```

follows from `sqrt(13)*10^13*10^500000 <10^500200`. The reserve in (12) is separate from the independently established coefficient conversion (7). The low/high argument applies unchanged to the disjoint two-copy sum: each copy's pointwise low norm is bounded by the global `S_n`, and the high `L2` terms sum in the retained norm. No finite-mode restriction appears.

## 5. Evaluated recursion and initial interval

Assume the initial pair has norm `d<=10^-6503020`. Set

```
L=10^500201,       C_J=10^6503000,
a_jets=10^-6507220.                                               (13)
```

Since `1+M<L`, (12) implies `S_(r+2)<=L S_(r+1)` whenever the preceding hierarchy satisfies (8). Starting from `S_1=d`, induction gives

```
S_n<=L^(n-1)d,                   1<=n<=14,
L^13=10^6502613 < C_J.                                          (14)
```

The induction is valid without circularity: every candidate right side in (14) is at most `C_J d<=10^-20`, which is strictly smaller than the assumed threshold `10^-12`. Thus all thirteen nonlinear rows are defined in the actual coefficient domain and obey the finite product estimates. The conclusion is

```
I14=S_14<=C_J d.                                                 (15)
```

This is a restriction of the original amplitude `e`, not a reparameterization or a replacement initial family. For the actual selected canonical family, (1) and (13) imply

```
a_jets<a_init,
10^4200 a_jets=10^-6503020,
sup_|e|<=a_jets I14(e)<=10^-20.                                  (16)
```

Every still smaller interval satisfies the sharper functional bound `sup I14(e)<=C_J d(a)`. The separate compatible-seed proof may require a further restriction using its evaluated lift norm. Equation (16) alone asserts no nonlinear solution on `[-T,0]`; it quantifies only the exact initial temporal hierarchy. Its enormous exponent is a conservative finite consequence of the raw coefficient envelope, not an optimal radius or a measured physical scale.

## 6. Corner, symmetry, and finite regularity completion

Every row of (5) is local. Because the initial fields vanish on fixed open end collars and `C(x,0)=0` identically, (6) gives `J_j=0` on these collars at every constructed order. Therefore all boundary corner traces that exist at the relevant finite Sobolev order vanish. No trace of a merely L2 top jet is claimed. Smooth approximants preserving the zero collars have the same conclusion at all smooth orders. The approximants are analytic tools; the canonical constraints are imposed on the original limiting datum, not newly assumed for arbitrary approximants.

All operations are the actual natural tensor operations of the retained equations and preserve Gamma invariance. Nonlinear products may create higher angular modes; the estimates and the recursion must allow them. A future frozen Picard successor with these same predecessor jets preserves them through order fourteen by (6). The separate regularized seed need only match the jets through order fourteen. Its fixed all-time zero endpoint collars, rather than any claim that it matches all higher nonlinear interior jets, provide smooth all-order corner compatibility. A later order-fifteen smoothing argument would require the corresponding additional jet/lift interface.

The present result gives only the initial hierarchy. A compatible regularized seed is a separate construction, and neither statement supplies the actual quasilinear evolution on the prescribed past slab or any proper-observer tube.

The map from the initial pair to `(J_0,...,J_14)` is also `C^2` on the relative open ball `d<2*10^-6503020` in the fixed collar-supported, Gamma-invariant `H14 x H13` pair space. This auxiliary analytic pair space does not assert that every pair solves the canonical constraints. On this slightly larger ball the same induction gives `S_14<2*10^-20<10^-12`, so the physical coefficient patch remains strict. Two Frechet variations replace field factors or add order-zero state-direction factors; the intrinsic derivative total remains at most fourteen. The row's twelve base/time differentiations, two variations, and optional mean-value extraction require at most fifteen coefficient derivatives, within the companion's sixteen-derivative bound. The same finite low/high argument applies to the directions without assuming that their norms are small, and one additional coefficient derivative controls continuity of the resulting second differential. Finite composition of the thirteen graded rows proves this `C^2` assertion.

Consequently the actual `C^2` initial parameter curve gives a `C^2` finite jet curve on the admitted interval. No numerical bound for its jet Hessian is claimed here. A separate bounded linear lift preserves this initial-curve regularity; neither fact asserts parameter regularity for an evolved PDE solution.

## 7. Source ledger

All retained paths below are relative to the repository `research/` directory, at the registered base `2592e90650cb440cab2ce74f726e3dcf9b9d1247`. They are read-only inputs; hashes pin the bytes, not an independent proof of every inherited premise.

| Retained source | SHA256 |
|---|---|
| `HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md` | `785b20a06ecac643d635a89f70a800a765d20164d22da7e2df3c777a31b3bfa9` |
| `HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md` | `7b28775e6c44b8ef1213e4789c4efb3c6b82727226cd18c5953c0d20038f9f39` |
| `HDBLAST_ORDERED_EVOLUTION_20261106/ORDERED_BACKGROUND.md` | `0b85d3f5dcbdd2b04c8732687ce7691f5649e250fda8a27cec209fafd741dd62` |
| `HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md` | `e78d8b4ece66335aabe2704b4eb7cec599a3bdb5686e6704d97e006d29d29a4c` |
| `HDBLAST_PROPER_SECOND_GEOMETRY_20261113/REVERSE_TIME.md` | `9bc3c93c9e35cea108d06726dd7d837cdf2e3f06c0da7733d2d0d10819fb8764` |
| `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md` | `d6f80e6158226b253c56f2f60c8d9c6cfa2ee2940e79bde42dd8fe8e48a5d6c4` |
| `HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md` (global tensor embedding, section 4) | `26871e9ed63a2345737b2f55fa855d0144ef8ebc85f985573e783c1dee923062` |

The previous accepted private `NEXT_NONLINEAR_EVOLUTION_PLAN.md` supplied the checkpoint question, not a numerical hypothesis. Its hash is `6c7156d67e920a15d31d7b5db8a292b8c07bf081441b969018e11cc593fb9a2f`. The new all-state coefficient companion is `MIXED_INITIAL_JET_PRODUCTS.md`, accepted original-proof SHA256 `777ece939f4e34b15f7537585135e1b832abf0b682a8252804aa23b9990b6d9c`.
