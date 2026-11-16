# Evaluated low-energy Picard convergence on the full retained slab

14 November 2026. Conditional mathematical checkpoint. Retain the authenticated smooth background, actual normalized Einstein–scalar coefficient map, and compatible smooth frozen construction used by the admitted Picard scheme. This note proves convergence of those iterates in the low Cauchy-pair topology, with an explicit factorial tail. Identifying the nonlinear PDE limit and its higher mixed regularity is a separate interface. No previous source bytes are changed.

## 1. Fixed norms and the admitted sequence

Let `ell=10^-6`, `Sigma=[0,6ell] x S3`, and `0<=s<=T=13ell/8`. The fixed positive product bundle is

```
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D.
```

Use full Frobenius norm in the symmetric block, every ordered covariant derivative slot, and measure `dz dmu_kappa`. Let `Dsp` denote the spatial product derivative and `Hsp=Dsp^2` the full spatial covariant Hessian. The full first jet has five direction slots:

```
Q_u=(u_s,Dsp u).
```

Define

```
X_m(u;s)=sum_(j=0)^m ||partial_s^j u(s)||H^(m-j),
X1(h)=||h||H1+||h_s||L2,
Z1(h)^2=||h||L2^2+||Dsp h||L2^2+||h_s||L2^2.
```

Thus `Z1<=X1<=sqrt(2) Z1`. All norms below apply both to one copy and to the retained identical two-copy Hilbert union. Each time level includes all angular modes.

The accepted admission result supplies, on the unchanged original parameter interval

```
|epsilon|<=a_iter=10^(-10^500209),
R=10^-12,              sup_s X14(U^n;s)<R/4,   n>=0,               (1)
```

a sequence whose first element `U^0` is the compatible auxiliary seed and whose successors are smooth frozen solutions. It also covers the admitted supported smooth approximants of the finite initial Cauchy pair. The common initial pair, homogeneous mixed boundary rows, and nonlinear initial jets through fourteen are retained. The same-state past convention is `U(s)=J_E V(-s)`, with the tensor time-slot reflection and opposite elapsed initial velocity; no new parameterization or future-state replacement is made.

The actual normalized operator and source are

```
P_v u=u_ss-A^i(x,v)D_i u_s-B^ij(x,v)H_ij u,
F(x,V,Q)=C(x,V)+E^mu(x,V)Q_mu+G^mu nu(x,V)(Q_mu,Q_nu),
P_(U^n) U^(n+1)=F(x,U^n,Q_(U^n)).                                 (2)
```

Here `P_v` denotes the differential operator and `Q_v=(v_s,Dsp v)` denotes the full first-jet array. The frozen field is `v=U^n` and the successor is `u=U^(n+1)`.

The coefficient definitions are `a=-q^ss`, `A=2q^si/a`, `B=q^ij/a`, `C=-cbar/a`, `E=-bbar/a`, `G=-dbar/a`. They contain the complete physical reference terms and the correct substitution from coordinate first jets to bundle first jets. The bilinear `G` is symmetric in its combined first-jet inputs, as in the retained coefficient definition. No scalar conformal truncation or omission of the time first-jet slot is used.

## 2. Exact difference source and its evaluated Lipschitz constant

Abbreviate the source at the same base point by `F(v,Q_v)`. Set

```
W^n=U^(n+1)-U^n.
```

For `n>=1`, put `v=U^n`, `w=U^(n-1)`, and `h=v-w=W^(n-1)`. Subtract the two actual frozen equations. The exact identity is

```
P_v W^n = F(v,Q_v)-F(w,Q_w)
 +(A(v)-A(w))^i D_i v_s
 +(B(v)-B(w))^ij H_ij v =: r_n.                                  (3)
```

Both principal-difference signs are positive. The differentiated field in these terms is the preceding frozen solution `v`, not the new successor or the older predecessor `w`. Every `W^n` has zero initial Cauchy pair and homogeneous mixed boundary rows. Equation (3) starts at `n=1`; it is not asserted for the auxiliary seed `U^0`.

Use the accepted all-state coefficient constant and tensor embedding

```
Cstar=10^500000,            E3=10^7,
p=E3 R=10^-5.                                                       (4)
```

We use the weaker upper bound `X14(v),X14(w)<=R`, so (1) has reserve. The full five-slot first jet obeys

```
||Q_v||infty
 <=E3 (||v_s||H13+||v||H14)<=p,
||Dsp v_s||infty<=E3||v_s||H4<=p,
||Hsp v||infty<=E3||v||H5<=p.                                    (5)
```

These are complete tensor arrays, with ranks plus the embedding reserve at most seven. No high-rank coefficient map is Sobolev-embedded. The same bounds hold for `w`. The segment `v_theta=w+theta h`, `Q_theta=Q_w+theta(Q_v-Q_w)` lies in the coefficient patch: its star state norm is at most `2E3R=2*10^-5<10^-3`, and its first jet has norm at most `p`.

The exact derivatives along this joint segment are

```
F_V[h]=C_V[h]+E_V[h]Q+G_V[h](Q,Q),
F_Q[Qdot]=E Qdot+G(Qdot,Q)+G(Q,Qdot)=E Qdot+2G(Q,Qdot).                            (6)
```

In particular both quadratic cross terms are retained. The full-array operator certificate gives

```
||F_V||op<=Cstar(1+p+p^2),
||F_Q||op<=Cstar(1+2p).
```

Integrating (6) over `0<=theta<=1`, taking L2, and using the fixed linear product derivative gives

```
||F(v,Q_v)-F(w,Q_w)||2
 <=Cstar[(1+p+p^2)||h||2+(1+2p)||Q_h||2].                         (7)
```

The principal coefficients depend only on the state. Their pointwise difference bounds are `|A(v)-A(w)|<=Cstar|h|` and `|B(v)-B(w)|<=Cstar|h|` in the full direction arrays. Contracting against (5) and adding the two terms in (3) costs at most `2p Cstar||h||2`; each contraction has norm one in these complete arrays. Hence

```
||r_n||2<=Cstar[(1+3p+p^2)||h||2+(1+2p)||Q_h||2].                 (8)
```

Since `p<=1/100`, both displayed coefficients are at most `11/10`. Cauchy–Schwarz gives the stronger intermediate bound

```
||r_n||2 <=(11 sqrt(2)/10) Cstar Z1(h)<=2 Cstar Z1(h).
```

We record the convenient evaluated interface

```
L_delta=4 Cstar=4*10^500000,
||r_n||L2<=L_delta X1(W^(n-1)),            n>=1.                  (9)
```

Every homogeneous inequality includes the zero difference. No second derivative of `h` occurs. The same source bound holds if the two principal coefficient differences multiply `D_i w_s,H_ij w` instead, since that field has the same high-norm bound; this is the form needed for comparison of two actual solutions. Norms and constants in (5)–(9) are unchanged on the two-copy union: use global L2 for the differing field and the maximum pointwise bound supplied by its global Sobolev norm.

## 3. Varying energy metrics and the common Volterra kernel

For each predecessor `v`, the accepted normalized energy is

```
mathcal E_v[h]=1/2 integral (|h_s|^2+B(v)^ij<D_i h,D_j h>+|h|^2),
Z1(h)^2/4<=mathcal E_v[h]<=Z1(h)^2.                               (10)
```

The shell coefficient traces `A^z=B^zA=0`, the preserved orthogonal Dirichlet/Neumann split, and the outer all-Dirichlet rows make its exact flux zero. For a smooth solution `P_v h=r`, the already evaluated energy calculation gives

```
d sqrt(mathcal E_v[h])/ds <=||r||2+G_en sqrt(mathcal E_v[h]),
G_en=16 Cstar.                                                    (11)
```

Use the regularized root at zero energy. The bilinear coefficient map `G` in (6) is distinct from the scalar growth constant `G_en`. From (10),

```
sqrt(mathcal E_v[h])<=X1(h),
X1(h)<=2 sqrt(2) sqrt(mathcal E_v[h])<=3 sqrt(mathcal E_v[h]).       (12)
```

The metric in the energy of `W^n` is the one built from `U^n`; it may vary with `n`. We never identify two such energies. The common norm `X1` and the uniform comparisons (10)–(12) link successive steps.

For `n>=1`, zero initial energy, (9), and (11) imply

```
X1(W^n;s)
 <=beta integral_0^s exp(G_en(s-r)) X1(W^(n-1);r) dr,
beta=3L_delta=12 Cstar.                                          (13)
```

There is no assertion `beta*T<1`. The accepted interval and amplitude are not shrunk to obtain a contraction.

For the starting difference, the admission estimate alone gives

```
sup_s X1(W^0;s)
 <=sup_s[X14(U^1;s)+X14(U^0;s)]<=R/2.                             (14)
```

No PDE equation or PDE residual bound is imposed on the seed in this step.

## 4. Evaluated factorial estimate and finite tails

Set `b_n(s)=exp(-G_en s)X1(W^n;s)`. Equations (13)–(14) give

```
b_0(s)<=R/2,
b_n(s)<=beta integral_0^s b_(n-1)(r)dr,           n>=1.
```

Induction by integration over the ordered time simplex proves

```
X1(W^n;s)<=(R/2) exp(G_en s)(beta s)^n/n!,       n>=0.             (15)
```

For integers `M>N>=0`, telescope and use the literal pair norm:

```
X1(U^M-U^N;s)
 <=(R/2)exp(G_en s) sum_(n=N)^(M-1)(beta s)^n/n!.                  (16)
```

For every `z>=0`, `(N+k)!>=N! k!` implies

```
sum_(n=N)^infinity z^n/n! <= (z^N/N!)exp(z).
```

Consequently a uniform finite tail bound is

```
sup_(0<=s<=T) X1(U^M-U^N;s) <=min(R/2,Tau_N),
Tau_N=(R/2)exp((G_en+beta)T)(beta T)^N/N!.                          (17)
```

The first entry of the minimum follows independently from admission. With the actual `T`, the constants in this expression are exactly

```
G_en T=26*10^499994,
beta T=(39/2)*10^499994,
(G_en+beta)T=(91/2)*10^499994.                                    (18)
```

These are finite numbers; no decimal expansion or floating-point evaluation of them is required. For fixed `beta T`, the factorial tail tends to zero. An explicit further tail certificate is available. Put

```
N_star=10^499997.
```

For every integer `N>=N_star`, write `z=beta T` and `g=G_en T`. Then `N>=2z`; the ratio of successive tail terms is at most one half, so (16) is at most `R exp(g) z^N/N!`. The elementary integral bound for `log(N!)` gives `N!>=(N/e)^N`. Using `e<3` and `log(2)>=1/2`, the needed finite comparisons reduce to

```
2z/N_star=39/1000<1,
3z/N_star=117/2000<1/2,
4g/N_star=104/1000<1.
```

Thus `z^N/N!<=2^-N` and `exp(g)<=2^(N/2)`, yielding

```
sup_s X1(U^M-U^N;s)<=R*2^(-N/2),       M>N>=N_star.                (19)
```

This enormous sufficient index is a mathematical tail certificate, not a practical iteration count or a physical scale. Equations (15)–(19) are uniform in the original parameter, the admitted supported smooth approximation, and angular content.

## 5. The low-topology limit and the stability handoff

The Banach space `C([0,T];H1(E)) x C([0,T];L2(E))` is complete. Equations (17)–(19) therefore yield a unique pair limit `(U,Z)` of `(U^n,partial_s U^n)` in this space. For every iterate,

```
U^n(s)=U^n(0)+integral_0^s partial_s U^n(r)dr
```

holds in L2. Passing to the uniform limits shows `U(s)=U(0)+integral_0^s Z(r)dr`; hence `Z=partial_s U` and

```
U in C([0,T];H1) intersect C^1([0,T];L2).
```

The same initial pair is retained. Letting `M` tend to infinity in (17) or (19) gives the same explicit tails for `X1(U-U^N)`. In particular the low limit is not merely a subsequential weak limit. This conclusion does not yet identify the nonlinear PDE or the Neumann trace from H1 convergence alone; the separate higher-mixed recovery supplies those facts.

For the companion uniqueness/stability calculation, here is the precise evaluated conditional input. If a smooth difference `h` obeys `P_v h=r`, the same homogeneous boundary rows, and `||r||2<=L_delta X1(h)`, then (11)–(12) give

```
X1(h;s)<=3 exp((G_en+beta)s)
                    [||h(0)||H1+||h_s(0)||L2],
G_en+beta=28 Cstar.                                               (20)
```

For two actual solutions the exact source identity uses second derivatives of the second solution, as stated after (9). Extending (20) to the later finite-regularity solution class is justified in the companion argument once that class and its energy identity are established. No existence or uniqueness assertion for a nonlinear finite-Sobolev solution is inferred solely from (20).

The present result makes no continuous-full-X14 assertion, no evolved C2 parameter assertion, and no observer or finite nonlinear observable claim. It uses the same admitted original amplitude and full prescribed past slab throughout.

## 6. Source pins and proof boundary

The source tree is the repository, registered base `9b74d57b6caef415775abc069838a4bcea303831`. The following paths are relative to its `research/` directory. Their byte hashes identify read-only inputs while retaining their analytic premises.

| Input | SHA256 |
|---|---|
| `HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/CONDITIONAL_PICARD_ADMISSION.md` | `9c6332825e6a0a503f300b5e0673b991d21a443409575efd506917b1646e22fd` |
| `HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/FROZEN_ENERGY_COMMUTATORS.md` | `32bae0deb33acda059e4acb55cfa7d511b8699cb5f59affd383d8e463aab746f` |
| `HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/NORMALIZED_METRIC_MARGINS.md` | `07f637fc7e32cc4ff3ebcef398687df43752e670825a5c4da85941509dc537a2` |
| `HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md` | `37ea4bbe555eba7aaae1140d778d7e2d6e86ee20f7a489e5e35825750411fcc9` |
| `HDBLAST_TEMPORAL_JETS_SEED_20261114/QUANTITATIVE_COMPATIBLE_SEED.md` | `94f0f75709c8f0ee385067971f431829c6881fa3e6e3c506a05cf6d9511fe3d3` |
| `HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md`, global tensor embedding | `26871e9ed63a2345737b2f55fa855d0144ef8ebc85f985573e783c1dee923062` |

The previously frozen `NEXT_DIFFERENCE_AND_SMOOTHING_PLAN.md` selected this calculation; its prospective constants are not assumptions. The pointwise segment argument above supplies the new finite Lipschitz constant. The already evaluated energy supplies the growth coefficient, and ordinary integration supplies the factorial bound. No qualitative Sobolev or energy constant has been promoted to a number. The actual finite-Sobolev PDE limit, weak top continuity, and nonlinear uniqueness use the separate lower-mixed and limit companions.
