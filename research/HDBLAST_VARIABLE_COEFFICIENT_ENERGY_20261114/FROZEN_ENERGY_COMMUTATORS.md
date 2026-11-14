# Evaluated frozen variable-coefficient energy and commutators

14 November 2026. Conditional mathematical checkpoint. The actual smooth background, reduced coefficient map, and fixed product geometry retain their authenticated source premises. For every smooth compatible predecessor in the stated ball and smooth solution of its frozen equation, this proof evaluates `C_Y=10^30`, `K_u=10^500101`, and `K_v=10^500100` in the integrated high-order energy inequality. It does not prove nonlinear existence, propagation of a Picard ball, or a finite-amplitude evolved family.

## 1. Spaces, hypotheses, and exact operator

Let `ell=10^-6`, `Sigma=[0,6ell] x S3`, `0<=s<=T=13ell/8`, and use the retained positive product bundle

```
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D.
```

The symmetric block has full Frobenius norm. The fixed product connection is trivial in `s,z`; its spatial Sobolev norm includes every ordered covariant derivative slot and measure `dz dmu_kappa`. Write

```
X_m(w;s)=sum_(j=0)^m ||partial_s^j w(s)||H^(m-j).
```

The statement applies to one copy and the specified identical two-copy Hilbert union with exactly these norms. No finite angular cutoff is used.

Assume a smooth predecessor `v` satisfies the homogeneous four-Dirichlet/twelve-Neumann shell rows, the outer all-Dirichlet rows, and

```
sup_(0<=s<=T) X14(v;s)<=R=10^-12.                                (1)
```

A smooth frozen solution `u` is assumed to satisfy the same boundary rows and the exact equation below. There is no smallness assumption on `u`. Thus the estimate may be used in a later first-exit argument without inserting its intended conclusion as a hypothesis.

The elapsed-time fields represent the same physical past by `U(s)=J_E V(-s)`, with the parallel tensor time-slot reflection and opposite elapsed initial velocity. The all-state coefficient certificate is invariant under this exact transformation. The normalized equation is

```
P_v u := u_ss-A^i(x,v)D_i u_s-B^ij(x,v)H_ij u=F(x,v,Dv),
F=C(x,v)+E^mu(x,v)D_mu v+G^mu nu(x,v)(D_mu v,D_nu v).              (2)
```

Spatial indices have four values, `mu,nu` have five, and `H_ij` is the spatial product Hessian. The coefficient definitions are

```
a=-q^ss, A^i=2q^si/a, B^ij=q^ij/a,
C=-cbar/a, E=-bbar/a, G=-dbar/a.
```

The exact bundle source is obtained from the coordinate source by the full substitution `partial V=DV-Omega V`, including the `-b Omega V` term and all physical reference-connection terms. Its zero-state constant is exactly zero, so `C(x,0)=0` identically. No coordinate first jet is silently relabeled as a covariant first jet.

The accepted temporal-seed coefficient certificate supplies

```
Cstar=10^500000
```

for all ordered horizontal/vertical derivatives of `A,B,C,E,G` through total order sixteen, in full output/base-array Frobenius norm and Hilbert state-input operator norm. The horizontal derivative includes `-K_V(Omega_i V)` and the actions on all output, direction and state-input slots. The global tensor embedding is `||T||infty<=10^7||T||H3` for at most sixteen slots. Hence (1) implies

```
||v||star,infty<=2*10^-5<10^-3,
||Dv||infty, ||v_s||infty<=10^-5.                                (3)
```

The coefficient evaluation stays in its actual all-state domain. High-rank coefficient arrays are controlled by `Cstar`, not by this tensor embedding.

The new accepted `NORMALIZED_METRIC_MARGINS.md`, accepted original-proof SHA256 `0a941cb3a994ef775edb3746234ee3d5556d467bd9c3317c9b81d58084d38c12`, proves the fixed-product-frame margins

```
|A|<=1/19,
(931/1280) I<=B<=(20/19)I,
19/20<=B^zz<=20/19.
```

In particular we may use `I/2<=B<=2I`. The ratios use cancellation of the same-point `rho^-2`, with the angular `sech` conversion retained. At the shell the predecessor's four Dirichlet fields give `A^z=B^zA=0` exactly. Bulk mixed components are not set to zero.

## 2. Exact energy identity and evaluated production

For a section `w`, define

```
E_v[w]=1/2 integral_Sigma (|w_s|^2+B^ij<D_iw,D_jw>+|w|^2).
```

Let `Z1(w)^2=||w_s||2^2+||Dw||2^2+||w||2^2`. The metric margins give

```
Z1(w)^2/4<=E_v[w]<=Z1(w)^2.                                      (4)
```

If `P_vw=g`, metric compatibility and one integration by parts give

```
dE_v[w]/ds = integral <g,w_s> + integral <w,w_s>
 -1/2 integral (D_i A^i)|w_s|^2
 +1/2 integral (partial_s B^ij)<D_iw,D_jw>
 -integral (D_i B^ij)<D_jw,w_s>
 +integral_boundary [(A^i n_i)|w_s|^2/2
                        +B^ij n_i<D_jw,w_s>].                  (5)
```

The fixed product Hessian includes its derivative-index connection. Equation (5) is not a componentwise integration that drops the base connection. At `z=0`, the coefficient traces and orthogonal split reduce the flux to `-B^zz<w_z,w_s>=0`: the Dirichlet part has `w_D,s=0`, and the Neumann part has `w_N,z=0`. At `z=6ell`, all components have `w_s=0`. These vanishings hold for every time/rotation-commuted section used below. No normal derivative is assigned an additional homogeneous boundary condition.

The all-state chain rule and (3) give

```
||DA||infty, ||DB||infty, ||B_s||infty <=2 Cstar.
```

Spatial divergence is a contraction over four directions, of operator norm two in the full arrays. Thus `||div A||infty,||div B||infty<=4 Cstar`. The four production contributions in (5) are bounded respectively by

```
E_v[w],       4 Cstar E_v[w],
4 Cstar E_v[w],       8 sqrt(2) Cstar E_v[w]<12 Cstar E_v[w].
```

The last bound uses `||Dw||2<=2sqrt(E)` and `||w_s||2<=sqrt(2E)` from (4) and the definition of energy. Consequently

```
dE_v[w]/ds <=sqrt(2E_v[w])||g||2+21 Cstar E_v[w].                 (6)
```

The homogeneous zero-energy case is included by regularization rather than division by zero.

## 3. Natural rotation words and their initial energy comparison

Let `R_1,...,R_6` be the natural lifts of the six rotations on every tensor slot, and let

```
w_(p,I)=partial_s^p R_I u,       p+|I|<=13,
Y14(u;s)^2=sum_(p+|I|<=13) E_v[w_(p,I)].                          (7)
```

All rotation words are ordered, with no commutation or finite-mode projection. The number of `(p,I)` words is

```
N13=sum_(k=0)^13(14-k)6^k=18,807,399,380.                         (8)
```

The natural lifts preserve the mixed split, commute with `D_z,partial_s`, and intertwine the full product covariant derivative including its direction slots. They act on coefficients as well as on the unknown. The estimates hold for arbitrary angular sections, hence also the actual Gamma-invariant fields; individual rotation words need not remain in that invariant subspace.

For clarity, an aggregate norm comparison can avoid charging the count (8) a second time. Put `A_k(T)^2=sum_|I|=k ||R_I T||2^2`. At a north-pole orthonormal frame, three rotation generators are translations with zero derivative and the other three are stabilizer rotations of operator norm one. Global isometry invariance therefore gives, for a rank-at-most-r tensor,

```
sum_a |R_a T|^2 <=|D_Omega T|^2+6r^2|T|^2.
```

More explicitly, apply this pointwise identity to each complete `D^j T` and sum `0<=j<=h`. Naturality gives

```
sum_a ||R_a T||Hh^2 <=[1+6(r+h)^2] ||T||H^(h+1)^2.
```

The angular derivative of `D^j T` is a component of the full next ordered derivative, with no derivative reordering. Iterating this inequality (or its equivalent aggregate-word recurrence) gives

```
A_k(T)<=40^k ||T||Hk,             r+k<=16,                       (9)
```

because `1+6*16^2=1537<40^2`. The same assertion holds for direct sums and the two-copy union. Apply it to `partial_s^p u`, `partial_s^(p+1)u`, and `D partial_s^p u`; their original ranks are at most two, two, and three, respectively. No derivative is commuted through another spatial derivative. For `p+k<=13`, (4) and (9) give

```
sum_|I|=k E_v[w_(p,I)]
 <=40^(2k)(2||partial_s^p u||H^(14-p)^2
                 +||partial_s^(p+1)u||H^(13-p)^2).
```

Summing the geometric series and then the time hierarchy yields

```
Y14(u;s)<=sqrt(6)*40^13 X14(u;s)<10^22 X14(u;s)
                                           <=10^30 X14(u;s).    (10)
```

If `u` has the common nonlinear initial jets `J0,...,J14`, then `X14(u;0)=I14`. In particular

```
Y14(u;0)<=C_Y I14,              C_Y=10^30.                       (11)
```

For a frozen successor, the common jets follow from the retained triangular recursion whenever the predecessor has those jets and the initial Cauchy pair is the prescribed one. No replacement of `I14` by the uncompleted Cauchy-pair norm is made.

## 4. Exact principal commutator

For `Z` equal to `partial_s` or a natural rotation lift, covariance gives

```
[P_v,Z]u=(Z A)^i D_i u_s+(Z B)^ij H_ij u,                         (12)
```

where the coefficient derivatives are the full natural derivatives, including all displayed tensor slots. In particular a rotation derivative of a coefficient is not its componentwise directional derivative. The spatial connection/Hessian is natural under the rotation action, and the fixed time connection is trivial; this is why no unrecorded Hessian commutator remains in (12).

For a length-n ordered word `W=Z_1...Z_n`, repeated ordered Leibniz gives

```
P_v(Wu)=W F + sum_(empty!=S subset {1,...,n})
 [ (Z_S A)^i D_i partial_s(Z_(S^c)u)
  +(Z_S B)^ij H_ij(Z_(S^c)u) ].                                (13)
```

Every subsequence retains its original order; repeated identical terms carry their actual multiplicity. This exact formula needs no commutation of two rotations. It has at most `2(2^n-1)` principal commutator terms. The signs are positive because `[P,W]=PW-WP` and the principal spatial/mixed terms in `P` carry minus signs.

## 5. Evaluated word-source and commutator norms

The independently accepted companion [ROTATION_MIXED_PRODUCTS.md](ROTATION_MIXED_PRODUCTS.md), accepted original-proof SHA256 `00ec8f8b2f0e4253fc8863d68adeb79769b04eb4445290a3eae368e771549715`, proves the following interfaces with the exact full word square sum from (7):

```
G_source := (sum_(p+|I|<=13)||partial_s^p R_I F(x,v,Dv)||2^2)^(1/2)
                      <=10^500100 X14(v),
G_comm := (sum_(p+|I|<=13)||[P_v,partial_s^p R_I]u||2^2)^(1/2)
                      <=10^500100 X14(u).                       (14)
```

Here is the finite ledger used in that companion. First compose the principal coefficients `A(x,v),B(x,v)` covariantly. At derivative degree `n<=13`, the horizontal/state chain rule has at most `(n+1)!` ordered leaves, counting multiplicities. Coefficient arrays use the all-state `Cstar` certificate. All but a maximal predecessor factor have intrinsic order at most six and use the low `H3` embedding; their pointwise factors are at most `10^-5`. The field-free coefficient leaf is bounded in L2 by `Cstar`, since the volume of one copy and also the two-copy union is less than one. The resulting coefficient estimate is

```
||partial_s^p K(x,v)||Hk
 <=sqrt(n+1)(n+1)! Cstar (1+X_n(v)),   p+k<=n<=13,
K=A or B.
```

Only after this composition is performed are rotations applied to the rank-one or rank-two coefficient output. The complete derivative-slot rotation estimate has `Drot=40^14<10^23`. Thus for a coefficient word of degree `q`,

```
q<=7:  ||T^J K(v)||infty <10^40 Cstar,
q<=13: ||T^J K(v)||2 <10^35 Cstar(1+X13(v)).
```

No high-rank vertical coefficient map is fed into the rank-sixteen Sobolev embedding or rotation comparison.

For the source, use the exact identity

```
C(x,v)=integral_0^1 C_V(x,theta v)[v] dtheta.
```

The three resulting source groups begin with at most two field factors and retain at least one after differentiation. Through thirteen derivatives their leaf count is

```
3*(4*5*...*16)=10,461,394,944,000 <10^14.
```

A source leaf has total intrinsic field order at most fifteen, individual order at most fourteen, and at most one factor of order at least eight. Put that factor in L2 and all others in L-infinity using three spatial derivatives and `E3 X10(v)<=10^-5`. The maximal low field rank with this reserve is `2+7+3=12<=16`. The source covariant Sobolev bound is `<10^15 Cstar X14(v)`. Applying the rotation bridge to the rank-two-or-less complete source costs less than `10^23`; taking the word square sum costs less than `10^6`. Therefore

```
G_source <=10^44 Cstar X14(v) <=10^500100 X14(v).
```

For a commutator term in (13), let `q>=1` derivatives hit the coefficient and let `d=n-q+2` be the unknown's intrinsic order. Since `n<=13`, `d<=14` and `q+d<=15`. If `q<=7`, the coefficient is bounded in L-infinity by `10^40 Cstar`, and the other factor is bounded in L2 by `Drot X14(u)<10^23 X14(u)`. If `q>=8`, then `d<=7`: the unknown factor is bounded in L-infinity by `E3 Drot X10(u)<10^30 X10(u)`, and the coefficient has the displayed L2 bound. Thus the high coefficient derivative is paired with a low unknown derivative; the unknown is not assumed small.

There are fewer than `10^5` ordered-subset terms per word. The per-word result is

```
||[P_v,T]u||2
 <=10^70 Cstar [X14(u)+X13(v) X10(u)].
```

Now use `sqrt(N13)<10^6`, `X13(v)<=R`, and `X10(u)<=X14(u)`. This gives

```
G_comm <=10^77 Cstar X14(u) <=10^500100 X14(u).
```

These calculations include the complete coefficient-slot action by composing covariantly first and then applying the natural rotation bridge to the low-rank outputs. They charge every array contraction, ordered derivative word and high/low factor; there is no unnamed product or commutator constant in (14).

## 6. Integrated inequality and final constants

With (14), summing (6) for the words and using Cauchy–Schwarz gives, in the regularized upper-derivative sense,

```
Y14'(s)<=G_source+G_comm+16 Cstar Y14(s).
```

Indeed the forcing coefficient before rounding is `1/sqrt(2)` and the production coefficient is `21 Cstar/2`. The regularization `sqrt(Y14^2+delta^2)` proves the same integrated inequality at zero, upon letting `delta` decrease to zero.

Using (10), (11), and (14) then gives

```
Y14(u;s)<=C_Y I14
  +integral_0^s [K_u X14(u;r)+K_v X14(v;r)]dr,
C_Y=10^30,       K_u=10^500101,       K_v=10^500100.              (15)
```

The last rounding is the exact comparison

```
10^500100+16*10^500000*10^30 <10^500101.
```

This is a conditional frozen estimate for smooth compatible inputs. It is independent of angular cutoff and allows arbitrary unknown norm. It applies unchanged to the retained two-copy union. It is not a proof that the seed, successive frozen solutions, or a nonlinear limit remains in the predecessor ball. Normal recovery, the complete Picard majorant, finite-Sobolev completion and evolved parameter estimates remain separate obligations.

## 7. Source pins and retained scope

Read-only source tree: the repository, registered base `d27cfd999aa2e57995c08b1228a73d3da7c973d1`. Paths below are relative to its `research/` directory. The byte pins identify the actual inputs; they do not remove their stated analytic premises.

| Input | SHA256 |
|---|---|
| `HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md` | `37ea4bbe555eba7aaae1140d778d7e2d6e86ee20f7a489e5e35825750411fcc9` |
| `HDBLAST_TEMPORAL_JETS_SEED_20261114/NONLINEAR_INITIAL_JETS.md` | `7c5366015f801c883d9e632c5cb3fdcf4774d0783bcd62d7cfdcab9e21f38d95` |
| `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md` | `d6f80e6158226b253c56f2f60c8d9c6cfa2ee2940e79bde42dd8fe8e48a5d6c4` |
| `HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/SYMMETRY.md` | `3d1f0fe3f94e7e60aa37048c27b51c694fe8664d2b18e24b0d0d56cb2e8425a0` |
| `HDBLAST_ORDERED_EVOLUTION_20261106/QUADRATIC_FORCING.md` | `26871e9ed63a2345737b2f55fa855d0144ef8ebc85f985573e783c1dee923062` |
| `HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md` | `e78d8b4ece66335aabe2704b4eb7cec599a3bdb5686e6704d97e006d29d29a4c` |

The new normalized-margin proof is pinned in Section 1 and the new product proof in Section 5. The earlier prospective plan selected this calculation; its formerly unevaluated coefficients are not used as numerical assumptions here. The preserved system is the actual half-action Einstein–scalar source, with the original reference connection and same-state past conventions. This proof makes no geometric or physical inference beyond its conditional frozen estimate.
