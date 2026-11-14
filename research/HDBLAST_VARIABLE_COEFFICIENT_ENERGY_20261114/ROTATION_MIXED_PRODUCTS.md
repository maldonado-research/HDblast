# Evaluated natural rotations and mixed products through X14

14 November 2026. This conditional calculation supplies explicit constants for the frozen energy and interior normal recovery. It is conditional on the authenticated fixed product equation and full-array coefficient certificate. It asserts no evolved solution, nonlinear ball closure, parameter-Hessian bound, or physical consequence.

The fixed spatial manifold is `Sigma=[0,Z] x S3`, `Z=6*10^-6`, with full measure and the product metric. The field bundle is

```
E = R_v + Sym^2(R dt + T*S3)_H + R_w + (R dt + T*S3)_D.
```

It has sixteen independent components and the positive direct-sum/full-Frobenius norm. Its maximal rotating tensor rank is two. All covariant derivative direction slots are retained. `D` is the product connection, including the parallel time slot. The same estimates apply to the retained two-copy convention: sum squared norms on the disjoint copies, use the identical pointwise constants on each, and note that their combined volume is still less than one. The literal mixed norm is

```
X_m(w;s) = sum_(p=0)^m ||partial_s^p w(s)||_{H^(m-p)(Sigma)}.
```

The accepted inputs are `Cstar=10^500000` for ordered horizontal/vertical derivatives of `A,B,C,E,G` through sixteen; `C(x,0)=0`; the full-tensor embedding `||T||infinity<=E3||T||H3`, `E3=10^7`, for the fields/derivative arrays used below; and `X14(v)<=R=10^-12`. The coefficient horizontal derivative is the certified derivative including `-K_U[Omega U]`. The predecessor hypothesis implies `|v|star<=2E3 R=2*10^-5<10^-3` and `E3 X10(v)<=10^-5<1`.

The coefficient certificate is `research/HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md`, sections 1, 3--5, in the repository at base `d27cfd999aa2e57995c08b1228a73d3da7c973d1`. The actual frozen equation is the one in `research/HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`, sections 2 and 4. No input file or public capsule is changed.

## 1. Six natural lifts, including every tensor slot

For `1<=a<b<=4`, let `J_ab` be the Euclidean skew matrix of operator norm one and `Omega_ab(x)=J_ab x`, `x in S3`. Choose the natural infinitesimal pullback/pushforward action `R_ab` on the tensor bundle. Up to the immaterial common sign convention,

```
R_ab T = D_(Omega_ab) T + spin_ab(T),
||spin_ab(T)|| <= r ||D Omega_ab||op |T| <= r |T|
```

on a tensor of rotating rank at most `r`. Parallel `dt,dz` slots contribute zero to the spin. Direct sums obey the same estimate with maximal rank. For covariant and contravariant slots the actions have opposite signs but identical bounds.

These are isometry lifts. Thus `R_ab D T = D(R_ab T)` as complete derivative tensors, and similarly for every iterated covariant derivative and the full Hessian. This is an identity of natural differential operators; omitting the direction-slot action would invalidate it. Also `[R_ab,D_z]=[R_ab,partial_s]=0`. The natural split of the boundary bundle is preserved. The estimates hold for all sections, including rotation words that need not retain a separate reflection-invariance condition.

For completeness, the rotation coefficients themselves have evaluated derivatives. `|Omega_ab|<=1`, `|D Omega_ab|F<=sqrt(2)`, and the Killing-field identity and parallel round curvature give

```
|D^m Omega_ab|F <= (sqrt(12))^floor(m/2) * sqrt(2),  0<=m<=17,
```

where the odd initial factor is `sqrt(2)` and the even initial factor may be replaced by one. Here `|Riem(S3)|F=sqrt(12)`. Consequently all these derivatives are `<10^5`; this is a direct finite bound, not a coordinate-cover constant. The argument below uses naturality and does not repeatedly differentiate local frame coefficients.

At `x=e4`, the three rotations `Omega_i4` form an orthonormal tangent frame and have `D Omega_i4=0`; the other three rotations vanish and their derivatives are the orthonormal-plane isotropy generators. SO(4) invariance transports the resulting identity to every point:

```
sum_(a<b) |R_ab T|^2
  = |D_S T|^2 + sum_(1<=i<j<=3) |rho(J_ij)T|^2,
|D_S T|^2 <= sum_(a<b)|R_ab T|^2
            <= |D_S T|^2 + 3 r^2 |T|^2.                    (1)
```

All norms in (1) are full tensor norms. In particular no trivialization, angular cutoff, or spectral-gap assumption occurs.

Write `A_k(T)^2=sum_(|I|=k)||R^I T||L2^2`, retaining all `6^k` ordered words and no reordering. Applying (1) to all derivative tensors, and using naturality, gives

```
||D_S^k T||L2 <= A_k(T),                                 (2)
A_k(T) <= 40^k ||T||Hk,       if r+k<=16.                  (3)
```

Here is a finite induction for the upper bound. At the next step,

```
A_k(T)^2 <= A_(k-1)(D_S T)^2 + 3 r^2 A_(k-1)(T)^2.
```

The induction hypotheses on `T` and `D_S T` apply because the latter has rank `r+1`. Since `1+3 r^2<=769<40^2` for `r<=16`, (3) follows. The same proof works with the weaker coefficient `6r^2`, since `1+6*16^2=1537<1600`. For an individual word, a mixed radial/time derivative, and additional spatial `h` derivatives, the consequence used below is

```
||D_z^a partial_s^p R^I T||Hh
 <= 40^|I| ||partial_s^p T||H^(a+|I|+h),                 (4)
```

provided the tensor ranks plus spatial derivative slots in this application do not exceed sixteen. To avoid a hidden multiplicity from nested Sobolev sums, the stronger one-step identity used here is

```
sum_(a<b)||R_ab T||Hh^2
 <= [1+3(r+h)^2] ||T||H^(h+1)^2.
```

Indeed sum (1) over the full spatial derivative tensors of degrees zero through `h`; the angular derivative term at degree `j` is bounded by the full spatial derivative of degree `j+1`. Iterating this inequality over the ordered rotation words proves the Hh form of (3), with `40^|I|`, directly. Product radial derivatives commute, and taking a fixed radial component is a contraction. This proves (4) without an extra nested-Sobolev multiplicity. All applications below obey the displayed slot restriction. Set once and for all

```
Drot=40^14=26843545600000000000000 <10^23.                 (5)
```

## 2. Word count and energy/covariant comparisons

For `Y14^2=sum_(p+|I|<=13) E_v[partial_s^p R^I u]`, the exact word count is

```
N13 = sum_(k=0)^13 (14-k)6^k = 18,807,399,380,
sqrt(N13)<10^6.                                         (6)
```

Using the evaluated metric margins `1/2<=B<=2`, one has `E_v[w]<=||w_s||2^2+||Dw||2^2+||w||2^2`. For each fixed pair `(p,k)`, (3), applied to `u`, `u_s`, and the complete `Du`, bounds the corresponding sum over words by `3*40^26*X14(u)^2`. There are `14*15/2=105` pairs. Consequently

```
Y14(u) <= sqrt(315)*40^13 X14(u) <10^23 X14(u),
Y14(u;0) <=10^23 I14.                                  (7)
```

The deliberately rounded interface `C_word=10^30` is therefore valid. No word-count factor is omitted: the angular `6^k` sums are already present in (3), and the remaining 105 time/length pairs are displayed.

For normal recovery define

```
W14(u)^2 = sum_(a+p+|I|<=14) ||D_z^a partial_s^p R^I u||2^2.
```

The exact product splitting gives

```
||D_sp^k T||2^2 = sum_(a=0)^k binom(k,a)||D_z^a D_S^(k-a)T||2^2.
```

Using (2), `binom(k,a)<=2^14`, and Cauchy--Schwarz over the fifteen temporal terms yields

```
X14(u) <= sqrt(15*2^14) W14(u) <500 W14(u).               (8)
```

Larger finite comparison factors, such as `1920` or `10^4`, remain valid. Equations (2) and (8) include every angular mode and all normal/angular placements.

## 3. Covariant composition along an arbitrary mixed predecessor

All coefficient differentiations in this section are performed covariantly before the rotation bridge is applied. A derivative of `K(x,v)` has the two exact chain branches

```
horizontal derivative of K,       K_U[Dv].
```

Horizontal differentiation has precisely the certified connection correction. Thus there is no extra ordinary-coordinate connection term and no missing state-input action. After `j` differentiations of a coefficient there are at most `j` field factors. The next differentiation has at most `j+2` branches. For `n` ordered mixed derivatives there are at most

```
2*3*...*(n+1)=(n+1)!                                   (9)
```

leaves, counting repeated leaves with their multiplicities. Every actual field factor in (9) has positive derivative order; their total order is at most `n`. For `n<=13`, select a largest-order factor for L2. Each remaining factor has order at most six and its infinity norm is at most `E3 X9(v)<=10^-5`. Hence each nonconstant leaf is bounded by `Cstar X_n(v)`. The one leaf without actual field factors has pointwise norm at most `Cstar`; its L2 norm is at most `Cstar`, since

```
Vol(Sigma)=6*10^-6 *2*pi^2 <1.
```

The same argument for each lower spatial derivative in a Sobolev norm proves, for `p+k<=n<=13`,

```
||partial_s^p K(x,v)||Hk
 <= sqrt(n+1)(n+1)! Cstar [1+X_n(v)].                    (10)
```

We use this only for the composed principal coefficients `K=A(v),B(v)`, whose rotating tensor ranks are one and two. When `n=0`, the bracket is a harmless enlargement of the direct zeroth-order bound.

Coefficient maps before composition may have many more than sixteen slots because their vertical inputs are full field spaces. Those maps are bounded solely by the given operator-to-full-array `Cstar` certificate. They are never fed into (3), (4), or the H3 embedding. This avoids imposing an unjustified rank-sixteen bound on high vertical coefficient jets.

For `n<=13`, `sqrt(n+1)(n+1)!<10^12`; for `n<=10`, the same quantity is `<10^9`. Thus (4), (5), and (10) yield the following explicit estimates for a single ordered derivative `T^J` from `partial_s,D_z,R_ab`, of degree `q`:

```
q<=7:   ||T^J K(v)||infinity
        <= E3 Drot *10^9 Cstar (1+R) <10^40 Cstar;      (11)
q<=13:  ||T^J K(v)||2
        <= Drot *10^12 Cstar [1+X13(v)]
        <10^35 Cstar [1+X13(v)].                        (12)
```

For (11), H3 embedding raises the covariant mixed order only to `q+3<=10`; all coefficient values still lie in the state patch. Lower-degree versions of (12) may replace `X13(v)` by the actual available degree. The array contractions and spatial/time selectors have operator norm at most one. No dimension factor is silently attached to a single index evaluation.

## 4. Actual lower source, with its zero-state factor

Use the exact radial extraction

```
C(x,v)=Ctilde(x,v)[v],
Ctilde(x,v)=integral_0^1 C_U(x,theta v)dtheta.
```

The certified connection commutes with radial scaling. Therefore `Ctilde` has mixed horizontal/vertical derivatives through fifteen bounded by `Cstar`; the available original coefficient order is sixteen. The exact source is

```
f(v)=Ctilde(x,v)[v]+E(x,v)Dv+G(x,v)(Dv,Dv).              (13)
```

It has three groups, initially at most two field factors, and at least one field factor in every leaf. After `j` differentiations there are at most `j+2` factors, so the next derivative has at most `j+4` branches. Consequently through `n<=13` differentiations there are at most

```
L_n = 3*(4*5*...*(n+3)) = 3*(n+3)!/3!,
L13 = 10,461,394,944,000 <10^14                         (14)
```

leaves. This includes the actual two-slot differentiation of `G(Dv,Dv)`. A leaf has total field derivative order at most `n+2<=15`, individual field order at most `n+1<=14`, and at most one factor of order at least eight. A maximal-order factor is placed in L2; every other factor has order at most seven and is bounded by `E3 X10(v)<=10^-5`. Thus the covariant source hierarchy satisfies

```
p+k<=n<=13:
||partial_s^p f(v)||Hk <= sqrt(14) L13 Cstar X_(n+1)(v)
                       <10^15 Cstar X_(n+1)(v).        (15)
```

Only the actual field and its derivative arrays use the tensor embedding. Their largest low rank, including the three-derivative reserve, is `2+7+3=12<=16`. Pure background-volume forcing is absent because of the exact zero-state extraction.

Applying (4) to the entire composed source, whose output rank is at most two, gives a per-word constant `<10^38 Cstar`. Multiplying by the square root of (6) gives

```
[sum_(p+|I|<=13)||partial_s^p R^I f(v)||2^2]^1/2
 <=10^44 Cstar X14(v).                                 (16)
```

For normal-recovery words with `a+p+|I|<=12`, the corresponding per-word estimate is

```
||D_z^a partial_s^p R^I f(v)||2
 <=10^38 Cstar X13(v).                                 (17)
```

No high coefficient-map rotation is used in (16)--(17): the covariant chain expansion is completed first, and then the output section is converted by naturality.

## 5. Principal commutator and the tame high-coefficient allocation

Let

```
P_v=partial_s^2-A(v)^i D_i partial_s-B(v)^ij H_ij.
```

Each natural lift commutes with the complete Hessian; time and radial derivatives also commute with this product Hessian. Therefore an ordered word `T` of degree `n` has exactly the ordered-subset expansion

```
[P_v,T]u = sum_(K=A,B) sum_(nonempty ordered subsets J of T)
              (T^J K(v)) * (remaining word acting on D u_s or H u).
```

Every displayed external sign is positive for `[P_v,T]=P_v T-T P_v`: the minus signs of the two principal terms are reversed in this commutator. The chosen generator sign convention is already inside each `T^J`. Subset words retain their original relative order. There are at most `2(2^n-1)<10^5` terms for `n<=13`; repeated terms retain their multiplicities. There are no curvature or direction-slot commutators missing from this expression. No rotation generators are reordered.

For a term let `q>=1` be the number of differentiations on the coefficient and `d=n-q+2` the total unknown order. For energy words, `n<=13`, so `d<=14` and `q+d<=15`.

* If `q<=7`, (11) puts the coefficient in infinity norm. By (4), the other factor is bounded by `Drot X14(u)<10^23 X14(u)`. Their product is `<10^63 Cstar X14(u)`.
* If `q>=8`, then `d<=7`. The unknown factor is put in infinity norm using (4) and H3 embedding: it is at most `E3 Drot X10(u)<10^30 X10(u)`. Equation (12) puts the coefficient in L2. Their product is `<10^65 Cstar [X10(u)+X13(v)X10(u)]`.

This is the explicit high coefficient jet times low unknown allocation. It gives the tame, per-energy-word bound

```
||[P_v,T]u||2
 <=10^70 Cstar [X14(u)+X13(v)X10(u)]
 <=10^70 Cstar [X14(u)+X14(v)X13(u)].                    (18)
```

No smallness or ball condition on `u` is used. Using (6), `X13(v)<=R`, and `X10(u)<=X14(u)`, one obtains

```
[sum_(p+|I|<=13)||[P_v,partial_s^p R^I]u||2^2]^1/2
 <=10^77 Cstar X14(u).                                 (19)
```

For normal recovery take `T=D_z^a partial_s^p R^I`, `n=a+p+|I|<=12`. A nonempty coefficient subset then leaves unknown total order at most thirteen. The same two allocations, with the lower available unknown norm, give

```
||[P_v,T]u||2 <=10^71 Cstar X13(u).                     (20)
```

Here high coefficient derivatives are still bounded by the predecessor ball `X14(v)<=R`. Keeping a more detailed tame term is possible, but (20) suffices for the normal recurrence. It involves only the unisolated operator; derivatives of `1/Bzz` are not required. It imposes no differentiated boundary condition.

## 6. Uniform rounded interfaces for the companion ledgers

The following single numerical product constant is supported by all the explicit inequalities above:

```
Mprod = 10^500100.
```

It can be used as each of:

```
energy source aggregate:        <= Mprod X14(v),
energy commutator aggregate:    <= Mprod X14(u),
normal source per ordered row: <= Mprod X13(v),
normal commutator per row:     <= Mprod X13(u).
```

The largest displayed expenditure is `10^77 Cstar=10^500077`, leaving a stated decimal reserve of twenty-three exponents. The energy word interface is independently `C_word=10^30`; the inverse angular interface is one; the full normal/word interface is `<500`. Exact-integer arithmetic checks are recorded in `evidence/rotation_product_arithmetic.json`; those checks certify the stated arithmetic, while the identities and allocations above provide the analytic argument. These constants are independent of angular cutoff, iteration number, and the unknown successor's size. They do not by themselves close the Picard ball or remove the finite-regularity and parameter gates in the prospective plan.
