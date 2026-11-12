# Evaluated high-order background first-wave estimate

This continuation follows merged PR #135. The immutable source snapshots
are authenticated in SOURCE_PINS.json. Read the retained analytic premises
in inputs/ORDERED_BACKGROUND_06I.md, inputs/INITIAL_TRANSPORT_06H.md and
inputs/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md before combining bounds.

The conclusion is conditional on the same actual smooth background and
the selected compatible initial first derivative as those sources. Put

```
T=13/8000000,  Z=6*10^-6,
M=[0,Z] x S3,  q_sp=dz^2+kappa,
u=V1,  N14(t)=||u(t)||H14(E)+||u_t(t)||H13(E).
```

All Sobolev norms use the fixed product connection, full sphere measure,
and full Frobenius norm on the symmetric block of the sixteen-dimensional
bundle E. The following explicit bound holds for the background linear
first wave on the retained mixed-wall problem:

```
sup_(0<=t<=T) N14(t) < 10^4000000 exp(10^152017).            (1)
```

This is a deliberately large analytic ceiling. It is not a bound for
a finite-amplitude nonlinear family, a physical energy, or a cosmological
origin conclusion. No spectral truncation or extra evolved
W13,infinity assumption is used.

## 1. Actual operator, coefficient interfaces, and boundary problem

The fixed bundle is

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D.
```

Its original tensor summands have at most two rotating covariant slots.
The radial and time factors are product-parallel. Let D denote the fixed
bundle/product connection and let
`Delta_E=trace_kappa D_Omega^2` be its angular rough Laplacian. In the
energy-sign convention the exact background operator is

```
L u = beta(z)^2 [u_tt-u_zz-s(t)^2 Delta_E u] + B D u + A u,
beta=1/rho,  s=sech(t),  L u=0.                            (2)
```

Here `B D u=B^t u_t+B^z u_z+B^Omega D_Omega u`; B is one map
from the full direct sum of five derivative directions into E. The
maps A and B are the actual zero-state fixed-bundle maps, after the
coordinate-to-bundle substitution. For example, in local components
with `D_i=partial_i+Omega_i`, they are

```
A=c_V|0-b^k|0 Omega_k
      +q0^{ij}(partial_i Omega_j-Omega_j Omega_i),
B^k=b^k|0+2q0^{ik}Omega_i-q0^{ij}gamma^k_ij I.              (3)
```

The c,b in (3) are the coordinate coefficients, not coefficients obtained
by merely relabeling partial derivatives as covariant derivatives. All
reference-connection terms already present in those coefficients remain
included. Formula (2) follows from the exact product-Hessian ordering
proved in ORDERED_BACKGROUND, including the cancellation of its
time-angular moving-frame and base-Hessian corrections.

At the shell impose `u_D=0`, `D_z u_N=0`, with the fixed orthogonal
4/12 split. At the outer wall impose the retained homogeneous Dirichlet
condition on every component. The initial fields and velocities vanish
on fixed open collars of both ends. These are the actual boundary and
corner hypotheses used below.

The input initial norm is

```
||u(0)||H14+||u_t(0)||H13 <10^5200.                         (4)
```

The retained energy E satisfies, for every section w,

```
||w||2^2 <=2 E[w],
||w_t||2^2, ||w_z||2^2 <=36000 E[w],
||D_Omega w||2^2 <=72000 E[w],
(||w||H1^2+||w_t||2^2)^(1/2) <=269 sqrt(E[w]).             (5)
```

It also obeys `E[w]<= (||w||H1^2+||w_t||2^2)/2`. For `Lw=f`
and the stated homogeneous boundary conditions, its proven differential
estimate is

```
E[w]' <= Kstar E[w]+190 sqrt(E[w])||f||2,
Kstar=10^150006.                                          (6)
```

This includes the time dependence of the moving frame and the actual
background coefficient derivatives. No unknown solution jet appears in
Kstar because the operator here is fixed at V=0.

## 2. Order-thirteen coefficient extension from the raw sixteen jets

This section supplies a new derivative allocation proof; it does not
reinterpret the preceding capsule's certified order-twelve statement.
Let `Craw=10^150000`. The raw coordinate coefficients c,b,d have the
certified entry-sum ceiling Craw for each ordinary mixed base/state
partial of total order at most sixteen.

The exact first-jet substitution, before differentiation, is

```
cbar(V)=c(V)-b^i(V)Omega_i V
 +d^{ij}(V)(Omega_iV,Omega_jV)
 +q^{ij}(V)(partial_iOmega_j-Omega_jOmega_i)V,

bbar(V)Q=b(V)Q-2d(V)(Omega V,Q)
 +q^{ij}(V)(Omega_iQ_j+Omega_jQ_i-gamma^k_ij Q_k),
dbar=d.                                                   (7)
```

In (2), `A=cbar_V|0` and `B=bbar|0`. Thirteen base derivatives
and at most two state derivatives of the maps in (7) require raw mixed
order at most fifteen. At most fourteen derivatives of Omega are needed,
because of the displayed `partial Omega`. Both counts are within the
pinned order sixteen bounds.

The predecessor's elementary complex-chart calculation gives ordinary
connection derivative entry sums below `10^60` through sixteen. Its
inverse-matrix normalized q-jet bound, with factorial and state-index
sums restored, is below `10^18000`. At the zero section, the finite
coefficient extraction has at most two connection factors. The same
`5^16<10^12` Leibniz allocation and `16^2` state-index sums therefore
give, through thirteen ordinary base derivatives, the bound

```
10^20 [Craw*(10^60)^2+10^18000*(10^60)^2] <10^150200.       (8)
```

For covariant conversion there are at most three E slots of geometric
rank two, at most two base coefficient slots, and thirteen derivative
slots: `6+2+13=21<32`. Every component index has dimension at most
sixteen. The enlarged stage count

```
5*(14+32*16)=2630 <10^4
```

charges differentiation of the existing factors and connection action on
all slots, including derivative slots. At stage j there are at most j
connection factors. Thus the final full component tensor bound, including
the `10^6` finite norm conversion, is

```
10^(150200+64*13+6)=10^151038.
```

Two optional moving-coframe conversion factors cost fewer than 127
further decades, including their `3^13` Leibniz count. Consequently the
common ceiling

```
Dstar=10^152000                                           (9)
```

bounds the full multilinear Hilbert norm of every mixed product
covariant base derivative through thirteen of
`cbar_V|0,cbar_VV|0,bbar|0,bbar_V|0,dbar|0,q_V|0`.
Only A and B are needed in this proof. In particular
`||partial_t^k A||infty,||partial_t^k B||infty<=Dstar` for
`k<=13`. These are fixed background coefficient derivatives.

The background rho jets give
`max_(j<=13)|partial_z^j rho^2|<10^210`, by the ordinary
two-factor Leibniz rule and the pinned rho derivatives below `10^100`.
Set

```
P=rho^2 B,  Q=rho^2 A,  H=10^152300.                       (10)
```

Since `2^13*10^(210+152000)<H`, their mixed coefficient derivatives
through thirteen have norm below H. We use only total order twelve of
P and Q below. Finally

```
|partial_t^j s^2|<10^60,    0<=j<=13.                      (11)
```

For example the same complex disc of radius 1/100 has
`|sech^2 t|<4`; `4*13!*100^13<10^60` proves (11).

## 3. Exact rotational naturality of the full operator

For `1<=a<b<=4`, let `K_ab=n_a partial_(n_b)-n_b partial_(n_a)`
be the six canonical Killing fields on the unit sphere in R4. Let
`R_ab` be their natural Lie derivative lifts to every summand of E.
They act on tensor slots as well as component functions; dt and dz
are invariant. Write R_alpha when the six pairs are enumerated.

The exact identities are

```
[L,R_alpha]=0,
[partial_t,R_alpha]=[partial_z,R_alpha]=0.                 (12)
```

To prove the first, let an element of SO(4) act by pullback on the
angular factor of spacetime and on all its tensor bundles. The actual
background `g0=rho^2(dz^2-dt^2+cosh^2(t)kappa)` and `phi0(z)` are
invariant. Ricci, scalar-wave, and the difference of physical and
reference connections in relative harmonic reduction are natural under
this pullback. The proper-y component root is a scalar, the radial
cross block is a covector, and the tangential metric block is a symmetric
tensor. Their reconstruction and the exact boundary flow commute with
the pullback: their pointwise operations are inverse, multiplication,
and contraction with invariant tangential tensors, with scalar t,z
coefficients. The fixed product connection is also preserved by the
pullback. Thus the entire transformed nonlinear equation is equivariant.
Differentiating this identity at the identically zero section shows that
its actual linearization commutes with the pullback. Differentiating the
group action gives (12). This includes all lower terms in (3), not
merely the principal symbol.

The split, shell, and outer wall are invariant, so every R_alpha
preserves their homogeneous boundary rows. These identities require
natural tensor lifts; they would generally be false with raw angular
component differentiation in their place. No assertion that the six
R_alpha commute with one another is needed.

## 4. Rotation words and covariant tensor derivatives

For an E-valued geometric tensor T, let

```
A_k(T)^2=sum_(|I|=k)||R^I T||2^2,
A_0(T)=||T||2.                                           (13)
```

All ordered words occur. The generators on derivative tensors always
act naturally on the added derivative slots. Isometry invariance of the
connection gives

```
R_alpha(D_sp^j T)=D_sp^j(R_alpha T).                       (14)
```

Write `R_alpha T=D_(K_alpha)T+S_alpha T`; on tensors of rank r,
`||S_alpha||op<=r`, since `||nabla K_ab||op<=1`. The following
pointwise identity is useful:

```
sum_alpha |R_alpha T|^2
 = |D_Omega T|^2 + sum_alpha |S_alpha T|^2.                (15)
```

At `n=e4`, the three generators K_a4 are `-e_a` and have zero
tangential covariant derivative. The other three generators vanish at
that point and act only through their tangent rotation matrices. Thus
the directional and spin terms have no cross term, proving (15) there.
An arbitrary SO(4) rotation moves any sphere point to e4 and transforms
the six generators by an orthogonal matrix in the adjoint representation.
The sums and fiber metrics are invariant, so (15) holds everywhere.
It remains valid on the direct sum and on derivative tensors.

Iterating (15) and (14) yields

```
||D_Omega^k T||2 <= A_k(T).                               (16)
```

No rank constant is needed in this direction. For the reverse direction
the deliberately weaker bound

```
sum_alpha ||R_alpha T||2^2
 <=2||D_Omega T||2^2+12r^2||T||2^2
```

gives

```
A_(k+1)(T)^2
 <=2 A_k(D_Omega T)^2+12r^2 A_k(T)^2.
```

For every k with `r+k<=16`, induction using
`2+12*16^2<60^2` implies

```
A_k(T)<=60^k ||T||Hk.                                    (17)
```

The Hk norm here can be the full product spatial norm, which dominates
the angular one. For an original E section and `j+k<=14`, the same
argument with its j derivative slots gives

```
sum_(|I|=k)||D_sp^j R^I v||2^2
 <=60^(2k)||D_sp^j v||Hk^2.                              (18)
```

In particular the sum of its H1 norms is at most
`2*60^(2k)||v||H^(k+1)^2`. These comparisons include the tensor
derivative slots and the full Frobenius metric.

Two convenient consequences, used with generous constants below, are

```
||D_Omega v||2 <=10 (sum_(|I|<=1)||R^I v||2^2)^(1/2),
||Delta_E v||2 <=100 (sum_(|I|<=2)||R^I v||2^2)^(1/2).      (19)
```

In fact (16) gives the smaller constants one and sqrt(3); (19) leaves
ample room without changing the later recurrence. The same inequalities
apply after fixed t,z derivatives, whose direction labels do not rotate.

## 5. Initial temporal jets: a fully evaluated finite recurrence

Put `U_j=partial_t^j u|_(t=0)` and
`a_j=||U_j||H^(14-j)` for `0<=j<=14`. Multiplying (2) by rho^2
and solving for the second time derivative gives

```
u_tt=u_zz+s^2 Delta_E u
     -P^t u_t-P^z u_z-P^Omega D_Omega u-Q u.               (20)
```

For `j=0,...,12`, its exact initial recursion is

```
U_(j+2)=partial_z^2 U_j
 +sum_(k=0)^j binom(j,k)(partial_t^k s^2)(0) Delta_E U_(j-k)
 -sum_(k=0)^j binom(j,k){
   (partial_t^k P^t)(0) U_(j-k+1)
  +(partial_t^k P^z)(0) partial_z U_(j-k)
  +(partial_t^k P^Omega)(0)D_Omega U_(j-k)
  +(partial_t^k Q)(0) U_(j-k)}.                            (21)
```

For any coefficient map with covariant derivatives through m bounded by
H, tensorial Leibniz and Cauchy-Schwarz give

```
||K v||Hm <=2^m sqrt(m+1)H||v||Hm <=10^5 H||v||Hm,
0<=m<=12.                                                (22)
```

This is a linear multiplier estimate. No high pointwise norm of v is
assumed. Fixed product traces also give

```
||partial_z^2 v||Hm<=||v||H^(m+2),
||Delta_E v||Hm<=sqrt(3)||v||H^(m+2).                     (23)
```

For the second inequality, commute no derivatives: the parallel angular
trace of the appropriate last two covariant slots is bounded by sqrt(3).
Every tensor on the right has the required ordered derivative slots.

Take `m=12-j` in (21). All its differentiated coefficient norms are
available: a coefficient with k time derivatives needs at most m spatial
derivatives, and `k+m<=j+12-j=12`. Each input field has enough of
the known H^(14-index) derivatives. From (11), (22), and (23),

```
a_(j+2) <= [1+2^j(2*10^60+4*10^5 H)]
                         max_(0<=l<=j+1) a_l
          <= R_init max_(0<=l<=j+1)a_l,
R_init=10^10 H=10^152310.                                 (24)
```

There are exactly thirteen steps, producing U2 through U14. Since
`max(a0,a1)<10^5200`,

```
max_(j<=14) a_j < R_init^13*10^5200 =10^1985230.           (25)
```

Every operation in (21) is local and homogeneous. The initial open zero
collars therefore persist in all these jets. They satisfy every used
shell and outer corner condition, rather than merely satisfying a
generic H14/H13 regularity requirement.

## 6. Tangential and time energy closes without derivative loss

Define

```
Y(t)^2=sum_(p+|I|<=13)E[partial_t^p R^I u](t).              (26)
```

There are
`sum_(k=0)^13(14-k)6^k<15*6^14<10^13` summands. Each section
obeys the same homogeneous boundary rows. The energy, including its
moving-frame term, is applied separately to these full E sections.

Using (18) and the upper energy inequality at t=0, for each `p+k<=13`
the sum over words of length k is bounded by
`60^(2k)(2a_p^2+a_(p+1)^2)/2`. There are 105 pairs (p,k), so

```
Y(0) <=sqrt(315/2)*60^13 max_(j<=14)a_j
      <=10^25 max_(j<=14)a_j <10^1985300.                 (27)
```

There is no additional factor 6^k in this line: the aggregate rotation
bound (18) already includes all words.

By (12), only time commutators remain. For `w=partial_t^p R^Iu`,
the exact equation is

```
L w=sum_(k=1)^p binom(p,k){
  beta^2(partial_t^k s^2)Delta_E partial_t^(p-k)R^Iu
 -(partial_t^k B)D(partial_t^(p-k)R^Iu)
 -(partial_t^k A)partial_t^(p-k)R^Iu}.                    (28)
```

The principal term has k>=1. For
`v=partial_t^(p-k)R^Iu`, one additional rotation still has total
word degree at most thirteen. This permits one angular derivative to
be placed inside the energy and the other in that additional rotation.
For an explicit bound, apply (15) to `D_Omega v` and use naturality:

```
||Delta_E v||2
 <=sqrt(3)(sum_alpha||D_Omega R_alpha v||2^2)^(1/2)
 <=sqrt(3)*269 Y <=4000 Y.                               (29)
```

Thus no fourteenth commuted word or H15 initial datum is required. In
particular the commutator does not need normal recovery to close its
own energy estimate. Its lower first-derivative terms are already
controlled by (5).

The k-th term of a fixed word's right side in (28), before its binomial
coefficient, has norm at most
`(4000*10^60+538*Dstar)Y`. Here beta^2<=1, and the two lower terms
use the full first-jet and zeroth-order norms from (5). Since

```
2^13(4000*10^60+538*Dstar)<10^7 Dstar,
```

the square sum over all words gives

```
F_comm(t):=(sum_(p+|I|<=13)||L(partial_t^pR^Iu)||2^2)^(1/2)
          <=C_comm Y(t),
C_comm=10^152020.                                        (30)
```

Indeed the word-count square root is below 10^7, so the smaller
`10^14 Dstar=10^152014` would already suffice.

Sum (6) over the words and apply Cauchy-Schwarz. Regularizing at zero
energy, if needed, gives

```
Y'(t) <=(Kstar/2+95 C_comm)Y(t),
Y(t) <=Y(0)exp[(Kstar/2+95 C_comm)t].                     (31)
```

The coefficient is below `10^152022`. On the whole retained interval,

```
T*(Kstar/2+95 C_comm) < T*10^152022 <10^152017.             (32)
```

Both derivative counts in this calculation are fixed: time derivatives
of the linear maps reach thirteen; no derivative of an evolving
coefficient field or unknown pointwise solution norm is used.

## 7. Normal recovery: all recurrence constants evaluated

For `0<=r<=14`, define

```
W_r(t)^2=
 sum_(j+p+|I|<=14, 0<=j<=r)||partial_z^j partial_t^p R^Iu||2^2.
                                                               (33)
```

First relate W1 to Y. Its terms with j=0 and total word degree at
most thirteen are bounded by `2Y^2`. Terms with j=1 are bounded by
`36000Y^2`. At total degree fourteen with a time derivative, use
the time-energy term of the word of degree thirteen, costing another
`36000Y^2`. Finally the pure rotation words of degree fourteen obey

```
sum_alpha||R_alpha w||2^2
 <=2||D_Omega w||2^2+48||w||2^2 <=144096 E[w]
```

for each word w of degree thirteen. It follows that

```
W1^2 <=216098 Y^2 <=10^6 Y^2,
W1 <=10^3 Y.                                             (34)
```

Solving the actual equation for the second radial derivative gives

```
u_zz=u_tt-s^2 Delta_E u
       +P^t u_t+P^z u_z+P^Omega D_Omega u+Q u.             (35)
```

Fix `2<=r<=14` and a new row
`partial_z^r partial_t^p R^Iu` with `r+p+|I|<=14`.
Because the full operator commutes with every rotation, first use
(35) for `R^Iu`, then apply
`partial_z^(r-2) partial_t^p`. This order of argument introduces
no angular derivative of a coefficient and does not replace a true
operator commutator by a componentwise one.

The leading time term has radial count r-2. The angular term also
has radial count r-2, and its two angular derivatives are controlled
by rotation words of length at most two using (19). The first-order
radial term has radial count at most r-1. All other terms have still
smaller radial count. Every solution word has total order at most
fourteen, and therefore occurs in W_(r-1).

The largest coefficient differentiation order is `r-2+p<=12`.
For the principal angular term its time allocation has binomial sum
at most `2^p`; for the lower terms the combined t,z allocation has
sum at most `2^(r-2+p)<=2^12`. Use H from (10), (11), and the
conservative constants ten and one hundred in (19). The four lower
terms cost `1+1+10+1=13` times H. Hence every new row obeys

```
||partial_z^r partial_t^p R^Iu||2
 <=[1+2^12(100*10^60+13 H)] W_(r-1)
 <=10^6 H W_(r-1).                                      (36)
```

These are L2 coefficient multiplier estimates at a fixed time. No
derivative of u is placed in L-infinity, and no normal boundary
condition is used.

The total number of terms in W14 is exactly

```
sum_(k=0)^14 binom(16-k,2)6^k =135413275530 <10^12.        (37)
```

Thus the number of new rows at any stage is smaller than 10^12.
Adding their squares to the previous rows gives

```
W_r <=sqrt[1+10^12(10^6H)^2] W_(r-1)
     <=M_norm W_(r-1),
M_norm=10^20 H=10^152320.                                (38)
```

There are thirteen stages r=2,...,14, so (34) yields

```
W14 <=10^(13*152320+3)Y=10^1980163 Y.                    (39)
```

This is a triangular recovery in radial count, with no derivative
loss. The procedure never applies the boundary energy theorem to
`u_z` or assigns a homogeneous Neumann or Dirichlet row to a normal
derivative. The only boundary-energy fields are those in (26).

## 8. Convert to the full spatial H14/H13 norm and conclude

The product connection commutes with the radial derivative. At a fixed
spatial derivative order j, its Hilbert tensor square is

```
||D_sp^j v||2^2
 =sum_(a+b=j)binom(j,a)||D_Omega^b partial_z^a v||2^2.     (40)
```

This is the multiplicity of normal/angular derivative-index placements;
no angular covariant derivatives are commuted past one another. By (16),
the terms are bounded by the corresponding rotation-word norms.
Since each binomial coefficient is at most `2^14`, summing j<=14
gives

```
||u||H14 <=2^7 W14,
||u_t||H13 <=2^(13/2) W14,
N14 <=256 W14.                                           (41)
```

Inserting (27), (31), and (39), and using `256<10^3`, gives the
more explicit finite-time estimate

```
N14(t) <10^3965466 exp[(10^150006/2+95*10^152020)t].        (42)
```

Indeed `1985300+1980163+3=3965466`. The rounded bound (1) follows
from (32). All exponents and all finite recurrence multipliers in this
calculation have now been evaluated; no unspecified product, commutator,
initial-jet, or normal-recovery constant remains in (1).

For justification at finite regularity, one may first use smooth initial
approximations vanishing on fixed slightly smaller open boundary collars.
Their recursive corner jets still vanish by (21). The retained smooth
linear solver supplies the solutions on the fixed slab. The preceding
estimates apply uniformly; applying them to differences yields convergence
at H14/H13 from convergence of the original data in those norms. The
limit agrees, by the retained energy uniqueness, with the linear solution
whose data are (4). It also agrees with the retained parameter derivative
V1, which solves this same linear initial-boundary problem. In particular
we do not infer X14 differentiability of the nonlinear family merely from
its previously stated C2 parameter dependence into X10.

The proof applies to the identical reflected union by using the summed
Hilbert energies and the initial bound already supplied for that union.
Every finite coefficient and norm bound is unchanged. No arbitrary number
of copies is asserted.

## 9. Immediate forcing consequence and limits

The predecessor proves for the complete second-wave force, including
the principal variation term and all factors two,

```
||R2[u](t)||H12 <=10^304240 N14(t)^2.                     (43)
```

Combining (1) and (43) is now legitimate throughout the background slab:

```
sup_t ||R2[u](t)||H12
 <10^8304240 exp(2*10^152017)
 <10^8400000 exp(10^152018).                              (44)
```

This closes the evolved whole-slice H12 forcing interface; it is not
the full collection of time-differentiated source and initial source-jet
norms required by a high-order second-wave estimate. Those mixed source
jets require a separate calculation. A quantitative finite-amplitude
metric patch, nonlinear evolution remainder, observer clock, actual area,
positive denominator, and common proper-time observable domain also
remain outside this theorem. The underlying conditional background and
linear existence premises are retained, not established by finite
arithmetic replay. No novelty or external peer-review claim is made.
