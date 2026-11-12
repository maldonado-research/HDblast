# Quantified initial temporal jets and triangular radial recovery

New manual proof for the retained background first wave.
The frozen initial-transport and ordered-evolution capsules are unchanged.
This note closes the initial-value and normal-recovery parts of the proposed
time/rotation energy argument. It imports the actual background linear
operator and coefficient certificates; it does not substitute a scalar
component Laplacian or a different boundary problem.

## 1. Exact operator, norms, and finite coefficient constants

Use the fixed sixteen-dimensional bundle E, positive Frobenius/direct-sum
fiber metric, product connection, spatial metric `dz^2+kappa`, and full
measure `dz dmu_kappa` of INITIAL_TRANSPORT and ORDERED_BACKGROUND.
Let `Delta_E=trace_kappa D_Omega^2` be its angular rough Laplacian.
The exact zero-state first wave satisfies

```
u_tt=u_zz+s(t)^2 Delta_E u
       -P^t u_t-P^z u_z-P^Omega D_Omega u-Q u,             (1)
s(t)=sech(t), P=rho^2 B, Q=rho^2 A.
```

A and B are the complete fixed-product lower maps from the actual
coordinate-to-bundle substitution, including the full varying reference
connection. They are not scalar estimates applied separately to components.

For every mixed time/spatial covariant derivative of total order at most
twelve, the preceding background coefficient proof supplies the A,B
operator bound `Dstar=10^152000`. Its finite derivative count applies to
mixed derivatives: the number of base derivatives is the total order, not
separately twelve time plus twelve spatial derivatives. The known background
rho jets give `|partial_z^j rho^2|<10^210` through order thirteen. Consequently

```
H=10^152300
```

is a bound for every required mixed derivative of each of P^t, P^z,
P^Omega, Q through total order twelve. Indeed

```
2^12 * 10^210 * 10^152000 < 10^152300.                    (2)
```

The scalar coefficients obey `|(s^2)^(k)|<10^60` for `k<=12` on the
retained time interval, by the existing time-frame Cauchy bound. A fixed
linear multiplier K whose first m spatial derivatives are bounded by H
satisfies

```
||K v||Hm <=2^m sqrt(m+1) H ||v||Hm <=10^5 H ||v||Hm,
0<=m<=12.                                                (3)
```

The last assertion follows from `2^24*13<10^10`. There is no solution
W13,infinity hypothesis. Only fixed background coefficient derivatives
appear in (2)-(3).

## 2. Literal initial temporal recurrence and its numerical ceiling

Let `U_j=partial_t^j u(0)` and

```
S_j=||U_j||H^(14-j), 0<=j<=14.
```

The verified initial interface is `S_0+S_1<10^5200`. For `j=0,...,12`,
differentiate (1) exactly j times in t and evaluate at zero:

```
U_(j+2)=partial_z^2 U_j
 +sum_(k=0)^j binom(j,k)(s^2)^(k)(0) Delta_E U_(j-k)
 -sum_(k=0)^j binom(j,k){
       (partial_t^k P^t)(0) U_(j-k+1)
      +(partial_t^k P^z)(0) partial_z U_(j-k)
      +(partial_t^k P^Omega)(0) D_Omega U_(j-k)
      +(partial_t^k Q)(0) U_(j-k)}.                       (4)
```

There are no independent higher Cauchy time derivatives to choose. The
induction has exactly the required regularity: set `m=12-j`; then
`partial_z^2 U_j` and `Delta_E U_j` use H^(m+2)=H^(14-j), while
every first-order term uses at most one further spatial derivative of an
already available U. In a coefficient term with k time derivatives and
at most m spatial derivatives, the coefficient order is at most
`k+m=k+12-j<=12`. No order-twenty-four coefficient bound is used.

The angular trace has norm at most sqrt(3)<2, hence its Hm estimate uses
`2||U||H^(m+2)` without commuting angular derivatives. Product connection
contractions are parallel, so this also includes all derivative slots.
Equation (4), with (2)-(3), gives the literal nonnegative triangular bound

```
S_(j+2) <= S_j
 +sum_(k=0)^j binom(j,k){
       2*10^60 S_(j-k)
       +10^5 H [S_(j-k+1)+3 S_(j-k)]}.                   (5)
```

For an exact integer replay, initialize `b_0=b_1=10^5200` and define
`b_(j+2)` by the right-hand side of (5) with S replaced by b. Then
`S_j<b_j` for all j. All constants in this finite recurrence are specified
integers. A convenient closed ceiling follows because, for `j<=12`,

```
1+2^j(2*10^60+4*10^5 H) <10^10 H.
Rinit=10^10 H=10^152310.
```

Starting from the maximum of S_0,S_1, the thirteen updates j=0,...,12 give

```
max_(0<=j<=14) S_j < Rinit^13 * 10^5200
                    =10^1985230.                       (6)
```

Every initial value and velocity is identically zero on fixed open collars
of both boundary components. Equation (4) is local and homogeneous, and
the coefficient maps are smooth. Induction shows that all U_j through
fourteen retain those open zero collars. Thus the shell and outer temporal
corner compatibility used here follows from the actual selected packet,
not from a claim that arbitrary H14/H13 data satisfy those corners.

## 3. From temporal jets to the initial tangential energy

Let R_alpha, alpha=1,...,6, be the natural lifts of the six canonical
SO(4) Killing fields on every E summand. Keep all ordered words R^I.
The separate symmetry proof establishes the exact identities

```
[L_bg,R_alpha]=0, [partial_t,R_alpha]=[partial_z,R_alpha]=0,
R_alpha(D_Omega^j v)=D_Omega^j(R_alpha v),                 (7)
```

where the natural action in the last line includes every derivative slot.
The original E tensor rank is at most two. Natural lifts preserve the
homogeneous mixed shell split and the homogeneous outer Dirichlet wall.

For a tensor T of rank r, let
`A_k(T)^2=sum_(|I|=k)||R^I T||L2^2`. The identities
`R_alpha T=D_(K_alpha)T+S_alpha T`, `||S_alpha||<=r`, and
`sum K_alpha tensor K_alpha=kappa^-1` imply

```
A_(k+1)(T)^2 <=2 A_k(D_Omega T)^2+12 r^2 A_k(T)^2.
```

Using the natural intertwining in (7), induction through the required
fourteen orders gives the deliberately enlarged comparison

```
A_k(T) <=60^k ||T||Hk
```

whenever the original rank plus the relevant derivative order is at most
sixteen. This follows from `2+12*16^2=3074<60^2`. Applying the same
estimate to the angular derivative slots, and commuting the parallel radial
derivative, yields

```
sum_(|I|=k)||R^I U_p||H1^2
             <=2*60^(2k)||U_p||H^(k+1)^2.                (8)
```

The factor two safely covers the value/radial/angular splitting. No scalar
component harmonic comparison is being substituted for the bundle action.

Define the actual commuted energy

```
Y13(t)^2=sum_(p+|I|<=13) E[partial_t^p R^Iu](t).
```

The preceding energy density is at most one half the fixed positive
first-jet square. At zero, use (8) for the field and the same rotation
estimate for its velocity U_(p+1). There are
`sum_(k=0)^13(14-k)=105` pairs (p,k), after each exact-length word family
has already been aggregated in A_k. Therefore

```
Y13(0)^2 <=315*60^26*(max_j S_j)^2,
Y13(0) <=10^30 max_j S_j <10^2000000.                    (9)
```

For the arithmetic it suffices to check `315*60^26<10^60` and
`1985230+30<2000000`. The estimate is substantially looser than the
calculated finite sum; there is no extra factor counting individual
rotation words once they are aggregated in (8).

## 4. Angular estimates needed for radial recovery

For any E section v, put `A_k(v)=(sum_(|I|=k)||R^I v||L2^2)^(1/2)`.
The pointwise natural-lift identities give

```
||D_Omega v||2 <= A_1(v)+5A_0(v),
||D_Omega^2 v||2 <= A_2(v)+13A_1(v)+40A_0(v).
```

For completeness, apply (7) to D_Omega v, whose tensor rank is at most
three. The resulting factor sqrt(6)*3 is less than eight. The stack of
first derivatives of R_alpha v is bounded by A_2+5A_1. Combining with
the first inequality yields `A_2+(5+8)A_1+8*5A_0` exactly as displayed.
Since `1+13^2+40^2=1770<50^2`, the trace inequality supplies

```
||D_Omega v||2 <=10*(sum_(|J|<=1)||R^Jv||2^2)^(1/2),
||Delta_E v||2 <=100*(sum_(|J|<=2)||R^Jv||2^2)^(1/2).     (10)
```

These estimates apply equally after any time or radial derivatives,
which are invariant/parallel directions. They require no fixed-isotypic
or Casimir-eight restriction; that optional sharper route is unnecessary
for the constants proved here.

## 5. Exact finite radial recurrence

At each time set

```
W_r(t)^2=sum_{0<=j<=r, j+p+|I|<=14}
                 ||partial_z^j partial_t^p R^Iu(t)||2^2,
0<=r<=14.                                               (11)
```

The index set is finite. Over all r<=14 its total number of summands is

```
Nall=sum_(k=0)^14 binom(16-k,2) 6^k
    =135413275530 <10^12.                               (12)
```

For `2<=r<=14`, fix a new row with `r+p+|I|<=14`; put a=r-2 and
v=R^Iu. By exact commutation with the complete operator, v satisfies (1)
with the same coefficient maps. Solve it for v_zz, then apply
`partial_z^a partial_t^p`. The literal identity is

```
partial_z^r partial_t^p v
 =partial_z^a partial_t^(p+2) v
  -sum_(k=0)^p binom(p,k)(s^2)^(k)
                         Delta_E partial_z^a partial_t^(p-k) v
  +sum_(b=0)^a sum_(k=0)^p binom(a,b)binom(p,k){
      (partial_z^b partial_t^k P^t)
                           partial_z^(a-b)partial_t^(p-k+1)v
     +(partial_z^b partial_t^k P^z)
                           partial_z^(a-b+1)partial_t^(p-k)v
     +(partial_z^b partial_t^k P^Omega)
                           D_Omega partial_z^(a-b)partial_t^(p-k)v
     +(partial_z^b partial_t^k Q)
                           partial_z^(a-b)partial_t^(p-k)v } . (13)
```

No angular derivative hits a coefficient in this identity: the full
rotated section v is first known to satisfy the same equation. Radial/time
coefficient differentiation has `b+k<=a+p<=12-|I|`. Thus (2) covers
all these terms. In the field factors the leading terms have radial
count a=r-2 and the first-order radial term has at most a+1=r-1.
After applying (10), every added rotation word has total differential
order at most fourteen and remains in W_(r-1). For example the angular
Laplacian term has total order
`a+(p-k)+|I|+2<=r+p+|I|<=14`. The recurrence is strictly triangular
in radial order.

The leading time term in (13) costs at most W_(r-1). The Laplacian term
costs `100*10^60 W_(r-1)` per time-binomial allocation. The four lower
terms cost at most `(1+1+10+1)H W_(r-1)=13H W_(r-1)` per allocation.
Both binomial sums are at most 2^12. Hence every new row has the
explicit bound

```
||partial_z^r partial_t^p R^Iu||2
 <=[1+2^12(100*10^60+13H)] W_(r-1)
 <=10^6 H W_(r-1).                                      (14)
```

Using (12) to sum the new rows gives a literal finite recurrence

```
W_r^2 <=[1+Nall*(10^6 H)^2] W_(r-1)^2,
W_r <=M W_(r-1),  r=2,...,14,
M=10^20 H=10^152320.                                    (15)
```

Indeed `1+10^12*(10^6H)^2<10^26 H^2<10^40 H^2`. The chosen
M leaves seven decades beyond even the intermediate square-root bound.
All inequalities multiplying a possibly zero field are non-strict.
Equation (15) is an interior L2 identity/estimate. It imposes no fictitious
homogeneous shell row on u_z or any higher normal derivative.

## 6. Starting radial recovery from Y13, and returning to N14

The previous energy controls all zero/one radial jets in (11), with the
only remaining degree-fourteen pure rotation row handled by one natural
lift estimate. Explicitly:

* zero radial derivatives of total time/rotation degree at most thirteen
  contribute at most `2Y13^2`;
* one radial derivative contributes at most `36000Y13^2`;
* zero radial derivatives of degree fourteen with at least one time
  derivative contribute at most `36000Y13^2`, by the energy velocity;
* pure rotation words of length fourteen contribute at most
  `144096Y13^2`.

For the last claim, write each word as R_alpha w with |I(w)|=13.
For the original rank-at-most-two bundle,
`sum_alpha||R_alpha w||2^2<=2||D_Omega w||2^2+48||w||2^2`.
The fixed energy has `||D_Omega w||2^2<=72000E[w]` and
`||w||2^2<=2E[w]`, yielding 144096. All these sums are subfamilies
of the same Y13 square, so no word-count multiplier is added. Thus

```
W_1^2 <=216098Y13^2 <=10^6Y13^2,
W_1 <=10^3 Y13.                                         (16)
```

The constant comparison is strict, while both field estimates remain valid
when Y13=0.

To recover full spatial Sobolev norms, the natural-lift identities give,
by induction through order fourteen,

```
||D_Omega^k v||2
 <=40^k*(sum_(|I|<=k)||R^I v||2^2)^(1/2).                (17)
```

At stage j<=13 the coefficient is at most
`1+sqrt(6)(j+2)<40`; derivative slots are explicitly included. In the
fixed product geometry radial and angular covariant derivatives commute.
The square norm of a full spatial derivative of order m therefore has
the sum over j radial slots with binomial multiplicity binom(m,j).
Combining with (17), for u through order fourteen and u_t through thirteen,

```
N14=||u||H14+||u_t||H13
 <=2*[sum_(m=0)^14 1601^m]^(1/2) W_14
 <=10^30 W_14.                                          (18)
```

The equality `1601=1+40^2` accounts for those mixed-slot multiplicities.
The constant follows, for example, from
`4*sum_(m=0)^14 1601^m<10^60`. It is independent of a coordinate
partition and uses the same fixed product norm as the forcing theorem.

Equations (15)-(18) give the fully evaluated recovery interface

```
N14(t) <=10^30 * M^13 * 10^3 Y13(t)
        =10^1980193 Y13(t).                              (19)
```

## 7. Completed interface and combination with the tangential propagator

The companion tangential proof establishes, from exact rotational
equivariance and the newly derived order-thirteen time-coefficient budget,

```
Y13(t) <=Y13(0) exp(10^152022 t),
0<=t<=13/(8*10^6),
10^152022*13/(8*10^6) <10^152017.
```

For negative t, apply that same proof in the variable tau=-t. The
coefficient P^t and odd time derivatives acquire signs, but every
coefficient envelope is absolute on |t|<=T. The initial velocity only
changes sign, the spatial boundary conditions are unchanged, and the
same energy/recurrence constants therefore apply in the reversed time
direction. This justifies the absolute-time version displayed here.

Combining that result with (9) and (19) yields

```
N14(t) <10^4000000 exp(10^152017)                         (20)
```

throughout the retained background slab. Indeed
`1980193+2000000=3980193<4000000`. This is the quantitative first-wave
N14 conclusion once the companion commutator/energy proof is bound to
the same sources. It is not a nonlinear family estimate. The prior
quadratic-force theorem can now receive this numerical first-wave input.

All calculations first apply to the smooth retained first packet and
smooth linear solution; alternatively use compatible smooth approximation
and the retained finite-regularity linear theorem. The estimates use only
the displayed H14/H13 initial norm. Energy is applied exclusively to
time/rotation derivatives preserving the actual boundary rows; interior
radial recovery does not alter those rows.

## 8. Exact replay ledger

A checker can evaluate the finite recurrences (5), (13)-(15) literally,
or verify their displayed domination by positive integer constants. The
following elementary checks cover every numerical reserve used here:

```
2^12*10^210*10^152000 <10^152300,
2^24*13 <10^10,
for j=0,...,12:
  1+2^j*(2*10^60+4*10^5*10^152300) <10^152310,
13*152310+5200 ==1985230,
2+12*16^2 <60^2,
315*60^26 <10^60,
1985230+30 <2000000,
1+13^2+40^2 <50^2,
sum_(k=0)^14 binom(16-k,2)*6^k ==135413275530 <10^12,
1+2^12*(100*10^60+13*10^152300) <10^6*10^152300,
1+10^12*(10^6*10^152300)^2 <(10^152320)^2,
2+36000+36000+144096 ==216098 <10^6,
for j=0,...,13: 6*(j+2)^2 <39^2,
4*sum_(m=0)^14 1601^m <10^60,
30+13*152320+3 ==1980193,
1980193+2000000 <4000000,
13*10^152022 <8*10^6*10^152017.
```

The derivative-allocation checks should additionally enumerate the ranges
in (4) and (13): `k+12-j<=12`; in radial recovery `b+k<=r-2+p<=12-|I|`,
all field radial counts at most r-1, and all field total counts after
the at-most-two added rotations at most fourteen. These are the actual
no-loss conditions, not a numerical test of a Sobolev theorem.

This note does not close the full high-order second-wave mixed source and
corner hierarchy, nonlinear nonzero-family commutators, an actual-family
metric patch, a finite amplitude interval, or the observer/area/denominator
remainder. The very large ceilings are in fixed model normalization and
carry no interpretation as a measured physical blast energy.
