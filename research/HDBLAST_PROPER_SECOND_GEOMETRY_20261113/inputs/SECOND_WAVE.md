# Evaluated forced second wave in the full X13 topology

This continuation follows merged PR #136. SOURCE_PINS.json authenticates
eight immutable public inputs at the base commit. The mixed source proof is
[MIXED_SOURCE.md](MIXED_SOURCE.md); the retained first-wave and energy
premises are preserved in the inputs directory.

Under its retained smooth-background, compatible-data and mixed-wave
existence premises, the actual second amplitude derivative at zero obeys

```
v=V2=partial_epsilon^2 V|0,
X13(v;t)=sum_(j=0)^13 ||partial_t^j v(t)||H^(13-j),

sup_(0<=t<=T) X13(v;t) <10^10200000 exp(10^152019),
T=13/8000000.                                            (1)
```

The norms are exactly the fixed positive product-bundle Sobolev norms
of the predecessors. In particular (1) bounds
`N13(v;t)=||v(t)||H13+||v_t(t)||H12`. It is a statement about the
second derivative at the zero state. It does not assert C2 dependence
of the nonlinear family into X13, a numerical amplitude interval, or a
physical observable remainder.

## 1. Actual forced equation, signs, and data

The spatial cylinder is `M=[0,Z]xS3`, `Z=6*10^-6`, with product
metric `q_sp=dz^2+kappa`, measure `dz dmu_kappa`, and the fixed bundle

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D.
```

The symmetric block uses the full Frobenius metric. Let D be the fixed
product connection and H its full spacetime Hessian, including the base
derivative-index connection. Let `u=V1` be the already bounded first
amplitude derivative. Write the actual fixed-bundle equation as

```
q(x,V):H V=c(x,V)+b(x,V)DV+d(x,V)(DV,DV),
c(x,0)=0,
A=c_V|0, B=b|0,
C2=c_VV|0, E2=b_V|0, D2=d|0, Q1=q_V|0.
```

These are the coefficients after the complete coordinate-to-bundle
substitution, including all reference and fixed-connection terms. In the
energy-sign convention the background operator and the exact second
force are

```
Lw=beta(z)^2[w_tt-w_zz-s(t)^2 Delta_E w]+BDw+Aw,
beta=1/rho, s=sech(t), Delta_E=trace_kappa D_Omega^2,

Lv=F,
F=2Q1[u]:H u-C2(u,u)-2E2[u]Du-2D2(Du,Du).                (2)
```

The amplitude convention is `V=epsilon u+epsilon^2 v/2+o(epsilon^2)`.
Every factor two in (2), including that of the principal variation
`Q1[u]:H u`, is retained. There is no assumption that first-order energy
alone bounds this source.

At the shell the exact straightened boundary rows are

```
v_D=0,  D_z v_N=0,   rank(E_D)=4, rank(E_N)=12.
```

At the retained auxiliary outer wall all sixteen components are
Dirichlet. The differentiated rows remain homogeneous because the
straightened boundary operator is linear and independent of epsilon.
Quadratic physical boundary effects are already incorporated in the
exact point transformation; they are not replaced by an uncharged
boundary forcing in this equation.

The second initial transport theorem supplies

```
||v(0)||H14+||v_t(0)||H13 <10^10000.                       (3)
```

Hence the same ceiling applies at H13/H12. Both second initial fields
vanish on fixed open collars of the shell and outer wall. The first
fields and their recursively determined initial time jets do also.

## 2. Complete mixed source interfaces used here

The companion [mixed-source proof](MIXED_SOURCE.md) uses the first-wave mixed derivative
bound W14, rather than only its displayed N14 consequence. Define the
following literal sums, not square sums:

```
S12(t)=sum_(p=0)^12 ||partial_t^p F(t)||H^(12-p),
I11=sum_(p=0)^11 ||partial_t^p F(0)||H^(11-p).
```

Its evaluated conclusions, for the complete force (2), are

```
sup_t S12(t) < Sstar=10^8200000 exp(10^152018),
I11 < Istar=10^4200000.                                  (4)
```

For clarity about the independent source interface, its inputs and
finite allocation are these. Put
`M14(u;t)=sum_(a=0)^14||partial_t^a u(t)||H^(14-a)`.
The first-wave mixed W14 proof gives

```
M14(u;t)<10^4000000 exp(10^152017),
M14(u;0)<10^1985300.
```

In every mixed derivative of (2) of total order at most twelve, the
two differentiated field factors have intrinsic orders drawn from
`(0,0),(0,1),(1,1),(0,2)`, and hence total field order at most fourteen.
The factor with lower total order has order at most seven. Putting that
factor in L-infinity uses only three further spatial derivatives and
the retained tensor H3 embedding constant `10^7`; the other factor
stays in L2 at order at most fourteen. All time and spatial allocations
are contained in M14. The coefficient mixed derivatives through twelve
have bound `10^152000`. The finite three-factor Leibniz sums, the
thirteen time slots, the Sobolev-order square sums, intrinsic-jet norm
factors, and the weights `1+2+2+2=7` yield the conservative estimate

```
S12(t)<=10^152020 M14(u;t)^2.                             (5)
```

Thus its unrounded exponents are 8152020 for the evolving source and
4122620 initially, implying (4). No initial H15/H14 first-wave
extension or evolved W13,infinity norm is used. The companion proof
records the full tensor product calculation; the forced evolution below
uses precisely the source norms stated in (4).

In particular the retained m=13 source costs are now evaluated:

```
sum_(p=0)^12 ||partial_t^p F||_(L1_t H^(12-p))
 =integral_0^T S12(t)dt <T*Sstar,
sum_(p=0)^11 ||partial_t^p F(0)||H^(11-p) <Istar.           (6)
```

The pointwise S12 ceiling also supplies the source terms needed in
normal recovery. A time/rotation source norm alone would not suffice
for that step; its radial derivatives are accounted for below.

## 3. Retained coefficient, energy, and rotation constants

The first-wave proof establishes the following fixed background inputs.
Set

```
Dstar=10^152000,
P=rho^2 B, Q=rho^2 A, Hcoef=10^152300,
Kstar=10^150006,
C_comm=10^152020,
R_init=10^152310, M_norm=10^152320.
```

Mixed coefficient derivatives through thirteen of A and B have full
tensor operator norm at most Dstar. This was derived from raw mixed
order sixteen, not assumed. The normalized maps P,Q have mixed
derivative norm below Hcoef through thirteen, since the rho-squared
jets are below `10^210` and
`2^13*10^(210+152000)<Hcoef`. Here only order twelve is needed
in time commutation and order eleven in initialization and normal
recovery. Also `|partial_t^j s^2|<10^60` through thirteen.

For any homogeneous mixed-boundary section w satisfying `Lw=f`, the
actual energy, including its moving-frame contribution, obeys

```
E[w]'<=Kstar E[w]+190 sqrt(E[w])||f||2,
||w||2^2<=2E[w],
||w_t||2^2,||w_z||2^2<=36000E[w],
||D_Omega w||2^2<=72000E[w],
(||w||H1^2+||w_t||2^2)^(1/2)<=269sqrt(E[w]).               (7)
```

The upper energy bound is half that fixed Hilbert square.

Let R_alpha be the six natural SO(4) lifts generated by
`K_ab=n_a partial_(n_b)-n_b partial_(n_a)` on the unit sphere.
They act on every tensor slot and preserve both boundary operators.
The full background operator satisfies

```
[L,R_alpha]=0,
[partial_t,R_alpha]=[partial_z,R_alpha]=0.                (8)
```

The predecessor proves (8) from equivariance of the actual reference
metric, reduced equations, reconstruction, and boundary flow. It also
proves the pointwise derivative/spin identity for arbitrary bundle
tensors and its two consequences

```
A_k(T)^2=sum_(|I|=k)||R^I T||2^2,
||D_Omega^k T||2<=A_k(T),
A_k(T)<=60^k||T||Hk     when rank(T)+k<=16.                (9)
```

Natural Lie differentiation commutes with the product connection when
its action on the new derivative slots is included. All sums retain
ordered rotation words; the six generators need not commute. No special
Casimir restriction on either u or v is used.

## 4. Normalize the force and charge its radial word array

Set `G=rho^2 F`. For a coefficient K whose spatial derivatives through
m have norm at most C, the finite tensorial multiplier estimate is

```
||K h||Hm<=2^m sqrt(m+1)C||h||Hm.                         (10)
```

For m<=11, the rho-squared derivative bound and (10) are smaller than
`10^220`. Since rho is time independent, the exact source comparisons
are

```
sum_(p=0)^11 ||partial_t^p G(t)||H^(11-p)<=10^220 S12(t),
sum_(p=0)^11 ||partial_t^p G(0)||H^(11-p)<=10^220 I11.      (11)
```

In particular no time derivative of rho is being silently added.

For the later normal recovery define the complete radial/time/rotation
array

```
Z11(G;t)^2=
 sum_(a+p+|I|<=11)||partial_z^a partial_t^p R^I G(t)||2^2.
```

Because radial differentiation commutes with the natural lifts, (9)
applies directly to each radial derivative of each time derivative.
The pure radial labels are product-parallel. For a fixed p there are
at most 78 pairs `(a,k)` with `a+k<=11-p`. Therefore

```
Z11(G;t)
 <=sqrt(78)*60^11
      sum_(p=0)^11||partial_t^p G(t)||H^(11-p)
 <=10^250 S12(t).                                       (12)
```

The numerical comparison uses
`sqrt(78)*60^11<10^22`, followed by (11); the extra eight decades
are reserve. This is a bound on radial source derivatives as well as
tangential ones. It does not infer them from the tangential forcing
energy alone.

## 5. Inhomogeneous initial time jets and corner compatibility

Let `J_j=partial_t^j v(0)` and `a_j=||J_j||H^(13-j)` for
`0<=j<=13`. Solving (2) for the second time derivative gives

```
v_tt=v_zz+s^2 Delta_E v-P^t v_t-P^z v_z
                              -P^Omega D_Omega v-Qv+G.   (13)
```

Thus for `0<=j<=11` the exact initial recurrence is

```
J_(j+2)=partial_z^2 J_j
 +sum_(k=0)^j binom(j,k)(partial_t^k s^2)(0)Delta_E J_(j-k)
 -sum_(k=0)^j binom(j,k){
   (partial_t^k P^t)(0)J_(j-k+1)
  +(partial_t^k P^z)(0)partial_z J_(j-k)
  +(partial_t^k P^Omega)(0)D_Omega J_(j-k)
  +(partial_t^k Q)(0)J_(j-k)}
 +partial_t^j G(0).                                     (14)
```

The last term is essential. It is not part of the homogeneous first-wave
recurrence and is not absorbed into the initial position or velocity.
Its sign is positive in this temporal recovery.

Take the H^(11-j) norm. The spatial trace estimates and coefficient
multiplier argument used for the first wave apply with one fewer
derivative. All coefficient allocations have total order at most eleven:
`k+(11-j)<=11`. The homogeneous row is bounded by R_init times
the largest earlier a_l. By (11), the inhomogeneous row satisfies

```
a_(j+2)<=R_init max_(l<=j+1)a_l +10^220 I11,
0<=j<=11.                                               (15)
```

There are twelve steps. Since `R_init>=2`, the geometric source sum
is at most `R_init^12`. The data bound (3) then gives

```
max_(j<=13)a_j
 <=R_init^12(10^10000+10^220 I11)
 <10^6027941.                                           (16)
```

For the last arithmetic,
`12*152310=1827720`,
`10^10000+10^220 Istar<10^4200221`, and
`1827720+4200221=6027941`.

The source jets used here vanish on the initial open collars. Indeed
F and each of its time derivatives through eleven is a finite local
quadratic expression in first-wave initial temporal/spatial jets, all
of which vanish there. The same is true of G. Starting from the
second initial fields in (3), locality of (14) proves that J2 through
J13 also vanish on those open collars. This establishes the initial
corner conditions used in the energy problem. It is not a statement
that arbitrary H13/H12 data and arbitrary source jets are compatible.

## 6. Forced time/rotation energy through twelve commutations

Define

```
Y(t)^2=sum_(p+|I|<=12)E[partial_t^p R^I v](t).              (17)
```

Every commuted section still obeys the exact shell 4/12 rows and outer
Dirichlet wall. There are fewer than `10^13` words. At t=0, the
aggregate rotation estimate, with derivative slots included, yields

```
Y(0)<=sqrt(273/2)*60^12 max_(j<=13)a_j
     <=10^25 max_(j<=13)a_j <10^6100000.                  (18)
```

Here 91 is the number of pairs `(p,k)` with `p+k<=12`, and
`273/2=3*91/2`. The aggregate rotation norm already counts all 6^k
words, so another such factor is not inserted in (18).

For `w=partial_t^p R^I v`, (8) gives its exact forced equation

```
Lw=partial_t^p R^I F
 +sum_(k=1)^p binom(p,k){
   beta^2(partial_t^k s^2)Delta_E partial_t^(p-k)R^I v
  -(partial_t^k B)D(partial_t^(p-k)R^I v)
  -(partial_t^k A)partial_t^(p-k)R^I v}.                 (19)
```

The commutator sum is bounded by the same C_comm Y as in the first
wave, now at a smaller order. The principal commutator has k>=1,
so one extra rotation fits in the word and its other angular derivative
is controlled by the energy. Precisely
`(p-k)+|I|+1<=12` for every term. No normal recovery is needed to
close the time/rotation energy, and no H14 datum for v is required
by this commutator count.

The external source aggregate has the independent bound

```
F_rot(t)^2=sum_(p+|I|<=12)||partial_t^p R^I F(t)||2^2,
F_rot(t)<=sqrt(13)*60^12 S12(t)<=10^25 S12(t).             (20)
```

For each p, sum the rotation comparison (9) over lengths
`0<=k<=12-p` and use `H^k` below `H^(12-p)`. The square sum
over p is no larger than the square of the literal sum S12.

Summing the energy inequalities (7) and using Cauchy-Schwarz gives,
with the standard regularization if Y=0,

```
Y'(t)<=Lambda Y(t)+95*10^25 S12(t),
Lambda=10^150006/2+95*10^152020 <10^152022.
```

Consequently

```
Y(t)<=exp(Lambda t)
 [Y(0)+95*10^25 integral_0^t S12(r)dr].                   (21)
```

Both the initial forcing hierarchy and the evolving forcing have now
been charged. Since `Lambda T<10^152017`, `T<1`, and (4),(18)
hold, a completely evaluated ceiling is

```
Y(t)<exp(10^152017)
       [10^6100000+10^27*10^8200000 exp(10^152018)]
     <10^8300000 exp(10^152019).                         (22)
```

For example the bracket is below
`10^8200028 exp(10^152018)`, and
`10^152018+10^152017<10^152019`. These generous roundings
avoid treating the inhomogeneous source as zero.

## 7. Forced radial recovery through total order thirteen

For `0<=r<=13` put

```
W_r(t)^2=
 sum_(a+p+|I|<=13, 0<=a<=r)
       ||partial_z^a partial_t^p R^I v(t)||2^2.           (23)
```

The same base comparison as for the first wave gives

```
W1(t)<=10^3 Y(t).                                       (24)
```

It counts separately word values through degree twelve, their first
normal derivatives, top degree thirteen time derivatives supplied by
the time-energy term, and pure rotation words of length thirteen
supplied by lifting one rotation onto an energy field. The previous
constant `216098` still bounds the squared ratio, so (24) follows
without a new boundary estimate.

The exact normal equation has the opposite source sign to (13):

```
v_zz=v_tt-s^2 Delta_E v+P^t v_t+P^z v_z
                               +P^Omega D_Omega v+Qv-G.  (25)
```

For a new row `partial_z^r partial_t^p R^I v`,
`r+p+|I|<=13` and `r>=2`, first use the full commutation (8)
to obtain (25) for R^I v with forcing R^I G. Then differentiate
`r-2` times in z and p times in t. Every homogeneous solution term
has strictly fewer than r normal derivatives and total order at most
thirteen. Its coefficient derivatives have total order at most
`r-2+p<=11`. The first-wave row bound therefore gives

```
||homogeneous part of this row||2<=10^6 Hcoef W_(r-1).
```

The source row is exactly
`-partial_z^(r-2)partial_t^p R^I G`, of total source order at
most eleven. This is a row of the complete array Z11 in (12), not
merely of the tangential source array (20).

The full number of terms in W13 is

```
sum_(k=0)^13 binom(15-k,2)6^k=22568879235 <10^12.          (26)
```

Apply the triangle inequality in the Hilbert direct sum consisting of
the old rows and the new rows. The homogeneous part is bounded by
`sqrt[1+10^12(10^6 Hcoef)^2] W_(r-1)<=M_norm W_(r-1)`.
The new forcing rows are a subset of Z11. Thus the complete forced
recurrence is

```
W_r(t)<=M_norm W_(r-1)(t)+10^250 S12(t),
2<=r<=13.                                               (27)
```

At each fixed t, S12(t) is a number independent of r. Since
`10^250+1<M_norm`, add it to both sides of (27) to obtain

```
W_r(t)+S12(t)<=M_norm[W_(r-1)(t)+S12(t)].
```

There are twelve stages. With (24), this gives

```
W13(t)<=10^1827840 [10^3Y(t)+S12(t)],
12*152320=1827840.                                      (28)
```

No time derivative of S12 is involved in this pointwise algebraic
recurrence. It is valid even when S12 or W1 is zero. No boundary
condition has been assigned to a normal derivative: every boundary
energy field is one of the boundary-preserving fields in (17).
The normal derivatives are recovered solely from the actual equation.

## 8. Full mixed X13 norm, numerical bound, and identification

For the product connection the radial and angular derivative-index
placements satisfy

```
||D_sp^m h||2^2=
 sum_(a+b=m)binom(m,a)||D_Omega^b partial_z^a h||2^2.
```

This commutes only radial derivatives with angular ones, not angular
derivatives with one another. By (9), the angular covariant tensor
norm is bounded by the corresponding natural rotation-word norm.
Each binomial coefficient is at most `2^13`. Summing also over the
time derivatives gives

```
sum_(p=0)^13 ||partial_t^p v||H^(13-p)^2 <=2^13 W13^2.
```

The sum of fourteen Hilbert norms is at most sqrt(14) times their
square sum. Hence

```
X13(v;t)<=sqrt(14)*2^(13/2) W13(t)<=512 W13(t).             (29)
```

Every term on the right side of (28) has now been bounded numerically.
Combining (4), (22), (28), and (29),

```
X13(v;t)
 <512*10^1827840 {
       10^8300003 exp(10^152019)
       +10^8200000 exp(10^152018)}
 <10^10200000 exp(10^152019),                             (30)
```

which proves (1). For an arithmetic check, the larger exponential
dominates both terms, `512<10^3`, and the braces are below
`10^8300004 exp(10^152019)`. The resulting finite prefactor
exponent is `1827840+8300004+3=10127847<10200000`.
There are no unspecified commutator, multiplier, initial-jet, or
normal-recovery constants left in this estimate.

The preceding derivation can first be made for smooth compatible
solutions. To justify its stated finite topology, approximate the first
and second initial data by smooth data vanishing on fixed slightly
smaller open collars. The retained first-wave estimates give convergence
of the first waves in the mixed order-fourteen topology. The mixed
bilinear source estimate, applied to differences, then gives convergence
of F and its initial source jets in the S12 and I11 interfaces.
Every approximating source has the compatible initial zero collars,
as does its second datum. The retained forced linear solver applies.
The same estimates above applied to differences show convergence of
the resulting second solutions in

```
X_13=intersection_(j=0)^13 C^j([0,T];H^(13-j)(M;E)).
```

The limit solves the same forced linear equation, initial data, and
homogeneous boundary problem. The retained actual parameter curve is
C2 into X10, so its second derivative V2 is already defined there
and satisfies (2) with the very same data and source. Energy uniqueness
identifies that derivative with the constructed X13 solution. Thus
(1) is a bound for the actual selected second derivative at epsilon=0;
it is not obtained by assuming the nonlinear curve is C2 into X13.

On the prescribed identical reflected union, use the summed Hilbert
energies and the transported initial/source interfaces for that union.
The same finite constants apply. No assertion is made for an arbitrary
number of copies.

The completed result bounds the full mixed second-wave X13 norm and
the retained second-wave source hierarchy. It does not prove the
second wave in X14, which would need one further source derivative.
It also does not supply a finite-amplitude nonlinear metric patch,
an admitted numerical amplitude radius, higher parameter remainder,
observer embedding and proper clock, actual cut-area, positive
denominator, or a common proper-time observable domain. Those are
distinct subsequent estimates. The enormous constants are analytic
ceilings in the fixed normalization, not physical blast energy or
evidence of a cosmological origin.
