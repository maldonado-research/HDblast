# Complete mixed source hierarchy for the second background wave

This conditional continuation uses the frozen first-wave
HIGH_ORDER theorem and the complete quadratic source from QUADRATIC_FORCING.
This note supplies the source norms required at second-wave order thirteen.
It does not assume first-wave order fifteen, an evolved W13,infinity bound,
or C2 dependence of the nonlinear family into a stronger topology.

## 1. Fixed-bundle equation and the actual quadratic source

All norms use the fixed product connection on

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D,
M=[0,Z]xS3, Z=6*10^-6,
dmu=dz dmu_kappa,
T=13/8000000.
```

The symmetric block has the full Frobenius norm. The time direction is
product-parallel and commutes with spatial covariant differentiation.
The symbol D denotes this fixed spacetime product connection, including
the derivative-index connection in its Hessian. Moving-frame coefficients
such as sech(t) are inside the previously derived coefficient tensors;
they are not silently treated as constants in a moving first-jet basis.

Write the exact transformed equation before changing its sign as

```
q(V):D^2V=c(V)+b(V)DV+d(V)(DV,DV),  c(0)=0.
```

These are the actual fixed-bundle coefficients after the complete
coordinate-to-bundle substitution P=DV-Omega V, not the old coordinate
coefficients with their first derivatives relabeled. Let u=V1 be the
first amplitude derivative and v=V2 the second amplitude derivative,
with Taylor convention V=e u+e^2 v/2+o(e^2). In the retained ordered
positive-energy sign convention,

```
L0 v=F2[u],
F2[u]=-C2(u,u)-2E2[u]Du-2D2(Du,Du)+2Q1[u]:D^2u,         (1)
C2=c_VV|0, E2=b_V|0, D2=d|0, Q1=q_V|0.
```

The full principal variation is present, including its time-time Hessian
part. Both coefficient-variation terms have their factor two. The d-term
also has factor two because d was defined as one half of the first-jet
Hessian. No q_VV term survives: it would multiply D^2V0=0. Likewise
d_V terms multiply a background first jet that vanishes.

The retained coefficient proof bounds all mixed time/spatial covariant
derivatives of C2,E2,D2,Q1 through total order twelve by

```
C=10^152000.                                            (2)
```

These are full Hilbert multilinear operator norms with every base and
bundle slot included. Twelve external derivatives plus at most two state
derivatives use raw mixed order at most fourteen, within the pinned raw
order sixteen certificate. The separate connection conversion also fits
its recorded derivative count. The full varying bulk reference connection
and all boundary-chart terms remain in these maps.

## 2. The first theorem controls the full mixed X14 norm

The source first-wave theorem proves the stronger intermediate quantity

```
W14(t)^2=sum_(a+p+|I|<=14)
                  ||partial_z^a partial_t^p R^Iu(t)||L2^2,
W14(t) <=10^1980163 Y13(t),
Y13(0)<10^1985300,
Y13(t)<=Y13(0) exp[(10^150006/2+95*10^152020)t].          (3)
```

R^I ranges over all ordered natural SO(4) rotation words, including their
action on tensor slots. The full product spatial derivative decomposition
and the proved natural-rotation inequality give, for each p<=14,

```
||partial_t^p u||H^(14-p) <=2^7 W14.
```

Thus define the actual mixed first-wave norm

```
M14(t)=sum_(p=0)^14 ||partial_t^p u(t)||H^(14-p).
```

It obeys

```
M14(t)<=1920 W14(t)
      <10^3965467 exp[(10^150006/2+95*10^152020)t]
      <B1:=10^4000000 exp(10^152017),  0<=t<=T.           (4)
```

This is a consequence of the mixed W14 theorem, not an inference from
the displayed N14=H14xH13 consequence alone.

At zero, the same source proves its recursive temporal jets satisfy
`max_(p<=14)||partial_t^p u(0)||H^(14-p)<10^1985230`.
Consequently

```
M14(0)<15*10^1985230<10^1985300.                          (5)
```

All these initial temporal jets vanish on the actual fixed open boundary
collars. The estimates refer to derivatives at the zero amplitude; they
do not produce a uniform nonlinear family radius.

## 3. Literal mixed derivative allocations

Let J_0u=u, J_1u=Du and J_2u=D^2u in the full fixed spacetime product
metric. Each term of (1) is a coefficient map K applied bilinearly to
`J_r u,J_s u`, with intrinsic derivative orders

```
(r,s)=(0,0), (0,1), (1,1), or (0,2).                    (6)
```

For a time derivative p and a spatial derivative tensor of order j,
with p+j<=12, its exact Leibniz expansion has allocations

```
p_K+p_1+p_2=p,
j_K+j_1+j_2=j,
```

and coefficients

```
multinomial(p;p_K,p_1,p_2)
 *multinomial(j;j_K,j_1,j_2).                            (7)
```

The spatial formula includes the usual derivative-slot shuffle maps,
which are isometries for the Hilbert tensor norms. Angular covariant
derivatives are not reordered past one another. Time derivatives commute
with spatial covariant derivatives because the fixed product connection
has no time/spatial curvature or time-dependent frame.

The coefficient derivative has total order p_K+j_K<=12. The two field
factors have total differential orders

```
d_1=p_1+j_1+r, d_2=p_2+j_2+s,
d_1+d_2<=p+j+r+s<=14.                                  (8)
```

It follows that one factor has order at most seven and the other has
order at most fourteen. This is the exact derivative allocation that
avoids order fifteen. In particular the principal term with all twelve
external derivatives on D^2u uses precisely order fourteen of u.

The full spacetime jet norms satisfy

```
||D_sp^b partial_t^a J_r u||L2 <=2 M14,
||D_sp^b partial_t^a J_r u||H3 <=2 M14
                         if a+b+r+3<=14,                (9)
```

for r=0,1,2 and a+b+r<=14 in the first line. To check the factor two,
decompose the r intrinsic derivative slots into l time slots and r-l
spatial slots. Their binomial multiplicities sum to 2^r<=4. For the H3
version, for each fixed l the four spatial derivative orders are a subset
of the H^(14-a-l) norm of partial_t^(a+l)u; they are already included in
that norm and do not cost an additional derivative-order factor.

The established fixed-cylinder tensor embedding is

```
||T||L-infinity <=10^7 ||T||H3.                          (10)
```

It is valid through tensor rank sixteen. The lower-order factor in (8)
has original plus derivative rank at most 2+7=9, and even after its three
additional spatial derivatives the rank is at most twelve. Thus (10)
applies to that factor. The higher-order factor is only put in L2. No
coefficient tensor is passed to (10), and no W13,infinity assumption on
the evolving wave is introduced.

Combining (8)-(10), each allocation in (7) has L2 norm at most

```
4*10^7 C M14(t)^2.                                     (11)
```

The sum of all products of the multinomial coefficients is
`3^p*3^j<=3^12`. This proves the same finite estimate for every mixed
time/spatial derivative needed below.

## 4. Complete mixed source norm and initial source jets

Use literal sums of Sobolev norms, not a hidden maximum or a different
quadrature convention:

```
S12(t)=sum_(p=0)^12 ||partial_t^p F2(t)||H^(12-p),
I11=sum_(p=0)^11 ||partial_t^p F2(0)||H^(11-p).             (12)
```

The sum of the absolute coefficients of the four terms in (1) is seven.
For each fixed p, the square sum over its at most thirteen spatial
derivative orders costs at most sqrt(13)<4. Summing the thirteen values
of p, (11) gives the explicit homogeneous bound

```
S12(t) <=[13*4*7*4*3^12*10^7] C M14(t)^2
        <=10^20 C M14(t)^2
         =10^152020 M14(t)^2.                            (13)
```

All estimates in (9)-(13) are non-strict, including when the field
vanishes. There are no discarded mixed-time allocations or scalar-only
component estimates in (13).

Substituting (4) produces the numerical mixed hierarchy throughout the slab:

```
S12(t)<10^8152020 exp(2*10^152017)
       <Bsrc:=10^8200000 exp(10^152018).                  (14)
```

At zero, use (5), which has no exponential evolution factor. Monotonicity
of the Sobolev norms gives I11<=S12(0), and therefore

```
I11<10^(152020+2*1985300)=10^4122620<10^4200000.           (15)
```

Only external order eleven is required in I11; bounding it by S12(0)
uses one extra already-certified order and introduces no new datum.
Because all first-wave temporal jets vanish on the initial open boundary
collars and (1) is local and quadratic, all these initial source jets
also vanish there. This is the actual corner support statement needed
by the inhomogeneous second-wave temporal recursion.

## 5. Natural rotation-word source norms for the forced energy

Define the source family required by time/rotation energy of order twelve:

```
G12(t)^2=sum_(p+|I|<=12)||partial_t^p R^I F2(t)||L2^2.    (16)
```

For any E section T the frozen natural-rotation comparison is
`A_k(T)<=60^k||T||Hk`, where A_k already sums all ordered six-letter
words of exact length k. Apply it to the completed E-valued source
T=partial_t^p F2. Its rank remains at most two. Since p+k<=12,

```
G12(t)^2<=13*60^24*S12(t)^2,
G12(t)<=10^25 S12(t).                                   (17)
```

The constant follows from `13*60^24<10^50`. There is no further 6^k
multiplier after A_k has aggregated the words. This route does not
differentiate the principal Hessian coefficient using raw component
angular derivatives. Its covariant derivatives are charged in (2) and
(7); the natural Lie derivatives of the complete source are then
controlled by the proved geometric norm comparison. No false assertion
that component derivatives commute with the operator is needed.

For normal recovery it is also useful to define, for m<=12,

```
Z_m(F)^2=sum_(a+p+|I|<=m)
                 ||partial_z^a partial_t^p R^I F||L2^2.
```

An individual row is bounded by
`60^12||partial_t^p F||H^(m-p)` because the radial derivative is
parallel and commutes with natural rotations. The total number of rows
is at most `sum_(k=0)^12 binom(14-k,2)6^k<10^12`. Therefore the
conservative full-array interface is

```
Z_m(F2;t)<=10^30 S12(t),  m<=12.                         (18)
```

This additionally includes arbitrary time/radial/rotation source words
needed in the triangular second-wave normal recurrence.

## 6. Normalized forcing and the exact retained X13 source cost

When solving the second-wave equation for its next time or radial
derivative, the source is multiplied by rho^2. The background factor is
time-independent, rotation-invariant, and has spatial derivatives through
twelve below 10^210. The linear multiplier estimate with
`2^m sqrt(m+1)<10^5`, m<=12, gives

```
sum_(p=0)^m ||partial_t^p(rho^2 F2)||H^(m-p)
 <=10^215 sum_(p=0)^m||partial_t^p F2||H^(m-p)
 <=10^220 sum_(p=0)^m||partial_t^p F2||H^(m-p),           (19)
```

for m<=12. The same bound applies to the initial hierarchy. Combining
(18)-(19), one may use

```
Z_11(rho^2 F2;t)<=10^250 S12(t),
sum_(p=0)^11||partial_t^p(rho^2 F2)(0)||H^(11-p)
                         <10^4200220<10^4200300.         (20)
```

These are enlarged numerical interfaces for the forced recurrence; they
are not additional unknown constants.

At m=13 the retained ENERGY_MAJORANT equation (20) requires exactly

```
Csource13 = sum_(p=0)^12
       ||partial_t^p F2||_(L1([0,T];H^(12-p)))
  +sum_(p=0)^11||partial_t^p F2(0)||H^(11-p).
```

The first term is the integral of the literal sum S12. Thus all of that
source cost, including initial source jets, is now evaluated:

```
Csource13 <= integral_0^T S12(t)dt+I11
          <T Bsrc+10^4200000
          <10^8300000 exp(10^152018).                    (21)
```

Here T<1. The same calculation gives
`||G12||L1_t<=10^25 T Bsrc` for the direct time/rotation energy route.
It is unnecessary to demand a first-wave H15/H14 extension to use these
order-thirteen source norms.

## 7. Identification with the parameter second derivative

The retained nonlinear theorem states C2 parameter dependence into X10,
where `X_s=intersection_(j=0)^s C^j_t H^(s-j)`. It does not state C2
dependence into X13 or X14. The present source calculation must preserve
that distinction.

The previously derived second differentiation identity identifies the
actual second jet v_actual in X10 as a solution of `L0v=F2[u]`, with
the registered second Cauchy data and the exactly differentiated homogeneous
boundary rows. Equation (1) is meaningful in a lower distributional/energy
topology under those retained premises. The first-wave theorem identifies
its stronger X14 solution with the actual first parameter derivative by
linear uniqueness; that supplies the u used in this source calculation.

If the forced linear construction using (12)-(21) gives v_strong in X13
with those same initial and boundary data, then v_actual-v_strong belongs
to the common lower energy class, solves the homogeneous background
problem with zero data, and vanishes by the retained energy uniqueness.
Consequently the particular second derivative at e=0 has the stronger
regularity and bound. This argument upgrades that jet by identification;
it does not upgrade the topology in which the whole parameter curve is C2.
The existence and quantitative forced-wave estimate are the companion
proof's task, not inferred from the size of the source ceiling alone.

## 8. Exact replay ledger and scope

The exact verifier should enumerate every tuple in (7) for p+j<=12 and
each intrinsic pair in (6), checking

```
p_K+j_K<=12,
d_1+d_2<=14,
min(d_1,d_2)<=7,
min(d_1,d_2)+3<=14,
max(d_1,d_2)<=14.
```

It should independently verify the time multinomial expansion of all four
quadratic source terms, preserving the principal Hessian term and all
factors two. Ordinary spatial tensor shuffles have the same multinomial
count; their analytic norm estimate is the proof above, not a sampled
identity. Useful exact arithmetic checks are

```
15*2^7==1920<10^4,
1980163+1985300+4==3965467<4000000,
15*10^1985230<10^1985300,
13*4*7*4*3^12*10^7<10^20,
152020+2*4000000==8152020<8200000,
2*10^152017<10^152018,
152020+2*1985300==4122620<4200000,
13*60^24<10^50,
sum_(k=0)^12 binom(14-k,2)*6^k<10^12,
10^12*60^24<10^60,
2^24*13<10^10,
4200000+220==4200220<4200300,
2*10^8200000<10^8300000.
```

All coefficient and first-wave inputs should be bound to their immutable
source snapshots. Replaying arithmetic and finite polynomial identities
does not replace the analytic product proof or retained existence theorem.

This completes the mixed source hierarchy at the derivative counts required
for the background second-wave X13 estimate. It does not claim the source
hierarchy for second-wave X14: its total external order thirteen would
allow a first-wave Hessian derivative of total order fifteen. It also does
not quantify a nonlinear-family amplitude interval, metric patch, observer
clock, actual cut area, denominator, or the common-domain observable remainder.
