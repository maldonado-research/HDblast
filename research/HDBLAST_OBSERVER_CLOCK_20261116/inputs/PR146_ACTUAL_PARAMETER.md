# Actual parameter derivatives by solution difference quotients

16 November 2026. New proof; no older proof or numerical
admission is changed. This note concerns the actual canonical one-parameter
family of sixteen-field solutions, with the fixed background, fixed density,
fixed shell chart and fixed elapsed-time convention of the retained source.
All norms below are ordinary bundle Sobolev norms with full derivative
arrays on the compact spatial collar. They are not domains of powers of a
mixed-boundary elliptic operator.

The proof constructs the derivatives from actual solution secants. It does
not obtain linearized existence by assuming a smooth nonlinear solution
flow, an X15 estimate, top-norm approximation, or top-norm continuity with
respect to general initial data. Its analytic inputs are stated precisely
in Section 1 and proved in the new companions
[DIFFERENTIATED_INITIAL_JETS.md](DIFFERENTIATED_INITIAL_JETS.md) and
[MIXED_DIFFERENCE_AND_LINEARIZED_ESTIMATES.md](MIXED_DIFFERENCE_AND_LINEARIZED_ESTIMATES.md).
The conclusion depends on those proofs and all retained hypotheses.

## 1. Inputs and resulting theorem

Put

```
Sigma=[0,6ell] x S3, ell=10^-6, T=13ell/8,
X_m(u;s)=sum_(j=0)^m ||partial_s^j u(s)||H^(m-j),
Y_m(J)=sum_(j=0)^m ||J_j||H^(m-j).
```

Precisely, C X_m is the closed Banach subspace of
product_(j=0)^m C([0,T];H^(m-j)) whose entries obey
partial_s u_j=u_(j+1) as distributions for j<m. Its norm is the sum
of the displayed uniform Sobolev norms. Thus all parameter derivatives
below are derivatives in a specified Banach space of compatible actual
temporal-jet tables.

The connection, all product slots, the four/twelve shell projections and
the all-Dirichlet outer projection are fixed. The actual equation is

```
P_v v=F(x,v,Dv),
P_v=partial_s^2-A(x,v)^i D_i partial_s-B(x,v)^ij Hsp_ij,
F(x,v,Q)=C(x,v)+E(x,v)Q+G(x,v)(Q,Q).                 (1)
```

Here Q has all five first-jet directions. The normalized maps include
their actual reference-connection terms, horizontal state-connection
corrections and division by the normalized positive time coefficient.
The original sixteen-order coefficient certificate and admitted principal
positivity bounds are retained.

The initial companion, Section 8, constructs an open normalized evolution
interval J containing the original closed interval
I=[-a_iter,a_iter], a_iter=10^(-10^500209), by the strict retained
initial-pair admission. Its closure lies inside that admission. Choose
a compact interval K with I contained in its interior and K contained
in J. The following four inputs hold uniformly on K. This open-domain
bridge is separate from the quotient argument; it supplies ordinary
two-sided derivatives also at the two endpoints of I.

1. For every a in K the actual family member v_a satisfies (1), the same
   homogeneous chart boundary rows, and v_a in C X14 with a common
   pointwise bound X14(v_a;s)<=M<R/4. No continuity of a↦v_a in C X14
   is assumed. Every temporal jet is an actual distributional derivative.
2. The actual compatible initial table Jcal(a)=(J_0,...,J_14) is C2 into
   Y14 on a neighborhood of K, with finite uniform bounds there. Write
   W_j(a)=partial_a J_j(a) and Z_j(a)=partial_a^2 J_j(a). The exact
   differentiated recurrences and homogeneous differentiated chart rows
   are those in the companion initial-jet proof.
3. The actual nonlinear difference estimate established in the energy
   companion is

   ```
   sup_s X13(v_b-v_a;s) <= K1 Y13(Jcal(b)-Jcal(a))
                          <= K2 |b-a|.              (2)
   ```

   Constants here and below are finite qualitative constants uniform on K.
   They impose no new numerical amplitude ceiling or shorter time slab.
4. The energy companion supplies the order-twelve comparison estimate
   for the already constructed first variations:

   ```
   sup_s X12(w_b-w_a;s)
     <= K3 [Y12(W(b)-W(a))
             +sup_s X13(v_b-v_a;s) sup_c,s X13(w_c;s)]
     <= K4 |b-a|,                                    (3)
   ```

   whenever each w_c solves the full linearized equation, belongs to
   C X12 with bounded X13, and has its inherited homogeneous chart rows.
   This application uses a spare derivative: commuted degree at most
   eleven sections have H2 fields, H1 velocities and L2 accelerations
   almost everywhere. It is not an assertion of an unproved generic
   weak linear existence theorem. Section 7 gives the exact residual
   to which that energy estimate is applied.

Low energy uniqueness for the full linearized equation, in the classes
constructed here, is part of the same energy input. Its direct
justification is also recorded in Section 5.

Under these inputs, the actual canonical evolution is twice continuously
differentiable into C([0,T];X10), with the exact first and second equations
in Sections 5 and 8. More precisely, the quotient construction gives

```
first derivative w_a in C X12,    sup_a,s X13(w_a;s)<infinity,
second derivative z_a in C X11,   sup_a,s X12(z_a;s)<infinity.   (4)
```

The first derivative depends continuously on a in C X12 and the second
in C X11. Thus the proof actually supplies C2 into C X11; C X10 is the
retained conservative public target. The two high bounds in (4) are
pointwise bounds on actual jet tables; their top jets are weakly
continuous, as proved below. Strong C X13 for w or C X12 for z is not
inferred.

The strict-admission argument in the initial companion permits K to
contain the original closed numerical interval in its interior. This
is a qualitative neighborhood for differentiation, not a newly
evaluated numerical amplitude interval. All conclusions asserted for
the physical canonical family retain the original numerical interval.

## 2. Ordinary compactness of the field/velocity pair

The compactness principle used twice below is elementary. Suppose
u_n has actual time derivatives and

```
sup_n,s X_m(u_n;s)<=L,    m>=12.                       (5)
```

Continuity of the individual high jet tables is unnecessary. It suffices
that u_n is continuous in H^(m-1), u_n,s is continuous in H^(m-2),
and that their distributional derivatives have the indicated bounded
high values. These properties hold for the first solution secants and
for the second first-variation secants below.

The derivatives in (5) imply

```
||u_n(t)-u_n(s)||H^(m-1) <= L |t-s|,
||u_n,s(t)-u_n,s(s)||H^(m-2) <= L |t-s|.              (6)
```

For the weaker-continuity version, the high-valued derivatives are
strongly measurable by separability and weak measurability. Integrating
in H^(m-1) and H^(m-2), respectively, proves (6); the resulting identity
agrees with the original Banach-valued distribution identity. This does
not assume strong continuity of a high jet.

At each time, the value sets are relatively compact by the ordinary
Rellich embeddings

```
H^m compactly embeds H^(m-1),
H^(m-1) compactly embeds H^(m-2).
```

Ascoli therefore gives a common subsequence with

```
u_n -> u_0 in C H^(m-1),
u_n,s -> u_1 in C H^(m-2).                            (7)
```

Pass the time integral identity to get partial_s u_0=u_1. Rellich and
Ascoli have only been used for these two actual fields. No unrelated
time jets are interpolated, and no high mixed-boundary spectral domain
is invoked. All qualitative compactness constants are outside the
evaluated numerical admission.

## 3. Finite product-continuity rule for the recovered rows

For later use, expand the actual differentiated rows by the ordered
Leibniz and coefficient chain rules. A typical leaf consists of a
certified coefficient array evaluated along one or more base states,
contracted with finitely many actual field factors and at most the
displayed first- or second-spatial/time derivatives of a variation.
Every spatial derivative is a full array; the product connection is
fixed and time independent, so it commutes with elapsed-time
differentiation. Spatial derivative slots need not be reordered.

Here is the precise continuity rule used below. Suppose the row after
its requested time and spatial differentiation has total field order
at most N, where N is 11, 12 or 13. Suppose all individual temporal
indices of the not-yet-fully-recovered variation are already available
in C H^(n-l), with total available order n>=11; and every leaf factor
of total order above n is instead a base/previous-variation factor
whose corresponding higher strong space has been established. Choose
one largest-order actual field factor for L2. Every other field factor
has total order at most floor(N/2)<=6. For each such factor the H3
embedding is applied to that very temporal jet: its available spatial
regularity has at least three unused derivatives, since 6+3<=11.
Consequently every low factor converges in L-infinity and the chosen
factor converges in L2, uniformly in time.

For example, convergence of D_sp^b u_l in L-infinity is justified only
by convergence of u_l in a space with l+b+3 within its own proven
mixed order. It is never justified by a bound on a different time jet.
The relevant low factor tensor rank, with the three embedding slots,
is at most 2+6+3=11. Coefficient arrays themselves are used in their
certified operator-to-full-array norms, without Sobolev embedding of
high-rank coefficient tensors.

Coefficient compositions converge by the same finite chain rule and
uniform state convergence. Their leaf coefficients converge uniformly,
using one additional vertical derivative and the segment of admissible
states. In the rows below at most eleven row derivatives, at most two
explicit state derivatives and one continuity derivative are needed;
even allowing the C(x,0)=0 extraction when it is used, the requirement
is at most fifteen derivatives, within the retained sixteen-order
certificate. Segment integrals are treated uniformly in their segment
parameter; the same estimates dominate their finite interval integrals.

Telescoping a difference of two leaves proves its convergence: a changed
L2 factor uses its strong L2 difference; a changed low factor uses its
strong infinity difference; a changed leaf coefficient uses its uniform
operator-array difference. All other factors are uniformly bounded in
their assigned norms. Finite sums and the finitely many requested
spatial derivative levels then give C H^k convergence of the row.
Repeated leaves retain their actual multiplicities. This rule supplies
all time, normal and angular derivatives of each source used below;
it is not a spatial-only source estimate.

## 4. Exact first secants and recovery of their temporal jets

Fix a parameter a. For b!=a put

```
h=b-a,       p_ab=(v_b-v_a)/h,
v_theta=(1-theta)v_a+theta v_b,   Q_theta=Dv_theta.
```

The coefficient state patch is convex at these admitted small states.
The segment is used only to evaluate smooth maps, not as a family of
solutions. Exact subtraction of (1) gives

```
P_(v_b) p_ab
 = integral_0^1 [F_V(v_theta,Q_theta)[p_ab]
                +F_Q(v_theta,Q_theta)[Dp_ab]] dtheta
   +integral_0^1 A_V(v_theta)^i[p_ab] dtheta D_i(v_a)_s
   +integral_0^1 B_V(v_theta)^ij[p_ab] dtheta Hsp_ij v_a.       (8)
```

In particular the second derivatives in the last two terms belong to
the base member v_a. The coefficients in P belong to v_b. No equation
has been assigned to an arbitrary auxiliary seed.

By (2) and the initial C2 bound, sup_b,s X13(p_ab;s)<=L1. Apply
Section 2 with m=13 to any sequence b_n→a. After a subsequence, the
pair converges in C H12 and C H11. Also v_b→v_a in C X13 by (2).

Recover the jets of p_ab in order j=2,...,12. Set r=j-2, request up
to k=12-j spatial derivatives, and differentiate (8) r times in s.
The zero-hit principal terms are

```
A(x,v_b)^i D_i (p_ab)_(j-1)
  +B(x,v_b)^ij Hsp_ij (p_ab)_(j-2).                  (9)
```

Every other principal coefficient-hit term contains a p_ab temporal
index at most j-2. In the right side of (8), differentiated first-jet
factors have temporal index at most r+1=j-1; the explicit variation
in A_V[p_ab] and B_V[p_ab] has time index at most r. All base factors
are already known in C X13. Thus the right side never calls on the
unproved temporal jet (p_ab)_j.

The total field derivative order of any expanded row, counting the
displayed principal second derivatives, is at most r+k+2=12.
The Section 3 rule therefore proves convergence of the entire right
side in C H^(12-j), using the already established lower temporal
indices. Define its limit to be w_j. The exact row proves convergence
of (p_ab)_j to w_j in that topology. Starting with the compact pair,
induction gives

```
(p_ab)_j -> w_j in C H^(12-j),    0<=j<=12.           (10)
```

Pass the time integral identities for j<=11 in the nested lower
spaces. They identify w_j=partial_s^j w_0. Hence (10) is convergence
of actual jet tables in C X12, not a collection of unrelated limits.

For weak retention of the thirteenth jet, differentiate (8) eleven
times in s, with no extra spatial derivative. Isolate (9) with j=13:

```
(p_ab)_13=A(x,v_b)^i D_i(p_ab)_12
                 +B(x,v_b)^ij Hsp_ij(p_ab)_11+R_ab.  (11)
```

All factors in R_ab have total order at most thirteen. Any p_ab
factor in R_ab has individual total order at most twelve, since a
principal coefficient has been hit or the term came from the source.
Base factors of order thirteen converge strongly by (2). The same
largest-factor rule, now using (10), gives R_ab→R in C L2. In (11)
the first two differentiated fields converge in C H^-1, because
(p_ab)_12→w_12 in C L2 and (p_ab)_11→w_11 in C H1. Here
H^-1=(H_0^1)^*. A and B converge in C W^(1,infinity), so they act
continuously on H^-1 by duality. It follows that (p_ab)_13 converges
in C H^-1 to an actual derivative w_13=partial_s w_12.

At every fixed s the uniform high bounds and the lower strong limits
identify a unique weak limit in H^(13-j), for each j<=13. The norm
bound survives by weak lower semicontinuity. Every subsequence has
that same weak limit, so no extra time-dependent extraction is needed.
Boundedness plus strong continuity in the lower injective spaces gives

```
w_j in C_weak H^(13-j), 0<=j<=12,
w_13 in C H^-1 intersect C_weak L2,
sup_s X13(w;s)<=L1.                                  (12)
```

For the weak-continuity assertion, test first against the dense dual
subspace on which lower strong continuity suffices, and approximate an
arbitrary high-space dual test using the uniform high bound. This also
gives strong measurability in the separable high spaces. These are the
representatives used in all later energy and source calculations.

## 5. First equation, initial values, rows and identification

Pass (8) using (10). The result is the full linearized equation

```
mathscrL_v w=0,
mathscrL_v u=P_v u-F_V[u]-F_Q[Du]
                  -A_V^i[u]D_i v_s-B_V^ij[u]Hsp_ij v,         (13)
```

at v=v_a, Q=Dv_a. The complete lower derivatives are

```
F_V[u]=C_V[u]+E_V[u]Q+G_V[u](Q,Q),
F_Q[Du]=E Du+G(Du,Q)+G(Q,Du).                        (14)
```

Evaluation at s=0 in (10), and in (11) in H^-1, gives exactly
w_j(0)=W_j(a), j<=13. The initial companion's differentiated
recurrence identifies these as the compatible jets of (13), rather
than a separately assigned high table. Equalities first obtained in
the lower injective spaces identify the same elements in their retained
high spaces. The initial companion defines W through order fourteen;
this evolution proof asserts the attained w_j time traces only through
j=13 and does not construct an evolved w_14.

Each p_ab has the same homogeneous four-Dirichlet/twelve-Neumann
shell rows and all-Dirichlet outer row. The strong C X12 passage
gives the field rows directly, and differentiated Dirichlet rows
through temporal degree eleven and Neumann rows through degree ten
by the ordinary trace maps. For the additional high trace needed by
the next energy argument, take an arbitrary time/rotation word T of
degree l. Its fixed spatial maps send (12) to weak H^(13-l).
For l<=11, this is H2 or better; the bounded H2→H^(1/2) normal trace
passes the homogeneous Neumann row weakly at every time. For l<=12,
the H1 Dirichlet trace passes weakly. The projections and normal
connection are fixed and commute with these time/natural angular
words. This is a passage of traces of actual high representatives.
It asserts no ordinary Neumann trace for a merely H1 top word.

The order-twelve energy actually used below is justified by those H2
Neumann traces and by H1 Dirichlet traces of its velocities, as detailed
in the energy companion, Section 6. Thus no unavailable top ordinary
normal trace or additional weak boundary-existence theorem is required.
The retained scalar principal structure, fixed D/N split, and zero wall
cross-principal traces give its homogeneous boundary flux. No normal
boundary condition is added to a normally differentiated field.

Low energy uniqueness for (13) is direct in these classes. For a
zero-data difference u, its inhomogeneity K_v u is bounded in L2 by
a finite coefficient times ||u||H1+||u_s||L2: v has bounded second
derivatives in infinity norm, and the full first-jet term is first
order in u. Apply the retained mixed energy identity and Gronwall.
The fixed rows and scalar principal wall structure make the boundary
flux zero. Our fields have far more than the H2/H1/L2 regularity needed
to justify this low-order identity almost everywhere. Thus u=0.

The fixed Gamma action, identical-copy relation and parallel tensor
time-reflection act by continuous linear maps in the spaces used here.
Every actual secant has those relations, so its strong lower limit does
as well; the same applies to the second secants below. Compactness is
global on the full compact collar and sphere, retaining every angular
mode. It needs neither a truncated angular basis nor a support claim
for the evolved variations.

Every subsequential limit constructed above has (13), the same initial
pair and rows, and lies in this uniqueness class. Consequently every
sequence b→a has the same limit w_a in C X12. This proves the actual
first derivative; it does not merely exhibit a formal linear solution.

## 6. Continuity of the first derivative with respect to the parameter

For any sequence a_n→a, the w_(a_n) have the uniform high bound (12).
The Section 2 argument applies to their field/velocity pairs, since
their actual derivatives and weak high representatives imply (6).
Take a compact pair subsequence and recover the other jets from (13).
Its base members v_(a_n) converge to v_a in C X13 by (2). Each
recovered row has the same order-twelve product structure as Section 4;
the full F_V/F_Q and base second-derivative terms converge by Section 3.
The initial pairs converge to W_0(a),W_1(a) by initial C2. The rows
pass by the same strong/weak trace argument. The limit therefore
equals w_a by low energy uniqueness. Sequential compactness and
uniqueness give

```
a -> w_a is continuous in C X12.                     (15)
```

This argument needs neither C X13 parameter continuity of w nor
C X14 parameter continuity of v.

## 7. Exact second secants and the order-twelve bound

Now every w_a has been constructed. For b!=a let

```
r_ab=(w_b-w_a)/(b-a),       p_ab=(v_b-v_a)/(b-a).
```

Subtract (13) at a,b without presuming a second derivative:

```
mathscrL_(v_b) r_ab = R_ab,
R_ab=((mathscrL_(v_a)-mathscrL_(v_b))w_a)/(b-a).       (16)
```

An exact integral expression makes every source factor explicit. Put
w=w_a and p=p_ab in the following formula; all unmarked coefficient
derivatives in its integrand are evaluated at (x,v_theta,Q_theta):

```
R_ab=integral_0^1 {
 F_VV[p,w]+F_VQ[w,Dp]+F_QV[Dw,p]+F_QQ[Dp,Dw]
 +A_VV^i[p,w]D_i(v_theta)_s
 +B_VV^ij[p,w]Hsp_ij v_theta
 +A_V^i[p]D_i w_s+B_V^ij[p]Hsp_ij w
 +A_V^i[w]D_i p_s+B_V^ij[w]Hsp_ij p } dtheta.         (17)
```

To verify (17), the principal difference P_a-P_b gives its third
line, whereas differentiating the base-dependent coefficients in
K_v(w)=F_V[w]+F_Q[Dw]+A_V[w]D v_s+B_V[w]Hsp v gives the other
lines. In particular the fourth line is present because the base
second derivatives in K_v also change. At p=w these two principal
cross lines add; neither is optional.

Every source row of total mixed derivative degree at most eleven
has total actual field order at most thirteen. The full-array product
estimate in the energy companion gives, uniformly on K,

```
S11(R_ab) <= K5 [sup_s X13(p_ab;s)]
                 [sup_s X13(w_a;s)],
S_k(f)=sum_(j=0)^k sup_s ||partial_s^j f||H^(k-j).    (18)
```

This controls both the degree-eleven time/rotation forcing and the
degree-ten time/normal/rotation forcing in that energy estimate.
The two high factors in (18) refer to different factors within each
leaf, not to placing two high derivatives in L2 simultaneously: one
largest factor is in L2 and the other factors have the three-slot
infinity reserve. The coefficient difference terms have been written
as segment integrals, so (18) uses no quotient of a merely continuous
coefficient map.

The data of r_ab are (W_j(b)-W_j(a))/(b-a), uniformly bounded by
initial C2. Equation (3), whose residual is exactly (16)-(18), gives

```
sup_b,s X12(r_ab;s)<=L2.                              (19)
```

For each fixed b!=a, r_ab belongs to C X12 because w_a,w_b do; the
additional X13 bounds of those two fields justify the commuted
order-twelve energy identity. Uniformity of those extra bounds after
division by b-a is not needed in that justification. Only the bound
(19) is used for the next compactness passage.

## 8. Second quotient limit and exact quadratic equation

Apply Section 2 with m=12 to any sequence b→a. A subsequence converges
for the pair in C H11 and C H10. Meanwhile p_ab→w_a in C X12 from
the first derivative construction, v_theta→v_a in C X13 uniformly in
theta, and w_a is fixed in C X12. Applying Section 3 to (17) shows

```
R_ab -> Q_v(w_a,w_a) in S10.                         (20)
```

Indeed a row of mixed degree at most ten in (17) has total field
order at most twelve. Every p_ab and w_a factor of that order
converges strongly in its own C X12 space; all other factors have
the infinity reserve. This verifies the full time/normal source
passage, including both principal cross lines and every mixed
state/first-jet term.

Recover r_ab temporal jets j=2,...,11 from (16), differentiating
r=j-2 times in time and up to 11-j times in space. Principal terms
use r_ab indices j-1,j-2. The lower linearized terms use indices at
most j-1. The differentiated source uses (20). Every row has total
field order at most eleven for its unknown variation. The Section 3
argument therefore gives convergence in C H^(11-j). Time integral
identities identify all limits as actual derivatives of a field z_a.
Thus

```
r_ab -> z_a in C X11 along the extracted subsequence. (21)
```

For the twelfth jet, differentiate (16) ten times in time and isolate
A(v_b)D(r_ab)_11+B(v_b)Hsp(r_ab)_10. These converge in C H^-1.
All remaining r_ab factors have individual total order at most eleven.
The base fields have their stronger C X13 convergence, and the
remaining source converges in C L2 by (20). Hence (r_ab)_12 converges
in C H^-1 to the actual derivative z_12. Weak lower semicontinuity
and (19) give the pointwise top representatives

```
z_j in C_weak H^(12-j), 0<=j<=11,
z_12 in C H^-1 intersect C_weak L2,
sup_s X12(z;s)<=L2.                                  (22)
```

The same density argument used for (12) proves these weak continuity
and measurability statements. The homogeneous field rows pass
strongly. Weak top traces give differentiated Dirichlet rows through
degree eleven and Neumann rows through degree ten; no nonexistent
ordinary trace is assigned above these thresholds.

Pass (16)-(17) to obtain

```
mathscrL_v z = Q_v(w,w),
Q_v(w,w)= F_VV[w,w]+F_VQ[w,Dw]+F_QV[Dw,w]+F_QQ[Dw,Dw]
 +A_VV^i[w,w]D_i v_s+B_VV^ij[w,w]Hsp_ij v
 +2A_V^i[w]D_i w_s+2B_V^ij[w]Hsp_ij w.               (23)
```

The terms linear in z remain in mathscrL_v. They are not duplicated
on the right side. For the literal quadratic source,

```
F_VV[w,w]=C_VV[w,w]+E_VV[w,w]Q+G_VV[w,w](Q,Q),
F_VQ[w,Dw]=E_V[w]Dw+G_V[w](Dw,Q)+G_V[w](Q,Dw),
F_QV[Dw,w]=F_VQ[w,Dw],
F_QQ[Dw,Dw]=2G(Dw,Dw).                               (24)
```

These identities use the full five-direction first jet. The actual
combined-slot symmetry of G permits the familiar combined factors,
but the displayed equations retain both slots explicitly.

Evaluation at time zero identifies z_j(0)=Z_j(a) in every lower
strong space through j=11 and at j=12 in H^-1; injectivity identifies
the same retained high-space elements. The differentiated initial
recurrence identifies the compatible table through the evolved orders
j<=12. Although the initial companion defines Z through fourteen,
this proof constructs neither evolved z_13,z_14 nor their time traces.
Low energy uniqueness
for the inhomogeneous equation (23) follows by applying Section 5 to
the difference of two candidates. Consequently every subsequential
limit in (21) is the same z_a, and the entire second quotient converges.
This identifies z_a as the actual derivative of a↦w_a in C X11.

## 9. Continuity, endpoints and Taylor remainder

For a_n→a, the z_(a_n) have the uniform high bound (22). Their pair
is compact in C H11/C H10. Their equations are (23). The base fields
converge in C X13 by (2), and the first variations converge in C X12
by (15). Every source term of mixed degree at most ten has total
order at most twelve, so (23)-(24) converge in S10 by the same product
rule. Recover the temporal jets exactly as in Section 8. Initial
data converge by initial C2 and the boundary rows pass at the stated
trace levels. The unique limit is z_a. Therefore

```
a -> z_a is continuous in C X11.                     (25)
```

It follows that the actual map a↦v_a is C2 into C X11, and in
particular into the requested C X10. On a compact admitted parameter
interval its Banach-valued Taylor formula has the uniform remainder

```
||v_(a+h)-v_a-h w_a-(h^2/2)z_a||C X10
 <= (h^2/2) sup_(|c-a|<=|h|) ||z_c-z_a||C X10
 = h^2 o(1),                                        (26)
```

where the segment stays in the admitted neighborhood. Formula (26)
comes from the two proved derivatives and their continuity; it is not
an assumed nonlinear-flow remainder. The proved open normalized
extension makes these ordinary two-sided derivatives at every point
of the original closed interval. Its restriction has the same continuous
one-sided endpoint derivatives if only secants within that interval
are considered. No third parameter derivative or numerical bound for
an observable remainder is added.

## 10. Scope

The construction concerns actual normalized fields at the same base
point and the same elapsed-time member. It also proves the actual
fixed-base reconstructed field curve is C2 in C X10. To verify this
composition step, use the exact physical point map X=Psi(x,v) in the
retained `inputs/COEFFICIENTS.md`, its literal chain rule, and its
all-point building-block envelope (Psi through eighteen derivatives).
The admitted uniform state margin gives an open neighborhood of each
v_a in C X10 on which this composition is defined. Expand at most ten
mixed derivatives of Psi and its first two state differentials. Every
leaf has total actual field/direction order at most ten. One largest
factor is placed in L2 and every other factor, of order at most five,
uses its own three-spatial-derivative infinity reserve within X10.
The full tensor, base and state-input slots are retained; coefficient
arrays use their operator envelope. At most ten mixed derivatives,
two state variations and one continuity derivative are required, a
total of thirteen within the retained eighteen. Finite coefficient
Taylor remainders and telescoping products therefore prove that the
Nemytskii map u↦Psi(x,u) is C2 on this open subset of the closed
compatible-jet Banach space C X10. Its ordinary Banach chain rule with
the actual constructed curve gives

```
X_a=Psi_V(x,v_a) w_a,
X_aa=Psi_V(x,v_a) z_a+Psi_VV(x,v_a)(w_a,w_a).         (27)
```

These are identified derivatives, with the same fixed-base meaning as
the normalized result. A parameter-dependent observer, worldline,
clock, source/detector/noise model or calibrated response is a separate
construction and has not been supplied here.

No X15 estimate, nonlinear C-infinity solution-flow theorem, strong
C X14 approximation theorem, top-norm data continuity, new numerical
constant, shorter time slab, changed shell normal or changed canonical
normalization enters the proof. The canonical initial curve and the
two lower comparison estimates are the substantive inputs; formal
coefficient identities alone do not replace them.
