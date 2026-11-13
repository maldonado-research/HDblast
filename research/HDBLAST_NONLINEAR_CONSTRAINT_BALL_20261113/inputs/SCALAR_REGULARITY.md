# Scalar weighted H16 and physical initial metric H14 bounds

6 November 2026. A conditional initial-data estimate for the explicitly
specified weight member. All numbers are in the retained model
normalization. Each field carrying subscript2 is the second amplitude
derivative, twice the coefficient of the amplitude squared.

## Inputs, conventions and conclusion

Use the retained exact-family construction and canonical conventions in
[EXACT_INITIAL_THEOREM.md](inputs/EXACT_INITIAL_THEOREM.md),
[REGISTERED_EXACT_FAMILY_BINDING.md](inputs/REGISTERED_EXACT_FAMILY_BINDING.md),
and [CANONICAL_CONSTRAINTS.md](inputs/CANONICAL_CONSTRAINTS.md).
The actual scalar source and prior weighted comparison are fixed in
[PRIOR_ACTUAL_SOURCE.md](inputs/PRIOR_ACTUAL_SOURCE.md) and
[PRIOR_SCALAR_OUTPUT.md](inputs/PRIOR_SCALAR_OUTPUT.md).
The background potential is the exact polynomial U(eta) in
[ACTUAL_POTENTIAL_JET_BOUNDS.json](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json),
eta=phi−1, and

```
r=(z−ell/2)(5ell/2−z)/ell², ell=10^-6, kappa=10^-8,
q=kappa²r², w=exp(−1/r), dy=−rho dz,
129<=rho<=131, |a|,|p|,|eta|<1.
```

The continuous proof in [SCALAR_COERCIVITY.md](SCALAR_COERCIVITY.md)
gives `||N||_w<9e128` for the actual selected scalar normal potential.
The source-jet calculation in [SOURCE_JETS.md](SOURCE_JETS.md) gives the
aggregate transformed scalar source bound G<3e25 used below.

Define the probability-sphere covariant norms

```
||N||_(H16_q,w)² = sum_{j=0}^16 integral w² q^(2j)|nabla_gamma0^j N|² dmu0,
||h||_(H14_gamma0)² = sum_{j=0}^14 integral |nabla_gamma0^j h|² dmu0.
```

For the three unit-harmonic coefficients define s by ds/dy=1/q and

```
M=sqrt(q rho³) w, v_l=M n_l,
g_l=sqrt(q rho³) w^-1 f_l,
G=(sum_l ||g_l||_(H12_s)^2)^(1/2),
||u||_(Hm_s)^2=sum_{j=0}^m ||D_s^j u||_L2(ds)^2.
```

The resulting explicit estimates are

```
(sum_l ||v_l||_(H16_s)^2)^(1/2) <=20 (G+9e128),              (1)
||N||_(H16_q,w) <=1e2000 (sum_l ||v_l||_(H16_s)^2)^(1/2),   (2)
||h2||_(H14_gamma0) <=1e2000 (sum_l ||v_l||_(H16_s)^2)^(1/2). (3)
```

Moreover h2 belongs to H14_0(A), and its extension by zero has the same
H14 integral norm on the ambient initial slice. This conclusion uses a
bounded map and compact approximation, not an assumed potential trace.

The actual source bound G<3e25, obtained from the three per-mode H12
derivative-sum bounds, and (1) give v-H16 less than2e130. Therefore
(2),(3) give

```
||N||_(H16_q,w)<2e2130,  ||h2||_(H14_gamma0)<2e2130.          (4)
```

The preceding estimates concern one bulk copy and probability sphere
measure. Full geometric sphere volume multiplies each norm by
`pi sqrt2<4.45`, so the following common per-copy ceilings hold in both
the probability and full geometric conventions:

```
||N||_(H16_q,w)<1e2131,  ||h2||_(H14_gamma0)<1e2131.
```

The transformed v-H16 ceiling2e130 uses unit harmonics in probability
sphere measure. A change to full sphere measure multiplies that norm
by the same `pi sqrt2`.

For two identical reflected bulk copies, squared integral norms add.
The union norm is therefore `sqrt2` times the corresponding per-copy
norm. Applied to the tighter probability bounds(4), the combined factor
`(pi sqrt2)*sqrt2=2pi<44/7` gives a full geometric two-copy ceiling
`2e2131` for both norms. This factor describes two bulk integrals; it is
not an additional estimate for an intrinsic shell field.

These enormous upper budgets describe the actual selected initial
second amplitude derivative. They are not a predicted magnitude or
physical energy. The numerical weights specify a new member; they do
not identify an earlier member whose defining function and weight ratio
were unspecified.

## Correct unitary transform and exact normal operator

The y interval is finite and q vanishes quadratically at both endpoints,
so s maps it diffeomorphically onto the whole line. The factor rho^(3/2)
in M is essential: the actual radial measure is rho³dy. In particular

```
||v_l||_L2(ds)=||n_l||_w.
```

Write b=1/rho, D=D_s, d=q_y, c=D_s log M and

```
c=d/2+(3/2)qa+kappa²r_y,
R_t=3b²+p²/3−3a².
```

Using q²D_y²=(qD_y−q_y)(qD_y), the four transformed scalar-adjoint
radial rows, without square-root angular normalization, are

```
A=−(3/2)aq(D−c)+q²[3a²+(lambda−3)b²]/2,
B=−D²/2+(c+d/2−aq)D−(c²+dc−D c)/2+aqc
  +q²[lambda b²/3−R_t/2],
C=qb(D−c−aq),
T=q²b².
```

For a unit scalar harmonic E with −Delta E=lambda E the normal equation
is exactly, in flat ds distributional pairing,

```
L v=g,
L=A* A+3 B* B+(lambda/2)C* C+[lambda(lambda−3)/6]T* T,       (5)
```

where star means flat formal adjoint and lambda=0,8,24. This follows
also by transforming the public energy identity: each q²P*g row is
conjugated by M, and the source pairing integral n f rho³dy becomes
integral v g ds. Its leading term is `(3/4)D^4`. The D³ coefficient
vanishes by self-adjointness; the estimates below need not exploit this.

For h2=w²q⁴P*gN, its normalized radial/angular components are

```
h2=H{(Av)E dy²+(Bv)E gamma_tan
     +(Cv)/2 [dy tensor dE+dE tensor dy]+(Tv)/2 T_E},
H=w q^(3/2)rho^(-3/2), T_E=Hess_S3 E+(lambda/3)E kappa.
```

Here dE and T_E are lifted to normalized angular components, so their
norms are those of the unit sphere; the displayed gamma_tan denotes the
identity angular block. Thus no rho or q Jacobian has been discarded.

## Explicit coefficient bounds, without numerical sampling

Set x=(z−3ell/2)/ell, so r=1−x², and use variables x,b,a,p,eta,
all of absolute value at most one. D_s acts on polynomials by the vector
field

```
D_s x=−10^-10(1−x²)² b,
D_s b=−q a b,
D_s a=q(−b²−p²/3),
D_s p=q[U'(eta)−4ap],
D_s eta=q p.
```

Each component has polynomial coefficient l1 norm at most1e−8 and
polynomial degree at most9. If a polynomial has coefficient norm A and
degree d0, its kth D_s derivative therefore has norm at most

```
A product_{i=0}^{k−1} [1e−8(d0+8i)].                       (6)
```

This uses the coefficient l1 norm of a product being at most the product
of the norms, and the total norm of all partial derivatives being at most
d0 times the original norm. It holds globally including arbitrarily near
either endpoint. Negative powers of r do not occur in c,d or the rows.

The exact replay constructs the four rows as rational
polynomials, applies (6) through order16, and expands the formal adjoint
composition by the exact identities

```
(a_i D^i)*=sum_{k=0}^i (−1)^i binom(i,k)(D^(i−k)a_i)D^k,
D^k(a_j D^j)=sum_{t=0}^k binom(k,t)(D^t a_j)D^(k−t+j).
```

Leibniz bounds then control all twelve derivatives of all lower
coefficients in `L=(3/4)D^4+sum_{j=0}^3 C_j D^j`. For each lambda,
the largest resulting coefficient-jet envelope is exactly
`15000003/5000000000000000=3.0000006e−9`, in particular

```
max_{j<=3, k<=12} ||D^k C_j||_infinity <1e−6.              (7)
```

Order16 row arrays exceed the order14 needed by the second-order
formal-adjoint composition and twelve coefficient derivatives. The public
reproduction instructions are in [REPRODUCE.md](REPRODUCE.md).

## Global regularity and numerical estimate

The selected normal equation gives (5) distributionally against compact
smooth s tests, since these correspond to compact smooth annular tests.
The low-order coercivity input ensures v is in L2(R). There is no missing
boundary condition at infinity.

Here is an elementary global bootstrap which avoids assuming the target
H16 estimate. Let L0=1+(3/4)D^4. Its Fourier inverse maps H^t(R) to
H^(t+4)(R) for every real t. If v is in H^r with integer0<=r<=15,
the lower part Rv=sum C_jD^jv belongs to H^(r−3). Bounded coefficient
jets through12 suffice: when r−3 is nonnegative, this is Leibniz;
when negative (only −3,−2,−1), it follows by duality from the positive
integer multiplier estimate. Therefore

```
v=L0^-1(g+v−Rv) is in H^min(16,r+1).
```

Iterating from r=0 proves H16 globally. The equality is valid in tempered
distributions; L0 has strictly positive Fourier symbol and no tempered
kernel. This argument supplies membership, with no numerical elliptic
constant used in (1).

Now differentiate (5) twelve times. With
`K=4*2^12*1e−6=.016384`, (7) gives

```
(3/4)||D16 v|| <=||D12 g||+K max_{0<=j<=15}||D^j v||.
```

Plancherel and |xi|^j<=1+|xi|^16 imply
`||D^jv||<=||v||+||D16v||`. Since K<3/4, absorption yields

```
||D16v||<=2||g||H12+||v||,
||v||H16<=20(||g||H12+||v||L2).
```

The final20 is loose; sqrt17 times the preceding bounds is smaller.
Apply the direct-sum triangle inequality to the three modes and the
potential L2 bound in [SCALAR_COERCIVITY.md](SCALAR_COERCIVITY.md) to prove (1).

## Radial physical and weighted-potential conversion

Put e=D_s log H=kappa²r_y+(3/2)d−(3/2)qa. For any row u,

```
D_y^j(Hu)=H q^-j F_j,
F_0=u, F_(j+1)=(D_s+e−j d)F_j.                            (8)
```

The coefficient-polynomial norm of a differential operator P with total
coefficient norm S and maximum coefficient degree dP increases under
this recurrence by at most the factor
`1+1e−8 dP+||e−jd||_poly`; degrees increase by at most8.
The exact replay starts with all four rows and proves that the sum of
their resulting coefficient norms is below10 through j=14. Consequently
the L2(ds) sums of F_j are bounded by a fixed number below100 times the
H16 norm of v. Crucially this accounts for differentiating the row itself.

Changing measure in (8) cancels rho exactly:

```
||D_y^j(Hu)||_L2(rho³dy)<=sup_A(w q^(2−j)) ||F_j||_L2(ds).
```

For j<=14, the supremum is at most1e226. For j>=3 this follows from
`sup_{t>=1} e^-t t^(2j−4) <=(2j−4)^(2j−4)`,
kappa^(4−2j)<=1e192, and24^24<1e34. For j<=2 the bound is immediate.
Thus, allowing the finite row and mode sums, every component radial
jet through14 has L2 norm below `1e232 ||v||H16`.

For the potential, the analogous identity is

```
M q^j D_y^j n=P_j v,
P_0=1, P_(j+1)=(D_s−c−j d)P_j.
```

The exact polynomial replay bounds each coefficient sum by2 through16.
This gives the radial weighted bridge without any unweighted potential
trace or unweighted potential integrability assumption.

## Covariant conversion with explicit angular and connection budgets

Radial component derivatives alone are not a physical H14 norm. The
following crude independent budgets account for the missing conversions.

1. Work in the bundle identified with R dy plus the unit-sphere tangent
   bundle, using normalized angular tensor components. Its reference
   connection is D_y radially and nabla_S3 angularly. The physical
   connection acts radially as D_y; in each normalized angular direction
   it is `b nabla_S3+a K`. K is the algebraic map exchanging a radial
   index and an angular index in the usual warped-product connection.
   At tensor rank R its direct-sum operator norm is at most2R sqrt3.
   The same description applies to derivative indices. K is parallel
   for the reference sphere connection. For all ranks used here R<=18,
   so a ceiling100 suffices, with any finite component sum covered by
   the term-count budget below.

2. In variables b,a,p,eta the unscaled D_y vector field is
   `(-ab,−b²−p²/3,U'−4ap,p)`. Its components have coefficient l1 norm
   at most100 and degree at most5. The same polynomial argument gives
   every D_y derivative through18 of a or b norm at most
   `product_{i=0}^{17}100(1+4i)<1e80`. This is a direct y-derivative
   bound; it is not inferred by dividing an s estimate by q.

3. All source harmonics are explicit ambient polynomial restrictions of
   degree at most4. Unit normalization uses E0=H0,
   E2=sqrt48 H2, E4=sqrt(720/19)H4. Their total ambient polynomial
   coefficient norms are below100. The sphere covariant derivative of
   a tangential rank R polynomial tensor is its ambient derivative
   projected with P=I−n tensor n on R+1 slots. The induced column-sum
   coefficient norm of each P is at most5. Hence one derivative multiplies
   the total coefficient norm by at most dP*5^(R+1), and increases its
   polynomial degree by at most2(R+1)−1. Iterating up to rank18 gives
   a norm below1e200. This also covers the bases dE and
   T_E=Hess E+(lambda/3)E kappa and their derivatives through16;
   the extra trace term costs less than100 and fits within1e200.
   Polynomial coefficient l1 bounds dominate pointwise tensor norm,
   hence also probability-sphere L2 norm.

4. Expand k physical covariant derivatives using item1 and Leibniz,
   leaving radial derivatives on scalar components, sphere derivatives
   on the fixed angular bases, and y derivatives on factors a,b.
   At each stage k<=16 radial Leibniz differentiation has at most17
   branches (at most16 scalar factors plus the primary coefficient).
   Intrinsic sphere differentiation is one tensor-valued branch. The
   connection-mixing branch has operator norm below100, already including
   the tensor slots and angular directions. Thus the total norm-weighted
   branching factor is below120; a deliberately loose1000 multiplier
   bounds it at every stage.
   Thus there are at most1000^k weighted terms. Each has at most k
   factors among the y jets in item2. Each angular factor is bounded
   by item3. The tensor rank and all derivative orders remain within
   the ranges already checked.

For h2, items1–4 and the radial bound give an upper multiplier

```
1000^14 * (1e80)^14 * 1e200 * 1e232 <1e1600.
```

The finite mode, block, and derivative-norm sums fit inside the remaining
margin up to1e2000, proving (3). For N, k<=16, multiply each physical
term by q^k. If i radial derivatives act on n then
`q^k D_y^i n=q^(k−i)(q^iD_y^i n)` with q^(k−i)<=1. The weighted
radial conversion and the same budgets give less than1e1600 before
finite sums, proving (2) with the same safe1e2000.

For zero extension, take compact smooth vector-valued v_m converging to
v in H16(R). The explicit formulas map these to smooth compact h_m on A.
The already proved bounded conversion implies h_m→h2 in H14(A).
Zero extension of compact h_m preserves the covariant integral norm;
completion therefore supplies the ambient H14 extension and proves its
zero boundary traces in the Sobolev sense. No trace assertion is needed
for N.

## Scope

The source derivation, low-order coercivity, global regularity and
covariant conversion are separate continuous arguments, with exact finite
arithmetic replay as described in [REPRODUCE.md](REPRODUCE.md). The result
remains conditional on the retained smooth background and exact-family
premises. Minimum energy supplies the low-order comparison; the
high-order bound follows from the quantified normal-operator estimate
and physical derivative conversion proved here.

No evolved second-field bound, common nonlinear proper-time cylinder,
physical energy or force, or uniform nonzero amplitude interval is
established by this initial estimate.
