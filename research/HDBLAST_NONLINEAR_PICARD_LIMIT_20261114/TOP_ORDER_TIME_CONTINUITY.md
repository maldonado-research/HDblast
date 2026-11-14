# Strong time continuity of the top jets from weak linear energy

14 November 2026 UTC. Conditional strengthening of the
accepted reduced-solution limit. This proves time continuity of each
solution's full X14 jet table by weak linear energy regularity and
interior normal recovery. It does not prove strong-X14 convergence of
the approximating sequence, top-norm continuity of the data-to-solution
map, an X15 tame estimate or evolved parameter differentiability.
All previous accepted files remain unchanged.

## 1. Inputs and the precise stronger conclusion

Use the fixed sixteen-component positive product bundle and the actual
normalized scalar-principal equation

```
P_U U=F(U),
P_U=partial_s^2-A(x,U)^i D_i partial_s-B(x,U)^ij H_ij,
F(U)=C(x,U)+E(x,U)DU+G(x,U)(DU,DU),
U_D=0, D_z U_N=0 at z=0;    U=0 at z=Z,
Z=6ell, ell=10^-6, 0<=s<=13ell/8.
```

The same-state elapsed reflection, complete covariant derivative slots
and full Frobenius norm are as in the accepted baseline
`NONLINEAR_LIMIT_AND_UNIQUENESS.md`, accepted original-proof SHA256
`44c54fe64f9287a350077f6ab43b7121481238f9f4727b6ad7d9e54f2dfa7862`.
It supplies the actual finite-data solution U, its exact initial jets,
the homogeneous mixed rows, and

```
U_j in C H^(13-j) intersect C_weak H^(14-j),    0<=j<=13,
U14 in C H^-1 intersect C_weak L2,
sup_s X14(U;s)<=M<2*10^-15<R/4, R=10^-12.              (1)
```

Every jet is an actual derivative of U. The desired addition is

```
U_j in C([0,T];H^(14-j)),                0<=j<=14.       (2)
```

The previously proved pointwise high bound and the admitted original
parameter interval will be retained unchanged.

Use the accepted same-root `WEAK_LINEAR_ENERGY_REGULARITY.md`, SHA256
`766ae3649fbf84dda9607c4e3dbbe587554e62001b9ba041cdfa10a597157eaa`.
Its closed-manifold wave corollary says: for fixed compatible positive
bundle connection, scalar principal coefficients Lipschitz in space
and time with B uniformly positive, a distributional solution
`w in L-infinity H1`, `w_s in L-infinity L2`, with source in
`L1 L2` and the stated weak initial traces, is in `C H1`, with
velocity in `C L2`. The proof uses the first-order symmetrizer
`diag(1,B,1)` and spatial Friedrichs commutators; no differentiability
of a high energy norm is being assumed. We verify below why it applies
to the mixed-wall problem here.

## 2. Commuted sections and their bulk equation

Use the canonical roster `T=partial_s^p R^I`, with I an ordered word
in the six natural SO(4) lifts and degree `n=p+|I|<=13`, and set
w=T U. Time commutes with these lifts, so this roster covers the
commuted sections without duplicating time placements. Natural lifts
act on every field and derivative slot. They commute with D_z and
intertwine the full product Hessian; no angular generators or
covariant derivatives are reordered. The accepted top bound gives

```
w in L-infinity H1,           w_s in L-infinity L2.     (3)
```

Their weak H1 and L2 continuity, including the original initial traces,
follow from (1) by bounded linear spatial word maps. For example a
word `partial_s^p R^I`, `p+|I|<=13`, maps U_p in
`H^(14-p)` to H1 and U_(p+1) in `H^(13-p)` to L2.
Original rank plus the spatial derivative slots is at most sixteen.

In the interior the distributional commuted equation is

```
P_U w = g_T := T F(U)+[P_U,T]U.                         (4)
```

The full ordered-subset principal commutator is the one in the accepted
product/energy proof; its coefficient-hit subsequences retain their
order and multiplicity. Their external signs are positive. Composed
coefficient arrays, including the horizontal state-connection term,
are used before the natural rotation bridge. The accepted aggregate
interfaces give, almost everywhere in time,

```
(sum_(p+|I|<=13)||g_(partial_s^p R^I)||L2^2)^(1/2)
 <=2 Mprod X14(U),               Mprod=10^500100.       (5)
```

These finite products are defined at (1): each leaf has at most one
top-order L2 field factor; all other actual field factors have the
three-spatial-derivative reserve for L-infinity. Products with a weakly
continuous top L2 factor and strongly continuous low factors are weakly
measurable L2 sections; separability and the displayed bound give strong
measurability. Thus each g_T is in L-infinity L2, in particular L1 L2.
Equation (4) is not obtained by applying a normal derivative to a
boundary energy estimate.

## 3. The actual weak conormal identity at degree thirteen

The missing Neumann trace of a degree-thirteen H1 section cannot be
assigned as an ordinary trace. What is needed instead is its weak
conormal identity. We prove it with the original mixed variational
domain

```
H_D^1={phi in H1: phi_D|_(z=0)=0, phi|_(z=Z)=0}.
```

For compactly time-supported smooth tests in this domain, the identity
for (4) is

```
-integral <w_s,phi_s>
+integral <w_s,D_i(A^i phi)>
+integral <D_j w,D_i(B^ij phi)>
 =integral <g_T,phi>.                                  (6)
```

All coefficients in (6) are evaluated along U. This is the
nondivergence equation written with one spatial integration by parts;
the derivatives of A and B have been retained. It encodes zero weak
Neumann conormal on the N block. It is not a claim that D_z w_N has
an H1 normal trace.

Here is why the commuted boundary identity has no extra forcing.
The actual state traces and scalar principal structure give

```
A^z=B^zA=0 at z=0.                                    (7)
```

At z=Z all fields U are Dirichlet zero. The auxiliary chart fixes
zero at every z, since its background-subtracted generator vanishes
at zero and hence `F_(plus/minus z)(0)=0`. Thus the physical metric
there is exactly the retained diagonal warped background, giving also

```
A^z=B^zA=0 at z=Z.                                    (8)
```

The mixed cross terms need not vanish in the interior. This use of
outer zero state follows from the actual zero-preserving chart, not
from a claimed causal support margin.

Time and natural angular derivatives preserve the boundary and the
fixed D/N split. Their derivatives of the zero traces (7),(8) are
zero distributions on either boundary. These are ordinary Sobolev
traces whenever they are used below: a composed coefficient word
of degree q<=13 has spatial regularity at least H^(14-q), hence at
least H1. Trace commutes with tangential differentiation and with
the distributional time derivatives in the weaker trace spaces.
For time derivatives this follows by applying the bounded H1 trace
operator to the corresponding Banach-valued distribution identity.
The coefficient regularity just used follows by the finite composition
rule through total order fourteen: a leaf has total actual field order
at most fourteen, and after selecting a maximal L2 factor every other
factor has order at most seven. Its three-derivative infinity reserve
is within X14. The scalar/rank-two coefficient output plus fourteen
spatial slots has rank at most sixteen, and the coefficient certificate
has order sixteen. This regularity assertion needs no new numerical
admission constant.

Commute the original weak equation by T, using the skew adjoints of
the natural rotations and the distributional time derivative. The
tests remain in H_D^1. When the coefficient-hit divergence terms
are rewritten as the L2 commutator source in (4), a possible shell
boundary term has the form

```
(T^J B)^(zj) D_j(T^(Jc) U_N),          J nonempty.       (9)
```

The remaining unknown word has degree l<=n-1<=12, so it is in H2.
For j=z, its ordinary Neumann trace is zero: commute its time and
natural angular word with the original normal row, at the available
top trace regularity. For j=A the coefficient trace in (9) is zero
by (7). These products are well-defined trace products, not formal
products of nonexistent traces. Indeed q=|J|>=1 and l=n-q give
coefficient trace regularity H^(13.5-q) and normal/tangential
derivative trace regularity at least H^(12.5-l); their summed
regularity is `26-q-l=26-n>=13`, well above the boundary dimension
three requirement. The relevant low factor can therefore multiply
the other distributional trace. The analogous first-order boundary
term from A has the zero coefficient trace `(T^J A)^z` at both ends;
the lower velocity word has H1 or better regularity. At z=Z every
test component is zero, so no B-conormal condition is needed there.

There are no field-matrix boundary couplings: A and B are scalar on
the sixteen-component fiber, and the fixed connection preserves
E_D and E_N. The original uncommuted boundary term vanishes by its
strong rows. Thus (6) follows. It states exactly the needed weak
conormal condition for w, while assigning no normal row to any
new normally differentiated field.

One may also justify this commutation by the smooth construction
used for U. Each time/rotation-commuted smooth successor satisfies
(6) with its own predecessor coefficients. Those coefficients converge
in C W1-infinity; all commutator products converge weakly in space-time
L2 because their single top factor is weak L2 and all low factors
converge strongly. The strong lower limit identifies the bulk source
as (4). Passing these weak identities yields (6), with no appeal to
a putative classical Neumann trace of w. The lower-word trace
calculation (9) explains why this weak limiting condition is exactly
the original homogeneous mixed condition rather than a new boundary
source.

## 4. A fourfold double with Lipschitz positive coefficients

Let

```
S0=P_N-P_D,                 SZ=-Id_E.
```

These are parallel fiber isometries commuting with the angular
connection and natural rotations. Reflect the field across z=0 by
`w_ext(-z)=S0 w(z)`: Neumann fields are even and Dirichlet fields
odd. Reflect across z=Z by `w_ext(2Z-z)=SZ w(z)`, making all fields
odd there. Repeating the two reflections gives a field on

```
S1_(circumference 4Z) x S3.
```

More explicitly, translation by 2Z acts by `SZ S0`, whose square is
the identity. Thus the second translation makes the field 4Z-periodic;
an unjustified 2Z-periodicity is not assumed. The bundle is the same
product bundle on this closed spatial double. Its connection remains
fixed and smooth. The H1 extension is legitimate: every odd component
has the required Dirichlet trace zero, and even components match their
values automatically. The velocity is extended in L2 by the same
fiber maps, and the source by `g_ext(reflected z)=S_boundary g(z)`.

Across either end extend the principal coefficients by

```
A^z odd,       A^A even,
B^zA odd,      B^zz and B^AB even.                      (10)
```

An even reflection of a Lipschitz function is Lipschitz. An odd
reflection is Lipschitz precisely with the continuity condition that
its boundary trace is zero; (7),(8) supply it here. All extended
coefficients are consequently W1-infinity in space and time, with
bounded one-sided derivatives and no delta in their first weak
derivatives. Higher normal smoothness across the seams is not needed.
For Rz=diag(-1,1,1,1) on spatial derivative directions, the reflection
of the principal coefficients is `A_ext=Rz A`,
`B_ext=Rz B Rz`. Hence their scalar-principal positivity and the
normalized energy margins are preserved on every reflected piece.
The coefficient extensions are defined by these base-direction
parities; w and g are extended by the fiber parities. This constructs
a linear commuted equation on the double. It requires no identification
of the extended coefficients with the nonlinear coefficient map
evaluated at a reflected state. The coefficient fields are already
2Z-periodic, while the field's fiber parities require the 4Z circle.

Identity (6) proves that the extended field solves

```
P_ext w_ext=g_ext                                     (11)
```

as a distribution on this closed double. To verify it for an arbitrary
smooth doubled test, add its reflected pullback to the original piece
with S0 or SZ as appropriate. At z=0 the resulting D test is odd and
has zero trace while the N test is unrestricted, exactly the test
domain of (6). At z=Z the resulting whole test is odd and zero.
The coefficients' parities (10), product connection, and parallel
fiber maps make all reflected volume terms match. This calculation
accounts for possible seam distributions: zero weak conormal in (6)
removes the even-Neumann seam term, and the odd-Dirichlet test
cancellation removes the Dirichlet seam term. No unavailable pointwise
normal derivative of w is invoked.

There are four identical-norm pieces. Thus each L2 or H1 squared norm
on the double is four times its original value, and the L2 source
norm is twice the original. In particular (3),(5) supply exactly the
finite energy-class/source assumptions of the closed-manifold lemma.
This temporary doubling has no effect on the physical one-copy or
identical-two-copy convention and changes no amplitude constant.

## 5. Strong continuity of every tangential energy section

Apply the accepted weak-linear lemma to (11). The doubled initial
data have weak H1/L2 traces given by the exact original nonlinear
jet table: (1) and bounded word maps already supply those traces,
and reflection is a bounded H1/L2 map. The lemma gives

```
w_ext in C H1,              partial_s w_ext in C L2.
```

Restriction to the original interval therefore gives, for each of
the finitely many ordered words of degree at most thirteen,

```
T U in C H1,                 partial_s(T U) in C L2.    (12)
```

The number of words is large but finite; no convergence of an infinite
angular harmonic sum is being assumed in this last step. Each section
contains all angular modes. Equivalently, its finite aggregate energy
is continuous. This is a strong continuity statement for one solution;
it has not shown that the energies of an approximating solution
sequence converge uniformly.

For normal recovery the useful base array is

```
W_r(s)^2=sum_(a+p+|I|<=14,a<=r)
                 ||D_z^a partial_s^p R^I U(s)||L2^2.
```

Interpret W_r here also as the corresponding finite L2-valued array.
Equation (12) makes W1 strongly continuous: terms of total degree
at most thirteen without normals are field parts; one normal uses
the H1 gradient; degree-fourteen terms with a time derivative use
the velocity in (12). Pure rotation words of length fourteen are
one natural rotation of a degree-thirteen section, a fixed bounded
map H1 to L2. These are exactly the four groups in the accepted
base estimate `W1<=16 Y14`; no energy is applied to D_z U.

## 6. Continuous interior normal recovery in all thirteen stages

For `j=2,...,14`, apply
`D_z^(j-2) partial_s^p R^I`, `j+p+|I|<=14`, to the unisolated
actual equation. The accepted interior proof differentiates before
dividing and obtains

```
Bzz D_z^j partial_s^p R^I U
 = zero-hit principal rows with normal count at most j-1
   + coefficient-hit commutator rows + source row.      (13)
```

The source and coefficient-hit words have degree at most twelve.
In their finite product expansion every individual actual field
factor has intrinsic order at most thirteen and total field order
is at most fourteen. Since a coefficient is hit at least once in
the commutator, the unknown factor there also has order at most
thirteen. The C X13 continuity from the baseline therefore makes
these remainders C L2: put one maximal factor in L2 and every other
factor, of order at most seven, in L-infinity with three spatial
derivatives in reserve. This conclusion does not require continuity
of an unproved top factor. Coefficient compositions are handled before
natural rotations, as in the accepted full-slot product proof.

Inductively assume the finite array W_(j-1) is strongly continuous.
Every zero-hit principal row in (13) is then C L2. These include
the interior A^z mixed time-normal term and both symmetric B^zA
terms, with their required factor two. Their vanishing boundary
traces do not remove these bulk terms. Angular derivatives of a
lower-normal row are continuous because the fixed inverse natural
rotation estimate bounds the derivative of the difference at two
times by the corresponding finite rotation-word difference array.
The full angular Hessian and its direction slots are retained.

The zeroth coefficient `1/Bzz` is bounded by `20/19` and continuous
in L-infinity since U is C H13 and the background is smooth.
Multiplying (13) by this continuous multiplier proves continuity of
the new L2 normal row. No normal derivative of that inverse and no
normal-commuted boundary condition is used. The same finite argument
at each of the thirteen stages gives

```
W14 in C([0,T]; finite direct sum of L2).               (14)
```

The accepted exact spatial-slot decomposition and inverse angular
bridge apply to the difference of the jet tables at any two times,
with the same harmless factor `1920<10^4`. Thus (14) implies (2).
In particular U14 is now strongly continuous in L2. At zero every
top jet has the same J_j trace already identified weakly in the
baseline; injectivity identifies the strong trace with it. The
pointwise estimate `X14(U;s)<=C0 d exp(2Hs)<=M` remains exactly the
old bound, because no new norm amplification enters (13) or this
continuity conclusion.

## 7. Scope of the strengthening

The constructed reduced solution has a time-continuous full X14
jet table on the unchanged original interval. This property also
holds for every solution of the baseline high class to which these
same weak-row and source hypotheses apply. The proof uses no X15
estimate and no order-fifteen seed.

It does not prove that the Picard or supported-data approximating
sequence converges strongly in C X14, nor top-norm continuity or
Lipschitz dependence on data. Those are distinct assertions; the
lower difference estimate and a smoothing comparison or substitute
still have to be supplied for them. No evolved C2 parameter,
subsidiary/Israel propagation, physical support/observer, or empirical
claim is added. All boundary and reflection arguments concern the
actual scalar-principal reduced system and its preserved orthogonal
fiber split.

Read-only retained source facts used specifically at the outer seam:
`research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md`,
Sections4--6, has `B_generator(0)=0` and `F_(plus/minus z)(0)=0`;
`METRIC_ENERGY_CERTIFICATE.md`, Sections1,4, identifies the exact
background `Gamma=J` and its diagonal warped principal tensor.
The accepted rotation/product, frozen-energy and interior-recovery
proofs are the source interfaces already authenticated by the
baseline capsule. The new weak-linear lemma and baseline hashes
above pin the only added analytic dependencies of this strengthening.
