# Evaluated interior normal recovery for the actual normalized system

14 November 2026 UTC. Conditional mathematical proof. This is a fixed-time
estimate for a smooth frozen solution under a predecessor state-ball
hypothesis. It does not construct that solution, show a Picard iterate
stays in the ball, complete a finite-Sobolev limit, or establish an evolved
parameter or observer theorem. No preceding capsule is changed.

## 1. Equation, norms and final constants

Use `Sigma=[0,6ell] x S3`, `ell=10^-6`, `0<=s<=13ell/8`, the fixed
positive product connection and full measure `dz dmu_kappa`. The field
bundle is

```
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D.
```

All symmetric off-diagonal entries are counted in the Frobenius norm.
The maximal original rotating tensor rank is two. All covariant
derivative slots are retained. Define

```
X_m(w;s)=sum_(p=0)^m ||partial_s^p w(s)||H^(m-p)(Sigma).
```

The smooth predecessor v is assumed to obey

```
sup_s X14(v;s)<=R=10^-12.                                (1)
```

The smooth unknown u need not lie in that ball. It solves the actual
normalized scalar-principal sixteen-component equation

```
P_v u=u_ss-A^i(x,v)D_i u_s-B^ij(x,v)H_ij u=f(v),
f(v)=C(x,v)+E^mu(x,v)D_mu v
                         +G^mu nu(x,v)(D_mu v,D_nu v).   (2)
```

Spatial indices have four values and mu,nu have five. H is the full
fixed-product Hessian, including its direction-index connection. The
normalized coefficients and the lower source are the actual maps, after
the complete first-jet connection substitution. In particular the
horizontal coefficient derivative includes `-K_V[Omega V]`, and
`C(x,0)=0` exactly. The coefficient certificate supplies
`Cstar=10^500000` for every ordered horizontal/vertical derivative
through total order sixteen in its declared full-array/operator norms.

The same-state past reflection is `U(s)=J_E V(-s)` with the retained
tensor time-slot action and opposite elapsed velocity. Its absolute
coefficient and norm bounds are the ones used here. The future field is
not relabelled as the original past field.

Let R_1,...,R_6 be the natural SO(4) lifts on E, and put

```
E_v[w]=(1/2)integral (|w_s|^2+B^ij<D_iw,D_jw>+|w|^2),
Y14(u;s)^2=sum_(p+|I|<=13) E_v[partial_s^p R^I u](s).     (3)
```

Every rotation word is ordered. Functions retain arbitrary angular
frequency. The final evaluated estimate is

```
X14(u;s) <=N Y14(u;s)+L_u X13(u;s)+L_v X13(v;s),
N=10^110,        L_u=L_v=10^500208.                     (4)
```

The same constants apply to the specified identical two-copy Hilbert
union. The estimate uses no smallness of u. Its lower unknown norm is
kept instantaneous; no invented short-time factor is applied to it.

## 2. Actual normalized metric margins and zeroth reciprocal

The retained metric patch gives, in its moving inverse frame,

```
||Gamma^-1-J||op<=1/39,
38/39<=-Gamma^-1,ss<=40/39,
Gamma^-1,spatial>=(38/39)I,  ||Gamma^-1||op<=40/39.
```

The common factor rho(z)^-2 cancels when normalizing at the same point.
The fixed angular derivative frame contributes sech(s). Since
`cosh(s)<=8/7`, the explicit fixed-product consequences are

```
|A|<=1/19<1,
(931/1280)I<=B<=(20/19)I,
19/20<=B^zz<=20/19,       1/B^zz<=20/19<2.                (5)
```

These are also derived in the companion
`NORMALIZED_METRIC_MARGINS.md`. For example the radial lower ratio is
`(38/39)/(40/39)=19/20`, and the fixed angular lower ratio is
`(19/20)*(49/64)=931/1280>1/2`. The mixed inverse time/spatial block
has norm at most 1/39, which gives `|A|<=2/38`.

Thus it is safe to use `I/2<=B<=2I`, and

```
(||w_s||2^2+||Dw||2^2+||w||2^2)/4 <=E_v[w].             (6)
```

The predecessor ball is a legitimate coefficient-domain hypothesis:
the actual cylinder tensor embedding has constant E3=10^7, so
`|v|star<=2E3 R=2*10^-5<10^-3`. It is not a conclusion about the
successor. No high derivative of `1/B^zz` is assumed or used below.

## 3. Rotation/covariant identities and evaluated lower rows

The companion `ROTATION_MIXED_PRODUCTS.md` proves the two exact
identities used for angular recovery. Natural rotations commute with
`D_z,partial_s` and intertwine the entire fixed covariant derivative
and Hessian, including their added slots. At a north pole the three
nonzero rotation vectors are an orthonormal tangent frame, while the
three isotropy rotations act only on tensor slots. Therefore

```
||D_S^k T||2^2<=sum_(|I|=k)||R^I T||2^2.                 (7)
```

There is no curvature-commutator remainder in this identity. The
connection and round curvature have already entered the natural lift's
direction-slot action. Replacing that lift by componentwise angular
directional differentiation would invalidate (7) and the exact principal
commutator below. No two rotation generators, and no two covariant
derivatives, are reordered.

The same finite proof gives the forward bridge
`(sum_|I|=k ||R^I T||2^2)^(1/2)<=40^k||T||Hk`
within the required field ranks and orders. Its common reserve is

```
Drot=40^14<10^23.                                       (8)
```

High-rank coefficient maps before state contraction are controlled by
Cstar, never by this field bridge or by the rank-sixteen embedding.

For clarity, here is the numerical lower-row interface needed in (4).
If `mathcal D=D_z^a partial_s^p R^I` has total degree at most twelve,
then

```
||mathcal D f(v)||2 <=Mprod X13(v),
||[P_v,mathcal D]u||2 <=Mprod X13(u),
Mprod=10^500100.                                        (9)
```

These are per ordered row, before any row-count square sum. Their finite
accounting is recorded fully in `ROTATION_MIXED_PRODUCTS.md` and can be
summarized directly from the accepted coefficient inputs as follows.

For composed rank-one/two principal coefficients K=A,B, covariant
chain differentiation has at most `(n+1)!` leaves at order n. Low
coefficient rows of degree q<=7 use H3 embedding and covariant order
q+3<=10, giving `<10^40 Cstar` in L-infinity after (8). High
coefficient rows through degree twelve have L2 bound
`<10^35 Cstar[1+X13(v)]`. No unknown field is hidden in these constants.

An unisolated principal commutator of degree n<=12 has at most
`2(2^n-1)<10^5` ordered-subset terms. A nonempty coefficient subset
of degree q>=1 leaves unknown total order

```
d=n-q+2<=13,               q+d<=14.                     (10)
```

For q<=7 use the coefficient L-infinity bound and the unknown L2
word/covariant bridge, at cost below `10^63 Cstar X13(u)` per
term. For q>=8, d<=6: the unknown has L-infinity bound
`E3 Drot X9(u)<=10^30 X13(u)`, while the coefficient uses the
displayed high L2 bound and `X13(v)<=R<1`. This costs below
`2*10^65 Cstar X13(u)` per term. Multiplying by the subset count
gives the per-row commutator bound `<10^71 Cstar X13(u)`, which is
smaller than the second line of (9). It is linear in the arbitrary u.

For the source, first use the exact zero-state extraction

```
C(x,v)=integral_0^1 C_V(x,theta v)[v]dtheta.
```

Thus every source leaf contains a predecessor field factor. The actual
three groups C, E Dv and G(Dv,Dv) have, through n<=13, at most
`3*(n+3)!/3!` leaves; at n=13 this is less than 10^14.
Through n<=12 their individual field orders are at most thirteen and
their total field derivative order is at most fourteen. At most one
factor has order at least eight. Put it in L2; every remaining factor
has order at most seven and L-infinity bound
`E3 X10(v)<=10^-5<1`. The covariant source bound is less than
`10^15 Cstar X13(v)`. Applying the output-section bridge (8) gives
`<10^38 Cstar X13(v)` for an ordered normal-recovery row, proving
the first line of (9). Derivatives falling on the two distinct Dv
slots in G are both included. The exact C(x,0)=0 prevents a false
constant forcing remainder.

The largest exponent just spent for normal rows is 500071, below
500100. The source and coefficient estimate is therefore fully
evaluated in (9), with no W14,infinity hypothesis and no ball condition
on u. It uses X10 low norms inside the predecessor's X14 ball.

## 4. Base normal count from time/rotation energies

Define the nested arrays

```
W_r(u)^2=sum_(a+p+|I|<=14, a<=r)
                     ||D_z^a partial_s^p R^I u||2^2,
0<=r<=14.                                               (11)
```

Only W_1 is initially compared with energy. Its terms split into four
disjoint groups:

1. No normal derivative and p+|I|<=13: the field part of (3) gives
   at most `2Y14^2`.
2. One normal derivative, p+|I|<=13: its spatial-gradient energy gives
   at most `4Y14^2` by B>=I/2.
3. No normal derivative, p+|I|=14 and p>=1: use the time derivative
   of the corresponding degree-thirteen word, costing `2Y14^2`.
4. Pure rotation words of length fourteen: write the first rotation
   followed by a word of length thirteen. For any original E section
   w, its rank is at most two and

```
sum_alpha||R_alpha w||2^2
 <=2||D_S w||2^2+48||w||2^2 <=104 E_v[w].                (12)
```

The first inequality follows from the six Killing lifts, each spin
operator having norm at most two, and the Parseval identity for their
directional parts. It is deliberately weaker than the exact isotropy
identity (7). The second uses (6) and the energy's zero-order term.
Summing over the degree-thirteen suffix words costs `104Y14^2`.

No factor 6^13 is added to these aggregate sums: those suffix words
already occur separately in Y14. It follows that

```
W_1^2 <=(2+4+2+104)Y14^2=112Y14^2<256Y14^2,
W_1 <=16Y14.                                            (13)
```

No energy is applied to D_z u in this step. The normal derivative is
only one component of the energy gradient of the original word.

## 5. Differentiate the unisolated equation before division

Fix `2<=j<=14` and a row `j+p+|I|<=14`. Apply

```
mathcal D=D_z^(j-2) partial_s^p R^I
```

to the **full equation (2)**. Its degree is at most twelve.
Expand the product derivatives of the principal coefficients before
performing any division. The unique term containing the desired highest
normal derivative is

```
B^zz D_z^j partial_s^p R^I u.
```

Terms in which at least one derivative hits B^zz belong to the second
line of (9) and contain unknown total derivative order at most thirteen,
as in (10). They are not derivatives of an inverse coefficient.

All zero-hit principal terms have strictly smaller normal count than j:

| Term | Normal count | L2 bound before division |
|---|---:|---:|
| `D_z^(j-2) partial_s^(p+2) R^I u` | j-2 | `W_(j-1)` |
| `A^z D_z^(j-1) partial_s^(p+1) R^I u` | j-1 | `W_(j-1)` |
| `A^A D_A D_z^(j-2) partial_s^(p+1) R^I u` | j-2 | `W_(j-1)` |
| `B^AB H_AB D_z^(j-2) partial_s^p R^I u` | j-2 | `4W_(j-1)` |
| `2B^zA H_zA D_z^(j-2) partial_s^p R^I u` | j-1 | `4W_(j-1)` |

The first three rows use (5), (7), and the corresponding fields in
(11). For the fourth, the angular block has Frobenius norm less than
`2sqrt(3)<4`, and the complete angular Hessian is bounded by the
two-rotation array in (7). For the last, the row-vector norm of B^zA
is less than two; the factor two is required by symmetry of B and
`H_zA=H_Az`. Its one angular derivative is bounded by the one-rotation
array. The product connection is flat in the z/ angular mixed plane,
so that last identity is exact.

In every row, total derivative order is at most fourteen and normal
count at most j-1. In particular the A^z and B^zA contributions
are retained **in the interior**. Their vanishing shell traces are
not used to remove them. The angular Hessian keeps both ordered
direction slots; its curvature terms are not discarded by replacing it
with raw component derivatives.

The zero-hit costs sum to `1+1+1+4+4=11`. Only after this expansion
multiply by the zeroth reciprocal in (5), bounded by two. With (9),

```
||D_z^j partial_s^p R^I u||2
 <=2[11W_(j-1)+Mprod X13(u)+Mprod X13(v)].               (14)
```

This procedure never differentiates `1/B^zz`. Solving first for
u_zz and then differentiating that isolated expression would require a
different inverse-jet ledger, which is not being assumed. Likewise
mathcal D u is not assigned any differentiated boundary condition.
The calculation is an interior equation identity integrated over the
cylinder, not a normal-commuted boundary-energy estimate.

## 6. All thirteen stages, their row multiplicities and finite sums

The full finite array in W14 has exactly

```
sum_(k=0)^14 binom(16-k,2)6^k=135413275530<10^12         (15)
```

rows. The binomial counts the allowed nonnegative time/normal orders
for a rotation word of length k. This counts ordered rotation words,
not commuting monomials. At stage j the new rows are

```
n_j=sum_(k=0)^(14-j)(15-j-k)6^k,
j=2,...,14:
(3134566561,522427758,87071291,14511880,2418645,403106,
 67183,11196,1865,310,51,8,1).                            (16)
```

Thus every square-root row multiplicity is below 10^6. Define the
instantaneous lower quantity

```
D_low=Mprod X13(u)+Mprod X13(v).
```

Adding squares of the new rows to W_(j-1), and using the triangle
inequality for the resulting nonnegative two-vector, (14) gives

```
W_j <=W_(j-1)+sqrt(n_j)[22W_(j-1)+2D_low]
 <=10^8 W_(j-1)+10^7 D_low.                             (17)
```

The finite comparisons are `1+22*10^6<10^8` and `2*10^6<10^7`.
There are exactly thirteen stages j=2,...,14. Consequently

```
W_14 <=10^(8*13) W_1
           +10^7 sum_(r=0)^12 10^(8r) D_low
 <=16*10^104 Y14+2*10^103 D_low
 <=10^106 Y14+10^104 D_low.                              (18)
```

For the geometric sum, factoring its largest term gives
`sum_(r=0)^12 10^(8r)<2*10^96`. All derivative-order and word
multiplicities are accounted for before this rounding. No high
coefficient bound is multiplied thirteen times: coefficient-hit rows
are already in the lower unknown/source norm D_low; the triangular
leading part uses only the evaluated zeroth metric bounds.

## 7. Recover the literal mixed X14 norm

The product connection commutes with radial differentiation and has
parallel radial/angular splitting. Therefore at each spatial order k,

```
||D_sp^k T||2^2
 =sum_(a=0)^k binom(k,a)||D_z^a D_S^(k-a)T||2^2.          (19)
```

This counts every ordered placement of normal/angular derivative
indices without commuting two angular derivatives. Use (7) on the
angular derivative tensors and `binom(k,a)<=2^14`. For each p<=14,

```
||partial_s^p u||H^(14-p) <=2^7 W_14.
```

Summing all fifteen temporal terms gives the convenient safe estimate

```
X14(u)<=15*2^7 W_14=1920W_14<10^4 W_14.                 (20)
```

The sharper aggregate constant from the companion is below 500; (20)
is retained here to display the literal time sum directly. Combining
(18) and (20) proves

```
X14(u)<=10^110 Y14+10^108 Mprod X13(u)
                             +10^108 Mprod X13(v),
```

which is (4), since `108+500100=500208`.

## 8. Copies, compatibility and stopping point

Every L2/energy identity above uses the full positive direct-sum fiber
norm. For the retained identical two-copy union all sums are Hilbert
sums over copies. The pointwise coefficient estimates use the maximum
over copies and mixed products use
`sum_i a_i^2 b_i^2<=(sum_i a_i^2)(sum_i b_i^2)`. The numerical
constants in (4) therefore remain valid for that specified union.

The boundary conditions need not be normally differentiated at any
point. In the separate energy proof, only time/rotation words preserve
the homogeneous four-Dirichlet/twelve-Neumann shell split and the outer
all-Dirichlet row. The scalar principal coefficient and actual traces
cancel their flux. Interior normal recovery neither adds a boundary
row nor supplies a general boundary-forcing theorem.

For a smooth compatible Picard successor the lower identities
`X13(w;s)<=I14+integral_0^s X14(w;r)dr` may subsequently be inserted
into (4). That is the retained finite-jet fundamental-theorem step;
it does not turn an instantaneous source into a small time multiple.
The constants (4) alone do not show the predecessor hypothesis (1),
close the Picard majorant, prove convergence or finite-Sobolev top
continuity, or establish evolved parameter differentiability.

Read-only retained sources are the accepted private next-variable plan
and, at base `d27cfd999aa2e57995c08b1228a73d3da7c973d1`,
`HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`
sections 2 and 4, `HDBLAST_BOUNDARY_FLOW_PATCH_20261106/METRIC_ENERGY_CERTIFICATE.md`,
and the full normalized coefficient certificate in
`HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md`.
The new companion rotation/product and normalized-margin proofs supply
their explicitly cited finite interfaces. Existence and regularity
premises remain conditional analytic inputs, not consequences of the
finite arithmetic ledger.
