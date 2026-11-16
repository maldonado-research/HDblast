# A quantitative all-mode compatible seed through order fourteen

14 November 2026 UTC. Conditional mathematical proof. The construction is a linear
extension of a finite initial jet table. It does not solve the nonlinear
wave equation and does not use a qualitative precontact Cauchy solution.
Every constant below is in the same fixed product-bundle norms as the
retained temporal hierarchy.

## 1. Input, target and result

Let `M=[0,Z] x S3`, `Z=6ell`, `ell=10^-6`, with product metric
`dz^2+kappa`, full sphere measure and the fixed positive bundle

```
E=R_v + Sym^2(R dt+T*S3)_H + R_w + (R dt+T*S3)_D.
```

The symmetric block has its full Frobenius norm. All covariant derivative
slots are included. Denote this fixed product connection by nabla.
The twelve Neumann-block fields are v,H,w and the four Dirichlet-block
fields are D. The actual straightened wave boundary conditions are
`nabla_z V_N=0`, `V_D=0` at z=0 and `V=0` at z=Z.

Suppose the given initial jets satisfy

```
J_j in H^(14-j)(E),             0<=j<=14,
I14=sum_(j=0)^14 ||J_j||_H^(14-j) <infinity,              (1)
```

are invariant under the retained reflection group Gamma, and vanish off
the common closed radial annulus
`A0=[ell/2,5ell/2] x S3`. In particular each jet has zero open endpoint
collars. The compact-support closure at its stated finite regularity is
part of this input; no new ordinary endpoint trace is imposed.

There is an explicit linear extension S(J), defined below, such that

```
partial_s^j S(J)(0)=J_j strongly in H^(14-j), j<=14,
S(J) in C^k(R;H^(14-k)),                            k<=14,
sup_s sum_(k=0)^14 ||partial_s^k S(J)(s)||_H^(14-k)
                                      <=10^255 I14.      (2)
```

It is reflection invariant and zero on the fixed open spatial collars
`z<=ell/4` and `z>=11ell/4` for every s. It therefore satisfies the
required mixed boundary conditions, and every available boundary-corner
jet is zero. For s not zero it is spatially and temporally smooth.
The same operator norm works on the specified identical two-copy Hilbert
union if each norm in (1) already uses that union; no additional copy
factor is needed for this linear extension.

Write `C_lift=10^255` for the jet-table coefficient in (2). The
coefficient multiplying the original Cauchy-pair size d(a), denoted
C_seed in the preceding plan, also includes the separate nonlinear
jet-map constant; that joined coefficient is recorded in section 9.

The claimed extension reproduces exactly J0 through J14. It does not
claim to reproduce every higher interior nonlinear temporal jet of smooth
input data.

## 2. The auxiliary positive operator and its exact spectrum

Fix a global unit left-quaternion frame `E_1,E_2,E_3` on S3, normalized by

```
[E_A,E_B]=2 epsilon_ABC E_C,
nabla_(E_A) E_B=epsilon_ABC E_C,
Delta_comp=sum_A E_A^2.
```

Components of E use an orthonormal fiber basis in this frame. In
particular an off-diagonal symmetric basis vector is
`(e_i tensor e_j+e_j tensor e_i)/sqrt(2)`, so independent matrix-entry
coordinates are not silently substituted for Frobenius norm.
Delta_comp acts on these scalar component functions, not by the geometric
tensor rough Laplacian.

Define the positive self-adjoint auxiliary operator

```
A=1-partial_z^2-Delta_comp                               (3)
```

componentwise, with ND boundary conditions for each Neumann field and
DD boundary conditions for each Dirichlet field. Its separated
orthonormal eigenbasis is the usual scalar spherical harmonic basis
combined with

```
ND: cos((n+1/2)pi z/Z), n>=0,
DD: sin(n pi z/Z),       n>=1,
```

with the usual L2 normalizations. If `xi_n` is the displayed radial
frequency and `lambda_l=l(l+2)`, its eigenvalue is

```
Lambda=1+xi_n^2+lambda_l >=1.                             (4)
```

This spectral definition also specifies its self-adjoint realization.
Put `H_A^m=D(A^(m/2))`, with norm `||A^(m/2)f||_2`.
The complete scalar harmonic tower is used, with no finite angular cutoff.

A is an auxiliary extension operator. It is not the nonlinear wave
principal operator, and it need not commute with the physical reflection
action on the chosen component frame.

## 3. Explicit all-slot graph/covariant norm comparison

For the scalar component array f let W_j(f) be the square-summed L2 norm
of all ordered words of length j in `partial_z,E_1,E_2,E_3`; include
W_0=||f||_2, and put `W_<=m^2=sum_(j<=m)W_j^2`.

Each E_A is divergence-free and Killing, and commutes with the scalar
Laplacian. The sum of squared angular words of length b on a spherical
eigenfunction is exactly `lambda_l^b` times its squared norm. This
follows inductively from `-sum E_A^2=-Delta_comp`; the E_A need not
commute with one another. Radial differentiation changes sine to cosine
or vice versa. The differentiated radial families remain orthogonal
with the same frequencies in both DD and ND cases. Radial and angular
derivatives commute. Thus for a finite spectral sum,

```
W_<=m(f)^2=sum_(n,l,components) sum_(j=0)^m
               (xi_n^2+lambda_l)^j |fhat_(n,l)|^2.
```

The binomial coefficients of `(1+t)^m` lie between one and 2^m, so

```
W_<=m(f) <=||f||H_A^m <=2^(m/2) W_<=m(f).                (5)
```

No collar-width factor is lost: the large normal frequency xi_n remains
in both norms.

Here is the full tensor/frame comparison. Let C_j be the raw components
of `nabla^j T`, for an original tensor of angular rank at most r<=2.
The original dt and z slots are parallel. In the left-quaternion frame,

```
(C_(j+1))_0=partial_z C_j,
(C_(j+1))_A=E_A C_j+Omega_A^(r+j) C_j.                   (6)
```

The connection matrices have constant raw frame entries. One covariant
or contravariant slot carries an elementary angular rotation of operator
norm at most one. On every original and derivative slot together their
collected three-direction norm is at most
`sqrt(3)(r+j)<=2(r+j)`. This includes the connections on all earlier
derivative indices. The same statement holds on the direct-sum bundle.

Expanding C_j in raw word derivatives gives constant slot maps. A further
raw derivative hits only the word, since those matrices have constant
entries, while (6) adds a connection term. The total coefficient-map norm
therefore grows by at most `1+2(r+j)<=5(j+1)` at stage j. Conversely,
solve (6) for the raw derivative of C_j and repeat the identical argument.
In a reverse expansion only C_a with a<=j occurs, whose connection rank
is at most r+j. Consequently both expansion costs through order m are
at most `5^m m!`. Summing the m+1 derivative norms yields

```
||T||Hcov^m <=sqrt(m+1)5^m m! W_<=m(f),
W_<=m(f) <=sqrt(m+1)5^m m! ||T||Hcov^m.                  (7)
```

These expansions retain the ordered noncommuting frame words; they never
replace them by commuting partial derivatives. Combining (5)–(7), for
m<=14 both directions are bounded by the single evaluated integer

```
C_B=512*5^14*14!
   =272432160000000000000000 <10^24,
C_B^-1 ||T||Hcov^m <=||f||H_A^m <=C_B ||T||Hcov^m.        (8)
```

Indeed `sqrt(15)<4` and `2^7*4=512`. The final tensor rank is at
most sixteen; every derivative slot was counted in (6).

For compact smooth tests vanishing near both endpoints, all integer
powers of the differential operator (3) preserve that compact support,
and every ND/DD boundary condition needed for its powers holds. Hence
these tests belong to every H_A^m. Their finite spectral projections
converge in H_A^m, and (7) shows convergence in Hcov^m. Applying (8)
to compact-smooth approximations proves that a jet with the stated
zero-collar/closure property belongs to H_A^m with the same bound.
Thus (8) applies to every J_j in (1), including the merely L2 jet J14.
No graph-domain statement is asserted for arbitrary Hcov^m data lacking
the required boundary compatibility.

## 4. Corrected heat coefficients and uniform heat derivative bounds

For j=0,...,14 define

```
K_j=sum_(r=0)^floor(j/2) j!/[r!(j-2r)!] A^r J_(j-2r),
U(s)=e^(-s^2 A) sum_(j=0)^14 s^j K_j/j!.                (9)
```

If h=j-2r then `A^r J_h in H_A^(14-j)`, with norm exactly
`||J_h||H_A^(14-h)`. Therefore (9) is well-defined from the finite
regularity data and

```
sum_(j=0)^14 ||K_j||H_A^(14-j)/j!
 <=sum_(h=0)^14 ||J_h||H_A^(14-h)/h! sum_(r>=0)1/r!
 <=3 sum_(h=0)^14 ||J_h||H_A^(14-h).                    (10)
```

The factor three follows directly from the exponential series e<3.

Let `P_(k,j)(x)e^(-x^2)=d_x^k[x^j e^(-x^2)]`.
Starting with `P_(0,j)=x^j`, the exact recurrence is

```
P_(k+1,j)=P_(k,j)' -2x P_(k,j).
```

Its degree is at most j+k<=28. Its coefficient l1 norm is multiplied
by at most j+k+2 at the next stage. Uniformly for j,k<=14 it is at most
`product_(i=0)^13(16+i)=29!/15!`. For every n<=28 and x>=0,

```
|x^n e^(-x^2)| <=1+n! <=2*28!.
```

For x<=1 use one; for x>=1 the nth term of the exponential series gives
`e^(x^2)>=x^(2n)/n!`, hence the bound n!. The same absolute bound
holds for negative x. Thus a common exact scalar multiplier bound is

```
C_H=2*28!*29!/15!
   =4122968597824599187678433233286039863296000000000
   <10^49.                                               (11)
```

For a spectral value Lambda and x=s sqrt(Lambda),

```
partial_s^k[s^j e^(-s^2 Lambda)]
          =Lambda^((k-j)/2) P_(k,j)(x)e^(-x^2).
```

The input and target weights cancel exactly:
`Lambda^((14-k)/2)*Lambda^((k-j)/2)=Lambda^((14-j)/2)`.
Hence the operator norm from H_A^(14-j) to H_A^(14-k) is at most
C_H, uniformly for every real s, every angular mode and all radial
frequencies. Combining with (10),

```
sup_s ||partial_s^k U(s)||H_A^(14-k)
 <=3 C_H sum_j ||J_j||H_A^(14-j),             k<=14.      (12)
```

## 5. Strong initial traces, not a formal unbounded Taylor expansion

For fixed Lambda, the derivative multiplier at s=0 vanishes unless
`k=j+2r`; in that case its value is
`k!(-1)^r Lambda^r/r!`. The uniform bound (11), with the exact
weight cancellation above, dominates the spectral L2 sum for every
term K_j. Dominated convergence therefore proves the strong limits in
H_A^(14-k). They are

```
partial_s^k U(0)
 =sum_(r=0)^floor(k/2) k!(-1)^r/[r!(k-2r)!] A^r K_(k-2r)
 =J_k.                                                   (13)
```

For the second equality substitute (9). The coefficient of
`A^h J_(k-2h)` is
`k!/(k-2h)! * sum_(r=0)^h (-1)^r/[r!(h-r)!]`.
The sum is one for h=0 and zero for h>0, by the binomial theorem.
All involved powers have the regularity required for their common
target H_A^(14-k). This explicitly includes the lower-j terms j<k.

Finite spectral projections of the given jet table give classical smooth
heat expressions. Inequality (12), applied to their differences, gives
uniform convergence of every derivative through k in the weaker fixed
space H_A^(14-k); lower derivatives take values in stronger spaces and
embed with norm at most one because Lambda>=1. Completeness of
`C^k([-T,T];H_A^(14-k))` for any finite T proves that these are the
actual derivatives and are continuous at zero. This is not a Taylor
expansion of an unbounded operator in its operator norm.

For s not zero, exponential spectral decay puts U(s) in every graph
domain, with all time derivatives. The same word/connection comparison
at any fixed finite order supplies ordinary spatial smoothness.

## 6. A fixed flat cutoff with an evaluated multiplier norm

Use the dimensionless flat bump

```
Q(t)=exp(4-1/[t(1-t)]) for 0<t<1, and Q=0 otherwise,
I_Q=integral_0^1 Q(t)dt,
S(x)=I_Q^-1 integral_(-infinity)^x Q(t)dt.
```

On [1/3,2/3], Q>=exp(-1/2)>1/2, so `I_Q>=1/6`.
Consequently `0<=S<=1`, S=0 for x<=0, S=1 for x>=1, and
`||S^(k)||infinity<=6||Q^(k-1)||infinity` for k>=1.

Here is a finite global derivative bound, without sampling the bump.
Put `B=1/[t(1-t)]>=4`. Its normalized derivatives obey
`|B^(k)|/k!<=B^(k+1)`, by writing B=1/t+1/(1-t).
The Bell-series majorant gives, for n>=1,

```
|Q^(n)|/n! <=exp(4-B) [u^n] exp(sum_(k>=1)B^(k+1)u^k)
 <=exp(4-B) 2^n B^(2n).
```

Indeed with r factors, the coefficient count is
`binom(n-1,r-1)/r!`, whose sum is at most 2^(n-1), and the
power is n+r<=2n. Using exp(4)<100 and the maximum of
`B^(2n)exp(-B)`, then discarding its decreasing exponential factor,
gives for n<=13

```
||Q^(n)||infinity <=B_Q:=100*13!*2^13*26^26
 =31403446004315944454719861034756724848550331023360000
 <10^53.                                                 (14)
```

The n=0 value is at most one and is covered as well. These estimates
also prove flat extension at the endpoints, since each fixed derivative
is bounded by a polynomial in B times exp(-B).

Choose the requested cutoff

```
chi(z)=S(8z/ell-2) S(22-8z/ell).                         (15)
```

It equals one on `[3ell/8,21ell/8]`, and vanishes outside
`(ell/4,11ell/4)`. Its two transition intervals are disjoint. On the
left transition the second factor is identically one; on the right
transition the first is identically one. Thus at positive derivative
orders only one differentiated factor is present, and for j<=14,

```
||chi^(j)||infinity <=M_chi:=6 B_Q (8/ell)^14 <10^151.    (16)
```

The product spatial connection has parallel dz, so these are the full
covariant cutoff derivative norms, with no hidden angular connection.
The complete tensor Leibniz estimate gives, for m<=14,

```
||chi T||Hcov^m <=2^m sqrt(m+1) M_chi ||T||Hcov^m
 <=C_chi ||T||Hcov^m,
C_chi=2^14*4*M_chi <10^156.                             (17)
```

Thus the cutoff's complete H14 multiplier cost has been charged; its
flat value alone was not treated as the multiplier constant.

## 7. Projection, final norm and exact jet preservation

Let Pi_Gamma be averaging under the natural physical pullback action of
the sixteen coordinate reflections on the bundle and all its derivative
slots. Each is an isometry of the fixed product metric, preserves the
D/N split, and commutes with covariant differentiation in its natural
tensor action. Therefore `||Pi_Gamma||Hcov^m->Hcov^m<=1`.
It commutes with time differentiation and with the radial cutoff.

Define the final extension by

```
S(J)(s)=Pi_Gamma[chi(z) U(s)].                            (18)
```

There is no assertion that A commutes with Gamma. The projection is
performed **after** the spectral lift. Since each original J_j is
invariant and supported where chi=1, (13) gives
`partial_s^j S(J)(0)=Pi_Gamma(chi J_j)=J_j` exactly.

Apply (8) once to the input table and once to each output derivative,
then (12) and (17). Summing the fifteen time-derivative norms gives

```
sup_s X14(S(J);s)
 <=45 C_B^2 C_H C_chi I14
 <=45*10^(48+49+156) I14 <=10^255 I14.                   (19)
```

The homogeneous final estimate is read with <=, including the zero
table. This proves (2). All operators act diagonally on the retained
two-copy Hilbert union. Its linear norm is unchanged when the input I14
already sums the global copy norms; the input must not silently be
switched from per-copy to union convention in this join.

## 8. Corners, smooth approximants and scope of compatibility

Equation (18) is zero on fixed open endpoint collars for every s,
before and after Gamma projection. Hence both the Neumann and Dirichlet
conditions, and every boundary-corner identity allowed by the finite
regularity, hold there with zero perturbation. The construction requires
no extension of the physical background through z=0.

Finite spectral truncations are used in section 5 only to prove graph
membership and strong trace convergence. Independently truncated J0
through J14 need not be the nonlinear temporal hierarchy of the
truncated first Cauchy pair. They are not asserted to be nonlinear
compatible initial data.

For a smooth compatible-data approximation in the parent Picard argument,
first approximate the original Cauchy pair by smooth Gamma-invariant
data with support in a fixed slightly larger subannulus strictly inside
the plateau of chi. Recompute its nonlinear J_j by the same local
finite recursion. These smooth jets remain supported in that subannulus.
Apply (9) and (18) to the recomputed table. They then match precisely
that nonlinear table through fourteen. Smooth compact jets belong to
every graph domain, so this seed is smooth also at s=0 and has all-order
zero boundary-corner compatibility. Continuity of the parent's finite
jet map gives convergence of its table in (1), and bounded linearity
of (18) gives convergence of the seeds in the claimed mixed norm.

The cutoff can accommodate this fixed support enlargement because the
original annulus lies strictly inside its plateau. For example
`[7ell/16,41ell/16]` contains A0 with room and is itself strictly inside
`[3ell/8,21ell/8]`. Compact approximation and finite averaging preserve
that fixed radial support. No quantitative smoothing constant is used
in (19), only density and the parent's continuous finite-jet recursion.

If the original J table obeys the differentiated nonlinear equation
through temporal order twelve, the seed shares those identities at
s=0, since they involve at most fourteen initial temporal jets. This
is an initial-jet statement, not an equation solved for s not zero.
For smooth input, the truncated lift is still not asserted to match
nonlinear interior jets of every order beyond fourteen.

Finally, the operator in (18) is bounded and linear from the jet space
in (1) to the mixed seed space. A C2 parameter curve of such jet tables
therefore gives a C2 parameter curve of seeds. No numerical finite-jet
parameter-Hessian constant, evolved solution bound, continuation time,
observer tube, or physical conclusion is supplied by this extension
lemma alone.

## 9. The selected initial curve: a seed on an evaluated subinterval

The accepted companion `NONLINEAR_INITIAL_JETS.md` (main jet
original-proof SHA256 `6e49654baf3bba785907b808b0fb95216cbc9f1a8f7babeda6cca2d5b4d07574`)
supplies, in these literal norms and with the same original parameter,

```
I14<=10^6503000 d(a),
d(a)<=min(10^-2200,10^4200 a),
d(a)<=10^-6503020 as its finite-jet smallness guard,
a_jets=10^-6507220.
```

Here d(a) is the supremum of the actual completed initial Cauchy-pair
norm on |e|<=a, and d(a)<=10^4200 a uses the accepted whole-interval
first derivative and the zero background data. It does not substitute
the older, looser first-jet ceiling 10^5200.

In the notation of the preceding seed plan, the complete joined constant
is therefore

```
C_seed=10^(6503000+255)=10^6503255,
sup_(|e|<=a,s) X14(S(J(e));s)<=C_seed d(a),               (20)
```

whenever the parent's jet-map hypotheses and smallness guard hold.
C_lift=10^255 and C_seed=10^6503255 have different input spaces and
must not be interchanged in a numerical interface.

Choose the smaller, explicitly evaluated seed subinterval

```
a_seed=10^-6507485,
6507485=6503000+255+4200+30.                              (21)
```

It lies strictly inside both the jet interval and the original admitted
initial interval `|e|<=10^-6402`. The exact exponent joins are

```
d(a_seed)<=10^(4200-6507485)=10^-6503285
                                      <10^-6503020,
sup_(|e|<=a_seed) I14(e)<=10^(6503000-6503285)=10^-285,
sup_(|e|<=a_seed,s) X14(S(J(e));s)<=10^(255-285)=10^-30. (22)
```

For each s the proven cylinder embedding gives the stronger initial
iterate-state guard

```
||S(J(e))(s)||_* <=|S(J(e))(s)|_E,infinity
                  <=10^7 X14(S(J(e));s)<=10^-23<10^-3.   (23)
```

The star block norm is bounded by the positive E Hilbert norm in the
same fixed bundle, so no extra moving-frame factor is needed here.
This is an explicit small state patch for the constructed seed, not for
an evolved solution or an entire Picard sequence.

The companion finite-jet map is C2 on the relative open space of
collar-supported Gamma-invariant Cauchy pairs with
`d<2*10^-6503020`. The closed actual curve in (21) lies strictly inside
that neighborhood. Bounded linearity of the present lift therefore gives
a C2 seed curve on a relative open parameter neighborhood of this closed
subinterval. No numerical seed parameter-Hessian norm is assigned.

All these joins use the parent's norm convention for the actual one-copy
or retained identical two-copy family. When the parent uses union norms,
I14 in (21) is already the union input and section 7 adds no copy factor.
The original parameter, nonzero first velocity, tensor reflection signs,
and physical zero collars remain those of the given jet table. The
construction neither replaces that table nor changes the selected curve.
In particular the seed is an explicit compatible extension of its initial
data; it is not the corresponding nonlinear spacetime solution.
