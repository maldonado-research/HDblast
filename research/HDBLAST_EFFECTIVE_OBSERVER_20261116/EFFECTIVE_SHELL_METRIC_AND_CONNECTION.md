# Effective actual shell metric, scalar and connection bounds

16 November 2026 (America/Los_Angeles; 17 November UTC). New proof at immutable public base
`b4f2f59ab3918264c9832d87be1e39ed4dbfc178`. The source ledger
`SOURCE_PINS.json` pins the actual
base blobs and their direct originating-manifest member joins. All
retained background, selected canonical-family, point-chart and smooth
frozen-linear construction premises remain in force. The prior
prospective admission map is motivation only; no prospective assertion
is used as an established input here.

## 1. Statement and exact conventions

Retain the same original canonical interval and slab,

```
I=[-a_iter,a_iter], a_iter=10^(-10^500209),
ell=10^-6, T=13ell/8, Sigma=[0,6ell] x S3.
```

Let U_e(s) be the actual elapsed normalized solution. Its original-time
representative is `V_e(t)=J_E U_e(-t)`, `-T<=t<=0`, where J_E is the
fixed parallel time-slot isometry. The accepted norm is the literal sum

```
X14(U;s)=sum_(j=0)^14 ||partial_s^j U(s)||H^(14-j)(Sigma;E),
E=R_v+Sym^2(R dt+T*S3)_H+R_w+(R dt+T*S3)_D.
```

The symmetric block uses full Frobenius norm, the spatial connection
is the fixed product connection, and the measure is the full
`dz dmu_kappa`. Use the chosen one-copy or specified identical-two-copy
Hilbert convention already built into this norm; do not switch from
its global norm to a per-copy bound without stating the restriction.

On the actual physical shell `iota(t,n)=(t,0,n)`, put

```
h_e=iota^*g_e, barh_e=h_e/rho_b^2,
barh_0=-dt^2+cosh(t)^2 kappa,
p=dt^2+kappa, D=D_p,
eta_e=barh_e-barh_0,
delta=sup_(e in I) sum_(m=0)^4 sup_([-T,0]xS3)|D^m eta_e|_p.
```

All derivative tensors in this definition are full ordered covariant
arrays, including every time/angular derivative and both symmetric
output slots. The positive p norm is not the Lorentzian metric barh.
The size delta here is not the fixed shell-model tension parameter.

With the already evaluated constants

```
C0=10^7003210, H=10^500212,
Dbar=10^4200 a_iter,
M=2 C0 Dbar exp(2HT),
```

the following bounds hold for the actual entire original I:

```
M<10^(-10^500208),
delta<=1936*10^7 M<10^11 M<10^-20.                  (1)
```

Let `A_e=nabla^(barh_e)-D` and `A0=nabla^(barh_0)-D` be the
intrinsic shell connection differences. Their full p-array bounds are

```
|A_e-A0|<=6delta, |D(A_e-A0)|<=48delta,
|A0|<1, |D A0|<14,
|A_e|<2, |D A_e|<15.                                (2)
```

Let `phi_b=phi0(0)` and let V_w denote the scalar component of the
normalized V, not a parameter variation. The exact scalar identity and
value bound are

```
phi_e(t,0,n)-phi_b=V_e,w(t,0,n),
sup_(e,t,n)|phi_e(t,0,n)-phi_b|<=10^7 M.              (3)
```

There is no unknown reconstruction multiplier in (1)--(3). These are
field-size bounds. They do not bound the numerical size of the first
or second parameter derivative, prove a nonzero response, or assign
SI units to the registered normalization.

## 2. Actual norm source and stronger finite exponent comparison

The nonlinear limit source, equation (9), defines the common admitted
majorant `M=2 C0 Dbar exp(2HT)` for all supported smooth approximants
and iterates, and its finite-data argument retains the actual-datum
bound by weak lower semicontinuity. The current unreduced theorem,
Section 1, identifies the same actual solution and states

```
partial_s^j U in C H^(14-j), j=0,...,14,
X14(U;s)<=C0 d exp(2Hs), d<=Dbar.
```

Thus the actual solution is bounded by the same M, with ample room;
the factor two is inherited from the existing global approximation
allowance. No extra copy factor or new amplitude is inserted here.
The physical-reconstruction source, Section 2, proves the same-member
tensor reflection and its parallel time-slot action. Therefore

```
X14(V;t)=X14(U;-t)<=M.                               (4)
```

For the exponent comparison set `P=10^500207`. The admitted constants
give exactly

```
2HT=(13/4)*10^500206=(13/40)P<P,
10^500209=100P,
M=2*10^(7007410-100P)*exp(2HT).
```

Using e<10 and 2<10, rather than computing a huge integer,

```
M<2*10^(7007410-99P)
 <10^(7007411-99P)
 <10^(-10P)=10^(-10^500208).                         (5)
```

The last inequality is equivalent to `7007411<89P`, witnessed by
`P>=10^6` and `7007411<89*10^6`. All comparisons are strict. The
constants and interval are precisely those of the retained admission;
(5) only strengthens its previously convenient rounding.

## 3. Exact F0 identity and global coframe conversion

The current physical reconstruction gives, in the tangential de Sitter
coframe `(dt,cosh(t)theta)` on a round orthonormal angular frame,

```
u=F_z(V)=(v,Htilde,w,Dcross),
G_ab=rho(z)^2 J4_ab+rho_b^2 Htilde_ab,
phi=phi0(z)+w, J4=diag(-1,1,1,1).
```

Its Section 6 proves the actual physical shell is z=0 and `F_0=Id`.
At that point `rho(0)=rho_b`; hence the tangential block and scalar
are exactly `Htilde=V_H`, `w=V_w`. In the positive product bundle,
define the global endomorphism

```
C(t)=Pi_time+cosh(t)Pi_angular,
L_t H=C(t)^T H C(t).
```

The two projections are p-orthogonal and D-parallel; in particular
their angular covariant derivatives vanish. The literal metric formula
therefore gives, at every actual nonzero amplitude,

```
eta_e(t,n)=L_t V_e,H(t,0,n).                          (6)
```

The global tensor statement (6) is independent of the chosen local
round orthonormal frame. It is not an assumption that local angular
frame components have zero connection: all angular connections are
already in D applied to the tensor V_H. The endomorphism scales the
entire angular subbundle and commutes with that connection.

On the Frobenius symmetric block, L_t has orthogonal blocks with
multipliers 1, cosh(t), cosh(t)^2. Off-diagonal symmetric components
are counted in the full Frobenius metric on both sides; no independent
entry norm or uncharged factor two is substituted. Because `2T<1`
and `cosh(2T)<2`,

```
||partial_t^a L_t||op<=16, 0<=a<=4.                  (7)
```

For a=0 the norm is at most cosh(T)^2<4. For a>=1, the angular
two-slot factor has derivative bounded by
`2^(a-1)cosh(2T)<2^a<=16`; the mixed time/angular and time/time
blocks are smaller. These are ordinary derivatives of a globally
defined bundle map, not coordinate-frame derivative estimates.

No bulk Psi derivative envelope is required for (6). The exact F0
identity removes the nonlinear flow at the shell, while rho_b^2
cancels in barh. The scalar identity in (3) follows simultaneously:
its coframe carries no tensor slots or scale factor, and phi0 is
independent of t.

## 4. Independent audit of the boundary-inclusive embedding constant

The primary `QUADRATIC_FORCING.md`, Section 4, proves for
`Z=6ell=6*10^-6`, using the full product spatial norm,

```
||T||infinity <=10^7 ||T||H3([0,Z]xS3).              (8)
```

It applies to tensor/direct-sum sections with at most sixteen rotating/
spatial slots and requires no endpoint boundary condition. The relevant
charges in its proof are as follows.

The complete round spherical spectrum has eigenvalues l(l+2),
multiplicity `(l+1)^2`, and full volume `2pi^2`. The Hilbert-valued
component evaluation constant is exactly `1/sqrt(12)<1`; it has no
component-count factor. In a global left-invariant round frame the
induced connection on r<=16 slots has elementary norms at most r.
The source's full tensor-to-component conversion bounds

```
||(1-Delta_comp)T||L2
 <=(1+3r^2)||T||L2+2sqrt(3)r||D_S T||L2
                             +sqrt(3)||D_S^2 T||L2
 <1000 ||T||H2(S3),                                 (9)
```

since the three coefficient ceilings are 769,64,2 and their sum is
below 1000. The Frobenius representation on symmetric tensors is
orthogonal, so these are the literal slotwise operator estimates.

The radial Hilbert-valued estimate is

```
sup_z ||T(z)||H2(S3)
 <=Z^-1/2 ||T||L2_z H2(S3)+sqrt(Z)||D_z T||L2_z H2(S3).
```

Here `(Z^-1/2+sqrt(Z))^2=1/Z+2+Z<10^6`, so this costs less than
1000. The needed radial/angular combinations are a subset of the
full H3 derivative square sum. The extra radial slot is parallel;
the argument does not reapply a rank-sixteen angular estimate to a
seventeenth rotating slot. Combining (9) with this radial bound costs
less than 10^6, so the retained rounded 10^7 in (8) is valid. This
audit uses the exact fixed measure and all spherical modes.

Apply (8) at each fixed original t to
`T=D_S^k partial_t^j V_H`, for j+k<=4. Projection onto H and then
onto angular derivative slots has norm at most one, and is parallel.
Its H3 norm is bounded by the same actual time jet in H^(k+3),
which is contained in H^(14-j) because k+3<=14-j. The tensor rank,
including all derivative and embedding reserve slots, is at most
`2+j+k+3<=9`, within sixteen. Consequently

```
sup_(t,n)|D_S^k partial_t^j V_H(t,0,n)|<=10^7 M,
j+k<=4.                                             (10)
```

This is boundary evaluation already covered by (8), not an additional
ordinary trace inequality. No normal derivative is assigned a new
homogeneous row. For an identical two-copy Hilbert union, restriction
to either copy has norm at most one. Its pointwise bound therefore
uses the same global M. The common physical shell metric is identical
on those two copies; delta is its norm, not a sum of two pointwise
norms. This closes the copy convention without a hidden sqrt(2).

## 5. Full C4 array count and the finite 10^-20 threshold

For fixed j time and k angular derivative slots with j+k<=4,
(6)--(7), Leibniz and (10) give

```
sup |D_S^k partial_t^j eta|
 <=16*sum_(a=0)^j binom(j,a)*10^7 M
 =16*2^j*10^7 M.                                    (11)
```

Time commutes with the angular product connection. In a full ordered
spacetime derivative tensor of length m, placing j time directions
in the word gives binom(m,j) isometric copies of that ordered angular
array; no two angular derivatives are reordered. Bounding the full
direct-sum norm by the sum of these block norms yields

```
sup |D^m eta|<=16*sum_(j=0)^m binom(m,j)2^j*10^7 M
              =16*3^m*10^7 M.
```

Sum m=0,...,4, then take the uniform parameter supremum:

```
delta<=16*(1+3+9+27+81)*10^7 M
      =1936*10^7 M<10^11 M.                          (12)
```

Let `N=10^500208` in this paragraph only. Equations (5), (12)
give `delta<10^(11-N)<10^-20`, since N>31. This needs only a
finite exponent comparison, for example `500208>=2` implies
`N>=100>31`; no huge integer is expanded. The original amplitude
interval has not changed.

The scalar component uses (8) with projection norm one. Thus its value
bound is exactly (3). One may also assemble its full product C4 array:
the row bound is 10^7 M without a coframe factor, so

```
sup_(e in I) sum_(m=0)^4 sup|D^m(phi_e|shell-phi_b)|
 <=(1+2+4+8+16)*10^7 M=31*10^7 M.                   (13)
```

In particular (3), (13) are derived from the actual scalar state
component. They do not rely on a generic unknown composition norm.

## 6. Inverse metric and connection perturbation ledger

The inverse and all norms in this section are measured using p. Since
`barh_0=-dt^2+cosh(t)^2 kappa`,

```
|barh_0^-1|op<=1,
|D barh_0|<8T<1/2,
|D^2 barh_0|<8.                                     (14)
```

Indeed `|kappa|p=sqrt(3)<2`, `D p=0`,
`sinh(2T)<4T` and `cosh(2T)<2`. All base/output slots in (14)
are retained as full arrays. The only nonzero background derivatives
are the time derivatives of cosh(t)^2 kappa.
The corresponding zeroth-order metric bounds are explicit as well:
`cosh(T)^2=(1+cosh(2T))/2<3/2` gives
`|barh_0|op<3/2`. Therefore `|barh|op<3/2+delta<2`, and in the
four-dimensional positive p norm `|barh|F<=2|barh|op<4`.

By (12), delta<10^-20<1/4. Neumann inversion of
`I+barh_0^-1 eta` therefore gives, uniformly in e,

```
|barh^-1|op<=1/(1-delta)<=4/3<2,
|barh^-1-barh_0^-1|op<=delta/(1-delta)<=2delta.        (15)
```

This inversion is valid for the Lorentzian metric: it uses the
positive auxiliary operator norm, not positive definiteness of barh.
The path `barh_0+theta eta` stays nondegenerate, so its inertia stays
Lorentzian. In addition `barh_tt<=-1+delta<0`, the angular block is
positive, and the inverse tt component remains negative by (15).
These are consequences of the same perturbation bound.

For brevity put `G=barh^-1`, `G0=barh_0^-1`, `K=G-G0`.
Differentiating inversion with the full connection gives

```
DG=-G(D barh)G, DG0=-G0(D barh_0)G0,
|DG|<4*(1+delta)<5.
```

Use the exact three-term decomposition

```
DK=-K(D barh_0)G-G(D eta)G-G0(D barh_0)K.             (16)
```

The three respective bounds are at most `2delta,4delta,2delta`:
the first uses `|D barh_0|<1/2`, while the last may use the looser
`<1`. Hence `|DK|<=8delta`. Matrix operator norms control their
actions on full Frobenius arrays, and contraction of the remaining
arrays is bounded by the product of their full norms.

The difference of Levi-Civita connections relative to D is

```
A^a_bc=(1/2)G^ad(D_b barh_dc+D_c barh_db-D_d barh_bc).
```

The three signed permutations of a derivative tensor have total norm
at most three times its full Frobenius norm. Their derivative-index
permutations are parallel and isometric. Thus

```
|A-A0|<=(3/2)(|K|op |D barh_0|+|G|op |D eta|)
        <=(3/2)(2delta+2delta)=6delta.                (17)
```

Differentiate this exact tensor formula before estimating. The four
terms in `D(A-A0)` are respectively `DK*D barh_0`,
`DG*D eta`, `K*D^2 barh_0`, and `G*D^2 eta`, with the same
three-permutation factor. Equations (14)--(16) give

```
|D(A-A0)|<=(3/2)(8+5+16+2)delta
          =(93/2)delta<=48delta.                    (18)
```

All estimates homogeneous in delta use non-strict inequalities, so
the zero-perturbation case delta=0 is included.

For the background itself,

```
|A0|<(3/2)*8T=12T<1,
|D A0|<(3/2)(1+8)=27/2<14.
```

Together with delta<10^-20 these yield `|A|<2` and `|D A|<15`.
Every bound is for the intrinsic connection of barh (and hence h,
since rho_b is constant). A general reduced-PDE coefficient envelope
is not substituted for a shell Christoffel norm.

## 7. Completed interface and limits

The delivered quantitative interface is the actual uniform metric
size delta<10^-20, its literal inverse/connection bounds (2), and
the actual scalar identity/value estimate (3), with the stronger
functional bounds in M retained. This uses original-time fields of
the same canonical member, the fixed physical shell and global copy
norms. No parameter is rescaled and no time slab is changed.

These bounds are inputs for the separate stopped intrinsic-geodesic
continuation and Jacobi/cut-metric proof. They do not by themselves
construct a longer observer domain, prove causal illumination or a
nonzero response, or calibrate a physical source/detector norm.
