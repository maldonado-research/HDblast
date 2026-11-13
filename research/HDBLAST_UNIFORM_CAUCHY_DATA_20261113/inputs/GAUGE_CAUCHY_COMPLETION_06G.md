# Quantitative initial relative-harmonic Cauchy completion

6 November 2026. Conditional mathematical checkpoint. This concerns initial
amplitude derivatives at zero for the selected weighted canonical curve.
It supplies no amplitude interval, evolved estimate or proper-clock claim.
The algebraic completion is retained work; the contribution here is its
application to the selected general second metric and a quantified norm bridge.

## Conclusion and sources

The generic gauge companion allows freely chosen initial lapse and shift.
The retained exact-family binding selects them completely:

```
alpha_e=rho(1+e f), beta_e=0, f=uY, Y=n1^2-n2^2.
```

Consequently `alpha2=beta2=0`. The actual physical second tensors together
with the specified first packet determine all fifteen metric and one scalar
wave fields and first time velocities. The single inequality `I2<10^2132`
does not specify those tensors; it gives a bound after applying this map.
There is no remaining second-order initial gauge obstruction for this
explicit binding. A generic physical-data statement without this choice
would have five freely specifiable initial second gauge fields.

The precise source pins, including SHA256 and Git blob hashes where retained,
are in `SOURCE_PINS.json`. Relevant source passages are:

* `REGISTERED_EXACT_FAMILY_BINDING.md` 64–66 (derivative convention), 74–113
  (first physical packet and canonical normalization), 167–180 (affine scalar
  sectors), 214–267 (selected lapse/shift, reference connection and full first
  gauge jet), and 317–337 (physical second conversion).
* `INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` 20–43 (normals and its alternate
  canonical normalization), 62–129 (freely selected gauge and exact local
  completion), and 138–152 (subsidiary normal derivative uses reduced equations).
* 06e `PHYSICAL_CAUCHY_SECOND_JET.md` 12–41 (fixed geometric norm and `I2`),
  45–76 (actual canonical conversion), 95–145 (known-product/background budgets).
* 06e `SOURCE_JETS.md` 145–209 and its pinned rational engine (global finite
  Taylor majorants, with sixteen stored orders).
* 06f `COORDINATE_DICTIONARY.md` 12–45 and 49–71 (coordinate/root distinction
  and retained varying bulk reference connection).

The 16 October `EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`, SHA256
`bc3ae85c3166575ed15b017468a3e9668615e15dfaca34bb2a13afe27ea28f92`,
already performs a specialized second completion for a different polynomial
physical second jet. Its formulas are consistent with the general map below;
its specialized second metric must not be substituted for this weighted curve.

## Conventions and exact initial identities

All subscripts 1,2 mean amplitude derivatives at zero, with the second term
`e^2/2`. Write `A=rho(y)`, `a=rho_y/rho`, and use proper background coordinate
`y`, with inward conformal `dy=-rho dz`. On the initial section `t=0`,

```
gamma0=dy^2+rho^2 kappa_S3, K0=Pi0=0,
h=gamma1=f(-4dy^2+2gamma_angular),
L=K1=S(-2dy^2+gamma_angular), S=ut Y/rho,
psi=phi1=vY, p=Pi1=vt Y/rho,
k=gamma2, M=K2, phi2=0, Pi2=-2fp.
```

The temporal ADM convention is `K=Lie_u gamma/2`,
`u=alpha^-1(partial_t-beta)`, and `Pi=u(phi)`. The actual half-action
canonical densities are

```
pi_act=sqrt(det gamma)(K^ij-gamma^ij tr K)/2,
varpi_act=sqrt(det gamma) Pi.
```

The generic companion uses twice these densities; its displayed inverse
canonical formula cannot be applied to `pi_act` literally. Its completion
identities in the physical variables `K,Pi` apply without a rescaling.
The 06e physical `K2` includes its factor two and quadratic conversion, and
`Pi2=-2fp` is retained. `alpha` is the temporal lapse, distinct from the
radial component root `N_y=sqrt(g_yy)` and from the boundary map `F_z`.

The gauge is exactly

```
C^A=g^BC(Gamma[g]^A_BC-GammaBar^A_BC), GammaBar=Gamma[g0],
R^A=g^BC GammaBar^A_BC.
```

There is no additional prescribed gauge-source vector. In particular the
metric dependence of `R` must be differentiated. At `t=0`, `beta_e=0`
for every amplitude. The reference `GammaBar^0_ij=GammaBar^0_00=0` for
spatial `i,j`; `GammaBar^0_0y=a` need not vanish, but `g_e^0y=0`.
Thus `R_e^0=0` identically on the initial section for every amplitude.

The retained harmonic completion reads

```
gamma_t=2alpha K, phi_t=alpha Pi,
alpha_t=alpha^2 tr_gamma K,
beta_t^i=alpha^2(Gamma_gamma^i-R^i)-alpha gamma^ij partial_j alpha.
```

No second time accelerations are independent Cauchy data. The physical
constraints and reduced equations imply the normal subsidiary derivative.

## All second wave fields and velocities

Use `D=nabla_gamma0`, `tr=tr_gamma0`, and set

```
e(k)_i=D^j k_ji-(1/2)D_i tr k.
```

In the fixed `(t,y,Omega)` identification the complete second metric values
are

```
(g_tt)2=-2rho^2 f^2, (g_ti)2=0, (g_ij)2=k_ij, phi2=0.       (1)
```

The spatial and scalar velocities are

```
(g_ij,t)2=2rho M_ij+4rho f L_ij,
(phi_t)2=rho Pi2+2rho f Pi1=0.                              (2)
```

The lapse and temporal metric velocities are

```
(alpha_t)2=rho^2(tr M-24fS),
(g_tt,t)2=-2rho^3 tr M+44rho^3 fS.                          (3)
```

Indeed `tr L=S`, `h:L=14fS`, and
`(tr_gamma K)2=tr M-28fS`. The scalar cancellation in (2) does not assert
that the physical normal `Pi2` vanishes.

The four mixed velocities are the covector `J` given by

```
(g_ty,t)2=rho^2[e(k)_y+a k_yy-26f f_y-2a f^2],
(g_ta,t)2=rho^2[e(k)_a+a k_ay+22f partial_a f].              (4)
```

Angular indices here denote coordinate covectors, so `partial_a f` is the
angular differential of `f`; this formula is tensorial over the initial
spatial manifold. The full first mixed velocity is also retained:
`(g_ty,t)1=-6rho^2(f_y+4af)` and the first angular mixed velocity is zero.

For a direct derivation retaining the reference source, put

```
C_i(gamma)=gamma^jk(D_j gamma_ki-(1/2)D_i gamma_jk).
```

Because `GammaBar^i_tt=rho D^i rho`, lowering the shift completion gives
the exact identity at initial `beta=0`

```
g_ti,t=alpha^2 C_i(gamma)+gamma_ij rho D^j rho-alpha partial_i alpha. (5)
```

Its second derivative is

```
rho^2{e(k)_i-2h^jk(D_j h_ki-(1/2)D_i h_jk)+4f e(h)_i}
 +rho k_i^j rho_j-2rho f partial_i(rho f).                  (6)
```

The reference contribution is the second term in (6). Omitting it loses
`rho^2 a k_iy` in (4) and even spoils the exact background cancellation
in (5). In the warped orthonormal splitting,

```
e(h)_y=-5f_y-18af, e(h)_a=partial_a f,
h^jk(D_j h_ky-(1/2)D_y h_jk)=2f f_y-36af^2,
h^jk(D_j h_ka-(1/2)D_a h_jk)=-10f partial_a f.
```

These produce (4) with signs fixed by the positive `K` convention.

## An explicit unstraightened wave Cauchy norm bridge

Extend the fixed spatial connection to `R dt + T*Sigma` by `D dt=0` and
give this bundle the positive fiber metric `dt^2+gamma0`. Use the actual
geometric measure `dvol_gamma0` and all spatial covariant derivatives in
the Sobolev norm, exactly as for 06e. Define the explicit unstraightened norm

```
W2=||g2||_H14(dt^2+gamma0)+||phi2||_H14
  +||g_t2||_H13(dt^2+gamma0)+||phi_t2||_H13.                 (7)
```

This is a specified positive analytic tensor norm with fifteen metric
components and one scalar, not an unspecified coordinate norm equivalence.
The claim is

```
W2 <10^3641 <10^4000 per copy.                              (8)
```

Here is a complete oversized budget for (8). The pinned 06e source engine
recomputes global rational proper-`y` derivative envelopes through fifteen:
`u<10^93`, `S/Y<10^91`. Through fourteen, its envelopes for
`rho,rho^2,rho^3,rho^2a,a` are respectively
`10^3,10^5,10^7,10^5,1`. The following radial profiles, including the
derivative needed in (4), all have derivative envelopes through fourteen
strictly below `10^100`:

```
rho^2 u^2, rho u^2, rho^3 u(ut/rho), rho u(ut/rho),
rho^2 u u_y, rho^2 a u^2.                                  (9)
```

Their largest decade ceilings are `86,84,86,82,93,85` respectively. These
are global finite Taylor majorants over the compact packet, not samples.
Stored order sixteen supplies the required order-fifteen margin for `u_y`.
In particular, the angular term `rho^2 u^2 Y dY` in (4) uses the radial
profile `rho u^2` with the reference sphere covector `Y Dbar Y`, because
`dY=rho^-1 Dbar Y` in the radially parallel orthonormal splitting. This
accounts explicitly for the angular physical frame scale.

To convert these radial bounds to full spatial derivative bounds, retain
the 06e background bounds `|partial_y^j a|,|partial_y^j rho^-1|<10^80`
through eighteen. In the warped orthonormal splitting

```
D_0=partial_y, D_a=rho^-1 Dbar_a+a K_a.
```

For rank `r+j<=16`, the full directional slot-mixing map has norm at most
`3(r+j)`. Expand each derivative into primitive radial coefficient jets,
reference angular derivatives, and parallel slot-insertion maps. After
`j` stages a term has at most `j` background factors. A radial derivative
branches into at most `j+1` differentiated factors; a reference angular
derivative adds one branch, and slot mixing adds norm at most `3(r+j)`.
For `r<=2`, `j<=13` their total is at most 60, hence the uniform recurrence
`E_0=1, E_(j+1)<=1000 E_j` is valid. The larger 1000 also covers fixed
directional decompositions. No derivative slot is omitted.

The base angular sections in (1)–(4) are made from `Y^2`, `Y Dbar Y`, the
radial covector and angular identity. The latter two are parallel under the
reference product connection. Their physical `D` derivatives include the `a K_a`
slot-mixing terms above. From
`Dbar^2 n_i=-n_i kappa` and `|Dbar^j n_i|<=2^j`, the finite product rule
bounds all these reference angular derivatives needed through order
fifteen by `10^200`, including fixed displayed coefficients. Therefore
for a base profile with radial derivative bound `B`, a valid stage-14
upper bound is

```
1000^14 (10^80)^14 B 10^200.                               (10)
```

For each fixed multiplier `rho`, `rho^2`, `rho^3`, or the covector
`rho^2 a dy`, use `B=10^80` in (10); all their `C14` derivative bounds
are below `10^1442`, hence below `Bfix=10^1500`. The combined
`rho^2 a dy` is treated as one multiplier. Its cost is not charged twice.

For known quadratic terms take `B=10^100` in (10). Their zero extensions
are flat. Their support lies in the annulus of width `2ell` in `z`, whose
full geometric volume is at most `40ell*131^4<10^6` using `pi^2<10`.
Thus its square root is below `10^3`; the derivative-sum square root is
below four. Allowing a factor 1000 for all finitely many displayed terms,
tensor coefficients and symmetric mixed slots gives the explicit budget

```
Bknown=4*10^6*1000^14*10^(14*80+100+200) <10^1600.           (11)
```

For any fixed tensor multiplier `b` and unknown tensor `T`, the covariant
Leibniz rule gives

```
||b T||_Hm <=2^m sqrt(m+1) max_{j<=m}||D^j b||infty ||T||_Hm.
```

For `m<=14` the first factor is below `10^6`. All contractions use the
fixed positive metric. Specifically,
`||e(k)||_H13<=3||k||_H14`, `||tr M||_H13<=2||M||_H13`.
The mixed symmetric tensor costs `sqrt(2)<2`. Equations (1)–(4) now yield

```
W2 <=(1+8*10^6*Bfix) I2+Bknown
   <10^1508*10^2132+10^1600 <10^3641.                      (12)
```

This uses the existing physical second data and zero-extension result;
it does not redo the weighted inverse or identify a historical weight.
For two identical reflected copies, the sum-of-integrals norm gains
`sqrt(2)` and the rounded ceiling `10^4000` still holds.

## Limits and generic gauge freedom

Without the explicit binding, replacing initial lapse by
`alpha_e+e^2 A2/2` and shift by `beta_e+e^2 B2/2`, for arbitrary smooth
annulus-supported `A2,B2`, preserves physical canonical data and the full
first tangent. It changes `(g_tt)2` by `-2rho A2` and `(g_ti)2` by
`gamma0_ij B2^j`. The harmonic equations fix their time velocities after
this choice; they do not select their initial values. Thus physical `I2`
alone could never bound those gauge additions uniformly. The selected
affine lapse and zero shift exclude this freedom.

Equation (8) concerns the unstraightened metric/scalar tensor norm (7).
It is not a bound for the nonlinear radial-root or `F_-z` straightened
`H14` data; those require their own differentiated chart budget. It gives
no bound for second time accelerations, higher corner jets, evolved
solutions, a nonzero admitted amplitude interval, observer clocks,
proper time, embeddings or actual-area observables.
