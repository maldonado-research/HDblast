# Actual first-scalar square source and its weighted Green integral

2026-11-06. Conditional expansion-point estimate for the registered actual first physical scalar.

For the actual first physical observer scalar `psi1(q,n)=S(q)Y(n)`, with
`Y=n1^2-n2^2`, `q=-tau/rho_b`, `ell=10^-6`, and `T=25ell/16`, this
conditional all-point calculation establishes

```
sup_(0<=q<=T) |S(q)|   < 51500,
sup_(0<=q<=T) |S_q(q)| < 3107000000000,

-766143/10^9 < P''(phi_b) < -766128/10^9,

-16000000000 < K_R(q) < 6520000000000000000000000,

K_R=(27/40)S_q^2-(1/5)sech^2(q)S^2+(9rho_b^2/20)P''_b S^2.
```

The bounds are in the fixed model units. The evaluator retains the sharper rational endpoints; their decimal displays are
`|S|<=51475.212388740125`, `|S_q|<=3106350839450.9443`,
`K_R in [-1.5950641521874395e10,6.513355487986372e24]`.
Only the exact rational endpoints establish the inequalities.

This is the known first-scalar part of the weighted second Ricci source.
The full identity remains

```
<R2>_w=K_R+rho_b^2 P'_b <psi2>_w-<E2>_w.
```

The second scalar and electric-Weyl terms are full pulled-back quantities,
including the physical event, observer, embedding, normal and tangent.
Neither is evaluated here, and this calculation does not give `Cr`, full
`R2`, a second solution, a uniform amplitude interval or a nonlinear remainder.

## Retained premises and source joins

Read `inputs/PAST_FUTURE_PROOF.md`, sections 1–4, and
`inputs/HIGHER_JETS.md`. We use their same smooth actual first solution,
branch existence and analytic background enclosures, all required power
domains of the positive dynamic-endpoint operator, local finite propagation,
physical reconstruction and linear uniqueness. These analytic/PDE premises
are conditional inputs; byte hashes and arithmetic do not prove them.
The higher-jet report's actual endpoint enclosures are accepted without
rerunning its producer or selecting or differentiating independent H cells.

`SOURCE_PINS.json` records eight previously public files copied byte for byte from immutable Git base `5ed540606f7ad5096a8baa9749cb13d753a551ef`, with source URLs, Git blob IDs, lengths and SHA256 values. The evaluator authenticates the retained bytes and these metadata. The copies contain no new private inputs. It verifies three
existing higher-report `input_pins`, the actual shell theorem status, model
parameters `delta=1/1000`, `c=96714321425692/161839258930731`, the same
waveform nonzero velocity `exp(-z)Q(z/ell)`, and the common `ell=10^-6`.
The current higher report's `mu` and `gamma` intervals equal the waveform
intervals; the analytic shell's actual `mu` interval is contained in that
reported interval. The report explicitly satisfies
`ell_over_mu=ell/mu_lower`, distinguishing the unscaled boundary `mu` from
the endpoint mass `Lambda=mu/ell`.

For `rho_b` and `p_b` we use the exact shell-box endpoints, contained in
`129<=rho_b<=131` and `-7/50000<=p_b<=-3/25000`; divisions by `p_b` are
therefore justified. The actual scalar interval is `phi_b=1+eta_b`, also
checked to equal `sigma_second_b/4` endpoint by endpoint. Approximately it
is `[0.9999127672122616,0.9999191592122617]`. No replacement by `phi_b=1`
is made. Only shell values of this older shell enclosure are used, so its
shorter radial collar is not used to claim a new evolution domain.

The model and first physical scalar identification match
`REGISTERED_EXACT_FAMILY_BINDING.md` and `PHYSICAL_SHELL_RESPONSE.md`.
The Gauss identity and angular coefficients are derived in the companion curvature theorem. These source/scale joins do not constitute an independent proof of the background.

The exact-family binding is used here for the fixed model, actual branch and initial tangent. Its stated nonlinear slab ends at `3ell/2`, which is shorter than the target `25ell/16`. It does not establish nonlinear continuation to the target cut. The retained smooth actual **linear** development, reconstruction and all-power operator premises through `3ell` independently cover this first-scalar calculation. Any conclusion involving `psi2`, `E2` or `D2` below additionally assumes that the full family and the original proper Gaussian observer cylinder exist through `T=25ell/16`. This target-cylinder existence is not proved here.

## The actual nonzero-velocity past development

Let `s=q/ell`, `r=z/ell`, and `U_plus(s,r)=w(ell*s,ell*r)`. On the common
auxiliary half-line the retained operator and its Hilbert norm are

```
H=L2(0,infinity) plus R,
||(u,y)||_H^2=int u^2 dr+Lambda y^2,
A0(u,y)=(-u_rr+V0*u,[-u_r(0)+b0*y]/Lambda), y=u(0),
Lambda=mu_b/ell>1000,  A0>=kappa I, kappa=5500ell^2,
omega(s)=-(35/4)ell^2 tanh^2(ell*s),
Z_ss+(A0+omega I)Z=0.
```

Both mass and spring are in the operator/form, with linked H1 trace and H2
operator domains and no extra static Robin boundary condition. The actual
past `U_minus(s,r)=w(-ell*s,ell*r)` has initial velocity `-g`, whereas the
future has `+g`, where `g=ell exp(-ell*r)Q(r)`. With
`Delta Z=Z_minus-Z_plus`,

```
Delta Z(0)=0,  Delta Z_s(0)=(-2g,0).
```

The factor two is retained. Flat support of `Q` on `[1,2]` puts these data
in every required power domain and makes all their initial endpoint rows
zero. It does not make the actual velocity zero. Both developments have
support in `[0,2+s]`, so for `s<=25/16` they stay inside `r<=57/16<6`,
where the auxiliary operator is the actual operator. No exterior extension
is used as an assertion about the global physical background. Endpoint
differences vanish for `s<=1` by finite propagation and smoothness.

On this shortened interval put

```
epsilon=(7/4400)(25/16)^2 ell^2<1.
```

Then `A0+omega I >= (1-epsilon)A0>0`. For `j=0,1/2,1`, define

```
Ej=1/2||A0^j Delta Z_s||_H^2
   +1/2||A0^(j+1/2)Delta Z||_H^2
   +omega/2||A0^j Delta Z||_H^2.
```

The retained regularity and scalar commutation give
`Ej_s=omega_s||A0^j Delta Z||_H^2/2<=0`. The complete endpoint terms are
included in this abstract spectral form identity. Half-power vectors need
not be compactly supported; no such support claim is made.

The accepted bump bounds are `||Q||_2<=1`, `||Q_r||_2<=16`,
`||Q_rr||_2<=256`, and `0<V0<103000ell^2` on the length-one initial
support. They give

```
||Delta Z_s||_H <= G0=2ell,
||A0 Delta Z||_H <= Ghalf
  =[4ell^2((16+ell)^2+103000ell^2)/(1-epsilon)]^(1/2)<33ell,
||A0 Delta Z_s||_H <= G1
  =2ell(256+32ell+103001ell^2)<513ell.
```

For the last inequality use `E1(0)=2||A0g||_H^2` and positivity to retain
its kinetic term. There is no coercivity loss in this kinetic estimate.
Since `g` is flat at the endpoint, `(A0g)_b=0` initially, and

```
g_rr=ell exp(-ell*r)[Q_rr-2ell Q_r+ell^2 Q],
||A0g||_H<=ell(256+32ell+103001ell^2).
```

The half-graph estimate retains `1-epsilon` in its potential norm bound.
The checker encloses its square root and `Lambda_min^(-1/2)` from above by
96-bit dyadic rationals, verifying both by exact squaring.

Let `e_b>=sqrt(ell/mu_lower)`, `d_j` bound `|partial_s^j delta y|`, and
`h=25/16`. The endpoint mass and the quiet interval give

```
d1=e_b G0,  d0=(h-1)d1,
Omega=(35/4)ell^4 h^2,   Omega1=(35/2)ell^4 h,
d2=e_b Ghalf+Omega d0,
d3=e_b G1+Omega d1+Omega1 d0.
```

Here `delta y_ss=-(A0 Delta Z)_b-omega delta y` and
`delta y_sss=-(A0 Delta Z_s)_b-omega delta y_s-omega_s delta y`.
Thus `d3` is an actual third-derivative transfer, approximately
`1.128133570640353e-5`. The signs of the physical past odd derivatives are
already represented by `U_minus`; no extra time-reflection sign is inserted.
Adding `d0,d1,d2,d3` to the retained actual future bounds supplies actual
past C3 bounds. No derivative is taken of an interval-cell selection.

## Scalar reconstruction and its derivative

Put `c0(q)=cosh(q)^(-3/2)` and
`a(q)=35sech(q)^2/4-25/4`. The retained physical reconstruction on either
development is

```
S=3mu_b c0/(rho_b^3 p_b)[y_ss/ell^2+a y],
S_q=3mu_b c0/(rho_b^3 p_b)
 [y_sss/ell^3+a y_s/ell+a_q y
  -(3/2)tanh(q)(y_ss/ell^2+a y)].
```

The scalar has the original physical clock, observer and shell embedding.
The constant shell background scalar removes first clock/label amplitude
terms, as in the retained reconstruction; this does not remove second
pullback terms. The first coefficient is distinct from the metric `u`.

Since `0<=q<=T`, we have `c0<=1`, `|tanh(q)|<=T`, `|a|<=5/2`, and
`|a_q|<=35T/2`. The last two follow directly from the hyperbolic formulas
on this short interval. With `Bj` the exact scaled future bounds read from
`HIGHER_JETS_NORMAL.json`, the checker sets

```
C=3mu_upper/(rho_lower^3 (-p_upper)),
Af=B2/ell^2+(5/2)B0,
Sf=C Af,
Sqf=C[B3/ell^3+(5/2)B1/ell+(35/2)T B0+(3/2)T Af].
```

All these are actual future bounds, not comparator-only values. In
particular `B3` is `exact.interior.ythird`; `B0,B1,B2` are the boundary
fields. The two occurrences of `ysecond` are checked equal.

For the difference there is an exact useful cancellation, valid because
the complete endpoint equation is retained:

```
delta S=3mu_b c0/(rho_b^3 p_b)
 [-(A0 Delta Z)_b/ell^2+(5/2)delta y],
delta S_q=3mu_b c0/(rho_b^3 p_b)
 [-(A0 Delta Z_s)_b/ell^3+(5/(2ell))delta y_s
  -(3/2)tanh(q)(-(A0 Delta Z)_b/ell^2+(5/2)delta y)].
```

Indeed `a=5/2+omega/ell^2`; the derivative cancels the `omega_s` term too.
There is no added spatial-reflection sign. Let
`D=e_b Ghalf/ell^2+(5/2)d0`. Then

```
|delta S|<=C D,
|delta S_q|<=C[e_b G1/ell^3+(5/(2ell))d1+(3/2)T D].
```

Approximately these bounds are `17.620859563959865` and
`281933770.70095325`. Adding them to `Sf,Sqf` gives the stated actual
past bounds. A second, uncanceled reconstruction using transferred past
`Bj+dj` independently encloses both computed values. The first scalar and
its derivative vanish through `q=ell`. This first-field quiet interval
does not extend the full `R2` or unknown-projection quiet interval past
the retained `ell/4`.

## Actual polynomial coefficient and weighted source

For the exact model

```
W=1-phi+phi^3/3,  V=(W')^2/2-(2/3)W^2,
sigma=2W+delta(1+c phi),  P=sigma'^2/16-V/2-sigma^2/12,
```

the checker expands the full potential and tension in exact polynomial
arithmetic and differentiates twice. It separately verifies

```
P=delta[c W'/4-(1+c phi)W/3]
   +delta^2[c^2/16-(1+c phi)^2/12],
P''=delta[7c/6-(4c/3)phi^2-(2/3)phi]-delta^2 c^2/6.
```

Since `c,delta,phi_b>0`, the derivative of this last quadratic is strictly
negative on the actual interval. Its endpoint evaluations give the exact
`[Plo,Phi]` in the report, approximately
`[-0.0007661426270170118,-0.0007661281803348285]`.

The coefficient of the unresolved second scalar is also nonzero. The
actual background laws `sigma=6a`, `sigma'=-2p`, `p_y=V'-4ap`, and the
retained physical response `B=p_y/p+sigma''/2`, `mu=(rho B)^(-1)`, give

```
P'=sigma' sigma''/8-V'/2-sigma sigma'/6
   =-(p/2)(p_y/p+sigma''/2)=-p/(2rho mu),
rho_b^2 P'_b=-rho_b p_b/(2mu_b).
```

Monotone exact interval arithmetic in the same actual shell/report boxes
therefore gives `3<rho_b^2 P'_b<5`. The report retains the sharper rational
endpoints. The unscaled `mu`, rather than `Lambda=mu/ell`, is used here.
This nonzero coefficient does not establish a value for the second scalar
projection or exclude its cancellation with the Weyl projection.

Using `1-T^2<=sech^2(q)<=1`, the entire coefficient of `S^2` obeys

```
Cmin=-1/5+(9/20)rho_upper^2 Plo,
Cmax=-(1-T^2)/5+(9/20)rho_lower^2 Phi,
Cmin<=-sech^2(q)/5+(9rho_b^2/20)P''_b<=Cmax<0.
```

The approximate interval is `[-6.019797195441258,-6.019637549940801]`.
If `Ms,Msq` denote the exact computed bounds, positivity of the square
and negativity of this coefficient prove the continuous enclosure

```
Cmin Ms^2 <= K_R(q) <= (27/40)Msq^2.
```

This deliberately coarse box uses no unproved correlation of endpoint
jets or background parameters and no temporal sampling. Its large upper
value is driven by the first physical scalar derivative square at the
registered short packet width. No focusing sign follows.

## A nonzero known scalar particular

Define the nonnegative zero-data Green weight

```
W_T(s)=cosh^2(s)[tanh(T)-tanh(s)], 0<=s<=T,
H=(1-T^2/2)^(-2)>cosh^2(T).
```

The bound on `cosh` follows term by term from `(2j)!>=2^j` in its
power series, then summing the geometric series. Since the derivative
of `tanh` is `sech^2`,

```
(T-s)/H <= W_T(s) <= H(T-s).
```

The actual future waveform's `sign_witnesses.rising` record gives an
actual interval for **ell times S** at `a=5ell/4`. Dividing its upper
endpoint by `ell` and adding the already established past/future scalar
transfer gives

```
S(a) <= -37378.588580306125... < -37376.
```

This uses the actual `ell_S` witness, including its actual-minus-comparator
enlargement, not the comparator-only `H0`. Smooth finite propagation gives
`S(ell)=0`. On `[ell,a]`, Cauchy–Schwarz therefore yields

```
int_ell^a S_q^2 ds >= S(a)^2/(a-ell) > 37376^2/(a-ell),
(27/40)int_0^T W_T S_q^2 ds
  > (27/40)(T-a)/(a-ell) *37376^2/H
  = (27/32)*37376^2/H.
```

The whole negative scalar-square coefficient has magnitude at most
`Cabs=-Cmin`. The first scalar vanishes through `ell`, so

```
int_0^T W_T Cabs S^2 ds <= H Cabs Ms^2 (T-ell)^2/2.
```

Combining these exact inequalities proves

```
int_0^T W_T K_R ds
 > (27/32)*37376^2/H-H Cabs Ms^2 (T-ell)^2/2
 > 1170000000.
```

The exact positive term is approximately `1178689535.9971223`, the
negative cost is at most `0.002523441334521446`, and the resulting lower
bound is approximately `1178689535.9945989`. The upper bound follows from
`K_R<=Khi` and the same first-field quiet interval:

```
int_0^T W_T K_R ds <= H Khi (T-ell)^2/2 < 1040000000000.
```

Thus the known scalar source has a nonzero weighted Green integral. Its
source-only Raychaudhuri particular is the negative of this integral,
strictly less than `-1.17e9`. It is a contribution to the second derivative
at amplitude zero, with no extra Taylor factor `1/2`. It establishes no
sign for the full Ricci source, complete denominator or finite-amplitude
observable.

## A necessary compensation if the second denominator is nonnegative

Assume additionally that the full family and target Gaussian cylinder
exist and that the complete pullbacks and exact identities hold there.
The completion shares the registered first jet, exact quiet collar, anchor
and original normal Gaussian observers, so the zero-data `K_D` identity
applies. An arbitrary second initial intrinsic metric is not admitted by
these hypotheses.
Define the actual unresolved combination

```
Z(q)=rho_b^2 P'_b <psi2(q)>_w-<E2(q)>_w.
```

The retained [quadratic-area identity](inputs/QUADRATIC_AREA_BOUND.md)
gives, in the same amplitude-derivative convention,

```
D2/D0=K_D-int_0^T W_T(s)[K_R(s)+Z(s)]ds,
K_D<5.357e-16.
```

The full nonlinear quiet interval is only `q<=ell/4`. Thus `Z` vanishes
there, while neither unknown is set to zero through `ell`. For
`Zmax=sup_(0<=q<=T)|Z(q)|`,

```
|int_0^T W_T Z ds| <= [H(T-ell/4)^2/2] Zmax <= (8.614e-13) Zmax.
```

Consequently **if** `D2/D0>=0`, its actual unknown projection must supply a
negative weighted integral large enough to offset the known scalar
particular. Exact arithmetic gives the necessary condition

```
Zmax > [known_integral_lower-5.357e-16]/(8.614e-13) > 1.35e21.
```

This is a conditional lower threshold needed for compensation, not an
established value or upper bound for `Zmax`. The actual second scalar and
Weyl projections are coupled outputs, rather than freely selectable
functions or independent preparation choices. The calculation neither
proves their compensation nor excludes it, and does not establish `Cr`,
net focusing, a force, a central zero or an origin model.

Positivity of the physical denominator `D_e` does **not** imply
`D2(T)>=0`. A negative second derivative, even one with large magnitude,
is compatible with a positive actual area for sufficiently small amplitude.
The compensation threshold applies only under the additional stated
second-derivative sign hypothesis.

## Reproduction and validation limits

From this directory run

```
python -I -B evaluate_curvature_source.py --root . --output /tmp/CURVATURE_SOURCE_NORMAL.json
python -I -O -B evaluate_curvature_source.py --root . --output /tmp/CURVATURE_SOURCE_OPTIMIZED.json
```

The evaluator uses the Python standard library and returns JSON-native exact
fraction strings. Floating point fields are explicitly display-only. All-point
polynomial identities and exact source/state/scale checks support the analytic
argument; no finite arithmetic proves the operator existence, regularity,
local uniqueness or target nonlinear-cylinder premises. The capsule check
wrapper separately verifies canonical reports and intended-diagnostic controls.
No old PDE/ODE campaign, second-field solve or finite-amplitude certification
is performed.
