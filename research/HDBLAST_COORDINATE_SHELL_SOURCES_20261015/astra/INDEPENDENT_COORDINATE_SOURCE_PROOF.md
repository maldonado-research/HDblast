# Independent derivation of the complete coordinate boundary sources

15 October 2026. Conditional mathematical checkpoint for the registered
Einstein–scalar shell model. This note independently derives the known quadratic
boundary sources in the specified first longitudinal representative. Its exact
arithmetic checks retain the prior analytic waveform certificate as an input;
they do not re-prove that upstream PDE estimate or construct the second-order
response.

The inputs are the actual public equation notes and waveform data copied under
`inputs/`, with their immutable source links in `inputs/INPUT_PINS.json`.
`check_coordinate_sources.py` independently pins ten of those inputs in its
source, verifies their byte hashes, and checks exact path/hash/size membership
in the retained physical, contact, and inertial manifests. It imports only the
Python standard library, uses `Fraction` for every mathematical decision, and
does not import any producer, prior checker, or provisional future-work code.

## 1. Geometry, normal orientation, and derivative coefficients

The original physical reconstruction note, equations (1), (2), (6), (9), and
(13), and the original scalar boundary note, equations (1)–(8), fix the
conventions. On both copies the increasing proper-normal coordinate points
outward toward the common shell. The copy-even one-copy junctions are

\[
K_{ab}=\frac{\sigma}{6}h_{ab},\qquad n\phi=-\frac{\sigma'}2,
\quad a_b=\frac{\sigma_b}6,\quad p_b=-\frac{\sigma'_b}2.
\]

The derivative convention is
`g=g0+epsilon*h1+epsilon^2*h2/2`, with the same convention for the scalar.
The fixed first representative has zero first shell graph and

\[
(h_1)_{yy}=-4uY,\quad (h_1)_{ya}=0,\quad
(h_1)_{ab}=2uYh_{0,ab},\quad \chi_1=vY,
\quad Y=n_1^2-n_2^2.
\]

For the linear first-variation path, the inverse unit lapse expands as

\[
\mathcal N^{-1}=(1-4\varepsilon uY)^{-1/2}
=1+2\varepsilon uY+6\varepsilon^2u^2Y^2+O(\varepsilon^3).
\]

Multiply this by the exact normal numerators
`a+epsilon*(u_y+2*a*u)*Y` and `p+epsilon*v_y*Y`. Twice their quadratic
coefficients give

\[
D^2K=(4uu_y+20au^2)Y^2h_0,\qquad
D^2(n\phi)=(4uv_y+12pu^2)Y^2.                 \tag{A1}
\]

These are unit-normal variations. The lapse terms cannot be removed by
using the background vector `partial_y` as the perturbed unit normal.
Differentiating the junction right sides gives the material quadratic terms

\[
\left(\frac{2\sigma'}3uv+\frac{\sigma''}6v^2\right)Y^2h_0,
\qquad -\frac{\sigma'''}2v^2Y^2=-2v^2Y^2,        \tag{A2}
\]

where the registered tension has `sigma'''=4`. Moving (A1) to the right of
the Eulerian second equations and using the actual first shell rows

\[
u_y=-2au-\frac p3v,\qquad v_y=-2pu-\frac{\sigma''}2v
\]

then yields, without radial numerical differentiation,

\[
\boxed{S_p=-4pu^2+2\sigma''uv-2v^2,\qquad
S_K=\frac{\sigma''}6v^2-12au^2.}              \tag{A3}
\]

In the metric row the `uv` contributions cancel using `sigma'=-2p`.
The checker derives these identities with an independent exact polynomial
ring, including the inverse-normal normalization, before comparing them
with (A3). A non-derivative Taylor convention would introduce a factor two;
an outward sum of the two scalar junctions would introduce another copy
normalization. In particular the previously published outward-sum material
source is `-4*v^2*Y^2`, whereas the material summand in (A3) is `-2*v^2`.

The second graph, lapse, shift, and gauge/constraint variables remain in the
linear second boundary operator of the original note. Formula (A3) is the
complete **known quadratic coordinate source** for the fixed first fields;
it is not the complete second boundary solution. A change of first gauge
requires the full associated bulk and boundary source conversion. A
second-only gauge transformation leaves the fixed first source unchanged.

## 2. Actual-H reconstruction and conservative global bounds

The public reconstruction has `u=V/rho^2` and

\[
V=\alpha(z)\cosh(t)^{-3/2}w,\qquad
\alpha=(\rho/\rho_b)^{3/2}p/p_b,\quad\alpha(0)=1.
\]

Consequently the shell amplitudes in precisely this representative are

\[
u_b(t)=\rho_b^{-2}\cosh(t)^{-3/2}w_b(t),\qquad v_b(t)=S(t).
\]

The first scalar amplitude equals the intrinsic-clock first trace because
the background shell scalar is constant. This does not make the metric
amplitude `u` gauge invariant.

From the public actual trace `H(x)=ell*(w_z-gamma*w_b)(ell*x)` and dynamic
boundary condition

\[
-w_z+\mu w_{tt}+\left[\gamma+\mu
(35\operatorname{sech}^2t/4-25/4)\right]w_b=0
\]

the `gamma` terms cancel exactly. Setting `y(x)=w_b(ell*x)` gives

\[
y''+A(x)y=\frac\ell\mu H(x),\qquad
A(x)=\ell^2[5/2-(35/4)\tanh^2(\ell x)],
\quad y(0)=y'(0)=0.                              \tag{A4}
\]

The zero initial shell jets follow from compact support away from the shell.
The force is `ell*H/mu`. It is the actual trace, not the free comparator.
This equation reconstructs the trace of the retained coupled solution.

For `ell=10^-6` and `0<=x<=L=3`,

\[
|A|\le M=\frac52\ell^2,\qquad
q=\frac{ML^2}2=\frac9{800000000000}<1.
\]

Twice integrating (A4) gives `y=F-K_A y`, where
`F=(ell/mu)*integral_0^x (x-s)H(s)ds`. The prior 4353 closed cells partition
`[0,3]`; their actual-H intervals are whole-cell bounds, not sample values.
With `M_j=max(abs(H_j,lo),abs(H_j,hi))`, define

\[
B_H=\sum_j M_j\left[3(b_j-a_j)-\frac{b_j^2-a_j^2}2\right],\quad
B_y=\frac{\ell B_H}{\mu_{\min}(1-q)},\quad E=qB_y.
\]

Then `||y||<=B_y`, `||y-F||<=E`, and `||u||<=U=H_squared_max*B_y`.
Here the retained JSON field `H_squared` is the background `rho_b^-2`, not
the waveform square. The exact checker independently reproduces

\[
B_H\simeq3.005488782829515,\quad B_y\simeq0.0014591371894211345,
\quad E\simeq1.6415293380987764\,10^{-14},
\quad U\simeq8.643987698350998\,10^{-8}.           \tag{A5}
\]

All displayed decimals in this note are explanatory. Reduced rational
endpoints in the audit receipts decide every inequality.

Differentiating the integral equation, without differentiating the data
cells, gives

\[
\|y'\|\le\frac\ell{\mu_{\min}}
\sum_jM_j(b_j-a_j)+MLB_y\simeq0.0009727646685349448.
\]

For background proper time `tau=rho_b*ell*x`,

\[
\partial_\tau u=\frac{\rho_b^{-2}\cosh(\ell x)^{-3/2}}
{\rho_b\ell}\,[y'-(3/2)\ell\tanh(\ell x)y].       \tag{A6}
\]

It follows conservatively that `||partial_tau u||<0.000443543`. The sharper
producer enclosures use signed convolution extrema. This independent bound
uses absolute history sums and does not rely on those extrema.

## 3. Two distinct zeros of each full coordinate coefficient

Let `A=1.49998`, `B=1.49999`. We need a positive `u` bound throughout the
entire bracket because the actual scalar zero is not a chosen sample point.

Split every old cell at `A` and `B`. For a portion `[l,r]` below `A`, its
convolution contribution is bounded below by

\[
\begin{cases}
H_{j,\mathrm{lo}}[A(r-l)-(r^2-l^2)/2],&H_{j,\mathrm{lo}}\ge0,\\
H_{j,\mathrm{lo}}[B(r-l)-(r^2-l^2)/2],&H_{j,\mathrm{lo}}<0.
\end{cases}
\]

For a portion in `[A,B]`, use
`min(0,H_j,lo)*[B(r-l)-(r^2-l^2)/2]`. This permits partial traversal or no
traversal of that portion. The sum `I_-` bounds the actual triangular
convolution below for **every** `x` in `[A,B]`. Exact sums give

\[
I_-\simeq0.19151209523222054,\qquad
y\ge\frac\ell{\mu_{\max}}I_--E
\simeq8.430484535909867\,10^{-5}>0.                 \tag{A7}
\]

The factor `mu_max` is necessary for this positive lower bound. The clock
factor satisfies `1-t^2 <= cosh(t)^(-3/2) <= 1` on the window: the positive
series gives `cosh(t) <= (1-t^2/2)^(-1)`, and convexity gives
`(1-t^2/2)^(3/2) >= 1-3t^2/4 >= 1-t^2`. Hence throughout the bracket

\[
u\ge u_->4.9942\,10^{-9}>0.                       \tag{A8}
\]

The pinned waveform report contains opposite-sign actual-H endpoint
witnesses, so the retained continuity theorem supplies at least one actual
scalar zero `x_*` in `(A,B)`. The gain in `v=(G/ell)*cosh^(-3/2)*H` is
strictly nonzero. At **every** such scalar zero,

\[
S_p(x_*)=-4p_bu(x_*)^2
\ge 1.2405064080821801\ldots\,10^{-20}>0,
\]
\[
S_K(x_*)=-12a_bu(x_*)^2
\le-3.333587486002262\ldots\,10^{-17}<0.           \tag{A9}
\]

For the endpoints themselves, divide the retained `ell_S` witness intervals
by the exact `ell`. They imply approximately

| Endpoint | `v` interval | upper bound on `Sp` | lower bound on `SK` |
|---|---|---|---|
| `A=1.49998` | `[-2.2100977466263743,-1.5506779229486782]` | `-4.809202513239786` | `1.6029281737218288` |
| `B=1.49999` | `[1.1881154154770994,1.8102564163266304]` | `-2.8232352292631044` | `0.9409967340666433` |

The exact witness intervals and exact margins are saved in
`central_endpoint_source_signs` in both receipts. The enclosure formula is,
for `m<=|v|<=V` at an endpoint,

\[
S_p\le4|p|_{\max}U^2+2\sigma''_{\max}UV-2m^2<0,
\quad S_K\ge\frac{\sigma''_{\min}}6m^2-12a_{\max}U^2>0.
\]

Thus **each of `Sp` and `SK` has at least two distinct zeros entirely inside
`(A,B)`**: one lies in `(A,x_*)` and another in `(x_*,B)`, by continuity and
the strict signs. The zero intervals for the two coefficients need not
coincide. No uniqueness, exact location, derivative at a zero, or count of
all source zeros is inferred. The sign statements concern scalar
coefficients; their angular fields vanish identically where `Y=0`.

The same independent conservative global bound (A5), combined with the
retained positive tail H interval and negative physical gain, gives on all
of `2<=x<=3`

\[
0.0003210861041795985\ldots\le|v|
\le0.0005816479561500038\ldots,
\]
\[
\boxed{S_p<-2.0579\,10^{-7},\qquad S_K>6.8724\,10^{-8}.} \tag{A10}
\]

These coefficient sign reversals also show why a pointwise inference from
the material square alone fails even when a global norm correction is tiny.
They establish no sign for the second scalar or metric response.

## 4. Harmonic factors and proper-volume norms

Put `J=n1^2+n2^2`, `F2=J-1/2`, and `H4=Y^2-J/2+1/12`. Direct polynomial
identities and sphere moments give

\[
Y^2=\frac16+\frac12F_2+H_4,
\quad \Delta F_2=-8F_2,\quad\Delta H_4=-24H_4.
\]

The checker verifies harmonicity using homogeneous extensions to `R^4`
and computes moments with
`E[prod n_i^(2k_i)]=prod(1/2)_(k_i)/(2)_(sum k_i)`. The three terms are
orthogonal, with mean-square factors

\[
\frac1{36},\qquad\frac1{48},\qquad\frac{19}{720},
\quad\text{sum}=\frac3{40}.                       \tag{A11}
\]

Their power shares are `10/27`, `5/18`, and `19/54`. Ordinary unit-sphere
volume is `2*pi^2`; thus `integral Y^4=3*pi^2/20` while
`integral H4^2=19*pi^2/360`.

The background shell proper-volume measure is
`rho_b^3*cosh(t)^3*dOmega*d_tau` and `d_tau=rho_b*ell*dx`. Consequently

\[
\|S_pY^2\|_{L^2(d\mu_0)}^2
=\frac{3\pi^2}{20}\rho_b^4\ell
\int_0^3\cosh^3(\ell x)S_p(x)^2dx,                \tag{A12}
\]

with the same scalar-coefficient norm for `SK*Y^2`. This is not a tensor
energy or a contraction of the full metric source tensor. The degree-four
Robin source `EN*H4` instead uses `19*pi^2/360`.

The independent checker validates these factors against every source and
material/correction norm in the new report. It encloses pi independently
using Machin's identity and alternating arctangent series with 70 and 20
terms. The separate Sol audit reconstructs the entire cell quadrature;
this audit checks the factor and ratio joins. In particular, it verifies
exactly

\[
\frac{\|\delta S_pY^2\|^2}{\|-2v^2Y^2\|^2}
\le7.406446493642206\ldots\,10^{-25}
<(9\,10^{-13})^2,
\]

so the **unsquared** norm ratio is less than `9e-13`. The numerator uses
its upper endpoint and the positive denominator its lower endpoint. This
does not bound the pointwise ratio near the scalar zeros, where the
material summand vanishes and the full coordinate coefficient is nonzero.

## 5. Robin joins and a degree-two distinction

For the retained degree-four scalar system, the original boundary note
gives `C4=16u^2-(4/3)v^2` and `EN4=C4_y-4SK`. The exact radial rows yield

\[
E_{N,4}=-16au^2-\frac{16}3puv+\frac23\sigma''v^2. \tag{A13}
\]

The same note gives the symmetrizer `D=diag(3,8)` and
`det(DM)=24a*sigma''-16p^2`. Independent endpoint boxes show its determinant
exceeds `10.69129` and its leading principal minor is positive on this
registered shell box. This is positivity of a boundary matrix, not a
proof that the amplitudes are unconstrained hyperbolic channels or an
Einstein–scalar energy estimate.

There is a useful distinction for degree two. The original bulk note,
equations (13), (21), and (23), gives the full derivative angular trace-free
source as

\[
12u^2\operatorname{Hess}(Y^2)^{\rm TF}
-(12u^2+2v^2)(dY\otimes dY)^{\rm TF}.
\]

Its degree-two projections use a factor `1/2` for both Hessian structures.
Their `u^2` terms cancel, leaving `sT2=-v^2`. Since the linear row is
`-C2/2=-sT2`, one obtains `C2=-2v^2`. Therefore

\[
E_{N,2}=(C_2)_y-4(S_K/2)
=24au^2+8puv+\frac53\sigma''v^2,                \tag{A14}
\]

in the `F2` basis. It is **not** `EN4/2`. The scalar physical boundary
sources (A3) still project as `1/6,1/2,1`; the distinction enters through
the bulk algebraic row used to form the Robin reduction. Equation (A14)
is an algebraic appendix here; no new degree-two numerical table, degree-zero
master reduction, or coupled second evolution is supplied.

## 6. Replay, controls, and scientific boundary

Run from the capsule root:

```sh
python -I -B astra/check_coordinate_sources.py --output astra/ASTRA_COORDINATE_NORMAL.json
python -I -O -B astra/check_coordinate_sources.py --output astra/ASTRA_COORDINATE_OPTIMIZED.json
python -I -B astra/check_coordinate_sources.py --run-controls --output astra/ASTRA_COORDINATE_CONTROLS.json
```

Normal and optimized positive results contain 83 explicitly named check
groups and are byte-identical. Seventeen mathematical controls run in both
modes: normal normalization, copy factor, Taylor factorial, metric source
sign, angular measure, proper-volume radius power, scaled force, Volterra
resolvent, inverse-mu endpoint, tail square, proper derivative, comparator
substitution, gamma cancellation, metric radius, central endpoint square,
central metric erasure, and norm-ratio denominator. Each computes the wrong
quantity and must fail its named mathematical gate with exit code one.
Python assertions are not used; `-O` does not disable any gate.

The final receipts distinguish immutable prior inputs from the newly checked
source table and norm report. All signs and norms refer to the fixed
one-copy coordinate source convention. The work supplies neither the bulk
quadratic source throughout the interior, the second initial correction,
the coupled second-order solution, an effective nonlinear remainder, nor
an invariant backreaction or cosmological-origin claim.
