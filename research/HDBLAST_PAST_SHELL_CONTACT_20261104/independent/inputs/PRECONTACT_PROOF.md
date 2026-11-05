# A short-past characteristic precursor of the registered first packet

5 November 2026 UTC. Conditional mathematical continuation. This
note concerns the actual registered **first** master packet and its coupled
linear reconstruction. It does not vary a second initial profile. No old
PDE/background/source program is imported or executed; the accompanying
new verifier performs rational arithmetic on selected retained inputs.

The main result holds on every closed slab `-tau<=x<=0`, `0<tau<1`, and
hence on the open preshell interval `-1<x<=0`. In particular the concrete
quarter slab has support in `3/4<=r<=9/4`. Both normalized-master principal
characteristic channels are nonzero on every such slice, with each share
strictly between `0.499999` and `0.500001`. A compatible pair of nonconstant
past null traces is constructed on a half-unit diamond.

## 1. Fixed actual state, equation and conventions

The retained initial state is

\[
w(0,z)=q_\ell(z),\qquad w_t(0,z)=e^{-z}q_\ell(z),
\qquad \ell=10^{-6},
\]

where `q_ell(z)=exp(4-ell^2/((z-ell)(2ell-z)))` on `(ell,2ell)` and is
zero elsewhere. It is smooth, flat at the endpoints, and nonzero, with
maximum one. This is the same state used by the current canonical and
prepared-family binding. The Oct31 native velocity change was a change
to a second-response approximation, not to this first state. The selected
proof copies in `inputs/` record that distinction.

Write `eta=phi0-1`, `p=phi0_y`, `a=rho_y/rho`, and `dy=-rho dz`. The
physical longitudinal master `V`, not the scalar perturbation itself,
satisfies the retained equation

\[
V_{tt}-V_{zz}+3\tanh(t)V_t-C(z)V_z+
[8\operatorname{sech}^2t-4+(2/3)\rho^2p^2]V=0,
\quad C=\rho(3a+2p_y/p).
\]

Set

\[
\alpha=(\rho/\rho_b)^{3/2}p/p_b>0,\quad
\alpha_z/\alpha=-C/2,\quad
V=\alpha\cosh^{-3/2}(t)w.
\]

Then, by direct differentiation,

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,
\quad \mathcal W=35\operatorname{sech}^2t/4-25/4
+C_z/2+C^2/4+(2/3)\rho^2p^2.                 \tag{1}
\]

Put `U(x,r)=w(ell*x,ell*r)` and `Q(r)=q_ell(ell*r)`. Thus

\[
U_{xx}-U_{rr}+\mathcal V(x,r)U=0,\quad
\mathcal V=\ell^2\mathcal W(\ell x,\ell r),\quad
U(0,r)=Q(r),\quad U_x(0,r)=\ell e^{-\ell r}Q(r). \tag{2}
\]

The first scaled velocity is not zero. Define
`Jplus=U_x+U_r` and `Jminus=U_x-U_r`. In increasing future x,
`Jplus` has principal speed `dr/dx=-1` toward the shell at r=0, and
`Jminus` has speed `+1` away from it. These definitions are retained when
describing a past slice; no orientation is silently reversed.

## 2. Fresh static collar and negative-time coefficient bound

The old narrow-collar prose proof has width `8e-10`, which does not cover
the packet. It is not used as the width premise here. The later selected
`WIDE_BACKGROUND.json` is the exact copy of
`HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/inputs/BACKGROUND_DERIVATIVE_BOUNDS.json`,
whose declared length is `Z=1/125000=8ell`. Its static outer box is

\[
129\le\rho\le131,\quad -9/10^5\le\eta\le-8/10^5,
\quad -7/50000\le p\le-3/25000,
\quad111/1000\le a\le14/125.                         \tag{3}
\]

The new verifier does not merely accept the stored flow boxes. Using the
raw selected `SHELL_BOX.json` and the actual model polynomial

\[
U_{\rm pot}(\eta)=-2/27+(14/9)\eta^2+(50/27)\eta^3
-\eta^4/6-(4/9)\eta^5-(2/27)\eta^6,
\]

it recomputes all four interval images of

\[
\rho_z=-\rho^2a,\quad\eta_z=-\rho p,\quad
p_z=-\rho(U_1-4ap),\quad a_z=\rho^{-1}+\rho p^2/3,
\]

and verifies `I_shell+[0,8ell]F(I_outer)` lies **strictly inside** each
outer component. The four positive lower and upper margins are reported.
Also `131*(8ell)<9`, the retained proper endpoint lower bound. The actual
regular branch and validity of the analytic shell enclosure remain
premises. The ordinary smooth ODE continuation/first-exit argument then
justifies the static box throughout `0<=z<=8ell`. No sampled background
orbit, new initial shell estimate or numerical integration is used.

For completeness the elementary coefficient algebra is explicit. Put
`h=U1/p`, with `Uj` the j-th eta derivative of the potential. The ODE gives

\[
C=\rho(2h-5a),\quad
C_z=\rho^2(2h^2-10ah+5a^2-2U_2-5p^2/3)-5.
\]

Substitution in (1), including the constant `-5`, yields

\[
\mathcal W=\rho^2 A-35\tanh^2(t)/4,\quad
A=2h^2-10ah+35a^2/4-U_2-p^2/6.                    \tag{4}
\]

The new scalar comparisons from the polynomial and (3) are

\[
7/4<h<7/3,\quad31/10<U_2<25/8,\quad1/3<A<6.
\]

For example with `e=-eta` and `d=9/100000`,
`-U1(-e)=28e/9-e^2(50+6e-20e^2+4e^3)/9`; the bracket is positive,
and bounding its negative monomials gives the stated lower h bound.
The U2 remainder from `28/9` is at most
`100d/9+2d^2+80d^3/9+20d^4/9`. Safe A majorants are

\[
A>1933333317/5000000000>1/3,
\qquad A<2680267/450000<6.
\]

These calculations do not import an older positive-time potential cap.
The time dependence in (4) is explicit. For all real t,
`|tanh(t)|<=|t|`, by integrating `0<sech(t)^2<=1` from zero. Therefore,
on the newly declared negative-time rectangle

\[
-\ell\le t\le0,\qquad0\le z\le8\ell,
\qquad5500<\mathcal W<103000,
\quad0<\mathcal V<\mu:=103/10^9.                    \tag{5}
\]

Indeed `129^2/3-35ell^2/4>5500` and
`131^2*6+35ell^2/4<103000`. The background itself is
`dy^2+rho^2[-dt^2+cosh(t)^2*kappa]`, so its smooth negative coordinate
times use the same static positive-depth radial branch. No physical
negative-z background, positive-time evolved-field error, source cap or
future operator norm is transferred to the past.

## 3. Constructive Cauchy development, support and uniform bounds

Fix any `0<tau<1` and put `s=-x`. The required packet development has
support inside

\[
1-s\le r\le2+s,\qquad0\le s\le\tau,              \tag{6}
\]

which is contained in `(0,3)` and hence in the verified static collar.
For a construction one may extend the smooth coefficient boundedly to an
auxiliary whole line, matching it on a neighborhood of this compact
support cone and keeping its absolute bound at most mu. Finite propagation
makes the resulting solution independent of that auxiliary extension.
This introduces no actual background outside the positive-depth collar.

Let `A(s,r)=U(-s,r)`. Its free solution from `(Q,-ell*exp(-ell*r)Q)` is

\[
A_f(s,r)=\tfrac12[Q(r-s)+Q(r+s)]
-\tfrac\ell2\int_{r-s}^{r+s}e^{-\ell y}Q(y)\,dy.
\]

The actual solution is the fixed point

\[
A(s,r)=A_f(s,r)-\tfrac12\int_0^s
\int_{r-s+v}^{r+s-v}\mathcal V(-v,y)A(v,y)\,dy\,dv. \tag{7}
\]

On the closed subspace of continuous functions supported as in (6), the
integral operator has sup-norm at most
`beta_tau=mu*tau^2/2<1`. Thus (7) is a constructive Volterra contraction.
It supplies uniqueness, continuous dependence and the asserted support.
The smooth compact data and smooth coefficient give the smooth linear
solution by the differentiated characteristic/Volterra equations. Higher
coefficient derivatives are bounded on this fixed compact set; they are
not needed in the following quantitative zeroth/first norm estimates.

Translations and convolution on L2 also give

\[
\sup_{0\le s\le\tau}\|A_f(s)\|_2
\le(1+\tau\ell)\|Q\|_2,
\quad
\sup\|A\|_2\le B_\tau\|Q\|_2,
\quad B_\tau=\frac{1+\tau\ell}{1-\mu\tau^2/2}.
\]

The same argument in C0, using `0<=Q<=1`, gives `sup|A|<=B_tau`.
Uniformly over every tau below one, use the limiting constants

\[
\beta=\mu/2<1,\qquad B=\frac{1+\ell}{1-\mu/2}<1.000002. \tag{8}
\]

Solutions for different closed tau agree by uniqueness; their union is
the actual Cauchy extension on `-1<x<=0`. Equation (6) shows no shell
contact there. On the quarter slab its support is `[3/4,9/4]`; on the
half slab it is `[1/2,5/2]`. The limiting front reaches r=0 only at x=-1.
This note makes no assertion of a nonzero jet or boundary interaction at
that endpoint, and no conclusion about earlier past times.

## 4. Both principal channels on every preshell slice

Write `D=||Q_r||_2>0`. The elementary two-endpoint Cauchy estimate on the
unit support interval gives `||Q||_2<=D/2`. Moreover `Q(3/2)=1` and zero
endpoint values give `D^2>=4`, by applying Cauchy separately to each
half interval. These are analytic bump facts, not new bump evaluations.

Since `||ell*exp(-ell*r)Q||_2<=ell D/2`, the initial characteristic norms
lie between `(1-ell/2)D` and `(1+ell/2)D`. The characteristic equations are

\[
(\partial_x-\partial_r)Jplus=-\mathcal VU,
\qquad(\partial_x+\partial_r)Jminus=-\mathcal VU.
\]

Their past characteristic integrals and (8) imply

\[
\|J_\pm(x)-\hbox{translated }J_\pm(0)\|_2
\le\mu|x|B\|Q\|_2\le\mu B D/2.
\]

Consequently throughout `-1<x<=0`,

\[
(1-\delta)D\le\|J_\pm(x)\|_2\le(1+\delta)D,
\quad\delta=\ell/2+\mu B/2<10^{-6}.                \tag{9}
\]

In particular both norms are strictly positive. For
`Eplus/minus(x)=||Jplus/minus(x)||_2^2/4`, monotonicity of a ratio in
its two nonnegative arguments gives

\[
\frac{(1-\delta)^2}{(1-\delta)^2+(1+\delta)^2}
\le\frac{E_\pm}{Eplus+Eminus}\le
\frac{(1+\delta)^2}{(1-\delta)^2+(1+\delta)^2}.
\]

The exact rational guards show both shares are strictly in
`(499999/1000000,500001/1000000)`. Thus a state satisfying either
`Jminus=0` or `Jplus=0` on an entire slice in this finite interval cannot
be the given actual packet's past development.

These are precisely the **normalized-master principal quadratic weights**.
The unscaled t,z weights differ by a common factor `1/ell`, so their
shares agree. A spatially varying conjugation does not preserve these
weights. They are not physical canonical/Noether energies or radiation
fluxes, and need not be conserved under the variable potential. The result
does not exclude more remote incoming histories followed by scattering,
earlier shell interaction, another physical notion of incoming data, or
an explicitly specified source.

## 5. Two nontrivial compatible past null traces

There is a concrete linear precursor, beyond merely an obstruction to a
pure channel. Take the diamond with past apex `(-1/2,3/2)` and present base
`x=0,1<=r<=2`. In null coordinates `u=x+r`, `v=x-r`, it is

\[
1\le u\le2,\quad-2\le v\le-1,\quad u+v\le0.
\]

Its two past faces are `u=1` and `v=-2`. Restrict the already constructed
actual solution to them, defining `f(v)=U(1,v)` and `g(u)=U(u,-2)` in
these coordinates. They obey the common-corner condition `f(-2)=g(1)`.

Put `beta_half=mu/8` and
`B_half=(1+ell/2)/(1-beta_half)`. Equation (7) gives
`|U-U_free|<=beta_half*B_half` throughout the half slab. At the face
midpoints `(-1/4,5/4)` and `(-1/4,7/4)`, the free value is at least
`1/2-ell/4`. Hence each actual midpoint exceeds `0.499999`.

At the common past apex the free displacement terms vanish, so

\[
U_f(-1/2,3/2)=-\tfrac\ell2\int_1^2e^{-\ell r}Q(r)\,dr.
\]

On `[5/4,7/4]`, `Q>=exp(-4/3)>1/4`. One exact Taylor bound is

\[
e^{4/3}\le1+4/3+(4/3)^2/2+
\frac{(4/3)^3/6}{1-(4/3)/4}=103/27<4.
\]

Also `exp(-ell*r)>=1-2ell`. Thus the exact guards give

\[
U(-1/2,3/2)\le-\ell(1-2\ell)/16+
\beta_{1/2}B_{1/2}<-1/25000000.                 \tag{10}
\]

Each actual face is therefore negative at the common apex and positive
at its midpoint: both are nonconstant and have an interior zero. Positive
master conjugations preserve these sign changes, unlike the principal
weight ratios. The corresponding face derivatives are nonzero somewhere;
their exact normalization is `2 f_v=Jminus` along `u=1` and
`2 g_u=Jplus` along `v=-2`.

These compatible traces construct the solution by a genuine characteristic
Volterra problem. Here `U_uv=-mathcal V U/4`, so

\[
U(u,v)=f(v)+g(u)-c-\tfrac14\int_1^u\int_{-2}^v
\mathcal V((a+b)/2,(a-b)/2)U(a,b)\,db\,da.
\]

The integration rectangle stays inside the diamond. Its area is at most
`1/4`, since `(u-1)+(v+2)<=1`. Hence the sup-norm contraction factor is
at most `mu/16<1`. Smoothness and uniqueness identify its solution with
the actual Cauchy development and reproduce the present data, including
the appropriate derivative. This is an existence construction of the
selected **linear** characteristic data by restriction and a converse
uniqueness check. It is not an independently selected physical preparation
rule, arbitrary characteristic-data surjectivity, or a full nonlinear
null/shell constraint solution.
This is an interior first-master diamond: its `v=-2` face is not a
shell-corner cut, and its nonzero corner is not a chosen second scalar datum.

## 6. Coupled physical reconstruction and exact-family limits

The scalar master is not an unrelated test field. The retained all-point
differential identities reconstruct

\[
P=VY,\quad\Psi=P/\rho^2,\quad
h_{yy}=-4\Psi,\ h_{ya}=0,\ h_{ab}=2\Psi h_{0,ab},\quad
\chi=-3P_y/(\rho^2p).
\]

These satisfy the full first linearized Einstein-scalar equations and
constraints wherever (1) holds. The identities use explicit smooth
hyperbolic functions and the same radial ODE, so apply to the newly
constructed negative-time rectangle. The first perturbation is zero in
an open shell collar for every time under consideration; its physical
linear junction rows are then satisfied there without shell forcing.
Its original canonical tangent and its nonzero momenta are unchanged.

Identification with the fixed common-harmonic representative retains the
full initial one-jet, including the mixed metric velocity, and the linear
gauge construction/uniqueness premises. No old future gauge approximation
or source-error cap is imported to obtain a past bound.

A conditional exact-family interpretation can be formulated by applying
the retained local development theorem to the momentum-reflected exact
constraint curve `(gamma,pi,phi,varpi)->(gamma,-pi,phi,-varpi)` and then
reflecting time. The constraints and pure-tension shell equations have
the required time-reversal covariance. This is a new conditional use of
analytic existence, not a numerical backward Einstein solve or an
effective uniform nonlinear parameter window. In particular the full
nonlinear constraint correction may have a larger support annulus than
the first packet, so its contact-free time must not be inferred from (6).
The whole open preshell interval asserted here is a **first linear**
packet statement, not a support theorem for all parameter orders.

No second constrained state, past physical incoming law, injected energy,
matter or radiation production, absolute proper-observable sign, physical
scale or higher-dimensional origin is supplied. The open interval result
is a finite past principal-purity obstruction, not an all-past no-go.

## 7. Exact arithmetic, local inputs and review

`verify.py` is self-contained standard-library rational arithmetic. Its
only runtime input closure is the five exact local copies listed in
`SOURCE_PINS.json` (105859 bytes altogether). Two numeric JSON inputs
enter the fresh static-collar arithmetic; three selected proof snapshots
record the model, reconstruction and current-packet identity premises.
No transitive old code is executed. Analytic statements in those snapshots
remain analytic premises rather than claims discharged by hash checks.

The verifier authenticates the selected input bytes before and after
arithmetic, reports all four fresh ODE-image inclusions and exact scalar
constants, and writes only to a fresh absolute output outside its source
directory. Normal and optimized Python must have identical JSON bytes:

```
python -I -S -B verify.py --output ABS_EXTERNAL/NORMAL.json
python -I -S -B -O verify.py --output ABS_EXTERNAL/OPTIMIZED.json
```

Activity is reported explicitly: zero project imports, old checker calls,
background/source/bump queries or numerical PDE solves; four exact
interval ODE images are recomputed per successful verifier process.
The independent coefficient audit separately rederived the fresh wide
collar and negative-time bounds. A separate independent past-characteristic
review accompanies the final delivery; no old full replay is implied.
