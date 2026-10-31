# A conditional TT screen coefficient at a shell-influenced point

At the registered point

\[
 P_\dagger=(x,r)=(5/4,1/8),\qquad \ell=10^{-6},\qquad t=\ell x,
\]

the complete scalar transverse-traceless (TT) response gives a strictly
negative coefficient for the selected screen magnetic tensor. The primary
outward rational enclosure is

\[
\boxed{\quad
h(P_\dagger)\in\left[
-\frac{1049651791613786574257886764672101}
       {2596148429267413814265248164610048},
-\frac{18954240643899648600527883239121}
       {1298074214633706907132624082305024}
\right]\subset[-0.404311163330,-0.014601815851].\quad}
\]

The displayed decimal interval is rounded outward. Exact rational endpoints,
inputs and error charges are in [the primary report](main/POSTCONTACT_NORMAL.json).
This calculation encloses the actual scalar formal response, including its
full bulk feedback and the retained inhomogeneous Neumann history. It also
includes the complete fixed first-field quadratic curvature term.

The earlier [screen-magnetic checkpoint](../HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/)
already enclosed a negative coefficient at `(x,r)=(1/4,5/4)`, whose backward
cone is separated from the shell. The present point is different:
`x-r=9/8>1`, so its registered backward problem samples shell history after
the initial support edge can reach the shell. This is a new bounded point
calculation using established Green-kernel, contraction and curvature
methods. It makes no claim to a new general mathematical principle or to
external novelty.

## Assumptions and provenance

The assertion is conditional on the retained Einstein–scalar model and
actual first-solution premises. In particular:

1. The fixed regular background branch at `delta=1/1000` exists on the
   required collar, with the original potential and tension, positive
   `rho`, and the stated background field equations and Hamilton constraint.
2. The retained master and gauge solutions exist with the regularity,
   physical shell conditions and compatible corner jets used by their
   differentiated integral identities. The global actual-minus-comparison
   error certificate encloses derivatives of these same functions through
   total scaled order three.
3. The complete common-gauge first metric and scalar are used in the
   registered bulk and shell source maps. At the physical shell, values and
   derivatives are physical-side one-sided traces. The zero initial collar
   supplies the required compatibility, including `G_T(0)=0`.
4. The TT quantity is the solution of the scalar integral equation below,
   with the registered zero TT Cauchy data. Identifying it with the TT
   projection of a complete physical second coefficient additionally retains
   the stated initial-constraint, harmonic-gauge and subsidiary-propagation
   premises.
5. The observer, screen, clock and angular orientation are the registered
   prescription stated below. No replacement by a transported observer is
   assumed.

The principal source chain is:

| Input | Retained source |
|---|---|
| Actual branch and collar equations | [Background field equations and Hamilton constraint](../HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015/background/inputs/PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md) |
| Actual global common-gauge first jets and differentiated errors | [Geometric-source numerical inputs](../HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/actual_first_jets.py), with its authenticated predecessor certificates |
| Complete 23-component operator and physical adapters | [Original linear operator](../HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py) and [source adapters](../HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/source_adapters.py) |
| Collected bulk TT forcing | [Original coefficient maps](../HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json) |
| Original selected-screen curvature map | [Screen derivation](../HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/DERIVATION.md) and [coefficient implementation](../HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py) |
| Current implementation-to-physical-equation joins | [Registered initial constraints](../HDBLAST_REGISTERED_INITIAL_CONSTRAINTS_20261030/), [natural first-field binding](../HDBLAST_NATURAL_FIELD_BINDING_20261030/), and [bulk/shell constraint binding](../HDBLAST_BULK_SHELL_CONSTRAINT_BINDING_20261030/) |

Source hashes and closed capsule inventories are checked before scientific
imports and after execution as specified by each runner. These checks
authenticate the calculation's inputs; they do not independently prove the
underlying analytic existence and regularity premises.

## The actual scalar TT subsystem

Both TT slots, indices 21 and 22 in the registered 23-component order, are
autonomous in every row and column of all 25 primitive matrices. Their
diagonal entries agree. Actual assembly gives, for either orientation,

\[
 U_{xx}-U_{rr}=A_TU_x+B_TU_r+D_TU+F_T,\qquad
 U_r(x,0)=G_T(x),\qquad U(0,r)=U_x(0,r)=0,
\]

where

\[
\begin{aligned}
A_T&=-3\ell\tanh(t),& B_T&=-3\ell\rho a,\\
D_T&=\ell^2\left[-\frac6{c^2}
+\rho^2\left(-2a^2+2h_0^2-2b-8d+\frac43V\right)\right],
\end{aligned}
\]

with `c=cosh(t)`, `h_0=tanh(t)/rho`, `a=rho_y/rho`,
`b=-(a_y+a^2)` and `d=rho^-2-a^2`. The TT Robin matrix is identically zero.
The symbol `h_0` here is a background coefficient, distinct from the final
screen coefficient `h`.

The actual background identities

\[
a_y=-\rho^{-2}-p^2/3,\qquad
V=p^2/2-6a^2+6\rho^{-2}
\]

reduce the expression exactly to `D_T=-8ell^2/cosh(t)^2`. This consequence
uses both the field equation and Hamilton constraint. The primary numerical
bound retains the full original expression for `D_T`; it does not infer
correlations between separate interval boxes by simplifying them after
evaluation.

The two complete `metric_forcing/ab/5,6` lists are equal and contain 78 terms
each. Exact collected-map algebra verifies their scaled units and sign,

\[
F_T=-2\ell^2\rho^2(\mathrm{QEred})_T.
\]

The generally nonzero physical shell formula is

\[
G_T=-\frac83\ell\rho
       \left(a h_s+\frac{\sigma'v_s}{3}\right)d_s
    =-\frac43\ell\rho
       \left(a A-\frac{2p\varphi_1}{3}\right)M.
\]

The equalities use `A=h_s`, `M=2d_s`, `varphi_1=v_s` and `sigma'=-2p` on
the shell. They agree with the original conversion `G_T=2ell rho QKab_T`.
Here the scalar amplitude `varphi_1` is not an integration-time variable.
The evaluated history and endpoint intervals include zero. Thus this
calculation retains a source which is generally nonzero as a formula; it
does **not** establish that its actual value at the endpoint is strictly
nonzero, nor does it establish a source-zero identity.

Equal operators, sources and zero TT Cauchy jets imply equal metric TT
orientations for the unique scalar solution. The screen curl acts with
opposite orientation eigenvalues, so this equality does not cancel the
selected magnetic coefficient.

## Reflected Green integration and the physical boundary sign

The registered scaled radius increases inward:

\[
\partial_y=-\frac1{\ell\rho}\partial_r.
\]

For a physical scalar forcing `f(s,v)` on `v>=0`, define

\[
(K_Nf)(x,r)=\frac12\int_0^x
 \int_{r-x+s}^{r+x-s} f(s,|v|)\,dv\,ds.
\]

This is the homogeneous-Neumann wave kernel. Its first derivatives are

\[
\begin{aligned}
(K_Nf)_x&=\frac12\int_0^x
 [f(s,r+x-s)+f(s,|r-x+s|)]\,ds,\\
(K_Nf)_r&=\frac12\int_0^x
 [f(s,r+x-s)-f(s,|r-x+s|)]\,ds.
\end{aligned}
\]

There is no extra sign on the folded endpoint in the time derivative.
For the convention `U_r=G_T`, the causal boundary lift is

\[
(LG_T)(x,r)=-\mathbf1_{x>r}\int_0^{x-r}G_T(s)\,ds,
\quad (LG_T)_x=-G_T(x-r),\quad (LG_T)_r=+G_T(x-r).
\]

These signs supply the correct positive-inward radial trace. The compatible
corner condition makes the lift C1 across `x=r`. Negative integration depth
is an image construction; it is not a second physical branch. The kernel
reflects the **complete physical composite forcing**. If factors are extended
separately, `B_T` and `U_r` must both be odd; reflecting only `U_r` and making
`B_T` even changes their product and changes the equation.

Set `T=5/4`, `R=1/8` and `L=T-R=9/8`. The backward physical set is

\[
\mathcal C_\dagger=
\{0\le s\le T,\ \max(0,s-L)\le v\le T+R-s\}.
\]

It reaches boundary history `[0,L]` and has initial radial footprint
`[0,T+R]=[0,11/8]`. It is backward closed: from `(t,r)` in this set, earlier
reflected source points obey

\[
\max(0,r-t+s)\le v\le r+t-s,
\]

and `r+t<=T+R`, `r-t>=-L` place them in the same set at time `s`.
Every boundary dependency has time at most `t-r<=L`. This closure is needed
for the feedback bound, not only for the first source integral.

The existing global first-field API covers every rectangle in
`0<=x<=3, r>=0, x+r<=8`. It therefore covers the cone, the shell history and
the primary rectangular source overcover `[0,T]x[0,T+R]`. Its full global
differentiated errors retain the reflected branches and dynamic shell lift.
The sharper old exterior API requires `r>=x` and is not evaluated at this
new point.

## A complete C1 feedback allowance

The scalar actual target solves

\[
U=LG_T+K_NF_T+K_N(A_TU_x+B_TU_r+D_TU).
\]

Let `V=LG_T+K_NF_T`. This is an integration seed for an enclosure of `U`.
It is distinct from the retained finite approximation
`W1_T=LG_T+K_N(La_T LG_T+F_T)`.

On `Cdagger`, use

\[
X(u)=\max\left(\frac2T\sup|u|,\sup|u_x|,\sup|u_r|\right).
\]

Direct absolute integration gives
`X(K_N f)<=T sup|f|` and `X(LG_T)<=2sup|G_T|`. If uniform actual bounds are

\[
|A_T|\le a_*,\quad |B_T|\le b_*,\quad |D_T|\le d_*,\quad
|F_T|\le f_*,\quad |G_T|\le g_*,
\]

then

\[
\sigma=a_*+b_*+Td_*/2,\qquad q=T\sigma<1
\]

is a contraction factor. Consequently the scalar integral problem has a
unique C1 solution with the prescribed data, and

\[
X(U)\le\frac{2g_*+Tf_*}{1-q},\qquad
\boxed{\ |U_x(P_\dagger)-V_x(P_\dagger)|
       \le\frac{q(2g_*+Tf_*)}{1-q}.\ }
\]

The denominator sums the complete scalar feedback series. This estimate
does not discard feedback after a finite iterate and does not promote a
finite `W1` value to the actual solution.

The frozen calculation gives the following explanatory decimal displays:

| Bound | Value |
|---|---:|
| Uniform bulk source cap `f_*` | `0.0826459738028` |
| Uniform shell source cap `g_*` | `1.32339556112e-21` |
| Full scalar contraction `q` | `5.42651760331e-5` |
| Complete derivative feedback allowance | `5.60630212302e-6` |
| Actual `U_Tx(Pdagger)` | `[-0.004549931453471,-0.000164323537266]` |

The exact endpoint derivative seed is

\[
V_x(P_\dagger)=-G_T(9/8)+\frac12\int_0^{5/4}
 [F_T(s,11/8-s)+F_T(s,|s-9/8|)]\,ds.
\]

Whole-box source ranges multiplied by exact segment lengths enclose these
integrals. They are not point quadrature values with an uncharged remainder.
The lower characteristic is split at `s=9/8`. No additional `ell` factor
is inserted: the source map already has the scaled PDE units.

## Complete curvature map and the observation convention

Write `Rphys=rho cosh(t)` and `theta=tanh(t)`. For the complete first metric
frame amplitudes `A,B,C,D,E,L,M`, define

\[
\begin{aligned}
F_{\rm quad}={}&
\frac{(D+4M-L)M_t-2ML_t}{4\rho}
-\frac{E(D+L)}{2R_{\rm phys}}\\
&+\left(\frac C{R_{\rm phys}}-\frac{M_y}2\right)
\left[\frac{R_{\rm phys}(C_t-\theta C)}{2\rho}
-\frac{R_{\rm phys}E_y}2+\frac B2\right].
\end{aligned}
\]

The derivative suffixes are ordinary physical coordinate derivatives. In
particular `u_t=u_x/ell`; the background observer derivative is
`e_0u=u_x/(ell rho)`. These quantities must not be interchanged.

For `g(lambda)=g0+lambda q+lambda^2 k/2`, the selected tensor is the second
parameter derivative itself:

\[
\Pi_{\rm pair}H''(0)=h(T_+-T_-),\qquad
h=\frac{3U_{Tx}}{2\ell\rho^2\cosh(t)}
  +\frac{4F_{\rm quad}}{R_{\rm phys}}.
\]

There is no additional factor of one half. `Fquad` is a kinematic first-field
quadratic curvature term, distinct from the bulk forcing `F_T`. It retains
the observer, radial-screen, connection and epsilon raising/volume
variations. At `Pdagger` its positive contribution is enclosed, in decimal
display, by `[3.48148648036e-18,1.42037292928e-17]` and is included in `h`.

The observer is the future unit normal to the registered `t` slices. Within
that slice the radial normal is the positive spatial `y` normal. The screen
is tangent to the coordinate `t,y` spheres, with the positive symmetric
inverse-square-root frame convention. The frozen outward ambient sphere
orientation gives curl eigenvalues `+3,-3`. Reversing it reverses `h`.
The complete first variation of this defined screen tensor vanishes for
the retained scalar first field; the quadratic map closes on the indicated
TT pair. Their normalized squared sphere norms are `1/6` each and their
cross inner product is zero. Thus the nonzero coefficient gives a nonzero
selected angular projection, without asserting a nonzero value at every
angular point.

## Executed verification and its independence limits

The primary experiment fixes its point and budget before evaluation:
64 bulk-cap boxes, 66 characteristic boxes, eight shell-history boxes,
one boundary endpoint and one apex, totalling 140 old actual-jet calls.
It uses no adaptive refinement and no additional TT iterate. Its normal
and optimized reports are byte identical.

The verification has several distinct parts:

- The [first-field audit](first/FIRST_NORMAL.json) checks the reflected
  comparator, the complete inherited global derivative-error replay,
  physical two-jet reconstruction and physical shell normalization.
- The [exact structural/kernel audit](independent/INDEPENDENT_NORMAL.json)
  checks all primitive TT closures and assembled coefficients, both full
  source rows and source units, the actual shell source and zero TT Cauchy
  data. Thirty exact polynomial kernel cases check images, endpoint signs,
  PDE and boundary traces.
- A separate screen audit verifies four unsimplified ADM/Codazzi
  contributions against an independently entered 14-term expression and
  32 exact rational arbitrary-jet cases.
- The [report arithmetic replay](independent/REPORT_NORMAL.json) imports no
  project code and reconstructs cover geometry and every rational endpoint
  from source-cell intervals through the full tail and final screen sum.
  Its source-cell intervals are producer premises; it does not repeat the
  source campaign.
- The [independent numerical evaluator](numeric/NUMERIC_NORMAL.json)
  owns integer interval arithmetic, literal evaluation of all 78 terms,
  shell and screen formulas, clipped cone-cover geometry and feedback
  aggregation. It makes 119 old actual-jet calls, independently of the new
  TT producer. Its exact endpoints differ by small outward-rounding amounts
  and confirm the same strict negative sign. The old actual-jet solution
  and its derivative-error theorem remain common inputs.

Real corrupted-source and corrupted-implementation controls exercise source
authentication, reflection and lift signs, source units, omitted actual
error channels, scalar feedback, fixed quadratic terms and physical clock
factors. Corrupted reports separately test consistency guards. These are
different verification categories; a deliberately changed implementation
fixture is not an alternate physical solution. Normal and optimized modes
agree within each successful calculation. Independent internal
AI-assisted checks and written analytic reviews are not external specialist
review or proof-assistant verification.

## Relation to exact nonlinear families and remaining work

The earlier [localized exact-data theorem](../HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/)
constructs a family at fixed high finite Sobolev regularity and a protected
short precontact evolution under its stated hypotheses. The later
[controlled-region shell-evolution theorem](../HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md)
also asserts conditional exact evolution for the selected aligned seed on
a prescribed finite interaction interval: \(H^{14}\times H^{13}\) data,
\(\mathcal X_{14}\) evolution, and C2 parameter dependence in the lower
\(\mathcal X_{10}\) topology. Its hypotheses
include the exact canonical family, a timelike shell and uniform state
margins, a fixed support annulus, and the strict actual causal margin
`Z>Lsupport+c_*Tphys`. Those prior statements are preserved.

That theorem explicitly distinguishes its second initial derivative from
the correction chosen in another finite-angular Taylor construction.
Matching a first quadrupole tangent does not match the second coefficient.
This numerical checkpoint does not establish the initial/gauge, second-jet
and causal-domain identification with an applicable exact nonlinear family.
For the complete second response, that identification must compare all
completed second Cauchy jets. For this selected TT coefficient, a narrower
route can use the proved TT autonomy: match the first solution and the
physical/gauge prescription, then prove the exact family's TT displacement
and velocity have the same zero initial projection, or charge their
homogeneous contribution. A difference in a scalar second initial component
alone does not settle this TT question. Nor has the registered slicing
observer been replaced by a geometric observer and clock transported with
such a family.

The [next target](NEXT_STEP.md) is to close those identifications and
calculate a response for a specified transported observation prescription.
A finite-amplitude sign would additionally require a quantitative
parameter remainder relative to its nonzero Taylor coefficient. The present
coefficient alone is not a detector response, physical radiation or energy
measurement, force, stability result, or cosmological-origin result.
