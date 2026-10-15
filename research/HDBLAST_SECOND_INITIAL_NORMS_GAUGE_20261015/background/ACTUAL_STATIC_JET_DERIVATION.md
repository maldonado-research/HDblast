# Actual static coefficient jets on the original packet annulus

Unpublished continuation, 15 October 2026. Exact rational coefficient bounds
on `t=0`, `ell<=z<=2ell`, `ell=10^-6`. No earlier frozen capsule is changed.
The mathematical input is the retained actual smooth regular background,
not a sampled or numerically integrated surrogate.

`ACTUAL_ANNULUS_STATIC_JETS.json` supplies all normalized derivatives through
order four in both z and y. Its `z_jets[name][j].exact` is the interval for
`D_z^j(name)/j!`; `y_jets` uses `D_y^j(name)/j!`. The names are
`rho,eta,p,a,rho_inv,h,alpha,C,B,p_inv,k,U_0,U_1,U_2,U_3,U_4,exp_minus_z`.
The `strict_abs_*_derivative_caps` are integer bounds for **unnormalized**
derivatives; the exact normalized interval entries are usually much sharper.

The packet reconstruction is

\[
u=Cq,\quad u_t=Ce^{-z}q,\qquad
v=B(q_z+kq),\quad v_t=Be^{-z}[q_z+(k-1)q].
\]

These data suffice for the explicit physical second initial constraint jet
and its first radial derivatives, after combining them with certified bump
derivatives. In particular the conservative derivative budget is

| Target | bump q order | C order | B order | k order |
|---|---:|---:|---:|---:|
| gamma2 values | 2 | 2 | 1 | 1 |
| gamma2_y, phi2 values | 3 | 3 | 2 | 2 |
| phi2_y | 4 | 4 | 3 | 3 |

K2 values use orders at most one in q, while Pi2 values use at most two.
The report provides order four for every coefficient for convenience.

## Authenticated actual input and whole-annulus value refinement

The three byte-identical input copies, their original paths and hashes are
recorded in `INPUT_PROVENANCE.json`; the checker hard-pins all three hashes.
The actual coefficient report is
`fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0`.
The copied analytic background note specifies the same exact model,
`eta=phi-1`, `delta=1/1000`, `c=96714321425692/161839258930731`, and

\[
U(1+\eta)=-2/27+(14/9)\eta^2+(50/27)\eta^3
             -\eta^4/6-(4/9)\eta^5-(2/27)\eta^6.
\]

The new checker requires exact equality of these polynomial coefficients
with the authenticated report. Its retained collar covers `[0,8ell]`, with

\[
\rho\in[129,131],\quad \eta\in[-9/10^5,-8/10^5],\quad
p\in[-14/10^5,-12/10^5],\quad a\in[111/1000,112/1000].
\]

The inward background equations are

\[
\rho_z=-\rho^2a,\quad \eta_z=-\rho p,\quad
p_z=-\rho(U_1-4ap),\quad a_z=\rho^{-1}+\rho p^2/3.
\]

The new checker re-evaluates their interval images on the retained rectangle
and requires exact endpoint equality with the authenticated audit. Since the
actual solution is already proved to remain in that rectangle, integration
from its authenticated shell box yields, for every `z in [ell,2ell]`,

\[
X(z)\in I_{X,b}+[\ell,2\ell]F_X(I_{\rm outer}).
\]

The initial interval for each annulus jet is the intersection of this image
with its retained outer box. This is a full interval enclosure, not endpoint
evaluation or a sampled trajectory. Every intersection and its ingredients
are recorded. In particular rho stays positive and p stays strictly negative.

Write `h=U_1/p`, `k=rho(5a/2-h)`. On the whole retained collar the checker
recomputes the original strictly negative k interval and requires exact
equality with the pinned report. If `-K<=k<=-m<0`, the identity
`alpha_z=k alpha`, `alpha(0)=1` gives throughout the original packet annulus

\[
1-2\ell K\le\alpha\le\frac1{1+\ell m}<1.
\]

This follows from `exp(x)>=1+x` and `exp(-x)<=1/(1+x)` for positive x.
All denominators are checked positive. No fractional power or numerical
transcendental value is evaluated. Likewise
`exp(-z) in [1-2ell,1/(1+ell)]`, and its normalized derivative j is exactly
`(-1)^j exp(-z)/j!`.

## Normalized jet recurrence and coefficient cancellations

At an arbitrary point `z0` in the annulus let
`X_j=D_z^j X(z0)/j!`. Multiplication is the Cauchy product
`(XY)_n=sum_(j=0)^n X_j Y_(n-j)`. Adjoin

\[
b=1/\rho,\quad h=U_1/p,\quad
\alpha=(\rho/\rho_b)^{3/2}p/p_b,\quad C=\alpha/\rho^2,
\quad B=3\alpha/(\rho^3p).
\]

Their exact additional ODEs are

\[
b_z=a,\qquad h_z=-\rho(U_2-h^2+4ah),
\]
\[
\alpha_z=k\alpha,\quad C_z=L C,\quad B_z=m_B B,
\qquad L=\rho(9a/2-h),\quad m_B=3\rho a/2.
\]

The B logarithmic derivative includes exact cancellation of the p factor;
its sign is positive even though B itself is negative. The initial values
of b,h,C,B are obtained from the actual annulus boxes by rational division.
Each next normalized coefficient is the corresponding RHS coefficient divided
by `n+1`. The same recurrence evaluates the exact polynomial potential at eta.

Induction on n proves each interval coefficient contains the actual derivative
at every possible `z0`. This argument requires only smoothness and the exact
ODE. It makes no claim about convergence of a Taylor series, and evaluates
neither an ODE step nor an approximation residual. The reciprocal p jet uses
`(1/p)_0=1/p_0` and
`(1/p)_n=-(1/p)_0 sum_(j=1)^n p_j(1/p)_(n-j)`.

A separate Laurent polynomial total-derivative algebra over the five formal
variables `rho,eta,p,a,alpha` checks every recurrence against repeated exact
chain-rule differentiation of its defining function. Sixteen coefficient
families through four positive derivative orders give **64 exact zero
polynomial defects**. The elementary exponential family is treated directly.
These are all-point algebraic identities, not tests at generic points.

All coefficients in both algebras use Fraction arithmetic; floating values
are explicitly forbidden in Taylor coefficients. Floats occur only in labeled
display fields after the exact calculations and comparisons.

## Radial derivative convention and remaining work

The repeated operator is `D_y=-b D_z`. Applying it to the finite z jets gives
the `y_jets` entries, including the product-rule terms. In particular

\[
D_y^2 f=b^2D_z^2f+abD_zf.
\]

This identity is independently checked by the formal algebra. These are
derivatives of scalar radial amplitudes. In a frame parallel along y they
are also the radial covariant derivatives of those amplitude components.
For coordinate angular metric components, include the metric factor:
`D_coordinate_y(rho^2 f_l kappa_ij)=rho^2(f_l,y+2af_l)kappa_ij`.
The report does not identify ordinary coordinate derivatives with covariant
ones.

The bump derivative enclosure, angular factors and complete quadratic source
algebra must still be combined with these coefficient intervals to certify
second initial-field norms. This report alone supplies no such norms, gauge
completion, second evolution, exact nonlinear initial family or remainder.

## Reproduction and controls

```sh
python -I -S -B bound_static_coefficient_jets.py --output /tmp/actual-annulus-static-jets.json
python -I -S -B run_static_jet_controls.py --output /tmp/static-jet-controls.json
```

Both programs accept an external `--output` and leave their source unchanged.
The controls runner uses system temporary storage, requires normal and
optimized complete reports to match byte for byte, and requires both modes
to reject a modified input byte, missing Taylor normalization, reversed rho
ODE sign, halved B logarithmic derivative, and a p box crossing zero.
Each failure must have the intended reason and leave no output report.
