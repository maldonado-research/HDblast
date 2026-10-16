# Independent common-gauge quadratic source derivation and bound audit

16 October 2026. **Accepted within the stated common-gauge bulk-source scope;
no blocking mathematical defect found.** This is an internal independent
mathematical review, not external specialist peer review. The formulas below
refer to the actual common first relative-harmonic pair, not the old
longitudinal spacetime representative. The source identities, frame bound
construction and actual continuous numerical caps have been independently
checked. Exact accepted file hashes are in `BULK_REVIEW_RECEIPT.json`.

## 1. Affine Hessians and the registered reduction

Let `q=h_long+L_xi g0`, `eta=chi_long+L_xi phi0`, using the same constructed
first vector xi in both fields. Along the affine path
`g(e)=g0+e q`, `phi(e)=phi0+e eta`, all Q below are **second parameter
derivatives**, twice the coefficient of e squared. Fix the reference
connection to that of g0. The natural and reduced residuals are

\[
E_{ab}=R_{ab}-\phi_a\phi_b-\tfrac23 V(\phi)g_{ab},\quad
S=\Box_g\phi-V'(\phi),\quad
E^{\rm red}_{ab}=E_{ab}-\nabla^g_{(a}C_{b)}.
\]

The parentheses include one half, and `C_b=g_ba C^a`.
Let every raised index and unmarked covariant derivative below use g0:

\[
P^a{}_{bc}=\tfrac12g_0^{ad}
(\nabla_bq_{cd}+\nabla_cq_{bd}-\nabla_dq_{bc}),\quad
J^a{}_{bc}=-2q^a{}_dP^d{}_{bc},\quad
c_1^a=g_0^{bc}P^a{}_{bc}.
\]

The inverse-metric expansion is
`g(e)^-1=g0^-1-e q^sharp+e^2 q^a_c q^cb+O(e^3)`.
The difference of connections satisfies

\[
\Gamma(g(e))-\Gamma(g_0)=eP+\tfrac12e^2J+O(e^3).
\]

The exact connection-difference Ricci identity, differentiated twice without
commuting covariant derivatives, gives

\[
\begin{split}
Q^E_{ab}={}&\nabla_cJ^c{}_{ab}-\nabla_bJ^c{}_{ac}
+2(P^c{}_{cd}P^d{}_{ab}-P^c{}_{bd}P^d{}_{ac})\\
&-2\eta_a\eta_b-\tfrac23[V''\eta^2(g_0)_{ab}+2V'\eta q_{ab}].
\end{split} \tag{1}
\]

In particular the scalar gradient contribution is `-2 eta_a eta_b`, and the
potential product has no `V q2` term on an affine metric path.
Expanding the inverse metric times the scalar Hessian independently gives

\[
\begin{split}
Q^S={}&2q^a{}_d q^{db}\nabla_a\nabla_b\phi_0
-2q^{ab}\nabla_a\nabla_b\eta
+2q^{ab}P^c{}_{ab}(\phi_0)_c-V'''\eta^2\\
&+2q^c{}_d c_1^d(\phi_0)_c-2c_1^c\eta_c.
\end{split} \tag{2}
\]

The vector's affine Hessian is

\[
Q_C^a=-2q^a{}_d c_1^d-2q^{bc}P^a{}_{bc}. \tag{3}
\]

For completeness the full gauge-subtraction Hessian, before imposing a
harmonic first field, is

\[
Q_{\nabla C,ab}
=\nabla_{(a}(Q_C)_{b)}
+2\nabla_{(a}(q_{b)c}c_1^c)-2P^c{}_{ab}(c_1)_c. \tag{4}
\]

This follows by differentiating both the lowering metric and the connection
in `nabla^g C`. Omitting these cross terms is incorrect for a general first
variation. The actual gauge construction gives `c1=0` as a field, hence
also `nabla c1=0`. Consequently (2) loses its last line, (3) becomes
`Q_C^a=-2q^bc P^a_bc`, and the registered common-gauge source is exactly

\[
Q^{E,\rm red}_{ab}=Q^E_{ab}-\nabla_{(a}(Q_C)_{b)}. \tag{5}
\]

The second equations use negatives of these affine Hessians as forcing.
Absolute bounds are unaffected. The simplification applies to the actual
correlated solution before interval evaluation; independent selections from
its enclosing boxes need not themselves be harmonic.

These formulas use q through two covariant derivatives and eta through two.
Thus the certified master and gauge vector C3 regularity suffice for the
direct calculation. No separate evaluation of `E'[2L_xi h+L_xi^2 Phi0]`
is used; that route can require a fourth xi derivative.

## 2. Independent coordinate two-jet verification

`check_common_gauge_sources.py` has two independently organized evaluations.
The direct evaluation expands the inverse metric, its coordinate derivative,
the Christoffel symbols and their coordinate derivatives through e squared.
It then forms coordinate Ricci, the scalar wave operator, the relative
harmonic vector, its lowering, and the full reduced residual. The recipe
evaluation first computes background-covariant derivatives of q from the
metric two-jets and applies (1)-(4).

Two witnesses use generic curved, nonnormal, nondiagonal Lorentzian metric
jets. Coordinate second derivatives of the background metric are nonzero.
The covariant recipe therefore includes the background connection
derivatives when differentiating `nabla q`. A targeted mutation removes
those derivatives and must fail. The check does not pretend that a
connection which vanishes at one point has zero derivative there.

A third witness has flat background and nontrivial time-dependent,
spatially trace-free derivative jets of q. Its value q may be a general
constant symmetric tensor. It has `c1=nabla c1=0` with nonzero P and nabla P,
so it tests the actual harmonic simplification in (5), not merely a zero
perturbation. Scalar value, gradient and Hessian jets remain independent.

These are exact rational finite-jet witnesses. They supplement the universal
connection-difference and Taylor arguments above; comparison counts are not
a claim of a physical trajectory, full all-point enclosure, or novelty.

## 3. Frame derivatives and the derivative-slot issue

The producer's baseline connection D is a valid way to avoid differentiating
a sphere-normal frame as if its connection vanished identically. It is
metric compatible, preserves the radial/time/angular splitting, retains the
intrinsic sphere connection, and is generally torsionful. Its difference
`K=nabla-D` contains the warped connection coefficients. Put
`s=(rho cosh t)^-1` and `beta=tanh(t)/rho`.

For an angular tensor T and scalar amplitude F(y,t), D differentiates every
angular tensor and derivative index intrinsically. In particular its ordered
second derivatives include

\[
D_yD_0(FT)=\rho^{-1}(F_{yt}-aF_t)T,\quad
D_0D_y(FT)=\rho^{-1}F_{ty}T,
\]
\[
D_yD_a(FT)=s(F_y-aF)D_aT,\quad
D_aD_y(FT)=sF_yD_aT,
\]
\[
D_0D_a(FT)=s\rho^{-1}(F_t-HF)D_aT,\quad
D_aD_0(FT)=s\rho^{-1}F_tD_aT,
\]

while the angular-angular term is `s^2 F D_aD_b T`.
These orders cannot be silently commuted: the baseline D is torsionful.
The producer's proof retains the displayed derivatives of the frame scales.
The intrinsic quadrupole identities provide the covariant sphere jet caps
`1,2,4,8,16`, including the derivatives of DY and Hess Y as tensors.

Let `K0=max_AB sum_C |K^C_AB|` and
`K1=max_ABC sum_E |(D_A K)^E_BC|`. Since the splitting tensors are D-parallel,
the angular derivatives of K vanish intrinsically; only the scalar
coefficient derivatives remain. The valid bounds are

\[
K_0\le |a|+|\beta|,\qquad
K_1\le\max\{|a_y|+|a\beta|,\ \operatorname{sech}^2t/\rho^2\}.
\]

If H0, D1, D2 bound the metric value and its baseline jets, then

\[
H_1\le D_1+2K_0H_0,\qquad
H_2\le D_2+5K_0D_1+(2K_1+6K_0^2)H_0. \tag{6}
\]

To see the coefficient five, first apply D to `nabla q=Dq-K*q-K*q`.
This gives `D2+2K1 H0+2K0 D1`. Then convert the rank-three covariant
tensor `nabla q` from D to nabla, subtracting K on its derivative slot as
well as its two metric slots. This adds `3K0(D1+2K0H0)`. The derivative
slot is essential. For the scalar, `E1=D1eta` and
`E2<=D2eta+K0 E1` include the Hessian's derivative-slot correction.

This argument handles intrinsic curvature through covariant sphere
differentiation and extrinsic warped connection terms through K and DK.
It never holds a sphere-normal connection equal to zero under differentiation.

## 4. Universal contractions

In the adapted Lorentz-orthonormal frame, raising an individual index changes
only a sign. For n=5 and component maxima as above,

\[
P_0\le\tfrac32H_1,\quad P_1\le\tfrac32H_2,\quad
J_1\le2n(H_1P_0+H_0P_1).
\]

Equation (1) has two single derivative contractions and two double product
contractions, the latter multiplied by the Taylor factor two. Hence

\[
|Q^E|\le2nJ_1+4n^2P_0^2+2E_1^2
+\tfrac23(|V''|E_0^2+2|V'|E_0H_0). \tag{7}
\]

Equations (3),(5) on the actual harmonic first field give

\[
|Q_C|\le2n^2H_0P_0,\quad
|\nabla Q_C|\le2n^2(H_1P_0+H_0P_1),\quad
|Q^{E,\rm red}|\le|Q^E|+|\nabla Q_C|. \tag{8}
\]

The coefficient on the last term is one because symmetrization is the
average of two components. Lowering uses the parallel background metric.
For the background scalar, the only gradient component is p, and its
Hessian has diagonal entries `p_y,-ap,ap,ap,ap`, with
`p_y=V'-4ap`. Put `Phi2=max(|p_y|,|ap|)`. The first contraction in (2)
then needs n squared summands, rather than n cubed. Thus

\[
|Q^S|\le2n^2(H_0^2\Phi_2+H_0E_2+H_0P_0|p|)
+|V'''|E_0^2. \tag{9}
\]

I independently re-counted these sums and obtained the producer's factors.
They cover all fifteen symmetric Einstein components and the complete
scalar source over S3, rather than a selected scalar angular projection.
The numerical audit instantiates H0,H1,H2,E0,E1,E2 from the actual first
fields and their certified continuous jets; those caps cannot be inferred
by differentiating a value-only comparator error.

## 5. Completed independent numerical audit and evidence

The accepted final producer report is `COMMON_GAUGE_BULK_NORMAL.json`,
SHA256 `4175c04c09bccf0a9c69814a129ee0a46d1b17403986933a77251b359cfb8ae4`.
I reviewed the final implementation and reproduced all six readable field
caps and all five source caps with `check_independent_bulk_bounds.py`.
This independent evaluator stores **actual** mixed derivatives and inserts
explicit Leibniz binomial coefficients, whereas the producer stores
normalized Taylor coefficients. It uses a distinct finer outward interval
grid of denominator 10^50 and imports no producer module. Both calculate
the actual transformed scalar including `+pf` and the actual master
normalization `cosh(t)^(-3/2)`.

The independent construction authenticates the frozen actual background,
master and gauge reports, reconstructs B,k and the potential jets globally,
and intersects only already-certified actual gauge C0/C1 bounds. It then
evaluates all six metric frame sectors, the scalar, the ordered baseline
derivatives, the full connection corrections, and (7)-(9). It never
differentiates the comparator's value error.

| Continuous component cap | Value | First covariant derivative | Second covariant derivative |
|---|---:|---:|---:|
| Common first metric | 0.01771 | 4497 | 8.792e8 |
| Common first scalar | 176500 | 2.173e10 | 4.355e15 |

| Quadratic affine Hessian | Strict component cap |
|---|---:|
| Natural Einstein Q^E | 9.444e20 |
| Registered reduced Einstein Q^(E,red) | 9.444e20 |
| Scalar Q^S | 3.855e15 |
| Relative harmonic vector Q_C | 5971 |
| Covariant derivative of Q_C | 2.684e9 |

The Einstein entries cover all fifteen independent symmetric components.
The scalar cap and five gauge-vector components are not angular projections.
The equal rounded natural/reduced caps are valid: their exact rational
bounds differ by the derivative-Q_C bound, which is much smaller than this
rounding interval. The additional normalized S3 Euclidean Frobenius bound
of five times the 5x5 component cap is valid, since there are twenty-five
squared matrix entries. It is not a Lorentzian contraction.

The source checker passes 2893 exact finite comparisons on its curved and
harmonic witnesses. The separate frame checker passes 81 sphere fourth-jet
and 625 rank-three connection comparisons. In the sphere chart
`n=(x1,x2,x3,sqrt(1-|x|^2))`, at its origin the coordinate fourth derivative
of `Y=x1^2-x2^2` is zero and the connection vanishes, but its connection
derivative produces `D_1^4Y=-8`. This explicitly detects the invalid step
of holding a normal-frame connection zero under differentiation.

All three independent checkers produce byte-identical normal/optimized
reports. Eighteen source-formula mutations and ten frame/bound mutations,
counting both interpreter modes, produce their named failures with no
output. The finite witnesses support the tensor proofs above; the continuous
enclosure follows from the authenticated actual jets and the interval
evaluation, not from the number of finite comparisons.

This acceptance retains the original branch/master assumptions and the
published first-gauge construction. It establishes continuous common-gauge
bulk source bounds, with no derivative of Q^E, Q^S or Q^(E,red), no completed
boundary/constraint source, no second evolution, no nonlinear remainder,
and no observational conclusion.
