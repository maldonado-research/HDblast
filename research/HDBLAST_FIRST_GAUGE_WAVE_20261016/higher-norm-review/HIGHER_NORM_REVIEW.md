# Independent audit of the continuous C3 raw gauge-amplitude bounds

16 October 2026. **Accepted within the stated coordinate-amplitude scope; no
blocking mathematical defect found.** This is an independent internal AI
calculation and review, not external specialist peer review. Acceptance is
restricted to the exact files pinned in `HIGHER_NORM_REVIEW_RECEIPT.json`.

The final candidate report is `higher-norms/HIGHER_GAUGE_NORMS_NORMAL.json`,
SHA256 `f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d`.
It bounds every mixed partial through total order three of the raw amplitudes
`U=(f,b,d)`, where `xi^y=fY`, `xi^t=bY`, `xi^i=d D_kappa^iY`, on
`0<=x<=3, r>=0, x+r<=8`. Here `x=t/ell`, `r=z/ell`, `ell=10^-6`.
No rho rescaling of the amplitudes is implicit.

## Analytic and algebraic audit

I checked the differentiated equations directly from

\[
U_{xx}-U_{rr}=AU_x+BU_r+DU+s,\qquad U(0,r)=U_x(0,r)=0.
\]

For `v=U_x`, the additional forcing is
`A_x U_x+B_x U_r+D_x U+s_x`. For `w=U_xx`, it is

\[
2A_xv_x+2B_xv_r+2D_xv+A_{xx}U_x+B_{xx}U_r+D_{xx}U+s_{xx}.
\]

All three first coefficient commutators carry the factor two. The candidate
retains them, including the field term `2D_xv`.

Initially `U_rr=U_xr=0`, so `v_x(0)=s(0)`. Differentiating the PDE once gives

\[
w(0,r)=s(0,r),\qquad
w_x(0,r)=A(0,r)s(0,r)+s_x(0,r).
\]

The `A_x` term multiplies the zero initial U_x. It cannot replace the
displayed `A(0,r)` term. The implementation uses the conservative global
absolute majorant of A to bound this initial value, retaining the induced
nonzero time-component velocity. It also retains `s_r(0,r)` as the radial
derivative of w's initial displacement.

For a compatible free mixed wave with data `(h0,h1)`, the half-line
d'Alembert formulas imply

\[
|V|\le\|h_0\|+X\|h_1\|,\qquad
|V_x|,|V_r|\le\|h_{0,r}\|+\|h_1\|.
\]

Thus, with `X=3`, the weighted norm is bounded by
`max(2||h0||/X+2||h1||,||h0_r||+||h1||)`.
The candidate applies this formula correctly to v and w, rather than using
the zero-data source bound for a field with nonzero differentiated data.
In particular `N_free(v)<=2||s0||` and the full w free-data maximum is present.

I checked the componentwise majorant matrix
`M=X(|A|+|B|+(X/2)|D|)`. Its row sum is below one, hence its inverse
`(I-M)^-1` is the nonnegative Neumann resolvent. The same matrix must be used
for U, v and w because they have the same principal and lower-order operator.
The candidate does so. The `X/2=3/2` field-value factor is retained in M,
in both differentiated coefficient matrices, and in the free-data norm.

The direct norms give `U`, `U_x`, `U_r`, `U_xx`, `U_xr`, `U_xxx` and
`U_xxr`. The remaining partials are recovered from the actual PDE:

\[
U_{rr}=U_{xx}-F,\quad U_{xrr}=U_{xxx}-F_x,\quad
U_{rrr}=U_{xxr}-F_r,\qquad F=AU_x+BU_r+DU+s.
\]

Both coefficient and field derivatives in F_x and F_r are included. In
particular F_r contains `B U_rr` and `s_r`; neither is omitted. These
identities account for all ten mixed derivative classes, not only pure time
or pure radial derivatives. The optional improvement
`|U_xx|<=min(n1,(3/2)n2)` is valid because it uses two independently valid
bounds for the same derivative. Conversion to t,z derivatives uses exactly
`ell^(-i-j)`.

The pinned regularity note proves convergence of U, U_x and U_xx in C1,
including the finite offset in the second time-derivative Picard data. It
then recovers the remaining radial derivatives as continuous distributions
with continuous right sides. This provides actual C3 regularity, not merely
uniform bounds on formal differentiated iterates. The initial flat wedge
supplies the needed corner compatibility. The proof retains the physical
reflection traces and does not require a smooth extension of the odd drift
or forcing across the shell. The claimed derivatives at the shell are
physical one-sided derivatives. Those are the appropriate conclusions.

## Independent numerical inequalities

The separate `check_independent_higher_majorants.py` imports no producer code
and computes no matrix inverse. It authenticates the final coefficient and
candidate reports, independently forms M and all commutator/source
majorants, and estimates each resolvent application by sixteen nonnegative
Neumann iterates plus the exact tail

\[
\frac{q^{16}}{1-q}\|\mathrm{rhs}\|_\infty,
\qquad q=\|M\|_\infty<1.
\]

That construction gives slightly larger valid bounds than the producer's
exact inverse. Every one of its thirty mixed component bounds remains
strictly below the corresponding twelve displayed total-order caps:

| Scaled order | f | b | d |
|---|---:|---:|---:|
| 0 | 3.32641e-6 | 2.22340e-12 | 3.33510e-18 |
| 1 | 2.21761e-6 | 1.48227e-12 | 2.22340e-18 |
| 2 | 3.76985e-5 | 2.51980e-11 | 3.77970e-17 |
| 3 | 9.59747e-4 | 6.42163e-10 | 9.50895e-16 |

It separately verifies the candidate's free-data and second-commutator
majorants, the third radial recovery, and all unscaled conversions. Exact
Fraction arithmetic is used throughout. Normal and optimized independent
runs are byte-identical. Three intentional wrong formulas, each exercised
in both modes, fail with the named diagnostic and create no output.

The producer's polynomial checks are useful exact witnesses for the
differentiation formulas. The general product-rule and regularity arguments
above, together with the pinned theorem, are the reason those formulas apply
to the actual solution; polynomial witnesses alone would not establish that.

## Scope and retained inputs

The global coefficient report uses the actual conjugated master with
`u=(alpha/rho^2)cosh(t)^(-3/2)w`. It supplies coefficient and source mixed
derivatives through order two, on the full `8ell` collar established by the
later first-exit proof. The original branch/shell and master existence
theorems remain conditional inputs; this audit does not re-prove them or
replace the full-collar certificate with the older short-collar prose.

These are continuous pointwise bounds for raw coordinate amplitudes and
their t,z derivatives. They do not bound all angular or covariant vector
derivatives, an orthonormal-frame norm, a Sobolev energy, transformed
second-order forcing, a second-order evolution, or a nonlinear remainder.
Within this stated scope the initial data, derivative multiplicities,
boundary regularity, resolvent argument, and numerical inequalities are
consistent.
