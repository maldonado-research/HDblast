# Continuous raw gauge-amplitude derivatives through order three

16 October 2026. This separate certificate instantiates the compatible
mixed-wave regularity theorem with the frozen whole-collar coefficient
report. The coefficient directory is unchanged.

The result applies to the actual constructed raw amplitudes
\(U=(f,b,d)\), on

\[
D=\{0\le x\le3,\ r\ge0,\ x+r\le8\},\qquad
x=t/\ell,\quad r=z/\ell,\quad\ell=10^{-6}.
\]

The vector is \(\xi^y=fY\), \(\xi^t=bY\),
\(\xi^i=dD_\kappa^iY\). There is no additional rho rescaling.
All derivatives below are ordinary coordinate partials of these amplitudes.
They are continuous up to the physical shell and the one-sided edges of D.

Every mixed partial of a given total *scaled* order has the following
strict upper bound, by channel:

| Total order | f | b | d |
| --- | ---: | ---: | ---: |
| 0 | 3.32641e-6 | 2.22340e-12 | 3.33510e-18 |
| 1 | 2.21761e-6 | 1.48227e-12 | 2.22340e-18 |
| 2 | 3.76985e-5 | 2.51980e-11 | 3.77970e-17 |
| 3 | 9.59747e-4 | 6.42163e-10 | 9.50895e-16 |

All entries are exact finite-decimal rational caps, rounded upward and
verified strictly against the Fraction calculation. The JSON retains each
of the ten mixed partials in each of the three channels, its exact rational
upper bound, and its conversion to the unscaled coordinates.

Across all channels, bounds for unscaled \((t,z)\) partials of orders
zero through three are respectively

\[
3.32641\,10^{-6},\qquad 2.21761,\qquad
3.76985\,10^7,\qquad9.59747\,10^{14}.
\]

Here “unscaled” means the coordinates t and z, with the factor
\(\ell^{-k}\) for total derivative order k. These are not proper-time,
proper-radial, orthonormal-frame, angular, covariant-vector or energy norms.

## 1. Inputs and analytic regularity basis

The pinned coefficient report has SHA256
`4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53`.
It supplies elementwise absolute bounds for A, B, D and all their mixed
scaled derivatives through order two, together with the actual source
S and all its mixed derivatives through order two. Its source is in the
first channel only. Those bounds come from the authenticated actual
background and master, rather than an arbitrary coefficient family.

The pinned `MIXED_WAVE_REGULARITY.md` has SHA256
`00f654647d3cd5578e3f9afcb00dc3ec81739036f40adbe43f53e84b87506e98`.
It constructs U in C1 and proves that Ux and Uxx converge in C1 by
summable Picard differences. The PDE then recovers the remaining continuous
radial derivatives. The initial flat wedge x+r<1 makes all corner jets
compatible; the differentiated fields obey the same time-independent
Dirichlet row for f and Neumann rows for b,d. No smooth parity extension
of the physical coefficients is assumed.

These analytic facts distinguish the result from formally differentiating
a C0 bound or merely obtaining uniform bounds on nonconvergent derivatives.
The retained actual branch/master assumptions remain conditional inputs.

## 2. The same nonnegative resolvent

The raw system is

\[
U_{xx}-U_{rr}=AU_x+BU_r+DU+s,\qquad s=(S,0,0),
\quad U|_{x=0}=U_x|_{x=0}=0.
\]

Set \(X=3\), \(\theta=X/2=3/2\), and use the componentwise norm

\[
N_i(V)=\max\{\theta^{-1}\|V_i\|_\infty,
                          \|V_{i,x}\|_\infty,\|V_{i,r}\|_\infty\}.
\]

Let A0, B0, D0 denote the nonnegative elementwise global coefficient
majorants, not values at the initial slice. The kernel bound is X, so

\[
M=X(A_0+B_0+\theta D_0),\qquad R=(I-M)^{-1}\ge0.
\tag{1}
\]

The checker reconstructs M, checks equality with the frozen report,
checks its row sum is below one, and verifies the exact rational inverse
product. The factor theta is necessary because field values are bounded
by theta times the weighted C1 norm. It is retained in every commutator.

For shorthand, all lower-case source majorants below are three-component
nonnegative vectors with zero second and third entries. Let
\(s_0^*,s_x^*,s_r^*,s_{xx}^*\) be the global source caps. Define

\[
L_x=A_x^*+B_x^*+\theta D_x^*,\qquad
L_{xx}=A_{xx}^*+B_{xx}^*+\theta D_{xx}^*.
\]

The first norm vector is
\(n_0=R Xs_0^*\). It equals the channel bound recomputed in the
coefficient certificate.

## 3. Differentiated equations and their initial data

Put v=Ux and w=Uxx. Direct differentiation gives

\[
v_{xx}-v_{rr}=Av_x+Bv_r+Dv
+A_xU_x+B_xU_r+D_xU+s_x,
\tag{2}
\]
\[
\begin{split}
w_{xx}-w_{rr}={}&Aw_x+Bw_r+Dw
+2A_xv_x+2B_xv_r+2D_xv\\
&+A_{xx}U_x+B_{xx}U_r+D_{xx}U+s_{xx}.
\end{split}\tag{3}
\]

The initial data are

\[
(v,v_x)|_{x=0}=(0,s|_{x=0}),\qquad
(w,w_x)|_{x=0}=(s|_{x=0},A(0,r)s(0,r)+s_x(0,r)).
\tag{4}
\]

In (4), **A(0,r) is the initial coefficient value, not A_x(0,r)**.
The latter multiplies the zero initial Ux and does not supply this term.
The coefficient A0 in the following inequalities is its conservative
global absolute majorant.

For a compatible free mixed wave with displacement h0 and velocity h1,
the half-line d'Alembert formulas give

\[
\|V\|\le\|h_0\|+X\|h_1\|,\qquad
\max(\|V_x\|,\|V_r\|)\le\|h_{0,r}\|+\|h_1\|.
\]

Thus its weighted norm is bounded componentwise by

\[
N_{\rm free}\le
\max\{2\|h_0\|/X+2\|h_1\|,\ \|h_{0,r}\|+\|h_1\|\}.
\tag{5}
\]

Apply (5) to (4), using the global source caps for the initial traces:

\[
g_v=2s_0^*,\quad h_0^*=s_0^*,\quad h_{0,r}^*=s_r^*,
\quad h_1^*=A_0s_0^*+s_x^*,
\]
\[
(g_w)_i=\max\{2(h_0^*)_i/X+2(h_1^*)_i,
                          (h_{0,r}^*)_i+(h_1^*)_i\}.
\]

The nonzero displacement derivative and the induced tangential initial
velocity are retained. The same resolvent now gives

\[
\boxed{n_1=R\{g_v+X(L_xn_0+s_x^*)\},\qquad
n_2=R\{g_w+X(2L_xn_1+L_{xx}n_0+s_{xx}^*)\}.}
\tag{6}
\]

The analytic theorem implies \(N(U)\le n_0\),
\(N(U_x)\le n_1\), and \(N(U_{xx})\le n_2\).
Every equation defining these nonnegative majorants is checked exactly.

## 4. Every remaining mixed partial

The direct consequences are

\[
|U|\le\theta n_0,\quad|U_x|,|U_r|\le n_0,\quad
|U_{xx}|\le\min(n_1,\theta n_2),\quad|U_{xr}|\le n_1,
\quad|U_{xxx}|,|U_{xxr}|\le n_2.
\]

The remaining three partials follow from the actual PDE:

\[
U_{rr}=U_{xx}-AU_x-BU_r-DU-s,
\]
\[
U_{xrr}=U_{xxx}-A_xU_x-AU_{xx}-B_xU_r-BU_{xr}-D_xU-DU_x-s_x,
\]
\[
U_{rrr}=U_{xxr}-A_rU_x-AU_{xr}-B_rU_r-BU_{rr}-D_rU-DU_r-s_r.
\tag{7}
\]

Taking absolute values componentwise uses the already computed bound for
Urr in the last line. The last source radial derivative cannot be omitted.
All right sides are continuous by the pinned C3 construction, so these
are continuous classical derivative bounds, including the appropriate
one-sided shell traces. This recovery does not remove the physical
Dirichlet/Neumann reflection terms from the actual derivatives.

## 5. Exact checks and portable replay

`check_gauge_higher_norms.py` imports only the Python standard library.
It checks the input hashes, uses exact Fraction matrix arithmetic, verifies
the inverse and majorant equations, and restores the ell inverse powers
for each mixed partial. It also verifies the differentiated equations,
initial-data identities, and radial recoveries as exact bivariate polynomial
identities. Those algebra checks do not replace the analytic regularity
proof; that proof is a pinned input.

After creating a new external output directory, run for example:

```sh
python -I -B check_gauge_higher_norms.py --output /tmp/new-higher-norm-replay/normal.json
python -I -B -O check_gauge_higher_norms.py --output /tmp/new-higher-norm-replay/optimized.json
python -I -B run_higher_norm_controls.py --output /tmp/new-higher-norm-replay/controls.json
```

All output paths must be new and outside this source/input tree. Normal
and optimized bounds are byte-identical. Five wrong formulas are tested
in both modes: replacing A(0) by A_x(0), losing the second-commutator
factor two, dropping the field-value factor 3/2, omitting s_r from the
third radial recovery, and using positive ell powers for unscaled
derivatives. Every control must raise its named exception.

This certificate supplies only the linear first gauge vector's raw
amplitude coordinate derivatives through order three. It supplies no
second-order metric/scalar gauge conversion, complete second-order
contact evolution, nonlinear remainder, energy norm, or observational
conclusion.
