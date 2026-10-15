# Continuous coordinate shell sources from the certified quadrupole

15 October 2026. Conditional mathematics for the registered five-dimensional
Einstein–scalar shell. Exact rational enclosures and internal AI-assisted
independent derivations are supplied; external peer review and novelty are not
claimed.

## What this checkpoint adds

The previous checkpoint bounded the explicit material term in the scalar
junction. This checkpoint reconstructs the actual first metric boundary trace
and its first derivatives, and evaluates the previously symbolic **complete
coordinate source coefficients** in the scalar and isotropic Israel rows.
It supplies 4353 continuous time cells, harmonic projections, proper-shell
source norms and source-sign certificates. These are inputs for the second-order
problem, rather than an evaluated second-order spacetime response.

An instructive distinction follows. The metric correction is less than
`9e-13` of the material scalar source in the specified L2 norm, yet determines
the sign of the complete coordinate source at a scalar zero, where the material
square vanishes. Smallness in an integrated norm does not give a relative
pointwise approximation through a zero. This concerns these coordinate source
coefficients, not a gauge-invariant force or the sign of the response.

## Registered state and conventions

The background has delta=1/1000 and the retained rational tension parameter;
all background boxes are copied exactly from the public analytic certificate.
The actual aligned first quadrupole has width ell=1e-6, angular polynomial
`Y=n1^2-n2^2` on the unit S3, and observation variable

\[
 x=t/\ell=\tau/(\rho_b\ell),\qquad 0\le x\le3.
\]

The initial master fields are `w(0,z)=q_ell(z)` and
`w_t(0,z)=exp(-z)q_ell(z)`. We retain their actual coupled solution and its
certified H trace, rather than replace it with the free comparator H0.
The physical-shell reconstruction note specifies the ambient clock conversion:
at initial ambient time, `partial_t=r partial_T` and `exp(-z)=r/b`.
Thus these data obey the proportional ambient profile relation `F_T=F/b`.
This resolves a coordinate-profile distinction; it does not supply a numerical
second-order initial constraint correction.

Increasing y is outward on each copy. The first representative is longitudinal,
with zero first shell displacement and fixed intrinsic identification:

\[
 h_{1,yy}=-4uY,\quad h_{1,ya}=0,\quad
 h_{1,ab}=2uYh_{0,ab},\qquad \chi_1=vY.
\]

Taylor coefficients are derivatives:
`U_epsilon=U0+epsilon U1+epsilon^2 U2/2`. Every Sp and SK below is a
**one-copy coordinate** coefficient. The prior material outward-normal-sum
source `Qphi=-4v^2Y^2` uses a different row and convention. Its identical-copy
material source is `-2v^2Y^2`; neither expression equals the full coordinate
source without the normal and pullback terms.

All 18 copied inputs are authenticated by [INPUT_PINS.json](inputs/INPUT_PINS.json)
against immutable public Git blobs. The inertial, physical reconstruction and
full-contact manifest memberships are replayable. No prior scientific program
is imported or executed by the new calculation.

## Reconstructing the actual metric trace

The retained master reconstruction gives

\[
 u_b(t)=\rho_b^{-2}\cosh(t)^{-3/2}w_b(t),\qquad v_b(t)=S(t).
\]

The actual trace and dynamic shell law are

\[
 H(x)=\ell(w_z-\gamma w_b)(\ell x,0),\qquad
 -w_z+\mu w_{tt}+\beta(t)w_b=0,
\]
\[
 \beta(t)=\gamma+\mu(35\operatorname{sech}^2t/4-25/4).
\]

For `y(x)=w_b(ell x)`, the gamma terms cancel, leaving the exact identity

\[
 y''+A(x)y=\frac\ell\mu H(x),\qquad
 A(x)=\ell^2\left(\frac52-\frac{35}4\tanh^2(\ell x)\right),
 \qquad y(0)=y'(0)=0.                                      \tag{1}
\]

The scaled force is ell H/mu. Compact support away from the initial shell
gives the zero initial jets. This equation reconstructs the boundary value of
the actual first solution from its already certified trace; it does not
solve the complete second-order bulk system.

Define

\[
 F(x)=\frac\ell\mu\int_0^x(x-s)H(s)\,ds,\qquad
 (K_A y)(x)=\int_0^x(x-s)A(s)y(s)\,ds.
\]

On the whole window, `|A|<=M=5ell^2/2` and

\[
 q=\|K_A\|_{C^0\to C^0}\text{ upper bound}
   =\frac{M3^2}{2}=\frac9{800000000000}<1.
\]

If `Fsup` bounds the exact convolution, set

\[
 B_y=\frac{F_{\rm sup}}{1-q},\quad
 E=qB_y,\quad E'=3MB_y.
\]

Then `|y-F|<=E` and
`|y'-(ell/mu)integral_0^x H|<=E'`. In particular,

\[
 \boxed{\|u_b\|_\infty<1.119\times10^{-8}.}                \tag{2}
\]

The independently checked coarser positive-kernel bound is
`||u_b||<8.644e-8`; signed convolution improves it. The global continuous
proper-time derivative enclosure lies inside

\[
 -8.519\times10^{-7}<\partial_\tau u_b
                       <2.218\times10^{-4}.
\]

Its computation uses

\[
 \partial_\tau u_b=
 \frac{\rho_b^{-2}\cosh(\ell x)^{-3/2}}{\rho_b\ell}
 \left(y'-\frac32\ell\tanh(\ell x)y\right).                \tag{3}
\]

A bound for y alone would not justify this derivative claim. Radial first
jets needed by these boundary rows follow from the actual first junctions,
without numerical differentiation:

\[
 u_y=-2a_bu_b-p_bv_b/3,\qquad
 v_y=-2p_bu_b-\sigma''_bv_b/2.                              \tag{4}
\]

## Complete coordinate coefficients and Robin reduction

Direct differentiation of the exact normalized normal gives quadratic jets
`D2p=4u v_y+12p u^2` and
`D2K_ab=(4u u_y+20a u^2)Y^2 h0_ab`. Subtracting these from the material
junction and using (4) yields the retained complete coordinate coefficients

\[
 \boxed{S_p=-4p_bu_b^2+2\sigma''_bu_bv_b-2v_b^2,\qquad
 S_K=\frac{\sigma''_b}{6}v_b^2-12a_bu_b^2.}                \tag{5}
\]

These multiply Y^2 in the scalar physical row and the isotropic metric
physical row, respectively. The full metric boundary equation is tensorial;
its time-angular and angular trace-free scalar rows still constrain the
unknown second metric and shell bending, even though their direct quadratic
boundary right-hand sides vanish in this first representative. The TT
boundary rows also remain part of the coupled problem. The prescribed gauge
and constraint boundary rows must be supplied in the representative used to
solve it; (5) is not a list of generalized-harmonic boundary data.

For `Z2=n1^2+n2^2-1/2` and
`H4=Y^2-(n1^2+n2^2)/2+1/12`,

\[
 Y^2=1/6+Z_2/2+H_4.
\]

Thus the full scalar Sp/SK sources project with factors 1/6, 1/2 and 1 into
these precise bases. Their squared angular L2 shares are 10/27, 5/18 and
19/54. A coefficient of a different harmonic normalization must be converted.

In the retained **degree-four** scalar reduction,

\[
 C_4=N_4-T_4+K_4=16u_b^2-4v_b^2/3,
\]
\[
 E_{N,4}=(C_4)_y-4S_K
        =-16a_bu_b^2-\frac{16}3p_bu_bv_b
                         +\frac23\sigma''_bv_b^2.          \tag{6}
\]

The table field EN denotes this degree-four coefficient. It is the first
entry paired with Sp in the retained Robin matrix

\[
 \partial_y\binom{N_4}{w_4}+
 \begin{pmatrix}2a&-4p/3\\-p/2&\sigma''/2\end{pmatrix}
 \binom{N_4}{w_4}=\binom{E_{N,4}}{S_p}.
\]

For the registered boxes, `D=diag(3,8)` gives a positive definite symmetric
matrix DM: its smallest eigenvalue is bounded below by 0.6677, and its
determinant exceeds 10.691. This is a boundary-matrix fact. It does not
identify these two amplitudes as independent physical wave channels or prove
an energy estimate for the entire Einstein–scalar problem.

One must not project C4 or EN as if it were a universal Y^2 scalar. The
degree-two angular trace-free equation gives instead
`C2=-2v^2` and `EN2=24a u^2+8p uv+(5/3)sigma''v^2` in the Z2 basis.
The independent notes derive this algebraic distinction. No numerical EN2
table or degree-zero master reduction is claimed here; the monopole has no
nondegenerate angular trace-free Hessian equation of this kind.

## Signs and norm comparison

The actual scalar has a zero in the certified bracket
`J=(1.49998,1.49999)`. The signed Volterra reconstruction gives `u_b>0`
throughout its closed bracket. At every scalar zero x* in J, therefore,

\[
 S_p(x_*)=-4p_bu_b(x_*)^2>0,\qquad
 S_K(x_*)=-12a_bu_b(x_*)^2<0.                              \tag{7}
\]

Independent exact endpoint certificates show `Sp<0<SK` at both endpoints
of J. Continuity therefore gives **at least two distinct zeros for each
complete coordinate coefficient inside J**, one on each side of any such
scalar zero. Neither a unique scalar zero nor a unique coefficient zero
is assumed or concluded. The second retained scalar-zero bracket supplies
another positive-u certificate, without a claim counting every source zero.

Uniformly on the continuous tail `2<=x<=3`, the primary enclosures give

\[
 S_p<-2.061\times10^{-7},\qquad S_K>6.872\times10^{-8},
 \qquad E_{N,4}>2.748\times10^{-7}.                         \tag{8}
\]

All numbers use registered model units. They do not provide an observational
scale. With background proper shell measure
`dmu0=rho_b^4 cosh^3(t) dt dOmega`, the continuous enclosures imply

\[
 1.444\times10^{21}<\|S_pY^2\|_{L^2(d\mu_0)}^2
                              <2.161\times10^{21},
\]
\[
 1.604\times10^{20}<\|S_KY^2\|_{L^2(d\mu_0)}^2
                              <2.401\times10^{20},
\]
\[
 9.035\times10^{20}<\|E_{N,4}H_4\|_{L^2(d\mu_0)}^2
                              <1.352\times10^{21}.         \tag{9}
\]

The angular factor for the first two squared norms is `3pi^2/20`, whereas
the degree-four Robin norm uses `19pi^2/360`. This distinction is explicit
in the report. These are source-amplitude norms in a stated background
measure, not an Einstein energy or a norm of the solved response.

For `dSp=-4p u^2+2sigma''uv` and `Sp_material=-2v^2`, exact norm intervals
give

\[
 \frac{\|dS_pY^2\|_{L^2}}{\|(-2v^2)Y^2\|_{L^2}}
                                      <9\times10^{-13}.   \tag{10}
\]

The denominator has a certified strictly positive lower bound. No quotient
by v(x*) is used; such a pointwise relative comparison would fail at its
zeros. Small metric corrections and the sign statement (7) are compatible.

## Whole-cell arithmetic and reproduction

Every inherited closed cell is used: one prefront cell, 4096 pulse cells
and 256 tail cells, with time denominator 4096. For each endpoint envelope
of H, prefix integrals reduce its triangular convolution to a quadratic
polynomial. Extrema are taken over both endpoints and any interior vertex.
Positive inverse-mu uncertainty uses all endpoint products, including
negative and sign-changing convolution intervals.

The clock bounds are `1-(ell b)^2<=cosh(ell x)^(-3/2)<=1` and
`0<=tanh(ell x)<=ell b` on a cell ending at b. For the norm quadrature,
`1<=cosh^3(ell x)<=1+2(ell b)^2`; the tiny registered domain lies inside
the explicitly verified polynomial inequality giving this upper bound.
Machin's arctangent identity with alternating rational sums encloses pi.
All final cell intervals are rounded outward to denominator 2^144 before
norm integration. Decimal displays do not decide any inequality.

The CSV fields are y, yprime, u, utau, uy, v, vy, Sp, SK, EN, dSp, dSK
and dEN. Values enclosed on each whole cell may have correlations; interval
products discard those correlations conservatively. The generator does not
claim optimal enclosures or pointwise uniqueness from them.

The source report, independent derivations, all-cell reconstruction checker,
error controls and exact reproduction instructions accompany this note.
The replay roster records the actual named groups and processes after the
reviewed outputs are frozen. Finite arithmetic verifies the stated
conditional bounds, not the analytic assumptions inherited from the prior
Einstein–scalar construction.

## A controlled next interior step

The independent bulk-source note also supplies an interior value certificate.
On the backward-closed triangle
`D={0<=t<=3ell, z>=0, z+t<=8ell}`, the copied background certificate gives
`|W|<102000`. The half-line Dirichlet retarded operator KD has row mass at
most `t^2/2`. With the actual boundary value `b(t)=y(t/ell)`, let w0 be the
free odd-extension initial solution plus the boundary lift
`1_(t>=z)b(t-z)`. Then

\[
 w=w_0-K_D(Ww),\qquad
 \|K_DW\|\le r=\frac{459}{10^9}<1.
\]

Every predecessor remains in D, so this uses no potential bound beyond the
controlled collar. The independently reconstructed coarse boundary bound
gives the strict consequences

\[
 \|w\|_{C(D)}<1.001463,\qquad
 \|w-w_0\|_{C(D)}<4.597\times10^{-7}.
\]

The latter uses the exact boundary lift. Replacing it by a represented
approximation adds its boundary error. A candidate p has the explicit
residual certificate

\[
 \|w-p\|_{C(D)}\le
 \frac{\|w_0-p-K_D(Wp)\|_{C(D)}}{1-r}.
\]

This is a first-order interior **value** construction and bound. The accompanying
note derives the full tensor Einstein/scalar Hessian source formulas and the
relative harmonic conversion, with independent exact generic-jet checks of
all components. It does not evaluate those sources throughout the interior.
The scalar reconstruction differentiates the master field, so interior
derivatives through sufficient order are still needed before a complete
quadratic bulk source can be numerically certified.

## What remains to compute the full response

The retained conditional finite-angular construction is preserved. Its full
numerical application still needs interior first-field derivatives and
background jets sufficient to evaluate the quadratic bulk source, a
particular localized second initial constraint correction, the transformed
physical/gauge/constraint boundary rows in one common representative, and
continuous residual and trace estimates for the actual coupled solution.
The independent source/scope note supplies a concrete interior Dirichlet
reconstruction route and universal Hessian source recipe. An enclosure of
interior values alone would not certify the derivatives in that recipe.

The [focused literature check](literature/LITERATURE_CHECK.md) identifies
recent related Einstein boundary results and possible validated fixed-point
and Volterra approximation methods. None inspected supplies the required
five-dimensional field-dependent scalar-shell theorem or the missing
model-specific numerical certificate.

No full quantitative U2, finite-amplitude exact nonlinear family, effective
epsilon/remainder bound, observational scale, higher-dimensional origin of
the Big Bang, or publication-readiness claim follows. No Zenodo version or
APS submission is part of this checkpoint.
