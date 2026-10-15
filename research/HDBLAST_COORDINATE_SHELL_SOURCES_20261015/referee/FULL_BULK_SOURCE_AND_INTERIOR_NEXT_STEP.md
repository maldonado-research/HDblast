# Full bulk source formula and a controlled next interior step

15 October 2026. Independent formulation for the current coordinate-source
continuation. The retained full finite-angular coefficient theorem remains
conditional. This note gives no numerical second-order solution.

## An interior value contraction using the newly recovered boundary value

The actual first master solves, on inward \(z\ge0\),

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,
\quad w(0,z)=q_\ell(z),\quad w_t(0,z)=e^{-z}q_\ell(z),
\quad b(t):=w(t,0)=y(t/\ell).
\]

Here \(\ell=10^{-6}\), \(T=3\ell\), and the retained actual collar audit
gives \(|\mathcal W|\le M=102000\) through the relevant causal support.
The new boundary ODE reconstructs this actual \(b\), not a free comparator.
Precisely, use the backward-closed compact triangle
\(D=\{0\le t\le3\ell,\ z\ge0,\ z+t\le8\ell\}\).
The pinned `BACKGROUND_DERIVATIVE_BOUNDS.json` gives a strict first-exit
rectangle, its \(\rho\) and \(Q\) boxes, and
\(\mathcal W=\rho^2Q-35\tanh^2t/4\). The checker independently combines
these boxes and \(|\tanh t|\le t\) to verify \(|\mathcal W|<102000\).

Let \(q_o,v_o\) be odd extensions of \(q_\ell,e^{-z}q_\ell\), and define

\[
w_0(t,z)=\tfrac12[q_o(z+t)+q_o(z-t)]
+\tfrac12\int_{z-t}^{z+t}v_o(r)\,dr
+\mathbf1_{t\ge z}b(t-z).
\]

The last term is smooth across \(t=z\) because all initial boundary jets
vanish. The Dirichlet retarded operator is

\[
(K_D f)(t,z)=\frac12\int_0^t\!\int_0^\infty
\left[\mathbf1_{|z-r|\le t-s}-\mathbf1_{z+r\le t-s}\right]
f(s,r)\,dr\,ds.
\]

Its kernel is nonnegative; the second indicator's support is a subset of the
first. Its row mass is at most \(t^2/2\). The exact actual master therefore
satisfies

\[
w=w_0-K_D(\mathcal Ww),\qquad
\|K_D\mathcal W\|_\infty\le r:=MT^2/2
=\frac{459}{10^9}=4.59\,10^{-7}<1.
\]

This is a contraction on \(C(D)\), over the complete observation time interval.
Every effective integration point obeys \(r+s\le z+t\le8\ell\); no potential
bound on an unaudited half-line is used. It supplies a convergent
Neumann series, a direct bound \(\|w\|_\infty\le\|w_0\|_\infty/(1-r)\),
and the useful residual certificate

\[
\boxed{\|w-p\|_\infty\le
\frac{\|w_0-p-K_D(\mathcal Wp)\|_\infty}{1-r}.}
\]

If a represented boundary function \(\widetilde b\) has a certified uniform
error \(e_b\), replace the numerator by the represented residual bound plus
\(e_b\); initial-data and potential approximation defects must also be included.
Because \(b\) came from the actual dynamic-boundary solution, this Dirichlet
reconstruction identifies that same interior solution by uniqueness. Arbitrary
choices inside the boundary interval boxes need not satisfy the original dynamic
law; they may be used only as approximations with explicit error accounting.

The bump obeys \(0\le q_\ell\le1\) and
\(|e^{-z}q_\ell|\le1\); therefore
\(\|w_0\|_{C(D)}\le1+T+B_b\), where \(B_b\ge\|b\|_\infty\).
An independent exact sum of the retained actual-H cells gives the coarse
\(B_b\approx0.0014591371894211345\). Consequently the receipts give explicit
rational bounds

\[
\|w\|_{C(D)}\le\frac{1+T+B_b}{1-r},\qquad
\|w-w_0\|_{C(D)}\le\frac{r(1+T+B_b)}{1-r}.
\]

The latter is an actual approximation error using the exact boundary lift.
It is a value bound only, and its construction does not evaluate the variable
potential integral at every point. A sharper certified boundary representation
can directly improve \(B_b\).

The initial support has \(z\le2\ell\), so the first solution's support has
\(z\le5\ell\) by time \(T\). The retained collar reaches \(8\ell\).
At any support point, the backward cone at time \(s\) reaches at most
\(z+t-s\le8\ell\); the stated coefficient tube suffices. One may solve
on an auxiliary half-line with bounded extensions and use finite propagation
to remove dependence on that extension. Alternatively the triangular domain
above closes every integral directly. This does not extend the physical
background past its controlled collar.

This advances certified interior **values**. It does not certify differentiated
jets merely by differentiating a bound. With \(V=\alpha(z)\cosh^{-3/2}t\,w\),
the first reconstruction is

\[
u=V/\rho^2,\qquad v=-3V_y/(\rho^2p)=3V_z/(\rho^3p),
\qquad \partial_y=-\rho^{-1}\partial_z.
\]

Thus a direct pointwise bulk Hessian evaluation requires \(u\) through two
derivatives and \(v\) through two derivatives, hence up to three derivatives
of \(w\) before any equation-based reduction. A future implementation can
differentiate the characteristic/Duhamel system and propagate smooth compatible
jet bounds; it must retain derivative forcing, background derivatives, and
boundary jets. The present actual-H cells contain no derivative enclosure for
\(H\), so they alone cannot produce all of those jets.

## A source recipe that keeps all tensor components

Let \(h=h_1\), \(\chi=\chi_1\), and let every covariant derivative and
index raising below use \(g_0\). Define

\[
A^A{}_{BC}=\frac12g_0^{AD}(\nabla_Bh_{CD}+\nabla_Ch_{BD}-\nabla_Dh_{BC}),
\qquad B^A{}_{BC}=-2h^A{}_D A^D{}_{BC}.
\]

These are \(D\Gamma[h]\) and \(D^2\Gamma[h,h]\) along the affine metric
path \(g_0+\epsilon h\). In particular \(B\) already includes the factorial
factor two. The complete Ricci Hessian is

\[
\mathcal R^{[2]}_{AB}=\nabla_C B^C{}_{AB}-\nabla_B B^C{}_{AC}
+2(A^C{}_{CD}A^D{}_{AB}-A^C{}_{BD}A^D{}_{AC}).
\]

For the registered trace-reversed five-dimensional Einstein residual,

\[
\boxed{S^E_{AB}=\mathcal R^{[2]}_{AB}-2\chi_A\chi_B
-\frac23[U''\chi^2g_{0,AB}+2U'\chi h_{AB}].}
\]

No \(Uh_2\) term belongs to this Hessian; that term belongs to the linear
operator acting on the unknown second metric coefficient. Define
\(J^{AB}=h^{AC}h_C{}^B\), \(p_C=\nabla_C\phi_0\). The full scalar Hessian is

\[
\boxed{\begin{aligned}
S^\phi={}&2J^{AB}\nabla_A\nabla_B\phi_0
-2h^{AB}\nabla_A\nabla_B\chi
+2h^{AB}A^C{}_{AB}p_C\\
&-g_0^{AB}B^C{}_{AB}p_C
-2g_0^{AB}A^C{}_{AB}\chi_C-U'''\chi^2.
\end{aligned}}
\]

The physical second equations are \(D\mathscr E[k,w_2]=-S^E\) and
\(D\mathscr S[k,w_2]=-S^\phi\). These formulas are coordinate tensor
identities, valid off shell as well as on shell. Source compatibility with
the linear constraints additionally requires the background and first fields
to solve their equations. Their tensor nature allows all degree-zero/two/four
scalar and degree-two tensor sectors to be projected after evaluation, without
assuming the unchanged scalar-master reconstruction at second order.

## Relative harmonic reduction and lower-index caution

For \(C^A(g)=g^{BC}(\Gamma[g]-\Gamma[g_0])^A{}_{BC}\),

\[
C_1^A=g_0^{BC}A^A{}_{BC},\qquad
Q_C^A=g_0^{BC}B^A{}_{BC}-2h^{BC}A^A{}_{BC},
\]

and, because the metric used to lower the index also varies,

\[
Q_{C,A}=g_{0,AB}Q_C^B+2h_{AB}C_1^B.
\]

The cross term vanishes only after the first gauge solve makes \(C_1=0\).
Using the transformed \(h_1^w,\chi_1^w\), the second reduced metric source is
\(-S^E+\nabla_{(A}Q_{C,B)}\), with the symmetrization factor \(1/2\).
The scalar source remains \(-S^\phi\). This is exactly the representative
used by the retained full contact boundary sources. Alternatively the covariance
conversion in the retained theorem can be used, with all
\(2\mathcal L_\xi h_1+\mathcal L_\xi^2g_0\) and scalar terms retained.

These equations provide a usable full-component implementation contract. The
current boundary source signs are not preserved by a first-gauge change in
isolation, and a material boundary source cannot be combined with an unrelated
longitudinal or harmonic bulk source.

## Independent checks and their limits

`check_full_source_formulation.py` imports only the Python standard library.
It builds generic rational five-dimensional Lorentzian metric jets, including
off-diagonal first perturbations and nonzero background first/second metric
derivatives. Direct coordinate Christoffel/Ricci/scalar path expansions are
compared with the covariant Hessian formulas above for all fifteen independent
Einstein components and the scalar; connection and lowered harmonic-vector
identities are checked separately. Two deterministic nonspecial jet sets are
used. These are local exact checks supplementing the analytic tensor derivation,
not an all-point computer proof or a numerical Einstein solution.

The checker also verifies source elimination using the first shell laws,
degree-four and degree-two Robin join algebra, the ambient velocity conversion,
the actual symmetrizer sign bounds, and the interior contraction constant.
Explicit corruptions remove the Ricci factorial, potential cross term, lower
wave-map cross term, or flip the Robin determinant sign. Ordinary and optimized
Python must reject each with an explicit exception.
