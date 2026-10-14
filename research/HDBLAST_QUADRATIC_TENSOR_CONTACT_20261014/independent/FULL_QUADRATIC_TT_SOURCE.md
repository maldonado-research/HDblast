# Full quadratic bulk TT source for the aligned quadrupole

14 October 2026. Independent continuation calculation confined to this
scratch directory. Published input capsules remain unchanged.

The complete geometric and scalar source is nonzero. With the conventions
specified below, its angular tensor coefficient is

\[
\boxed{\big(D^2\mathscr E_0[(h_1,\chi_1),(h_1,\chi_1)]\big)_{ij}^{\rm TT}
 =\left(16u^2+\frac83v^2\right)\mathcal T_{ij},\qquad
 (-D^2\mathscr E_0)_{ij}^{\rm TT}
 =-\left(16u^2+\frac83v^2\right)\mathcal T_{ij}.}
\tag{1}
\]

Here the longitudinal metric amplitude is \(\Psi=u(y,t)Y(n)\), the scalar
amplitude is \(\chi_1=v(y,t)Y(n)\), and
\(Y=n_1^2-n_2^2\). The angular tensor \(\mathcal T\) has exactly the Weyl
normalization used in the published precontact calculation. Geometry
contributes \(16u^2\); scalar stress contributes \(8v^2/3\). They reinforce
each other in this gauge. No first-order master equation, background
Einstein equation, or radial equation has been used to obtain this local
quadratic identity; those equations ensure that the specified first-order
fields are actual solutions.

## 1. Input, first-order gauge, and Taylor convention

The actual spacetime ansatz is taken from the published
`HDBLAST_REGULAR_SHELL_CONTACT_20261013` manuscript, equations (3)--(6),
and its `SMOOTH_AMBIENT_RECONSTRUCTION_AND_LOCAL_FLUX.md`, sections 2--3.
It is not the initial ADM constraint solver. Its exterior-cone form is

\[
g_0=dy^2+\rho(y)^2\gamma,\quad
\gamma=-dt^2+\cosh^2t\,q_{S^3},\quad \phi_0=\phi_0(y),\quad
p=\phi_{0,y},\quad a=\rho_y/\rho,
\]
\[
(h_1)_{yy}=-4\Psi,\quad(h_1)_{y\mu}=0,\quad
(h_1)_{\mu\nu}=2\rho^2\Psi\gamma_{\mu\nu},\qquad
\chi_1=-\frac{3P_y}{\rho^2p},\quad\Psi=\frac P{\rho^2}.
\tag{2}
\]

Use the registered bulk equations

\[
\mathscr E_{ab}=R_{ab}-\phi_a\phi_b-\frac23U(\phi)g_{ab},\qquad
\mathscr S=\Box_g\phi-U_{,\phi},
\tag{3}
\]

and the second-variation convention

\[
g(\varepsilon)=g_0+\varepsilon h_1+\frac{\varepsilon^2}{2}h_2,
\quad
\phi(\varepsilon)=\phi_0+\varepsilon\chi_1
                   +\frac{\varepsilon^2}{2}\chi_2.
\tag{4}
\]

Then the second-order Einstein equation is
\(D\mathscr E_0[h_2,\chi_2]=-D^2\mathscr E_0[(h_1,\chi_1)]^2\).
The Hessian in (1) is computed along the *linear* path with \(h_2=\chi_2=0\).
The radial coordinate \(y\) is spacelike; its foliation lapse is not the
temporal Cauchy lapse. The operator TT projects only the angular \(ij\)
block on unit \(S^3\), at fixed \((y,t)\). Spatial sphere connection terms
in the four-dimensional Hessian differ from the sphere Hessian by a pure
angular trace; multiplying the angular metric by \(\cosh^2t\) does not
change its TT subspace.

For \(A=B=\operatorname{diag}(1,-1,0,0)\), linear evolution preserves the
one angular harmonic \(Y\), so (2) has the amplitudes in (1). These
amplitudes can include both displacement and velocity data.

## 2. Exact geometric anisotropy from the radial foliation

Consider the general zero-shift metric

\[
g=n^2dy^2+k_{\mu\nu}dx^\mu dx^\nu,\qquad
k_{\mu\nu}=\rho^2e^{2\omega}\gamma_{\mu\nu}.
\]

Its second fundamental form is pointwise proportional to \(k\):
\(K_{\mu\nu}=(a+\omega_y)k_{\mu\nu}/n\).
All radial derivative and extrinsic curvature terms in \(R_{\mu\nu}[g]\)
therefore have a factor \(k_{\mu\nu}\). Although their coefficients depend
on every spacetime coordinate, their angular blocks are pure trace.
The remaining anisotropic part is exactly

\[
R_{\mu\nu}[k]-\frac1n\nabla^{(k)}_\mu\nabla^{(k)}_\nu n.
\tag{5}
\]

In four slice dimensions, after discarding multiples of \(\gamma_{\mu\nu}\),
this becomes

\[
-2D_\mu D_\nu\omega+2D_\mu\omega D_\nu\omega
-\frac1nD_\mu D_\nu n
+\frac1n(D_\mu\omega D_\nu n+D_\nu\omega D_\mu n).
\tag{6}
\]

The linear Taylor path (4) has

\[
n=\sqrt{1-4\varepsilon\Psi}
  =1-2\varepsilon\Psi-2\varepsilon^2\Psi^2+O(\varepsilon^3),
\quad
\omega=\tfrac12\log(1+2\varepsilon\Psi)
       =\varepsilon\Psi-\varepsilon^2\Psi^2+O(\varepsilon^3).
\]

Substitution in (6) cancels its linear Hessians and yields the quadratic
coefficient

\[
\varepsilon^2[12\Psi D_\mu D_\nu\Psi
                  +6D_\mu\Psi D_\nu\Psi].
\tag{7}
\]

Thus \(D^2R\), as opposed to the \(\varepsilon^2\) coefficient, is
\(24\Psi D^2\Psi+12d\Psi\otimes d\Psi\), up to angular trace.
For every smooth scalar \(f\) on closed \(S^3\),
\(\Pi_{\rm TT}(D_iD_jf)=0\): integrating twice against a transverse
tensor makes the pairing zero. Since

\[
D_iD_j(\Psi^2)=2\Psi D_iD_j\Psi+2D_i\Psi D_j\Psi,
\]

equation (7) gives

\[
\boxed{(D^2R)_{ij}^{\rm TT}
             =-12\Pi_{\rm TT}(D_i\Psi D_j\Psi).}
\tag{8}
\]

There are no omitted radial or temporal derivative contributions to this
angular TT block: in this ansatz they occur in the terms already proved
to be angular trace. Equations (5)--(8) compute spacetime Ricci, not an
initial constraint residual.

## 3. Exponential-path cross-check and its precise limitation

The nonlinear scalar completion

\[
g_{\exp}(\varepsilon)
 =e^{-4\varepsilon\Psi}dy^2+\rho^2e^{2\varepsilon\Psi}\gamma
\tag{9}
\]

has the same \(h_1\), but nonzero second derivative

\[
b_{yy}=16\Psi^2,\quad b_{\mu\nu}=4\rho^2\Psi^2\gamma_{\mu\nu}.
\]

For (9), \(n=e^{-2\omega}\), \(\omega=\varepsilon\Psi\), and (6) is
exactly \(-6\varepsilon^2d\Psi\otimes d\Psi\). This verifies (8)
provided the completion is handled: its linear Einstein anisotropy is

\[
(D\mathscr E_0[b,0])_{\mu\nu}^{\rm anisotropic}
               =-12D_\mu D_\nu(\Psi^2),
\tag{10}
\]

which has zero angular TT projection. The chain rule gives
\(\partial_\varepsilon^2\mathscr E[g_{\exp},\phi_0+
\varepsilon\chi_1]|_0=D^2\mathscr E_0[(h_1,\chi_1)]^2+
D\mathscr E_0[b,0]\). It is therefore legitimate to use (9) for *this*
projection, but not to identify its unprojected second derivative with the
linear-path Hessian. This invariance concerns a scalar second-order path
completion with the fixed first-order longitudinal gauge. It is not a
claim that the bare quadratic source is invariant under arbitrary
first-order diffeomorphisms.

## 4. Scalar source, full angular decomposition, and radial coefficient

The background scalar has no angular gradient. Hence

\[
\big(D^2[-\phi_i\phi_j]\big)=-2D_i\chi_1D_j\chi_1.
\]

Every term in \(D^2[-2U(\phi)g_{ij}/3]\) is a multiple of \(q_{ij}\)
for the linear path. Its TT projection vanishes. Combining this with (8)
gives

\[
(D^2\mathscr E_0)_{ij}^{\rm TT}
=-(12u^2+2v^2)\Pi_{\rm TT}(D_iY D_jY).
\tag{11}
\]

Let \(B=\operatorname{diag}(1,-1,0,0)\), form the Weyl part \(C\) of
\(R_{ikjl}=B_{ij}B_{kl}-B_{il}B_{kj}\), and define
\(\mathcal T_{ij}=C_{ikjl}n^kn^l\), pulled back to sphere coordinates.
Its sectional normalization is

\[
C_{1212}=C_{3434}=-\tfrac13,\qquad
C_{1313}=C_{1414}=C_{2323}=C_{2424}=\tfrac16.
\]

With \(F_2=n_1^2+n_2^2-1/2\) and
\(F_4=Y^2-(n_1^2+n_2^2)/2+1/12\), the exact published identity is

\[
(DY\otimes DY)^{\rm TF}
=\tfrac12(\operatorname{Hess}F_2)^{\rm TF}
 +\tfrac13(\operatorname{Hess}F_4)^{\rm TF}-\tfrac43\mathcal T.
\tag{12}
\]

The tensor is transverse and trace-free, so
\(\Pi_{\rm TT}(DY\otimes DY)=-4\mathcal T/3\). Equations (11)--(12)
prove (1), including the total geometric contribution.

For the actual master write \(P=\mathcal P(y,t)Y\). Then

\[
u=\mathcal P/\rho^2,\quad v=-3\mathcal P_y/(\rho^2p),\qquad
\boxed{\mathcal C(y,t)=16u^2+\tfrac83v^2
=\frac1{\rho^4}\left(16\mathcal P^2
                   +\frac{24}{p^2}\mathcal P_y^2\right).}
\tag{13}
\]

The divided expression is used only where the exterior chart has
\(p\ne0\). The smooth ambient reconstruction gives the undivided form.
Put \(r=\sqrt{s}>0\), \(dy=\Omega dr\), \(\rho=r\Omega\),
\(p=2rf/\Omega\), \(f=\Phi_s\), and
\(F=\mathcal F(r,t)Y\), \(P=r^4F\). Since the ambient Euler field is
\(E=r\partial_r\) at fixed \(t\),

\[
u=\frac{r^2}{\Omega^2}\mathcal F,\quad
v=-\frac3{2\Omega^2f}(r\mathcal F_r+4\mathcal F),\qquad
\boxed{\mathcal C=\Omega^{-4}\left[
16r^4\mathcal F^2+
\frac6{f^2}(r\mathcal F_r+4\mathcal F)^2\right].}
\tag{14}
\]

For real fields, \(\mathcal C\ge0\), and it is strictly positive wherever
\(u\) or \(v\) is nonzero. In particular, a nonzero master amplitude
has a nonzero full quadratic TT source somewhere in its positive-cone
support. The scalar stress witness has therefore been upgraded here to
a full-source statement for this aligned quadrupole and gauge.

## 5. Implication for continuation and the boundary of the calculation

If the second-order correction is kept entirely in scalar angular sectors
in this fixed first-order gauge, its linearized Einstein angular TT block
is zero. Equation (1) rules out such a scalar-only second-order ansatz for
a nonzero aligned quadrupole. A tensor correction is required.

This statement supplies a bulk source. It does not solve the sourced
tensor evolution, material-shell conditions, or the through-contact
initial-boundary problem. The published exact theorem remains a strict
precontact theorem. This calculation supplies no novelty claim, global
stability claim, or higher-dimensional Big Bang conclusion.

## 6. Verification contract

`check_full_quadratic_tt.py` independently constructs five-dimensional
Christoffels and Ricci from both (9) and the linear Taylor metric in a
flat dS4 chart. Its \(\Psi\) is an arbitrary function of all five
coordinates and its \(\rho\) is arbitrary. It checks angular off-diagonal
components and diagonal differences; hence it verifies the entire
angular trace-free Ricci structure without assuming a separated master
equation or background equations. It then independently checks (12),
the trace, and every divergence component as rational identities in
stereographic coordinates on \(S^3\), and checks a nonzero component.
The local Ricci and global sphere stages together verify the projection
without conflating a pointwise trace-free part with the nonlocal TT
projection. Algebraic success is not a PDE existence certificate.

Both of these commands passed on 14 October 2026 with Sympy 1.14.0:

```sh
python check_full_quadratic_tt.py
python -O check_full_quadratic_tt.py --output full_quadratic_tt_results_optimized.json
```

The receipts are `full_quadratic_tt_results.json` (optimization level 0,
62.020 seconds) and `full_quadratic_tt_results_optimized.json` (level 1,
35.053 seconds). They record the same checker SHA256,
`3526d3667a95071e190f90dc5f8fcfc69ecb78fe92cb49e1671052a59721c34b`.
Explicit guards rejected omission of geometry, the wrong geometric sign,
the missing Taylor factor two, and the wrong TT projection sign. The
guards use `raise` rather than Python's removable `assert` statements.
