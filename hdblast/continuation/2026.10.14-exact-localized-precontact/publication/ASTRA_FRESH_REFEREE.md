# Fresh independent mathematical review of the exact finite-regularity initial-data argument

14 October 2026. Internal mathematical agent review; not external peer review or a novelty assessment.

Reviewed proof: `publish-exact-precontact-20261014/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`.

Reviewed SHA-256: `95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b`.

**Judgment.** I find no substantive new gap in the fixed-finite-regularity exact initial-data theorem, interpreted with the checkpoint's explicitly imported local adjoint classifications. The Banach mapping, mixed differential orders, density conventions, weighted inverse and physical zero extension support an actual parameter-smooth solution curve. This judgment was obtained by recalculating those steps, without taking the prior internal acceptance as evidence. A focused independent subreview separately checked the scalar adjoint/API input and its background scope. The precontact Cauchy argument is assigned to a separate fresh reviewer and is not certified by this report.

There is one meaningful scope clarification: generic radial time-symmetric Einstein–scalar data with a nonconstant radial scalar do not by themselves imply the imported kernel classification. The proof's explicit conditional-input paragraph prevents this from being a new proof gap, but the theorem must retain that conditional reading. The specific checkpoint background or the two required adjoint kernels should be named explicitly when stating the result independently.

## 1. Linear operator and weighted isomorphism

The checkpoint uses

\[
H=\tfrac12 R+\frac{2}{\det\gamma}\left(\frac{(\operatorname{tr}\pi)^2}{3}-|\pi|^2\right)-\frac{\varpi^2}{2\det\gamma}-\tfrac12|d\phi|_\gamma^2-U(\phi),
\qquad M_i=-2D_j\pi^j{}_i+\varpi\partial_i\phi,
\]

with contravariant metric momentum density and scalar momentum density. At the time-symmetric background, fixing the scalar sectors leaves independent metric and metric-momentum blocks. The metric-only adjoint is

\[
P_g^*N=\tfrac12\bigl(\nabla^2N-\gamma_0\Delta N-N\operatorname{Ric}(\gamma_0)+N\,d\phi_0\otimes d\phi_0\bigr).
\]

Thus the matter term is present. This is the metric-only adjoint whose kernel the checkpoint classifies, not merely the kernel of a larger coupled operator. Its highest derivative map is the invertible map \(B\mapsto(B-\gamma_0\operatorname{tr}B)/2\). Delay Lemma 8.1 covers this operator with its smooth lower-order terms; the proof of that lemma explicitly absorbs those terms near the boundary. This does not require treating every injective-symbol system as if it had full-derivative injectivity.

For momentum, writing \(\pi=\rho_0\tau\), \(\rho_0=\sqrt{\det\gamma_0}\), and using \(\widehat M=M/\rho_0\) gives

\[
(P_m\tau)_i=-2\gamma_{0,ik}\nabla_j^0\tau^{jk},\qquad
(P_m^*Z)^{ij}=\nabla^i Z^j+\nabla^j Z^i.
\]

Integration against the fixed background one-form pairing gives precisely this adjoint and its raised tensor indices. The required API is the symmetric-divergence/Killing-operator estimate identified in Delay §9.3, not the inapplicable full-derivative hypothesis of Lemma 8.1. Chruściel–Delay's exponential boundary estimate is consistent with this use.

Delay Theorem 3.8 really states the dotted-space isomorphism

\[
\pi_{\mathcal K^\perp}\,w^{-2}Pw^2q^{2m}P^*:
\mathcal K^\perp\cap\mathring H^{k+2m}_{q,w}
\longrightarrow\mathcal K^\perp\cap\mathring H^k_{q,w}.
\]

The weights in the draft coincide with Delay's \(\phi=d^2\), \(\psi=d^{2(a-n/2)}e^{-\beta/d}\) for \(n=4\). The orders are \(m=2\) and \(m=1\), giving exactly \(k+4\) and \(k+2\) potential regularity. The radial weights and fixed background pairing are invariant under the reflection group. The four boost lapses transform by the nontrivial characters \(\sigma_a\), and the six rotations by \(\sigma_a\sigma_b\). Averaging therefore annihilates the entire imported adjoint kernel. Invariant sources are orthogonal to it, and uniqueness of the complementary inverse makes the solution invariant. The projected inverse is consequently an inverse of the full operator on the invariant spaces. No perturbed-kernel classification is needed.

## 2. Mixed-order nonlinear map

On a ball of radius comparable to \(q_i\), the four-dimensional volume factor in an \(L^2\) norm is \(q_i^2\). Thus the local weighted norm factor is \(w_iq_i^2\). With the draft's \(V_{H,i}\) and \(V_{M,i}\), direct scaling yields

\[
\|h\|_{H^{k+2}_{\rm sc}}
\lesssim w_i^2q_i^{4-2}\frac{V_{H,i}}{w_iq_i^2}
=w_iV_{H,i},
\]

\[
\|\Pi\|_{H^{k+1}_{\rm sc}}
\lesssim w_i^2q_i^{2-1}\frac{V_{M,i}}{w_iq_i^2}
=w_iq_i^{-1}V_{M,i}.
\]

The multiplier derivatives have the same bounds by the exponential-weight scaling estimates. Hence curvature remainders \(hD^2h\) and \((Dh)^2\), momentum squares \(\Pi^2\), and momentum-constraint remainders \(hD\Pi\), \(\Pi Dh\) all have local size at most \(Cw_i^2q_i^{-2}(V_{H,i}+V_{M,i})^2\). Lower derivative background terms fit this estimate. All four compactly supported seed sectors can be included with finite local constants; the prescribed scalar terms create no boundary singularity.

The exact momentum-density identity in the draft is correct:

\[
M_i=-2\gamma_{ik}\partial_j\Pi^{jk}
-2\Pi^{jk}\partial_j\gamma_{ki}
+\Pi^{jk}\partial_i\gamma_{jk}
+\varpi\partial_i\phi.
\]

For a weight-one density, the contracted connection terms cancel. Contracting the remaining Christoffel term with \(\gamma\) gives the displayed expression, so there is no extra derivative or inverse-metric factor hidden in this block.

Multiplying the remainder by \(w^{-2}\) and applying the local range norm factor \(w_iq_i^2\) leaves \(Cw_i(V_{H,i}+V_{M,i})^2\). Summability follows from \(\ell^2\subset\ell^\infty\) and \(\sum V_i^4\le(\sum V_i^2)^2\). The same reasoning gives the difference estimate. At each higher derivative order, extra metric directions contribute further factors \(w_iV_i\); there are still at most two momentum directions or two total physical spatial derivative powers in the nonlinear constraint expression. The resulting multilinear bounds are summable. Uniform metric positivity supplies the smooth inverse-matrix and determinant compositions. This establishes the claimed smooth map on an open neighborhood in the fixed Banach spaces, rather than only a formal Taylor expansion.

The dotted-range argument is valid: compactly supported smooth potentials produce compactly supported constraint differences, and the estimates extend them in the closed dotted range. Multiplication by the rapidly growing \(w^{-2}\) is performed with the displayed compensated estimates; it is not being presumed bounded on arbitrary physical Sobolev data.

## 3. Physical support and zero extension

Passing from scaled to physical \(H^{k+2}\) metric norms multiplies the highest derivative by \(q_i^{2-(k+2)}=q_i^{-k}\). For momentum, the corresponding factor is \(q_i^{2-(k+1)}\), combined with the already present \(q_i^{-1}\). Both local bounds are therefore exactly \(Cw_iq_i^{-k}V_i\).

For every fixed finite \(k\), \(\sup_i w_iq_i^{-k}<\infty\). Bounded overlap then gives continuous physical maps into \(H^{k+2}\) and \(H^{k+1}\). Approximation by compactly supported potentials proves that the outputs lie in the compact-support closures \(H_0^{k+2}\) and \(H_0^{k+1}\). Zero extension is consequently legitimate and introduces no boundary delta terms in the constraints. The fixed closed annulus stays away from shell and vertex, giving common open background collars.

The Banach implicit function theorem now applies with an actual bounded inverse. Since \(Pv=0\), the first correction derivative vanishes, and the asserted tangent and second-derivative formula follow. Nothing in this reasoning supplies one common parameter interval for all spatial Sobolev orders; the draft correctly avoids that claim.

## 4. Background scope: a concrete reason to preserve the conditional input

The checkpoint classification uses its special Lorentz-invariant background equations, including

\[
(1/\Omega)''=\tfrac13(\Phi')^2(1/\Omega).
\]

Radiality, time symmetry and \(\Phi'\ne0\) alone do not replace this information. As a local example on any closed annulus in \(s=r^2>1\), set

\[
\gamma=\Omega(s)^2\delta,\qquad
\Omega=\frac{\sqrt{s^2-1}}s,\qquad
\Phi=\frac{\sqrt6}{2}\log\frac{s-1}{s+1},\qquad U=0.
\]

Writing \(f=\log\Omega\), one has

\[
f'=\frac1{s(s^2-1)},\quad
f''=-\frac{3s^2-1}{s^2(s^2-1)^2},\quad
\Phi'=\frac{\sqrt6}{s^2-1}.
\]

The conformal Ricci formula in four dimensions gives

\[
\operatorname{Ric}_{ij}
=(-12f'-4sf''-8s(f')^2)\delta_{ij}
 +(-8f''+8(f')^2)x_ix_j
=\frac{24x_ix_j}{(s^2-1)^2}
=\partial_i\Phi\partial_j\Phi.
\]

Also \(\Delta_\gamma\Phi=\Omega^{-2}(8\Phi'+4s\Phi''+8sf'\Phi')=0\). These are valid smooth time-symmetric Einstein–scalar data on the annulus, indeed induced by the static product metric \(-dt^2+\gamma\), yet \(P_g^*1=0\). The constant lapse is reflection invariant. This example does not satisfy the imported checkpoint kernel assumption, and so does not refute the stated conditional argument. It does show why an abbreviated statement must not silently omit that assumption.

Suggested clarification: name the fixed checkpoint background \(g_0=\Omega(R^2-T^2)^2\eta\), \(\phi_0=\Phi(R^2-T^2)\), or explicitly require \(\ker P_g^*=\operatorname{span}\{\Omega x_a\}\) and \(\ker P_m^*=\mathfrak{so}(4)\). This is hypothesis precision, not a request to redo the new weighted analysis.

Primary passages directly checked: Delay §§2, 3 (Theorem 3.8), 8 (Lemma 8.1), 9.3, and Appendix 10; the checkpoint manuscript §§1–3 and reflection continuation; Chruściel–Delay Appendix D exponential boundary estimate.
