# Targeted independent review of the explicit background

15 October 2026. Analytical review of the newly inserted subsection `sec:explicit-background`, its stated witness, and the title/abstract scope. No manuscript edits were made by this reviewer. This is an internal mathematical review, not external peer review or a novelty certification.

**Reviewed final source:**
`manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex`

**SHA-256:** `f23e98635c8e0e1cf0f9e73e93c5149e3186cb5bb597f96b3f9cf896094dd845`

**Verdict:** the explicit background solves the stated five-dimensional bulk and shell equations. The displayed radial positivity and the specific pre-contact witness are valid. No blocking mathematical or title/abstract scope defect was found. The example removes the earlier need to rely on an externally registered background for nonemptiness of the hypotheses. It deliberately uses different fixed functions \(U,\sigma\), as the manuscript now states.

## 1. Independent bulk calculation

For \(g=e^{2u}\eta\) in five dimensions and \(u=u(s)\), the conformal Ricci formula gives
\[
R_{AB}=-3\partial_A\partial_Bu-\eta_{AB}\Box_\eta u
+3(\partial_Au\partial_Bu-\eta_{AB}(\partial u)^2).
\]
Here \(\partial_As=2X_A\), \(\Box_\eta s=10\), and \((\partial s)^2=4s\). Thus the coefficients of \(X_AX_B\) and \(\eta_{AB}\) are respectively
\[
12(u_s^2-u_{ss}),\qquad -16u_s-4su_{ss}-12su_s^2.
\]
For \(z=1+s\), \(u=-2\log z\), one has \(u_s=-2/z\), \(u_{ss}=2/z^2\), and therefore
\[
R_{AB}=\frac{24X_AX_B+(32-24s)\eta_{AB}}{z^2}.
\]
Meanwhile \(\Phi_A=2\sqrt6X_A/z\), and the fixed potential in source lines 423–439 obeys
\[
U(\Phi)=84z^2-36z^3=(48-36s)z^2.
\]
Consequently \(R_{AB}=\Phi_A\Phi_B+\frac23U(\Phi)g_{AB}\), which is the correct trace-reversed Einstein equation for the action normalization \(R/2\) in five dimensions.

The scalar equation also checks independently. The conformal transformation of the scalar wave operator is
\[
\Box_g\Phi=e^{-2u}\left[\Box_\eta\Phi+3\,\eta^{AB}(\partial_Au)(\partial_B\Phi)\right].
\]
Since \(\Phi_s=\sqrt6/z\), \(\Phi_{ss}=-\sqrt6/z^2\), this is
\[
z^4\left(4s\Phi_{ss}+10\Phi_s+12su_s\Phi_s\right)
=\sqrt6(10-18s)z^2.
\]
Direct differentiation of the prescribed potential gives the same value:
\[
U'(\Phi)=\frac{168}{\sqrt6}z^2-\frac{108}{\sqrt6}z^3
=\sqrt6(10-18s)z^2.
\]
Thus source lines 429–445 have the correct coefficients and signs. The scalar equation is not merely inferred from a Bianchi identity in this check.

## 2. Outward shell signs and regularity

For the interior \(s\le b^2\), the normal \(n^A=X^A/(b\Omega_b)\) is spacelike unit and outward: \(n(s)=2b/\Omega_b>0\). Under the conformal rescaling, the hyperboloid curvature is
\[
K_{ab}=\Omega_b^{-1}\left(\frac1b+2bu_s(b^2)\right)h_{ab}
=\frac{(1+b^2)(1-3b^2)}b\,h_{ab}.
\]
The outward scalar derivative is
\[
n\Phi=\frac{2bf_b}{\Omega_b}=2\sqrt6\,b(1+b^2).
\]
At \(b=1/2\), these give \(a_b=5/8\), \(p_b=5\sqrt6/4\). A four-dimensional shell has \(K_i=4a_b\), so each side contributes \(-3a_bh_{ab}\) to \(K_{i,ab}-K_ih_{ab}\). The two-side Israel residual is therefore \((\sigma_b-6a_b)h_{ab}=0\) for \(\sigma_b=15/4\). The scalar residual is \(2p_b+\sigma'_b=0\) for \(\sigma'_b=-5\sqrt6/2\). This verifies the outward-normal sums and single-tension normalization in source lines 447–461.

On the full physical slab \(|T|\le1/2\), \(-1/4\le s\le1/4\), hence \(3/4\le z\le5/4\). All displayed Cartesian fields are smooth on each interior, with strictly positive conformal factor and nonzero scalar derivative. At \(T=0\), the time derivatives of the bulk fields vanish and the shell has zero velocity; the initial canonical momenta and rapidities vanish. Rotational invariance is immediate.

The nonzero tension gives a derivative jump at the material shell, so the glued metric is piecewise smooth. My initial wording recommendation has been incorporated: the abstract now says “an exact background with smooth interiors” (line 44), and lines 463–466 expressly distinguish smooth center/cone behavior from the junction-fixed jump across the shell. This clarification is resolved in the pinned source.

## 3. Exact positivity and the specified bump

At \(T=0\), the actual weights are
\[
L=-\frac4{1+r^2},\qquad
m_r=\frac38r^3(1+r^2)^4,\qquad
n_r=6r^7(1+r^2)^2.
\]
Their radial zeroth-order combination is
\[
\frac{m_r(2+9r^2L)+n_r}{m_r}
=2-\frac{36r^2}{1+r^2}+\frac{16r^4}{(1+r^2)^2}
=\frac{2-32r^2-18r^4}{(1+r^2)^2}.
\]
The numerator has derivative \(-64r-72r^3<0\) for \(r>0\), and at \(r=1/5\) it equals \(432/625>0\). Thus the radial integrand is the sum of a nonnegative derivative term and a strictly positive multiple of \(Q^2\) throughout the stipulated interval. Any nonzero real \(Q\in C_c^\infty((1/10,1/5))\) has \(\mathcal B_{\rm rad}(Q,Q)>0\). This direct argument correctly eliminates the need for the thin-annulus estimate in this particular example (source lines 469–483).

The displayed exponential bump is smooth and flat at both endpoints, positive between \(11/100\) and \(19/100\), and has support compactly contained in \((1/10,1/5)\). Its very small amplitude is mathematically nonzero; no floating-point evaluation of the bump or its integral is needed for positivity.

## 4. Causal support at the stated time

Let \(a=11/100\), \(d=19/100\), and \(t_*=1/32\). The two geometric restrictions are strict:
\[
t_*<a/2=11/200,\qquad t_*<(b-d)/2=31/200.
\]
For \(|T|\le t_*\), finite propagation gives
\[
\frac{63}{800}=a-t_*\le R\le d+t_*=\frac{177}{800}<b.
\]
More sharply, throughout that causal support,
\[
\frac{209}{40000}=a^2-2at_*
\le s=R^2-T^2
\le d^2+2dt_* =\frac{1919}{40000}.
\]
Hence the support stays strictly in \(s>0\), away from the cone, center, and physical shell, and wholly inside the actual smooth \(f\ne0\) background region. These rational endpoints and the positivity numerator were independently checked with exact `fractions.Fraction` arithmetic. No shell-contact boundary theorem is needed. With \(\omega=1\), the previously derived coefficient \(16\pi^2\mathcal B_{\rm rad}(Q,Q)/3\) is strictly positive.

## 5. Title/abstract and remaining publication questions

The title describes the result accurately. The abstract's explicit-background claim is now supported by a fully specified solution and fixed potential/tension. Its distinction between a necessary spacetime obstruction and the finite commuting-matrix geometry is maintained. It does not claim a sufficient PDE continuation criterion, a new general Noether formalism, or dynamical instability.

No publication-level mathematical blocker was identified in this targeted review. Whether the explicit application provides enough new significance for a particular journal remains a literature and editorial question. This calculation supplies no acceptance assurance and makes no priority claim for the exact background. It neither validates nor replaces the separate registered HDBLAST potential/tension branch; source lines 491–494 state that distinction explicitly.
