# Fournodavlos–Smulevici primary boundary check

Source: Grigorios Fournodavlos and Jacques Smulevici, **The initial boundary value problem in General Relativity: the umbilic case**, arXiv:2104.08851v1, 18 April 2021. [Primary HTML](https://arxiv.org/html/2104.08851); [version record](https://arxiv.org/abs/2104.08851v1). The received full HTML identifies v1 in its header; canonical metadata independently verifies title and authors. This is the umbilic paper, not their earlier totally geodesic paper arXiv:2006.01498.

Theorem 1 is **four-dimensional vacuum**, with a three-dimensional initial hypersurface, constant \(\lambda\), and the stated corner compatibility conditions.

Section 2.2, equation (16):
\[
g_{xt}=g_{xA}=0\quad\text{on }\mathcal T.
\]
Equation (17):
\[
\partial_xg_{ij}=2\lambda(g_{xx})^{1/2}g_{ij}\quad\text{on }\mathcal T.
\]
Equation (18), derived from their normal wave-coordinate condition:
\[
\partial_xg_{xx}=6\lambda(g_{xx})^{3/2}\quad\text{on }\mathcal T.
\]
These are boundary coordinate/Robin conditions, not Gaussian normal gauge throughout the collar.

Section 4.1 Lemma 1 writes the contracted Codazzi relation in their conventions:
\[
\operatorname{Ric}_{\alpha i}N^\alpha
=\partial_i\operatorname{tr}\chi-\nabla^j\chi_{ji}.
\]
Their assumption is quoted exactly as “for some constant \(\lambda\in\mathbb R\).” It makes this vacuum mixed curvature vanish.

Section 4.2 states: “Since we are not allowed to impose any further boundary conditions, these will have to be induced from the umbilic condition.”

Proposition 2 derives the subsidiary conditions
\[
\Gamma_x=0,\qquad \partial_x\Gamma_i=0,\qquad
\Gamma_\mu|_\Sigma=\partial_0\Gamma_\mu|_\Sigma=0.
\]
Its equation (24) uses the paper's definition of \(\Gamma_\mu\):
\[
\operatorname{Ric}_{\mu\nu}=\partial_\mu\Gamma_\nu+\partial_\nu\Gamma_\mu
\]
for a solution of their reduced equations. The normalization must not be substituted unexamined for a different generalized harmonic vector.

Their outward normal is \(N=-(g_{xx})^{-1/2}\partial_x\), with interior \(x\ge0\). HDBLAST instead uses increasing proper radius outward. The intrinsic dimension is four in HDBLAST rather than three in this theorem, and \(\lambda=\sigma(\varphi)/6\) is a field-dependent function. Thus coefficient 6, signs and literal coordinate Neumann conditions are not transferable without derivation.

In HDBLAST's own conventions, its physical even-copy laws imply the separate exact identity
\[
R_{na}=D_bK^b{}_a-D_aK
=-3D_a(\sigma/6)
=-\frac{\sigma'}2D_a\varphi
=(n\phi)D_a\varphi.
\]
It is the **Einstein–scalar mixed residual** which vanishes; neither term need vanish alone. This model identity suggests the same boundary-completion mechanism, but it is not a theorem supplied by the vacuum paper.

Appendix A Proposition 3/(29) supplies scalar Robin energy bookkeeping: integrate the actual Robin flux by parts in time, use the trace inequality and a sufficiently large bulk \(L^2\) coefficient, then Gronwall. A coupled matrix boundary estimate, source/corner compatibility and matter/gauge subsidiary propagation remain separate proof obligations for the current calculation.
