# Fresh adversarial review of the integrated draft

Reviewed and refreshed against the final integrated TeX on 15 October 2026. This is an AI-assisted mathematical review of the actual integrated sources, not external peer review or an acceptance prediction. I did not use earlier PASS summaries as evidence. No manuscript or repository files were modified by this review.

## Conclusion

**No blocking mathematical defect was identified in the final reviewed version. No unresolved actionable mathematical or citation-scope finding remains from this review.** The central result is a necessary obstruction for parameter-\(C^2\), classically regular exact families through a specified actual radial background. The final draft now gives such an exact background, fixed potential/tension, and specified obstructed first-order state inside the paper. Its assumptions and its distinction from nonlinear sufficiency are materially narrower than a general linearization-stability or shell-contact theorem, and the draft respects those distinctions.

The finite commuting-matrix material also received a fresh delegated check. That check independently confirmed the \(C^1\) velocity criterion, global real Cayley curve, Frobenius inverse norm, six charge constants, and radial integration by parts. This delegation is another AI-assisted check, not human peer review.

The three clarifications raised in the first review were nonblocking. All three are now resolved in the actual main source as described below. General background hypotheses remain conditional for other potentials/tensions; the explicit example establishes nonemptiness and supplies a complete instance.

## Resolved findings

### 1. Exact background and witness are now exhibited inside the paper — resolved

**Location:** `manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex`, lines 269–282, equation `eq:background`; subsection lines 420–494, `sec:explicit-background`, equations `eq:explicit-background` through `eq:explicit-positive`; abstract lines 41–51.

The original review requested either explicit qualification of the assumed nonlinear background or an actual example in the standalone article. The final source retains the correctly conditional general result and adds the example independently derived in `review/fresh/EXPLICIT_BACKGROUND.md`:


\[
 \Omega(s)=(1+s)^{-2},\qquad
 \Phi(s)=\sqrt6\log(1+s),\qquad
 U(\phi)=84e^{2\phi/\sqrt6}-36e^{3\phi/\sqrt6}
\]

The displayed Ricci and scalar calculations check directly and solve the five-dimensional bulk equations for \(s>-1\). At a doubled shell \(s=b^2\), put \(\phi_b=\Phi(b^2)\). Its outward principal curvature is

\[
 a_b=\frac{(1-3b^2)(1+b^2)}b,
\]

and its scalar outward derivative is \(p_b=2\sqrt6\,b(1+b^2)\). Thus the two identical-copy junction equations require

\[
 \sigma(\phi_b)=6a_b,\qquad
 \sigma'(\phi_b)=-4\sqrt6\,b(1+b^2)
\]

The final text sets \(b=1/2\), yielding \(a_b=5/8\), \(p_b=5\sqrt6/4\), \(\phi_b=\sqrt6\log(5/4)\), and the fixed affine \(\sigma=15/4-(5\sqrt6/2)(\phi-\phi_b)\). Both residuals vanish with the stated outward-normal conventions. On every full interior ball for \(|T|\le1/2\), \(-1/4\le s\le1/4\), so \(3/4\le1+s\le5/4\). This verifies smoothness through the center/cone within each interior, uniform metric nondegeneracy, timelike shell, \(f\ne0\), and time-symmetric rotationally invariant initial data. The final abstract and subsection now correctly distinguish these smooth interiors from the piecewise smooth metric across the material shell, where the nonzero extrinsic-curvature jump is fixed by the junction equations.

The inserted radial coefficients also check: \(L=-4/(1+s)\), \(m_r=(3/8)r^3(1+r^2)^4\), \(n_r=6r^7(1+r^2)^2\). Hence the zeroth-order coefficient equals \(m_r(2-32r^2-18r^4)/(1+r^2)^2\); its numerator decreases and is \(432/625>0\) at \(r=1/5\). The displayed bump has support \([11/100,19/100]\), compactly contained in \((1/10,1/5)\), and is smooth/flat at both endpoints. For the stated \(T_*=1/32\), taking the coarse envelope \(r_1=1/10\), \(r_2=1/5\) gives

\[
 T_*<r_1/2=1/20,\qquad T_*<(b-r_2)/2=3/20,
\]

and on its entire two-sided causal support

\[
 R\ge11/160>|T|,\quad R\le37/160<b,\quad
 3/800\le s\le21/400.
\]

Thus this specified time slab is valid with a positive margin from both cone and shell and stays in the actual coefficient domain. With \(\omega=1\), the calculated obstruction is strictly positive. The final paragraph explicitly distinguishes this illustrative potential/tension from the registered HDBLAST model and makes no exact-solution priority or long-time stability claim. No action remains on this finding.

Using the actual tighter support endpoints gives the independently checked bounds \(209/40000\le s\le1919/40000\); the coarser bounds above already suffice.

### 2. “Admissible pure gauge” is now defined — resolved

**Location:** `manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex`, lines 412–418, immediately after `eq:explicit-obstruction`.

The final text now defines the joint diffeomorphism action on bulk fields and shell embeddings, requiring a local flow near the compact initial section. This supplies the local exact diffeomorphic family needed by the argument without asserting completeness of an arbitrary smooth generator on an entire slab. The nonzero obstruction excludes such a pure gauge tangent. No action remains on this finding.

### 3. Historical-source limitation is now linked from the standalone article — resolved

**Location:** `manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex`, lines 79–94 and 555–562; `literature/PRIMARY_COMPARISON.md`, “Classical sources and access limits”; `REFERENCES.bib`.

The draft correctly frames the work as a physical-interface specialization of classical Noether/linearization-obstruction methods. It now links the primary-source comparison and identifies the remaining historical access gaps directly. The package/supplement is linked in the scope section. The comparison expressly says that complete original proofs of several classical sources and the Hájíček–Kijowski erratum were not inspected in this reassessment. The manuscript neither imports a numbered unverified theorem nor claims nonlinear sufficiency; I found no citation-scope error or misleading novelty claim. The full-primary gaps remain open and must not be represented as completed by this AI-assisted review. No action remains on this finding within the present deliberately bounded claim scope.

## Direct mathematical checks and their limits

### Exact cut identity and arbitrary second corrections

**Locations:** main source lines 148–196 (`eq:boost`, `eq:canonical`, `eq:tilt-identity`, `eq:exact-functional`), lines 198–264 (`eq:invariant-background`, `eq:first-zero`, `eq:hessian`).

The sign in the tilt identity checks directly: differentiating \(n=\cosh\eta\,s+\sinh\eta\,u\) along \(J\) and contracting with \(\tau\) gives \(K(\tau,J)=-J\eta-k(s,J)\). With \(\pi^{ab}=\sqrt\gamma(k^{ab}-k\gamma^{ab})/2\), integration by parts gives exactly the boundary term \(-\sqrt c\,k(s,J)\). Its cancellation against the mixed Israel residual leaves \(+\sqrt c\,J(\eta_++\eta_-)\). The mixed trace and tension projections vanish because \(h(\tau,J)=0\). No dimension-specific trace factor enters this calculation.

The density invariance assumptions are sufficient to annihilate the remaining first variation: \(\int\pi_0\mathcal L_J\delta\gamma\), \(\int\varpi_0J\delta\phi\), and \(\int\rho_{S0}J\delta\vartheta\) become background density Lie derivatives by integration by parts. There is no flux because \(J\) is tangent to the common cut and the cut is closed. Thus \(D\mathcal C_J=0\) holds on arbitrary exactly matching variations, independently of first-order support. This is precisely the point required to allow every compatible second-order bulk/shell/angular correction. The twice-the-quadratic-coefficient Hessian convention and chart independence at a critical point are consistent.

The assumed parameter-\(C^2\), piecewise-\(C^k\), \(k\ge4\) classical regularity is more than enough for these local products, normalizations, traces and two parameter differentiations after a \(C^2\) family of collar pullbacks. The theorem does not exclude continuations having weaker parameter regularity; the draft explicitly limits its claim accordingly.

### Genuine precontact solution and reconstruction sufficiency

**Locations:** main source lines 269–349 (`eq:master`, `eq:reconstruction`, `eq:residual-factorization`, `eq:precontact-time`); Appendix A, lines 577–621, `app:reconstruction`.

The principal symbol of \(\mathscr L\) is that of \(\Box_\eta\); its first/zeroth order coefficients are smooth where the actual background has \(f\ne0\). A smooth coefficient extension to Minkowski space supports the cited normally hyperbolic Cauchy theorem and causal-support estimate. The bound \(T_*<r_1/2\) gives \(R>|T|\) on the entire causal support; \(T_*<(b-r_2)/2\) keeps it strictly below \(R=b<\sqrt{b^2+T^2}\). The additional explicit shrinking to the actual \(f\ne0\) neighborhood closes the local-coefficient requirement. These facts permit a smooth zero extension across a neighborhood of the cone, center and shell. Zero first-order embedding and field jets near the shell satisfy the full linearized interface equations there without an interface evolution theorem.

The appendix is sufficient to check the residual factorization, rather than merely a necessary master equation. In particular, \(\Psi=P/\rho^2\), \(\chi=-3P_y/(\rho^2p_b)\) make the mixed residual vanish, and the other Einstein residuals are the stated multiples of \(C=-3\mathcal M[P]/\rho^2\). Converting \(\mathcal M[r^4F]=\Omega^{-2}r^4\mathscr LF\) back to ambient components gives the tensor factor in `eq:residual-factorization`. The scalar factor is also consistent: with \(C=-3r^2\Omega^{-4}\mathscr LF\) and \(p_b=2rf/\Omega\), \(p_b\delta\mathscr S=(\partial_y+4a)C\) gives exactly \(-3(\mathcal E+6)\mathscr LF/(2\Omega^4f)\). Thus \(\mathscr LF=0\) implies all linearized bulk equations on the support. This establishes sufficiency for the displayed reconstruction, not a classification of all scalar perturbations.

### Charge normalization, positivity, and non-gauge inference

**Locations:** main source lines 350–418 and 496–525, `eq:initial-reconstruction` through `eq:matrix-charge`.

At \(T=0\), the metric reconstruction gives \(\delta k_{ij}=p(r^2\delta_{ij}-3x_ix_j)/\Omega\), whence \(\delta\pi^{ij}=-3px_ix_j/(2\Omega)\). The scalar time derivative has the shift \(D+5\) from \([\partial_T,\mathcal E]=\partial_T\). Contracting the Lie derivative of the metric produces \(6r^7qp/\Omega\) after the radial volume factor, and the scalar term gives \(m_r=9r^3/(4\Omega f^2)\). Both copies and the Hessian factor together give \(4\mathcal B(Jq,p)\).

The radial identity follows from \(rm_r'/m_r=3-2r^2L\). Its zeroth order coefficient is \(20m_r-(9/2)(rm_r)'=m_r(2+9r^2L)\). Since \(m_rr^2>0\) on an interior regular annulus and the negative potential coefficient stays bounded on shrinking annuli, the Dirichlet Poincaré estimate yields a strictly positive diagonal radial form for sufficiently thin intervals. This proves the required nonzero test value without claiming positivity of a gravitational energy.

The angular constants also check: \(\int_{S^3}n_1^2n_2^2=\pi^2/12\), so \(\int(JH)^2=4\pi^2/3\). For general quadrupoles, \(\int H_AH_B=\operatorname{vol}(S^3)\operatorname{tr}(AB)/12\), \(J_MH_A=H_{[A,M]}\), and the stated \(J_{ij}\) gives \(\operatorname{tr}(B[A,M_{ij}])=-2[A,B]_{ij}\). The six rotations therefore test all six components of the skew commutator, provided the fixed radial prefactor is nonzero. No statement about the full adjoint constraint kernel follows, and the draft says so.

### Finite geometry

**Locations:** `manuscript/quadrupole_geometry.tex`, Proposition `prop:inverse`, Theorem `thm:velocities`, Corollary `cor:regular`, and the dimension-four seeds.

With Frobenius-normalized symmetric/skew edge units, each nontrivial edge contributes the row \((c_{ij},d_{ij})\); its sole positive singular value is \(g_{ij}\). The inverse formula is the Moore–Penrose inverse, including on the nonsurjective strata. Its norm is exactly the reciprocal of the smallest positive joint gap; no extra \(\sqrt2\) belongs in that formula. At the replacement seed the joint gap is \(2/\sqrt5\), giving \(\|R\|=\sqrt5/2\).

The \(C^1\) necessity argument does not silently require a second derivative: projecting the exact commutator to a maximal joint block gives \(t^2P_s[a,b]P_s+o(t^2)=0\). The linear kernel relation cancels the off-block contributions to that diagonal block, leaving \([a_{ss},b_{ss}]=0\). The stated skew \(K\) reconstructs every off-block velocity. Its Cayley transform exists for every real parameter, is orthogonal, and has derivative \(K\); conjugating the commuting block-diagonal pencils gives the asserted global real-analytic commuting curve and acceleration. The equivalence with the second-order equation is therefore valid in this finite problem.

The repeated-block velocity set is not a vector space, correctly ruling out a local \(C^1\) embedded manifold there. Distinct joint eigenpairs admit a simple real pencil after excluding finitely many parameters; the centered polynomial powers give the claimed regular chart. The draft consistently prevents these facts from being read as a PDE right inverse or nonlinear spacetime sufficiency theorem.

## Review snapshot and limitations

Sources reviewed directly and reread on the final refresh: `manuscript/ROTATIONAL_SHELL_OBSTRUCTIONS.tex`, `manuscript/quadrupole_geometry.tex`, `manuscript/REFERENCES.bib`, and `literature/PRIMARY_COMPARISON.md`. Line references above refer to the main-source version containing the fully specified background/witness subsection, local-flow gauge definition, source-comparison link, coefficient-neighborhood shrinking, first-pullback-jet paragraph, and final smooth-interior/piecewise-shell wording. Equation labels are the more stable identifiers.

SHA-256 snapshot:

- Main TeX: `f23e98635c8e0e1cf0f9e73e93c5149e3186cb5bb597f96b3f9cf896094dd845`
- Matrix TeX: `d932ca4fb9311912bbeacc7dfd236ab28f6205cc62e36c6f7f04d9c4bc69f93c`
- Bibliography: `87aa5f33f42db35a1c83a15d2e10b35bf6c8eb7dcee42b34da757bf747dd59da`

I did not freshly obtain/read the inaccessible complete classical historical primary proofs or erratum, validate the registered earlier branch/contact construction, establish a shell-contact well-posedness theorem, or determine journal significance. The present source comparison honestly leaves the historical gaps open. No acceptance guarantee, novelty certificate, full nonlinear continuation theorem, or human specialist endorsement follows from this review.
