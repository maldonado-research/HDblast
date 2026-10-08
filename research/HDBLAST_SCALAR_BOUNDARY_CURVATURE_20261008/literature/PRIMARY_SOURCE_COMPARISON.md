# What the proposed scalar-curvature result must add

8 October 2026. This is a bounded, fresh primary-source comparison of seven papers, not an exhaustive priority search or an explanation of the APS decision. The uploaded assessment was treated as reference material. Its proposed scalar-response note and its claimed closest precedents were checked against newly fetched author PDFs over verified HTTPS. Original source bodies are held only in `private-source-cache/`; they should not be copied into a publication packet.

**Finding:** a focused Einstein–scalar boundary-value paper is a more defensible research target than a new bulk-scattering mechanism. The exact shell constraint, coupled scalar/metric junctions, regular curved interiors, and perturbative changes of brane position all have close published precedents. The potentially useful increment is a *quantitative, uniformly controlled detuning theorem for the precisely specified doubled-interior branch*. The proposed coefficient alone has not yet been demonstrated to meet an APS significance threshold.

## Primary papers actually read

1. O. DeWolfe, D. Z. Freedman, S. S. Gubser and A. Karch, *Modeling the fifth dimension with scalars and gravity*, **Phys. Rev. D 62, 046008 (2000)**, [doi:10.1103/PhysRevD.62.046008](https://doi.org/10.1103/PhysRevD.62.046008), [hep-th/9909134](https://arxiv.org/abs/hep-th/9909134). Read §§3.1–3.2, especially Eqs. (6)–(10), (18)–(20), and the parameter-counting/existence qualifications. They construct scalar backreaction and detuned de Sitter/anti-de Sitter branes; the curved problem is not just the flat superpotential flow.
2. C. Barceló and M. Visser, *Moduli fields and brane tensions: generalizing the junction conditions*, **Phys. Rev. D 63, 024004 (2001)**, [doi:10.1103/PhysRevD.63.024004](https://doi.org/10.1103/PhysRevD.63.024004), [gr-qc/0008008](https://arxiv.org/abs/gr-qc/0008008). Read minimal-coupling Eqs. (2.15)–(2.23) and general-frame Eqs. (2.40)–(2.45). Scalar and geometric jumps are established variational conditions; the Einstein-frame specialization does not supply a new junction law.
3. J. K. Ghosh, E. Kiritsis, F. Nitti and L. T. Witkowski, *De Sitter and Anti-de Sitter branes in self-tuning models*, **JHEP 11 (2018) 128**, [1807.09794](https://arxiv.org/abs/1807.09794). Read §§2.2–2.3, Eqs. (2.18)–(2.25), (2.41)–(2.55), and **Appendix D, Eqs. (D.1)–(D.7)**. Appendix D explicitly expands the equilibrium scalar/brane position in curvature and solves its first correction. This comparison is necessary; a generic claim that scalar adjustment shifts brane curvature or position is already anticipated there.
4. The same four authors, *Holographic RG flows on curved manifolds and quantum phase transitions*, **JHEP 05 (2018) 034**, [1711.08462v3](https://arxiv.org/abs/1711.08462). Read §4.1, Eqs. (4.9)–(4.18), and Appendix E, Eqs. (E.1)–(E.11); inspected the small-curvature expansion organization in Appendix D. They derive regular positive-curvature interior endpoints, distinguish a Lorentzian coordinate horizon from a Euclidean smooth cap, and analyze the scalar linearization near an AdS minimum. These are direct regular-branch precedents, beyond merely citing a curved-brane title.
5. M. Mintchev and L. Pilo, *Localization of Quantum Fields on Branes*, [hep-th/0007002](https://arxiv.org/abs/hep-th/0007002), Eqs. (32)–(34) and their interpretation. The Robin boundary spectral density and induced generalized free field coincide with the HDBLAST bath under the stated positive-Robin identification.
6. A. George, *The Massive Klein-Gordon Field Coupled to a Harmonic Oscillator at the Boundary*, [hep-th/0412067](https://arxiv.org/abs/hep-th/0412067), §2, Eqs. (14), (19)–(21), (26)–(27), (41)–(43). The fixed-spatial-momentum HDBLAST oscillator/bulk model, bound-state cubic, stability boundary and reflection amplitude map exactly to this model.
7. C. Charmousis, E. Kiritsis and F. Nitti, *Holographic self-tuning of the cosmological constant*, **JHEP 09 (2017) 031**, [1704.05075v3](https://arxiv.org/abs/1704.05075). The targeted follow-up read Appendix C, Eqs. (C.9)–(C.13), and §5.3, Eqs. (5.20)–(5.25). They integrate out both bulk regions to obtain an effective interface action for scalar and warp factor, recover junctions by extremizing it, and analyze scalar stability through an effective brane mass. This is the directly relevant effective-action precedent for assessing whether a susceptibility sign is a generic relaxation consequence.

Fresh access, byte counts and SHA256 hashes are recorded in `PRIMARY_FETCH_RECEIPTS_ALLOWED.json`, `PRIMARY_FOLLOWUP_RECEIPT.json`, `EFFECTIVE_POTENTIAL_FOLLOWUP_RECEIPT.json` and `EXTRACTION_RECEIPTS.json`. The two follow-up receipts include their text-extraction hashes. Initial sandbox network failures remain recorded separately; no insecure TLS fallback was used. The comparison does not claim to have reviewed all later citing papers.

## Exact translation to the older Einstein–scalar construction

Use HDBLAST subscript H and DeWolfe et al. subscript D. After accounting for their opposite signature/curvature conventions, the scalar normalization and functions are related by

\[
\phi_H=\sqrt2\phi_D,\quad U_H=2V_D,\quad
\sigma_H=2\lambda_D,\quad W_H(\phi_H)=W_D(\phi_H/\sqrt2).
\]

DFGK Eq. (8), \(V_D=W_{D,\phi_D}^2/8-W_D^2/3\), becomes exactly \(U_H=W_{H,\phi_H}^2/2-2W_H^2/3\). Their curved constraint (18) becomes

\[
A_y^2-H^2=-U_H/6+\phi_{H,y}^2/12.
\]

For the HDBLAST doubled interiors, the two outward normals meet at the shell. Imposing the reflected version of both DFGK jumps gives the HDBLAST one-sided conditions \(A_y=\sigma_H/6\) and \(\phi_{H,y}=-\sigma_{H,\phi_H}/2\). Consequently

\[
H^2=\frac{\sigma_H^2}{36}-\frac{\sigma_{H,\phi_H}^2}{48}+\frac{U_H}{6}
\]

is an algebraic consequence of the established constraint and junction system. It should be introduced as that consequence, not as a new field equation or new scalar-gravity mechanism. The flat tuning \(\sigma_H=2W_H\) maps to the corresponding signed DFGK tuning \(\lambda_D=W_D\). Different ends of their interval require their stated different orientation signs.

For Ghosh et al. (1807.09794), set their scalar \(\varphi_G=\sqrt2\phi_H\), bulk potential \(V_G=2U_H\), brane potential \(W_{B,G}=2\sigma_H\), and induced scalar curvature \(R_{B,G}=12H_H^2\) in four brane dimensions. Their symbol \(U_G\) is an **induced brane Einstein term**, not HDBLAST's bulk potential; set \(U_G=0\) for this comparison. Their flow function \(W_G=-6A_y\) is also not automatically the chosen HDBLAST generating function \(W_H\) away from a flat first-order branch.

Their matching problem glues a UV side to a regular IR side and may include an induced gravity term. The proposed HDBLAST problem glues two identical regular interiors with no UV exterior. This is a real difference in the boundary-value problem. It does **not** justify declaring the coefficient new, and it means that Appendix D's coefficient cannot simply be substituted without translating the global conditions and the control parameter. They perturb imposed UV curvature; the proposed study perturbs tension \(\delta f\).

## Claim-by-claim assessment

| Proposed claim | Established comparison | Defensible status and required increment |
|---|---|---|
| Both scalar and metric junction conditions must hold | DFGK (7); Barceló–Visser (2.15)–(2.23) | Established physics; repairing an omitted scalar condition is a necessary correction. |
| A detuned scalar brane can have de Sitter curvature | DFGK §3.2; Ghosh et al. 1807.09794 §§2–3 | Established construction. A new parameter choice alone does not establish significance. |
| Exact shell curvature identity above | Directly follows from DFGK (7),(18) after the displayed normalization | Established consequence, useful as a consistency check. |
| A regular curved interior has vanishing scalar derivative at its cap | Ghosh et al. 1711.08462 §4.1 and Appendix E | Established local endpoint behavior. A finite truncated cone series is not a global existence proof. |
| The shell scalar shifts under a small change of boundary data | Ghosh et al. 1807.09794 Appendix D has an explicit perturbative brane-position response | Established broad method; compare the different control parameter and doubled-interior geometry explicitly. |
| \(c_2=f_0^2/36-kf_1^2/[24(w-2k)]\) for \(w>4k\), \(f_0>0\) | Exact identity plus regular growing-mode matching; earlier HDBLAST special case already in DOI 23111008 | A candidate general expression in the uploaded note. No identical expression was identified in the specifically read sections, which is insufficient for a priority claim. Its branch assumption and uniform remainder are decisive gaps. |
| Independence from \(W'''_0\) and \(f''_0\) at second order | Follows directly from the coefficient and perturbative order counting | Useful universality within this expansion, not proof of universality beyond its hypotheses. |
| Forty-eight numerical configurations and nine algebra identities | Internal reproduction/consistency results | Supporting evidence, not a nonlinear existence theorem, outward-rounded error bound, peer review, or APS eligibility certificate. |
| Robin spectral density and a four-dimensional continuum realization | Mintchev–Pilo (33)–(34) | Direct published antecedent. Do not use this as the main original contribution. |
| Coupled boundary oscillator scattering and stability threshold | George (14),(19),(21),(26),(27),(43) | Exact fixed-momentum antecedent, with the translation below. Static continuum/bound occupations are already diagonal in (41)–(42). |
| A higher-dimensional blast caused our Big Bang | No result in this static comparison establishes it | Unestablished. A controlled static mathematical result would remain useful without making this claim. |

## Exact comparison of the bulk model

George's coordinate/parameters map as

\[
x=-y,\quad\lambda_G=c,\quad\beta_G=-g,\quad
m_G^2=M^2+|\mathbf k|^2,\quad\mu_G^2=m^2+|\mathbf k|^2.
\]

Equation (26) becomes \((s+c)(s^2+m^2-M^2)-g^2=0\). Equation (27) at \(\mathbf k=0\) is the same nonnegative stability boundary \(g^2\le m^2(c+M)\), while HDBLAST chooses its strict interior. Equation (43) becomes

\[
R(p)=\frac{(ip+c)D-g^2}{(ip-c)D+g^2},\qquad D=m^2-M^2-p^2.
\]

Mintchev–Pilo's parameter \(\eta=c>0\) in Eq. (34) gives exactly

\[
\rho_b(z)=\frac{\sqrt{z-M^2}}{\pi(c^2+z-M^2)}\Theta(z-M^2).
\]

These are substantive equalities, not merely similar terminology. `verify_prior_art_translations.py` checks eight convention/algebra identities independently with SymPy. Its pass status certifies these algebraic translations only. A time-dependent coupling with external work can change capture, but quench physics has its own prior-art burden and still does not distinguish a bulk geometry from an equivalent continuum bath.

## What would make a focused paper materially stronger

1. **State one mathematical problem and one new result.** A suitable target is existence and local uniqueness of a regular doubled-interior Einstein–scalar branch for a specified class of \(W,f\), together with the signed second-order curvature susceptibility and a uniform remainder. Avoid making unrelated endpoint certificates and oscillator scattering carry the significance claim for this paper.
2. **Control the singular limit.** At zero detuning the shell radius diverges. Ordinary smooth dependence of ODE solutions on a fixed finite interval, a small shooting residual, or a nonzero Jacobian at one positive detuning does not alone prove a uniform \(\delta\to0^+\) branch. Establish estimates that remain valid as the interval grows, including the regular cap and both boundary conditions.
3. **Make the error statement quantitative.** Prefer a theorem or validated enclosure with explicit \(\delta_0,C\) and hypotheses ensuring \(|H^2-a_1\delta-a_2\delta^2|\le C\delta^3\) on \(0<\delta\le\delta_0\). If only an asymptotic remainder is proved, label it accurately. State any compact parameter range, lower bounds on \(k,f_0,w-4k\), and regularity/norm bounds on \(W,f\). The denominator \(w-2k\) has no resonance in the assumed \(w>4k\) domain.
4. **Compare like physical branches.** For \(f_1\ne0\), the constant-scalar metric-only comparator generally violates the scalar junction. A negative adjustment relative to that comparator is a mathematical correction, not a demonstrated reduction relative to another admissible universe. Show whether it changes an allowed curvature range, parameter inference or stability conclusion in an actual family of complete solutions.
5. **Separate existence from stability and cosmology.** A regular radial solution need not be dynamically stable. A Lorentzian radial coordinate horizon is not by itself a globally completed cosmological spacetime. The problem contains no demonstrated blast, lasting energy deposition, reheating, or observational prediction.
6. **Write the contribution against the exact baselines.** The opening should acknowledge DFGK's curved scalar branes, the standard scalar junctions, Ghosh et al.'s regular endpoints and perturbative matching, then identify precisely which uniform theorem or useful quantitative restriction they do not supply for this stated problem. An expert comparison should precede any first-result claim.

No manuscript submission, appeal, email, GitHub write or Zenodo publication was performed by this literature task. No editorial reason beyond the supplied letter is inferred. Meeting these research targets would strengthen a new manuscript; journal consideration and acceptance remain editorial judgments.

## Assessment of the stronger finite-detuning result derived in this session

The independently developed `../independent-math/INDEPENDENT_MATHEMATICAL_REVIEW.md` changes the perturbation parameter. It fixes any \(\delta>0\) and expands in a separate scalar-coupling parameter \(c\) in \(\sigma=2W+\delta[f_0+c g(\phi)]\), with \(g(\phi_0)=0\). This is a defensible way to make a stronger, more focused claim while respecting the long-interval difficulty of the original detuning expansion.

Its proposed result is: a locally unique analytic regular branch about the exact constant-scalar solution at each fixed positive detuning; an exact quadratic curvature susceptibility

\[
H^2(c)=H_0^2+C_2c^2+O(c^3),\qquad
C_2=-\frac{\delta^2g_1^2}{48(D+w)^2}
\bigl[4A(D+w)-D'(y_0)\bigr]<0
\]

for \(g_1\ne0\), where \(A=k+\delta f_0/6\), \(D=F'/F\) is the regular linear scalar Dirichlet-to-Neumann ratio at the unperturbed shell, and \(y_0=k^{-1}\operatorname{arccoth}(A/k)\). Its barrier proof of the positive bracket and its nonzero two-residual shooting Jacobian are substantive additions to the uploaded purely assumed-branch algebra. Mathematical review of those steps is a separate task; this subsection assesses their relation to the primary texts.

**Exact permissible comparison:** DFGK supplies the underlying field equations and curved-brane construction; Ghosh et al. supplies regular curved interiors and an explicit perturbative matching expansion for a different UV–IR problem. The selected sections do not present the displayed all-positive-detuning, weak-coupling sign formula for this doubled-interior problem. Thus the manuscript may present a *specified susceptibility theorem and sign proof*, with the older constructions acknowledged. It may not infer from this seven-paper comparison that the formula is the first such result in all literature or that it is sufficiently significant for PRD.

The phrase “all positive detunings” quantifies the **sign of the local coefficient at each detuning**. It does not mean all coupling strengths, a globally unique nonlinear branch, one common coupling radius down to \(\delta=0\), or a finite error bound at the registered coupling. The analytic Taylor remainder on compact positive-detuning intervals is an existential bound unless its constant and radius are evaluated. The limit of \(C_2/\delta^2\) recovers the uploaded coefficient, but does not justify exchanging detuning and coupling limits.

A careful contribution sentence would be: “Within an Einstein–scalar model formed by gluing two identical regular interiors, we establish a local scalar-coupling branch at every fixed positive tension detuning and derive a strictly negative quadratic curvature susceptibility in terms of the regular linear scalar response. We identify its small-detuning limit and the limitations of the associated remainder.” Replace “establish” with “derive under the stated analytic hypotheses” if unresolved gaps remain after independent mathematical review. A numerical radius or a new observable consequence would strengthen its practical value; reproduction counts alone would not.

## Is the sign merely a familiar relaxation effect?

The targeted primary follow-up to Charmousis–Kiritsis–Nitti is important because a sign can look novel when only its formula is searched. Their flat-brane effective action, schematically Eq. (C.12), contains \(e^{4A}[W_{IR}(\phi)-W_{UV}(\phi;C_{UV})-W_B(\phi)]\), together with a UV boundary contribution. The scalar and metric matching conditions follow by extremizing over the interface variables. The UV constant and fixed UV data are treated carefully and are not ordinary variables to minimize. The sign of a particular curved response cannot be read off by treating every one of these terms as a positive energy.

For comparison, an elementary positive-Hessian relaxation problem

\[
E(q,c)=E_0+\tfrac12q^{T}Kq+c\,j^{T}q+\tfrac12r c^2,
\qquad K>0,
\]

has stationary \(q=-cK^{-1}j\) and effective quadratic coefficient \((r-j^{T}K^{-1}j)/2\). This is a generic negative Schur-complement correction, not a new physical principle. This displayed algebra is our comparison, **not a theorem attributed to the cited paper about HDBLAST's coefficient**. Gravity includes constrained metric directions; \(H^2\) is not automatically the minimized scalar energy. A reduction from the proposed curved doubled-interior problem to a positive-Hessian effective potential has not been demonstrated here.

Accordingly, neither extreme is supported: the negative sign alone is not convincing originality, and this effective-action precedent does not by itself prove that the exact finite-detuning sign theorem has already been established for the stated boundary-value problem. The manuscript's value would rest on the precise theorem, rigorous domain, controlled error and a useful consequence, compared explicitly with these frameworks.

## Relation of the uniform detuning proof candidate to the cited expansions

The new `../branch-theorem/UNIFORM_BRANCH_THEOREM.md` proposes a different and stronger mathematical increment: for \(W\in C^4\), \(f\in C^3\), \(w>4k>0\), \(f_0>0\), it constructs a small-field regular branch on an interval whose shell endpoint tends to infinity as \(\delta\to0^+\). It gives a weighted linear inverse, a metric Volterra estimate, contraction inequalities, a unique shell in that small neighborhood, and formulas for sufficient \(\delta_*\) and \(K_H\) in a uniform \(K_H\delta^3\) curvature remainder. **At the time of this comparison, that proof is a candidate undergoing separate adversarial mathematical review.**

These explicit expanding-interval inequalities and uniform nonlinear remainder are more specific than the local cap series, parameter-counting argument, perturbative matching coefficient and effective-action formulas identified in the selected primary passages. If the proof survives review, a defensible statement is that this work supplies such a quantitative sufficient theorem for the **specified** doubled-interior problem. This comparison has not established that no equivalent theorem exists in other domain-wall, boundary-value or mathematical-relativity literature. Classical maximum principles, Green operators, Volterra equations, contraction mapping and the implicit-function theorem remain established tools; combining them here is an application whose significance must be assessed.

The constants are currently formulas rather than evaluated certificates for the 48 numerical configurations. An explicit conservative radius can be useful even if small, provided the manuscript distinguishes it from the tested sample range and explains what the proven bound allows a reader to conclude. A substantive physics consequence or a sharply useful theorem remains a more persuasive significance argument than the number of symbolic identities or decimal digits.
