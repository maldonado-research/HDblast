# Primary brane-paper reading notes for the cone phase-space checkpoint
Public adaptation: executor paths and cache-location details are removed; scientific statements, source attribution, reading scopes and limitations are preserved. Original source bodies are not redistributed.

Reading date: 2026-10-13. Sources are the previously retrieved original arXiv papers identified in CACHED_BRANE_SOURCE_RECEIPTS.json; no network retrieval or secondary summaries were used. Page references below are printed paper pages (also the PDF page indices for the passages cited). These are reading notes and cautions, not novelty claims or a new regularity theorem.

## 1. Garriga and Sasaki, *Brane-World Creation and Black Holes*

Source: **hep-th/9912118v2**, arXiv version dated 7 March 2000, previously retrieved original PDF and its text conversion. Read §§I–III, printed pp. 2–11, §IV relevant conical-periodicity passage on p. 12, and §V, p. 14. The bulk-black-hole thermodynamic calculation was not reviewed in full.

- **The de Sitter slicing ends at a coordinate horizon of the background.** §II, p. 3, Eqs. (2.1) and (2.4), uses (a(r)=\ell\sinh(r/\ell)). The Lorentzian chart covers only the exterior of the Rindler horizon (r=0). On p. 4 the authors identify this horizon as the light cone from the center of symmetry and continue inside it by (r\to it_c), (Ht\to r_c-i\pi/2), obtaining an open-FRW description of the AdS interior. This passage makes the background extension explicit; it does not establish that arbitrary separated perturbations extend with it. The Euclidean metric (2.1), with a unit four-sphere and (a\sim r), provides the smooth background center; this local smoothness observation follows directly from the displayed metric.
- **Radial normalizability is explicitly compatible with a singular perturbed horizon.** §II, pp. 6–7, RS gauge (2.12), separated radial equation (2.13), and de Sitter spin-two equation (2.14), give the normalizable graviton zero mode (\widehat h_{\mu\nu}\propto a^{3/2}) and KK threshold (m=3H/2). §III, pp. 8–9, considers the factorized Einstein-slice family (3.1)–(3.3). On p. 9, immediately before (3.8), the unnumbered identity is
  
  \[
  C_{(5)}^2=\frac{C_{(\gamma)}^2}{a^4(r)}.
  \]
  
  The authors state that nonzero numerator gives divergent curvature at (r=0), including generic cosmological gravity-wave zero-mode perturbations; they explicitly call the preceding zero mode singular on the Rindler horizon. Their discussion also mentions null plane-wave exceptions in the flat-slice case. Accordingly this curvature contraction detects the stated generic singularity; its vanishing should not be promoted to a complete regularity criterion.
- **A separated mode can be a representation device rather than an individually regular physical perturbation.** The same p. 9 paragraph says mode expansions are useful instruments for constructing Green functions. It cites the flat Randall–Sundrum example in which horizon-singular individual KK modes yield a nonsingular field of isolated brane sources after use of the Green function. This is a cited flat-slice example and the authors' explanation of the de Sitter concern, not a demonstrated completeness/regularity theorem for arbitrary de Sitter superpositions.
- **Finite action and a vanishing small-radius boundary term do not establish curvature smoothness.** §III, p. 10, following (3.10), describes the spatial horizon section as closing off smoothly while explicitly retaining the divergent (C^2) at the cigar tip. The authors also state that small zero-mode fluctuations have finite action despite this curvature divergence. The unnumbered Gibbons–Hawking boundary expression after (3.11) scales as (a^3\partial_r a) and vanishes as (r\to0). They regard the singular instanton as a legitimate saddle, pp. 10–11; possible smoothing by quantum fluctuations is an argument/speculation, not a smooth classical-extension proof. §V, p. 14, reiterates that non-conformally-flat Einstein slices give singular bulk centers without changing the action value.
- **Global data matter.** §I, p. 2, discusses bulk boundary conditions and incoming flux; §II, p. 7, argues that their creation scenario supplies a compact initial hypersurface and avoids a past-Cauchy-horizon boundary-condition problem. Detailed five-dimensional quantum-fluctuation evaluation is explicitly left to future work there. Neither a covariant symplectic form nor a global gauge/charge quotient is constructed in the passages read.

For the present checkpoint: this source directly prevents equating a radial (L^2) graviton mode, finite quadratic/Euclidean action, or vanishing background small-radius boundary contribution with a globally smooth cone perturbation or a conserved physical symplectic norm.

## 2. Gen and Sasaki, *Radion on the de Sitter brane*

Source: **gr-qc/0011078v2**, arXiv version dated 4 December 2000, previously retrieved original PDF and its text conversion. Read §§I–IV, printed pp. 1–7, and Appendix A, pp. 7–10, including the complete thin-wall surface-action derivation. Reference pages were not substantively reviewed.

- **The radial orthogonality norm belongs to the projected-Weyl variable.** §II, p. 3, defines (\widehat E_{\mu\nu}=a^{n-2}E_{\mu\nu}), Eq. (2.8), with bulk equation (2.9), junction conditions (2.10), and a Green-function construction (2.13)–(2.14). On p. 4, Eqs. (2.15)–(2.17) separate (\widehat E=\psi_mY^{(m)}). The authors explicitly determine the radial orthogonality (2.19)–(2.20) by requiring self-adjointness: its weight is (a^{3-n}\,dz). This is an operator/Green-function normalization; it is not presented as the reduced metric/scalar action's kinetic or covariant symplectic norm.
- **Its special lowest mode is distinct from the graviton zero mode in source 1.** Here (\psi=\text{constant}) has (m^2=(n-2)K), with normalization (N_{(n-2)K}), Eq. (2.21), p. 4. For (n=4,K=1), (N=\cosh z_- -\cosh z_+\) diverges as the second brane approaches (z_-\to\infty). §III, p. 5, calls retaining this mode the zero-mode truncation and notes that (m^2=0) is absent from the projected-(E\) spectrum for (n>3). These mass labels and norms cannot be substituted directly for those of the metric tensor equation in Garriga–Sasaki.
- **Local displacement and physical relative displacement are distinguished.** §III, p. 6, Eqs. (3.8)–(3.11), identifies the scalar with brane displacement in RS gauge and intrinsic brane scalar curvature in Gaussian normal coordinates. The coupled displacement equation is (3.9); the Gaussian-normal curvature equation (3.11) is the trace of the effective Einstein equation, not an additional scalar degree of freedom. §IV, pp. 6–7, states that displacement alone cannot be detected because local coordinates can set it to zero; physical effects require bulk Weyl perturbations that propagate back to a brane. Relative displacement of two branes is physical. This local-coordinate argument does not classify globally smooth cone diffeomorphisms, their extension, or possible boundary charges.
- **The one-brane result has specific assumptions and a precise limit.** For a single positive-tension brane in vacuum AdS, Eq. (3.8) and the divergent (N\) make the special-mode Weyl response vanish, pp. 6–7. The authors conclude that no radion correction remains at linear order, apart from KK effects. They identify (z\to\infty) with the center of the spatial sphere, p. 6. This result is about this vacuum-AdS response problem; it is not a radion-absence theorem for arbitrary Einstein–scalar backgrounds or a cone phase-space construction.
- **Brane action terms must be varied together.** Appendix A starts from bulk Einstein–Hilbert terms, both extrinsic-curvature terms, tension, and matter, p. 7. The variation (A2), p. 8, includes variation of the integration domain and displacement of the boundary. Eqs. (A9)–(A16), pp. 8–10, retain these contributions before specializing to the (Z_2) system and obtaining the coupled surface action (A18), displacement equation (A19), and junction condition (A22). The ϕ-kinetic-looking term inside (A18) is not by itself a reduced gauge-invariant radion norm. §IV, p. 7, explicitly says that an effective action for the de Sitter radion still needs to be derived.
- **Negative mass squared is a separate issue from kinetic sign.** In the two-brane system the scalar equation has (m_s^2=-nK), §§III–IV, pp. 5–7. The source's tachyonic mass statement should not be rewritten as a ghost statement. It also warns on p. 7 that a flat-brane effective metric ansatz cannot simply be transferred to curved backgrounds.

The causal Green-function condition is left unspecified but must be common to all Green functions, §II, p. 3; obtaining (3.7), p. 5, also invokes an appropriate common causal condition. Homogeneous additions and global solution selection therefore should not be silently erased when reusing the displayed response formulas.

## 3. Kofman, Martin and Peloso, *Exact identification of the radion and its coupling to the observable sector*

Source: **hep-ph/0401189v1**, arXiv version dated 26 January 2004, previously retrieved original PDF and its text conversion. Directly cross-checked printed pp. 2–7 and Appendix A, p. 16. A delegated primary-reader also read §§I–III, pp. 1–7, §IV opening, p. 8, and §VI, pp. 13–16, and visually checked PDF Eqs. (26)–(27). Later model spectra, simulations, and Appendices B–C were not reviewed comprehensively.

- **Eigen-equations fix masses and shape ratios; the action fixes physical amplitude.** §III opening, pp. 4–5, and footnote 4, p. 4, explicitly distinguish solving the Schrödinger equation and linear junction conditions from identifying/normalizing the physical mode. The constraint/dynamical system is (9)–(12), p. 4. Canonical normalization in (2)–(4), p. 2, and (27), p. 7, determines brane couplings; it cannot be extracted from the linear equations alone.
- **Boundary derivatives enter the normalization.** §III B, p. 6, explicitly says total derivatives in (y) cannot be dropped because there are two branes. Reduced action (19) retains them. On p. 7, the antisymmetric bulk pairing (22) is a Wronskian boundary contribution; it cancels with the difference of the boundary terms in (19). Footnote 7 states that the junction conditions (12) are used in that calculation. The coefficient (26) first has the form
  
  \[
  C_{mn}=\frac12\int_0^{y_0}dy\,\widetilde v_m\widetilde v_n
  -\left[\frac{3M^3A^4}{4A'}\widetilde\Phi_m\widetilde\Phi_n\right]_{\rm endpoints},
  \]
  
  and then becomes the diagonal positive bilinear expression. The printed upper endpoint of the bracket is (r_0), apparently inconsistent with the surrounding (y_0) notation; the endpoint symbol is left generic here.
- **Physical kinetic coefficient, with PDF-verified measure.** Defining (Y_n=A^2\widetilde\Phi_n), Eq. (27), p. 7, is
  
  \[
  S_n=C_n\int d^4x\,Q_n(\Box-m_n^2)Q_n,
  \qquad
  C_n=\frac{3M^3}{2}\int_0^{y_0}\frac{dy}{A}
  \left[\frac{3M^3}{2\phi'^2}(Y_n')^2+Y_n^2\right],
  \qquad C_n\equiv\frac12.
  \]
  
  The outer measure is **(dy/A), not (dy/A^2)**; plain-text layout can obscure this. Positivity of this coefficient concerns the kinetic term and does not independently prove nonnegative (m_n^2).
- **The earlier Gibbons–Hawking cancellation has a different role.** Appendix A, p. 16, Eq. (A2), explains that the gravitational (y)-total derivative from the ADM decomposition cancels the starting action's Gibbons–Hawking term (5). This cancellation does not permit discarding the later reduced-action boundary terms in (19).
- **The setting is not a de Sitter cone endpoint.** §II, p. 3, Eqs. (5)–(8), specifies static warped Minkowski slices, a finite (S^1/Z_2) interval, two branes, a bulk scalar, and brane potentials. Footnote 3 sets the brane displacements to zero in the absence of fluids or scalar kinetic terms on the branes. General §III normalization permits arbitrary scalar backreaction and brane-potential (U''); large (U'') is an additional assumption only in later examples, p. 4.

For transfer to a cone, denominators involving (A') and (\phi') in (11), (19), (20), and (27) require scrutiny wherever these vanish or limiting coordinate variables degenerate. This is a scope caution inferred from the formulas, not a failure established by the paper. The scalar-free problem has a separate treatment in §III A, p. 5, Eqs. (13)–(17), and is not obtained by substituting (\phi'=0) in (27). The source supplies no cone-extension criterion or covariant phase-space theorem.

## Reading-based checkpoint boundary

The originals support three distinct questions: radial operator domain/orthogonality, action-derived physical kinetic normalization including boundaries, and global regularity/gauge admissibility across the cone. None of these papers proves their equivalence. A local positive kinetic coefficient and a self-adjoint radial problem therefore leave the cone-extension domain, regularity of reconstructed metric/scalar perturbations, allowed global bending diffeomorphisms, and boundary symplectic flux/charges to be checked in the actual background. These notes do not establish that such a completion succeeds or fails.
