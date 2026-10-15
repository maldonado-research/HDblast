# Harmonic initialization: primary comparison and scope

15 October 2026 UTC. This bounded review addresses the planned quantitative second-initial spatial frame-C1 bounds and local relative-harmonic initial completion. It used five targeted queries and selected passages from five primary papers, after reading the prior shell/core review, initial-constraint comparison and saved local gauge derivation. No first gauge wave, evolution, shell IBVP, exhaustive search, priority or publication-readiness claim is made.

**The literature directly supports the initialization architecture.** Lapse/shift values can be selected as initial coordinate data; the harmonic condition then determines their time derivatives. For a solution of the reduced equations, the physical Cauchy constraints and a vanishing initial harmonic constraint imply its vanishing initial normal derivative. For this checkpoint those conclusions must be applied through amplitude order two, with the complete quadratic terms. No cited paper supplies the actual numerical frame bounds or an evolution estimate for the five-dimensional scalar shell.

## Five primary sources actually inspected

**Explicit lapse/shift equations and normal gauge constraint.** Lee Lindblom, Mark A. Scheel, Lawrence E. Kidder, Robert Owen and Oliver Rinne, [A New Generalized Harmonic Evolution System, gr-qc/0512093v3](https://arxiv.org/abs/gr-qc/0512093v3), revised 4 December 2007; CQG 23 (2006), S447–S462. Primary §2, equations (9)–(10), expresses the time derivatives of lapse and shift in terms of the specified algebraic gauge source. Equations (11)–(16) derive the subsidiary equation and relate its normal derivative to the Hamiltonian/momentum constraints and tangential derivatives. The displayed Einstein system is four-dimensional vacuum; its boundary discussion explicitly distinguishes initial constraint control from incoming boundary constraint violations. Its generated HTML title-page date is not a revision date. The dimension/matter extension used below is an explicit local algebraic inference, not an imported shell theorem.

**Matter conservation and a useful compact-support caution.** Frans Pretorius, [Numerical Relativity Using a Generalized Harmonic Decomposition, gr-qc/0407110v2](https://arxiv.org/abs/gr-qc/0407110v2), revised 23 December 2004. Primary §II, equations (11)–(14), derives the homogeneous gauge-constraint wave equation using Bianchi and stress conservation and identifies the usual initial Einstein constraints. Primary §§IV.4–IV.5 treats a massless scalar and time-symmetric, conformally flat initial data. Compact scalar support does not imply compact metric correction: the discussed conformal factor has an asymptotic tail. This numerical four-dimensional construction neither certifies the current potential nor proves the finite scalar-shell IBVP; the paper expressly limits its general gauge-evolution well-posedness discussion.

**The scalar Cauchy system is genuinely coupled.** Yvonne Choquet-Bruhat, James Isenberg and Daniel Pollack, [Applications of theorems of Jean Leray to the Einstein-scalar field equations, gr-qc/0611009v1](https://arxiv.org/abs/gr-qc/0611009v1), submitted 2 November 2006. Primary §§2–3 displays the scalar stress tensor, scalar equation and n+1 Cauchy decomposition, including unit-normal scalar velocity. Theorem 4.1 concerns a different causal hyperbolic formulation, not this relative-harmonic shell reduction. The scalar equation implies stress conservation; that implication is the one needed here. Conservation alone should not be used to infer the scalar equation where the scalar gradient vanishes. No theorem from its alternate evolution gauge is imported.

**Physical scalar initial variables and constraint freedom.** The same authors, [The constraint equations for the Einstein-scalar field system on compact manifolds, gr-qc/0610045v2](https://arxiv.org/abs/gr-qc/0610045v2), revised 11 December 2006. Primary §2, equations (4)–(13), gives the physical constraints and conformal split; §4.1 and Theorem 2 retain CMC, sign-class and conformal-Killing compatibility conditions. “Compact manifold” does not mean a compactly supported perturbation. This paper supports the physical-variable conventions after translating the extrinsic-curvature sign. It does not establish the capsule's exact support, its explicit second correction or its frame-C1 constants.

**Localized inverses retain hypotheses.** Erwann Delay, [Smooth compactly supported solutions of some underdetermined elliptic PDE, with gluing applications, 1003.0535v4](https://arxiv.org/abs/1003.0535v4), revised 17 November 2011. Primary Definitions 1.1–1.2 and Theorem 1.3 require the stated asymptotic Poincaré/kernel regularity conditions and orthogonality to the formal-adjoint kernel. They construct a smooth supported solution in a chosen region. Theorem 1.5 adds flux compatibility for gluing. This remains a comparison for localization, not a verified Einstein–scalar right inverse here. The registered explicit formulas themselves establish their support and source compatibility; the general paper supplies no numerical norm bound for them.

## Precise use at the second-initial checkpoint

Use the future unit normal u to the Cauchy surface, distinct from the outward spacelike normal to the timelike shell. With the retained convention K=½ Lᵤγ, write

\[
g=-\alpha^2dt^2+\gamma_{ij}(dx^i+\beta^idt)(dx^j+\beta^jdt),\quad
D_t=\partial_t-\beta^i\partial_i,
\]
\[
C^A=g^{BC}(\Gamma[g]^A_{BC}-\Gamma[g_0]^A_{BC}),\qquad
R^A=g^{BC}\Gamma[g_0]^A_{BC}.
\]

The fixed reference connection is contracted with the **perturbed** inverse metric. Rewriting the contracted-Christoffel identity gives the local formulas already derived in the saved gauge note:

\[
D_t\alpha=\alpha^2\operatorname{tr}_\gamma K-\alpha^3R^0,
\quad
D_t\beta^i=\alpha^2(\Gamma_\gamma^i-R^i-\beta^iR^0)
-\alpha\gamma^{ij}\partial_j\alpha.
\]

These coordinate identities work in four spatial dimensions with its corresponding trace. They determine the five remaining initial gauge velocities after the physical data and initial lapse/shift values are chosen. They are local; they preserve an open background collar when all perturbations vanish there. This is initialization, not the integration of a gauge-vector wave. Preserve the selected full first spacetime tangent and include all quadratic conversion terms: zero second lapse variation need not give zero second g_tt variation.

For the model residual E_AB=Ric_AB−φ_Aφ_B−(2/3)Ug_AB and reduced equation E_AB=∇_(A C_B), set H_AB=E_AB−½g_AB tr(E). The physical constraints are H(u,u)=H(u,e_i)=0. If C vanishes initially, its tangential derivatives vanish. Let X_A=(∇ᵤC)_A there. Direct contraction of the reduced equations gives

\[
H(u,u)=X_u/2,\qquad H(u,e_i)=X_i/2.
\]

Hence the physical constraints imply the normal gauge derivative vanishes, **conditional on the reduced equations at initial time**. It is not an extra independently prescribed second metric velocity. In a second-jet calculation, apply this identity coefficientwise only through the verified amplitude order, with the lower-order identities and quadratic sources retained. A second jet satisfying the constraints does not establish an exact all-amplitude constraint family.

Bulk propagation additionally uses the actual scalar equation: the divergence of the scalar stress tensor is the scalar-equation residual times dφ. Inserting that equation and the correctly reduced Einstein equations produces the homogeneous subsidiary system. Initial gauge completion alone supplies no boundary data for its incoming fields. Vanishing on an open initial shell collar aids initial corner compatibility but does not establish later timelike-shell constraint preservation.

## Quantitative boundary of the claim

Full spatial frame-C1 bounds must include angular derivatives, radial derivatives, frame connection terms and the declared weights. Published initialization formulas contain no constants for this packet; those must come from the authenticated background and packet calculations. Distinguish values, physical velocities and gauge velocities in the bound table.

There is also an extra derivative to account for: the gauge-velocity formulas contain first spatial derivatives of the metric/lapse. A spatial C1 bound for those **resulting gauge velocities** generally needs second spatial input derivatives, as well as first derivatives of K. A C1 bound for physical initial values alone does not supply this. If that stronger bound has not been enclosed, state the local gauge completion algebraically and keep its quantitative norm claim separate.

Neither finite frame-C1 bounds nor point-jet identities provide the high Sobolev estimates, lifespan, differentiable solution map, shell derivative trace estimate or validated evolution response required later. The next checkpoint can defensibly report quantitative initial norms and a verified formal harmonic initialization through order two, with the conditional normal-gauge identity stated explicitly.

The [ledger](HARMONIC_INITIALIZATION_SOURCE_LEDGER.json) records versions, access depth and inherited-input hashes. Only bounded paraphrases and mathematical comparisons are saved; no raw primary body, private data, remote mutation or frozen-file edit is included.
