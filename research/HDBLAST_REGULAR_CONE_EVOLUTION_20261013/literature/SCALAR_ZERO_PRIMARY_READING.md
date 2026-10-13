# Scalar-master variables at background zeros: bounded original reading

Date: 13 October 2026. Scope: four established primary originals and two recent primary preprints, with selected sections actually read from their cached PDFs/text. The companion `PUBLIC_PRIMARY_READING_REGISTER.json` identifies the PDF/text bytes and reading scopes. No new network request was made for this subtask. Page numbers below are **PDF page numbers**; the two recent manuscripts' printed page numbers are one less than these PDF page numbers.

The useful precedent is **change the divided variable and inspect the original constraints when a reduced equation looks singular**. None of these six sources proves that the current `rho ~ y`, `phi' ~ kappa y` de Sitter cone has a regular ambient master, or that its reconstruction yields cone-crossing Einstein–scalar solutions. An interior scalar turning point, a stationary positive-warp endpoint, and a vanishing-warp coordinate cone are different problems.

## 1. Giovannini: complementary master variables need different regularity checks

Massimo Giovannini, *Gauge-invariant fluctuations of scalar branes*, [hep-th/0106041v1](https://arxiv.org/abs/hep-th/0106041v1), 5 June 2001. Read pp. 4–5, 12–18, 20–23, 26–27.

The background is the conformal, flat-sliced metric `a(w)^2 eta_5`, with `H=a'/a`; Eqs. (2.4)–(2.10), p. 4. The explicit regular thick wall has `a=(1+b^2 w^2)^(-1/2)` and `phi=sqrt(6) arctan(bw)`, Eqs. (2.11)–(2.13), p. 5. At its center `a` is positive, `H=0`, and `phi'` is nonzero. It is not a cone with vanishing warp factor.

The undivided scalar constraint is Eq. (4.9), p. 12; solving it for the scalar perturbation produces Eq. (4.11), which divides by `phi'`. The Bardeen master `Phi=a^(3/2) Psi/phi'` and background function `z=a^(3/2) phi'/H` occur in Eqs. (4.13)–(4.15), p. 13. The reciprocal potential is `z(1/z)''`. The canonical combination `G=a^(3/2)X-z Psi`, Eqs. (4.17), (4.22), pp. 13–14, instead obeys a potential `z''/z`, Eq. (4.18). Gauge invariance is verified through Eq. (4.21); it does not remove the division by `H` in `G` or the division by `phi'` in `Phi`.

The scalar bulk quadratic action is Eq. (4.41), p. 17. The paper obtains it after eliminating total derivatives by using the connection form of the Einstein–Hilbert action, pp. 16–17. This is a bulk thick-wall calculation; it is not a justification for dropping endpoint terms in the present shell problem. See Kofman–Martin–Peloso below.

The zero-mode normalization discussions make the domain issue visible. Equations (5.14)–(5.16), p. 20, require global convergence for the reciprocal variable. Equations (5.20)–(5.21), p. 21, show the canonical zero-mode integrand diverging at `w=0` for the explicit geometrically regular thick wall. The partner factorizations, Eq. (5.25), p. 22, are therefore structural identities; importing them as a stability or completeness theorem still requires an operator domain and boundary terms. Appendix B, Eqs. (B.3)–(B.6), pp. 26–27, explicitly demonstrates cancellations using the background equations even when an intermediate representation contains `1/H`.

**Transfer:** constrained reduction, alternative variables, background-identity cancellation, and the separation of geometry regularity from master-variable regularity. **Do not transfer:** flat-slice potentials, normalization measures, zero-mode nonlocalization, or a statement that gauge invariance alone makes a variable smooth at all background zeros.

## 2. Lesgourgues–Sorbo: scalar turning points and the stationary endpoint

Julien Lesgourgues and Lorenzo Sorbo, *Goldberger-Wise variations: stabilizing brane models with a bulk scalar*, [hep-th/0310007v2](https://arxiv.org/abs/hep-th/0310007v2), 21 June 2004. Read pp. 1–6 and p. 8 notes [15]–[16]; p. 3 visually checked.

The setup is compact and Minkowski-sliced, with branes at orbifold fixed points; Eq. (1), p. 1. Equations (7)–(9), p. 2, contain the original scalar constraints and dynamical equation. Eliminating the scalar perturbation gives Eq. (10), with coefficients containing `phi''/phi'`. The divided variable is `u=a^(3/2) Psi/phi'`, Eq. (11).

The alternative variable is explicitly

`v=a^(3/2) [Delta - (phi'/H) Psi]`, Eq. (12), p. 2.

Its potential `z''/z`, expanded in Eq. (14), is regular at `phi'=0` **as long as `H != 0`**, a condition stated directly beneath that equation. The reconstruction identity is Eq. (15); it determines `Box Psi`, not an unconditional inverse map for all massless or exceptional modes. Equation (16) supplies eigenvalue-dependent brane conditions. These are finite radial interval boundary conditions, not an ambient Cauchy theorem.

The one-brane construction has a second stationary, positive-warp fixed point: `phi'(0)=H(0)=0`, Eq. (18). Near it, `phi-phi(0)=O(y^2)` and `H=O(y^3)`, Eq. (19), p. 3. The leading alternative-variable equation contains `6/y^2`, Eq. (25), rather than a regular cone operator. Footnote [15], p. 8, explicitly distinguishes a locally admissible `Psi~y^3` branch from the even orbifold branch selected by their geometry. Thus a locally singular divided ODE alone does not choose the physical domain.

The two-brane stability criterion, Eq. (44), p. 6, assumes a monotonic positive warp factor and endpoint signs `g_- > 0`, `g_+ < 0` after choosing `H<0`. Its discussion of an interior zero of `phi'` and its one-brane instability conclusion use this compact flat-sliced setup. They do not imply instability at the present cone, where the warp factor itself vanishes and the slice is de Sitter.

**Transfer:** regularity at a scalar turning point is conditional on the chosen variable and background; exceptional modes must be handled when reconstruction contains a slice wave operator; a physical reflection condition can select among local ODE branches. **Do not transfer:** the stationary endpoint expansion, compact one-brane tachyon result, or the finite-interval sign criterion to the current cone.

## 3. Konikowska–Olechowski–Schmidt: explicit warning about divided Sturm–Liouville coefficients

D. Konikowska, M. Olechowski and M. G. Schmidt, *Radion stabilization with(out) Gauss–Bonnet interactions and inflation*, [hep-th/0603014v2](https://arxiv.org/abs/hep-th/0603014v2), 31 May 2006. Read p. 4 and pp. 6–14; p. 10 visually checked.

The background is conformal and has de Sitter or flat slices between two branes, Eqs. (3)–(8), p. 4. The undivided perturbation system and scalar boundary condition are Eqs. (14)–(17), p. 6. The reduced equation contains the shifted eigenvalue `lambda=m^2+4H^2`; the eigenparameter boundary conditions are Eqs. (23), p. 7. Equations (24)–(25), p. 8, are the Sturm–Liouville form with `p=3/(2a phi'^2)`. Equations (27)–(29), p. 9, retain endpoint terms in the radial norm.

At p. 10 the authors explicitly warn that Eqs. (20), (24), (32) have coefficients diverging at interior `phi'=0`. They introduce an alternative Mukhanov variable and expand its potential in Eq. (35), whose coefficients are regular at a scalar turning point when the other denominators are controlled. **Footnote 3 on this page restricts the simple Eqs. (32), (34) to `H=0`**; those equations must not be quoted as the paper's de Sitter regular-master theorem.

There is a source-fidelity defect that should be carried openly. Equation (33) in this v2 PDF prints `v=a^(3/2)[F3-(a/a')F1]`, without a `phi'` multiplying `a/a'`. This was checked in a rendered PDF page, so it is not a text extraction error. Under the standard radial gauge transformations, the missing factor prevents the expected cancellation; compare the explicit gauge-invariant combinations in Lesgourgues–Sorbo Eq. (12), Giovannini Eq. (4.22), and Kofman–Martin–Peloso Eq. (20). Treat this as an **apparent printed-variable defect**, not a formula to copy or silently repair. The warning about the divided equation and the displayed potential cancellation remain useful evidence; a present derivation should check its own variable and identities independently.

Their regularity/stability discussion is carefully bounded in the original. Page 11 states that positivity of the Rayleigh expression is justified only when the coefficient functions are finite, and that a singular Sturm–Liouville system needs more care. The necessary-and-sufficient criterion (44), pp. 12–14, is in the `H=0` part of the analysis and includes two brane endpoint signs and the graviton kinetic-sign condition. The interior scalar-zero tachyon argument, Eqs. (45)–(47), p. 13, uses a point strictly between those branes. Page 14 separately explains that `a'=0` does **not** make `Q` singular: Eq. (43) cancels the apparent denominator in `q` provided `phi'` and `xi` are nonzero and `a` stays positive. A zero of `a'` is not a zero of `a`.

**Transfer:** retain endpoint contributions; inspect the undivided system; use background identities before declaring a coefficient singular; distinguish `a'=0` from `a=0`; check the literal printed variable. **Do not transfer:** flat two-brane stability/tachyon conclusions, or the simple Mukhanov equation at nonzero `H`, or an unqualified statement that simultaneous scalar/warp zeros are covered.

## 4. Kofman–Martin–Peloso: bulk master is not the physical normalization or endpoint action

Lev Kofman, Johannes Martin and Marco Peloso, *Exact identification of the radion and its coupling to the observable sector*, [hep-ph/0401189v1](https://arxiv.org/abs/hep-ph/0401189v1), 26 January 2004. Read pp. 3–7 and p. 16; p. 7 visually checked.

The action includes bulk, two-brane potentials and extrinsic-curvature terms, Eq. (5), p. 3. The background is static with flat branes, Eq. (6). Footnote 3 excludes fluids and scalar kinetic terms on the branes when eliminating brane displacements; it does not cover the later dynamical kinetic shell boundary without additional work.

The master `U=A^(3/2)Psi/phi'` and `z=sqrt(2) A^(5/2)phi'/A'` are Eq. (11), p. 4. The authors state here that the bulk differential equation plus homogeneous boundary conditions determines the mass and relative wave-function shape, but **does not identify the physical normalization**. The scalar-free sector is treated separately on p. 5, Eqs. (13)–(17); setting `phi'=0` in their divided master is not that treatment.

The quadratic action, Eq. (19), p. 6, contains a radial total derivative. The canonical variable is Eq. (20), and its bulk equation is Eq. (21), p. 7. The boundary term cancels the antisymmetric bulk Lagrange form in Eq. (22). Equations (26)–(27) then give the action-derived normalization, including endpoint terms before rewriting it in the physical fields. The outer measure in Eq. (27) is **`dy/A`**, visually checked; it must not be mistranscribed as `dy/A^2`. Appendix A, p. 16, Eq. (A2), identifies the Einstein–Hilbert radial total derivative and its cancellation against the original Gibbons–Hawking term.

**Transfer:** derive the norm/kinetic coefficient from the full quadratic action and keep radial boundary terms; separate scalar-free/degenerate branches. **Do not transfer:** flat-brane measures or canonical normalization coefficients to the present de Sitter master; the stated displacement elimination to a boundary with additional kinetic dynamics; a regular ambient reconstruction theorem.

## 5. Recent originals: explicit reduced evolution, different geometry and domain

Weike Deng, Sheng Long, Qin Tan, Zu-Cheng Chen and Jiliang Jing, *Scalar-gravitational quasinormal modes and echoes in a five dimensional thick brane*, [2508.20937v1](https://arxiv.org/abs/2508.20937v1), 28 August 2025. Read PDF pp. 4–12, concentrating on §2 and the boundary/evolution statements in §3; PDF pp. 4 and 7 visually checked. This is the **v1 preprint** actually read, not a claim to have read the later JHEP published text.

The model has static Poincare slices, Eq. (2.4), and conformal radial coordinates, Eqs. (2.12)–(2.13). Its reconstruction divides by `partial_z phi`, Eq. (2.20), PDF p. 7. The rescaled master, Eq. (2.22), has the same scalar-divided structure as the older reciprocal Bardeen variables; Eqs. (2.23)–(2.24) give a `1+1` wave equation and a potential containing `phi'''/phi'` and `(phi''/phi')^2`. QNM conditions, Eq. (3.2), PDF p. 9, are outgoing conditions at radial infinities. The null coordinates of Eq. (3.6), PDF p. 12, are characteristics of this reduced `1+1` wave equation. They do not constitute a smooth ambient Einstein–scalar cone construction.

Do not treat its formulas or split-brane numerics as an independently certified regularity benchmark. In the v1 PDF, Eq. (2.9) literally prints `phi = constant * arctan(ky/delta)`, without an `s` power inside the arctangent, while Eq. (2.8) and the adjacent discussion/figures describe an `s>1` double-kink/split-brane family. This is a PDF-level source inconsistency, not a text extraction issue. The printed linearized equations also need checking before reuse: for example, Eq. (2.15) contains `varphi (partial_z delta phi)^2`, and Eq. (2.18) prints a first potential derivative multiplying `delta phi`. The quoted master structure is familiar, but this selective read does not repair or validate these source equations or reproduce the QNM calculations.

Gahan Chattopadhyay and Soumitra Sengupta, *Relic of the Kink Shape Mode in a Topological Thick Brane - some implications*, [2609.37346v1](https://arxiv.org/abs/2609.37346v1), 29 September 2026. Read PDF pp. 8–11, p. 13 reconstruction paragraph, and pp. 15–16 §§6.1–6.3.

Here `phi=v tanh(ky)`, Eq. (3.6), and the warp factor is Eq. (3.10), PDF p. 8: positive at every finite radial point, flat-sliced and centered at `A(0)=0`. The scalar derivative is nonzero at every finite point for the nondegenerate kink. Thus the divided reconstruction, Eqs. (4.14)–(4.17), PDF p. 10, and master potential, Eqs. (4.18)–(4.21), PDF p. 11, do not confront the present simultaneous vanishing-warp/scalar-gradient cone. The selected resonance profiles are reconstructed via Eqs. (4.17), (4.15), (4.14), as explicitly stated on PDF p. 13. This is useful precedent for checking the scalar field together with the metric rather than treating the master as the physical field itself.

The time evolution, Eq. (6.1), PDF p. 15, restricts to zero brane spatial momentum and evolves the scalar master as a `1+1` numerical wave system. Equations (6.6)–(6.10), PDF p. 16, add a sponge, numerical dissipation and tapered scattering-profile initial data. The original explicitly calls the sponge an approximate absorbing treatment. These calculations study a resonance/leakage problem; they do not prove a full Einstein–scalar timelike boundary theorem, general constrained Cauchy existence, or regular ambient cone crossing.

**Assessment:** the two recent readings add explicit reduced evolution/reconstruction examples and current resonance context, but neither supplies a theorem or geometry that materially closes the current ambient-cone obligation. Their presence does not support a priority or novelty declaration.

## Actionable checks for the current cone-evolution calculation

The proposed current variable `v=P/(rho^3 p)`, where `P=rho^2 Psi` and `p=phi'`, is a new model-specific reduction. Its name `v` does not identify it with any of the Mukhanov variables above. The readings support its motivation, not its correctness.

1. Derive its exact PDE from the current **undivided** Einstein–scalar constraints and background equations. Exhibit smooth ambient coefficients rather than only a cancellation in a separated radial potential. Check the principal symbol and the actual potential through the cone.
2. Reconstruct the fields locally, beginning with the identities `Psi=rho p v` and `X=-3[(3rho'+rho p'/p)v+rho partial_y v]` on the punctured patch. Verify in ambient coordinates that the apparent `p'/p` and radial coordinate tensor factors cancel using the even analytic background. The limiting leading expression for `rho~y`, `p~kappa y` is `X~-12v-3y partial_y v`; this is a structural check, not an existence proof.
3. Verify the reconstructed metric, scalar and all undivided field-equation residuals across the cone. A smooth master solution is useful only when this map produces a smooth tensor/scalar solution and the constraints hold. Do not infer this from a radial norm or gauge invariance.
4. State a one-way construction theorem first if that is what can be proved: admissible smooth ambient Cauchy data for the regular master yield constrained Einstein–scalar perturbations in a local cone-crossing patch. A bijection onto the full scalar physical phase space, including exceptional slice modes and gauge kernels, requires a separate analysis.
5. Keep the local cone-crossing result separate from the doubled-shell boundary problem. The positive radial endpoint norm is evidence about the reduced spectral problem; the complete shell kinetic sign and a constraint-preserving timelike Einstein–scalar IBVP remain separate obligations.

This reading makes no exhaustive-prior-art, novelty, global-stability, or completed physical-domain claim.
