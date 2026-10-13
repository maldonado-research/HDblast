# Independent mathematical review of the integrated cone-evolution draft

Date: 13 October 2026. Internal mathematical review by a separate agent process; this is not external peer review or a novelty certification.

## Reviewed artifacts and conclusion

Exact reviewed bytes were saved beside this report:

- `manuscript-authoring/REGULAR_CONE_EVOLUTION.md`: SHA-256 `8b0305dffde1c29d0b9a0f8b8ce6a95d8237b2ea879ebce0f7ddaac1dfa643ff`.
- `root-audit/check_canonical_jets.py`: SHA-256 `cc277580673b1bfc5df99d6c496122711373783f01d8b3f87ea20f148a15685b`.
- Supporting `physical-reconstruction/check_ambient_reconstruction.py`: SHA-256 `7ca7ff8b6af686d60ceed05a88d5bd4cd9d6d7f42147d4fa34fcc9611145779a`.

No substantive mathematical error was found in equations (8)–(19), the local sufficiency theorem, the energy estimate, or the explicit support construction in this snapshot. Within the stated actual-analytic-background assumptions, a rigorous **local theorem** is warranted: there exist smooth compact-data solutions of the full linearized coupled equations, including the remote homogeneous shell junctions before contact, with nonzero finite-segment covariant null symplectic flux and nonzero local gauge-invariant witnesses. This is sufficient mathematical substance to state a conditional local theorem in a research publication. It does not establish editorial acceptability, independent human validation, novelty, global evolution, or a completed physical phase space.

Two short wording clarifications below are recommended before final release. They do not require changing the mathematical result.

## 1. Full residual identities and signs

I rederived the master transformation and reconstruction independently in the earlier note in this directory. For the integrated off-shell residuals, an independent consistency chain is particularly useful. Write `D=Box_gamma+4` and

\[
R_\chi=\chi'-(p'/p)\chi-[3\rho^{-2}D/p-2p]\Psi.
\]

The momentum reconstruction `P'=-rho²p chi/3` gives, for arbitrary `P`,

\[
R_\chi=-{3M[P]\over\rho^2p}.
\]

The untraced ADM Hamiltonian residual is `−p R_chi=3M/rho²`. The manuscript's trace-reversed residuals have trace `−2M/rho²`; they therefore reproduce precisely that ADM constraint:

\[
\delta\mathscr E_{yy}-{1\over2}g_{yy}
 g^{AB}\delta\mathscr E_{AB}=3M/\rho^2.
\]

Direct differentiation, with the background equations and the exact momentum reconstruction, gives

\[
\delta\mathscr S=(\partial_y+p'/p+4a)R_\chi
=-{3\over\rho^2p}(\partial_y+2a)M[P].
\]

Thus both signs and the `2a` coefficient in (10) are correct. Substituting `M=Omega⁻²r⁴ L F` gives

\[
(\partial_y+2a)(\Omega^{-2}r^4 Q)
=\Omega^{-3}r^3(E+6)Q,
\]

which independently verifies the coefficient `−3/(2Omega⁴ f)` and the operator `E+6` in (11). The Einstein tensorial transformation is also correct:

\[
{M\over\rho^2}(3dy^2-g)
=\Omega^{-2}(3X_AX_B-s\eta_{AB})\mathcal LF.
\]

There is no lapse, trace-reversal, or radial-to-Cartesian sign mismatch. The field equations vanish because `L F=0`, even at the vertex where the Einstein prefactor itself vanishes. The scalar residual is also zero there; no extra vertex equation is being discarded.

The direct-index checker computes the connection variation and all 15 independent Ricci components for arbitrary `Psi(y,t,x1,x2,x3)`, rather than only a separated harmonic. Its subsequent algebra connects these Ricci identities to the off-shell residual factors. Together with the independent ADM derivation this supports the advertised factorization; it is not a completeness proof for the reconstruction.

## 2. Analytic continuation, smooth solutions, and support

The extension argument is valid in the stated analytic setting. The relevant claims are identities of differential-operator coefficients in Cartesian coordinates, whose coefficients are analytic functions of `s` multiplied by Cartesian polynomial tensors. Analytic continuation is applied to those coefficients, not to an arbitrary smooth solution. This distinction correctly allows smooth nonanalytic compact-data solutions through the cone. The actual background residuals extend by the same analytic identity principle. Nonvanishing `f(0)=kappa/2` gives a neighborhood where all reconstruction denominators remain regular.

After this local coefficient extension the Cauchy problem has exactly Minkowski principal symbol. Smooth extension of the lower-order coefficients outside a compact causal neighborhood is legitimate, and domain of dependence makes the physical restriction independent of that auxiliary extension. Neither analyticity of the perturbation nor a characteristic initial-value theorem is needed.

The snapshot's explicit support criterion is sound. Initial support `t0/2<R<3t0/2` and `−t0<T<3t0` imply

\[
R\le7t_0/2,\quad |T|<3t_0,\quad
|s|=|R^2-T^2|\le\max(R^2,T^2)<13t_0^2.
\]

The strict inequalities `13t0²<sigma` and `4t0<r_b` therefore keep the support inside the actual analytic coefficient neighborhood and separated from the shell. The slab contains the vertex. The annular initial harmonic is smooth in Cartesian coordinates; later support reaching `R=0` causes no singularity because the evolution is carried out in Cartesian coordinates. Finite propagation applies equally forward and backward from `T=t0`.

Every local derivative in the reconstruction has the same support inclusion. Thus both perturbations and all derivatives entering any local linearized junction condition vanish on a common open neighborhood of the shell. Zero shell displacement is compatible with either reflection parity. This proves the finite-time shell statement without assuming an unproved post-contact Israel-shell evolution theorem.

Minor clarification C2 below concerns only the more general “usable criterion”; the explicit `sigma` criterion already removes the potential ambiguity.

## 3. Canonical potential and energy bound

The positive scalar weight is correctly `A=1/(Omega f²)`. Its logarithmic derivative is `−L`, and `psi=sqrt(A)F` gives

\[
(\Box_\eta-V_c)\psi=0,\qquad
V_c=3L-2sL_s+sL^2+\tfrac83sf^2.
\]

For the signature in the draft, `psi_TT=Delta psi−V_c psi`. With the displayed unit-mass auxiliary energy, compact support permits integration by parts without spatial boundary terms and gives exactly

\[
\mathcal E'(T)=\int (1-V_c)\psi\psi_T\,d^4X.
\]

Consequently `|E'|≤sup|1−V_c| E`, by `2|psi psi_T|≤psi²+psi_T²`. No `V_{c,T}` term occurs because the chosen energy uses the fixed coefficient one in front of `psi²`, not the time-dependent potential. The exponent `K|T−t0|` and the absence of an extra factor two are correct. The supremum need only cover the compact causal support over the time interval. It need not be a bound on the physical coefficients at arbitrarily large spatial radius.

This is a regular master control norm, not the covariant gravitational norm or a positive physical Hamiltonian. The manuscript expressly preserves that distinction. Smooth compact reconstructed metric and scalar Cauchy data have finite local Sobolev norms, with sufficient derivatives of `F`; the statement does not require a claimed equality of the corresponding energy spaces.

The root checker also correctly obtains the regular background coefficients `k=U1/5`, `rho3=−U/36`, `phi4=k(U2/56+U/126)`, giving `L0=U2/7−5U/56`. Its conformal-radius correction is necessary and correctly included.

## 4. Null current, gauge witness, and normalization

The exact initial-cone jets in (16) follow from the master equation and the local constancy of the radial cutoff. They give

\[
W=\chi_1\chi_{2v}-\chi_2\chi_{1v}
=C^2ab(10+\lambda/2)H_l^2.
\]

For normalized `l=2`, integration over `S³` gives `14C²ab`. The doubled per-`dv` density at `v0=2t0` is indeed

\[
{(2t_0)^3\over4}\,14C^2ab=28t_0^3C^2ab.
\]

The per-copy factor `v³/8` is consistent with the inherited covariant current and the exact cone coordinates. Smoothness of the reconstructed fields gives the differentiated asymptotics required by that earlier formula. No scalar canonical norm was substituted for the Einstein–scalar current. Continuity makes the integral strictly nonzero on a short interval; positivity here pertains only to this ordered antisymmetric pair and orientation.

The relational scalar witness is valid. Under a background Lie derivative, `delta K_g` and `(K_g)_s chi/f` change identically. On the cone `h^{AB}phi_,A phi_,B=0`, so

\[
Z=4f(0)(v\chi_v-\chi).
\]

For `l=2` and `|t0²L0|<1`, the bracket in `Z1` is strictly greater than four, and the bracket in `Z2` is strictly less than minus two. Since `f0 C=3/2`, `ab>0`, and `H_l` is a nonzero function, each witness is nonzero. The solutions therefore cannot be pure infinitesimal diffeomorphisms. Zeros of a spherical harmonic at individual angles do not undermine this functional statement.

The `l≥2` choice also avoids the constant and Hessian-degenerate de Sitter subspaces, whose spatial angular content is `l=0,1`. This exclusion is consistent with rotational invariance of the master equation. It is not being used to invert a whole de Sitter mass eigenvalue.

## 5. Recommended wording corrections

**C1 — specify common parity when doubling the flux.** In Section 4 replace “Identical or opposite copying gives twice (14)” with: “For two solutions with the same reflection parity, doubling gives twice (14); mixed-parity bulk pairings cancel.” The intended construction already takes common parity, so this is a clarification of scope rather than a change in its normalization.

**C2 — state the domain for the general support/energy criteria.** The general “usable criterion” combines negative-`s` continuation with the original exterior branch. Say explicitly that the usable background region is the union of `−epsilon<s≤0` and the already available `0<s≤r_b²`. Alternatively add `(R0+Delta)²<s_+` when using only the displayed local interval `−epsilon<s<s_+`. The later `|s|<sigma` construction already supplies a complete sufficient criterion and is correct. In the energy paragraph, “on the compact causal support over the slab” is more precise than “on a compact time slab.”

No modification of equations (8)–(19) or either proposition is required by this review. Describing the scalar branch as “nonconstant-scalar” would also be more precise than “nonstationary,” but that terminology is not a mathematical obstruction.

## 6. Executed checks and portable adaptation

`INTEGRATED_REVIEW_CHECKS.json` records the reviewed hashes and actual commands. The root canonical checker passes 18 positive identities and five wrong-formula controls in both normal and optimized Python. Normal and optimized JSON output is byte-identical. Its injected-failure mode exits one in both modes. These are executed checks, not inferred success from source inspection.

The reviewed direct-index reconstruction checker was additionally executed in normal and optimized mode: both exited zero, with 32 positive checks and 6 wrong-formula controls; their JSON output is byte-identical. The result hashes are recorded in `MATHEMATICAL_REVIEW.json`.

The independent master checker in this directory previously passed 22 checks in both modes and rejected three actual coefficient mutations in both modes. That code uses no `assert`, so optimization cannot remove its validations.

A later portable adaptation may remove the `optimized` field from that checker's output and add output-path/CLI handling without changing the mathematical checks. Preserve the source hash and original receipts, record the precise adaptation, rerun both modes and every mutation, and put runtime-mode evidence in the runner receipt rather than pretending it remains in deterministic algebra output. Such an adaptation is a packaging change to the same independent derivation, not another independent mathematical verification. No portable adaptation has been made as part of this review.

This review does not recertify the background branch, the cited paper metadata, or the earlier action-derived symplectic-current derivation. It checks the local continuation theorem and algebra against those clearly stated inputs.
