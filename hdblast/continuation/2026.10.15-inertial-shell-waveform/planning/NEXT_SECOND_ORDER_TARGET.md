# Next quantitative second-order target for the selected S3 quadrupole

Status: planning candidate for independent review, outside the frozen waveform replay. The existing `HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014` checkpoint already constructively specifies a conditional **complete finite-angular second Taylor coefficient through contact**, including scalar degrees zero, two and four, tensor degree two, physical junctions and constraint propagation. `HDBLAST_ANGULAR_UNIFORM_CONTACT_20261015` adds all-angular frozen-linear control and an exact local nonlinear boundary reduction. These results are retained, not rediscovered here.

Earlier conditional exact initial-data/precontact families retain their stated finite spatial regularity and parameter hypotheses. Any subsequent conditional quasilinear/exact-family result also retains its own gauge, inverse, support, regularity and lifespan assumptions; this planning note neither withdraws nor strengthens it or identifies its second derivative with a particular earlier correction.

The new target is **explicit validated quadratic forcing and a quantitative second-order shell stress/clock observable for the wider `ell=10^-6` state**, followed by a quantitative nonlinear remainder under the applicable exact-family and stability hypotheses. The angular identities below are known coefficient-computation support, with no novelty or new conditional-existence claim. Their receipt adds no result to the frozen waveform replay or theorem register.

## 1. Fix the derivative convention and complete unknown

Keep the registered five-dimensional canonical Einstein–scalar action, potential and common shell tension. Denote the registered bulk potential by `V` here to distinguish it from the configuration `U`. The physical bulk residuals are

\[
\mathscr E_{AB}=R_{AB}-\phi_A\phi_B-\frac23 V(\phi)g_{AB},
\qquad \mathscr S=\Box_g\phi-V'(\phi).
\]

Let `U` contain every required metric/scalar component on both sides, the shell identification and embedding data before material pullback, and prescribed gauge/clock information. Fix these choices and the first-order reconstruction before evaluating a quadratic source. Use the derivative convention

\[
U_\varepsilon=U_0+\varepsilon U_1+\frac{\varepsilon^2}{2}U_2,
\qquad U_j=\left.\partial_\varepsilon^jU_\varepsilon\right|_0.
\]

For the same residual map `F` at all orders,

\[
F(U_\varepsilon)=F(U_0)+\varepsilon DF[U_1]
+\frac{\varepsilon^2}{2}\bigl(DF[U_2]+D^2F[U_1,U_1]\bigr)
+O(\varepsilon^3).
\]

The equations to instantiate quantitatively for the wider state are therefore

\[
\boxed{L U_2=-D^2F[U_1,U_1],\qquad L=DF|_{U_0},}
\]

with `F=(E,S)` for the physical bulk equations. Reduced-gauge equations instead use derivatives of the **same reduced residual**, including its quadratic gauge terms. An extra factor `1/2` on the right-hand side is inconsistent with this expansion. A first-order scalar master determines `U1` only after its full reconstruction and gauge transformation are specified. Its scalar waveform operator alone does not solve these complete second-order equations.

## 2. Exact angular preparation for Y=n1^2-n2^2

Use unit round `S3` inside `R4`, with `sum n_i^2=1`. Set

\[
s=r^2=\sum_{i=1}^4x_i^2,\quad A=x_1^2+x_2^2,
\quad B=x_1^2-x_2^2,\quad P_4=B^2,
\]

\[
Z_2=A-\frac{s}{2},\qquad H_4=P_4-\frac{sA}{2}+\frac{s^2}{12}.
\]

The subscripts indicate homogeneous degree. Direct Euclidean differentiation gives

\[
\Delta_{\mathbb R^4}B=\Delta_{\mathbb R^4}Z_2
=\Delta_{\mathbb R^4}H_4=0.
\]

Indeed, `Delta P4=8A`, `Delta(sA)=16A+4s` and `Delta(s^2)=24s`; these independently check the coefficient `1/12`. The exact ambient decomposition is

\[
\boxed{P_4=H_4+\frac{s}{2}Z_2+\frac{s^2}{6}.}
\]

Restricting to the sphere, and using `Z2,H4` also for their restrictions,

\[
\boxed{Y^2=H_4+\frac12Z_2+\frac16.}
\]

For the Laplace convention with nonpositive spectrum,

\[
\Delta_{S^3}Y=-8Y,\quad \Delta_{S^3}Z_2=-8Z_2,
\quad \Delta_{S^3}H_4=-24H_4,
\]

because homogeneous harmonic degree `l` gives eigenvalue `-l(l+2)` on unit `S3`. Thus

\[
|\nabla_{S^3}Y|^2=\tfrac12\Delta_{S^3}(Y^2)+8Y^2
=\boxed{\frac43+2Z_2-4H_4.}
\]

A separate geometric check projects `grad_R4 B` onto the tangent space:

\[
|\nabla_{S^3}Y|^2=4(n_1^2+n_2^2)-4Y^2.
\]

It yields the same decomposition without the Laplace product formula. With normalized spherical average, `mean(Y^2)=1/6`, `mean(Z2)=mean(H4)=0`, and `mean(|grad Y|^2)=4/3`. The selected `Y` does not have unit L2 norm. For physical angular metric `R^2 kappa_S3`, the contracted angular gradient and Laplacian have the additional factor `R^-2`.

For this single aligned first-order `l=2` mode on a rotationally invariant background, scalar quadratic angular sources have degrees `0,2,4` only. The identities explicitly cover the basic scalar contractions; rotationally invariant scalar derivative contractions preserve that scalar content. This presumes that the fixed first-order reconstruction/gauge introduces no additional independent angular data.

This scalar calculation does not determine the vector and tensor projections of `D_iY D_jY` or `Y D_iD_jY`. The prior full-contact construction uses

\[
\operatorname{Sym}^2(1,1)=(0,0)\oplus(1,1)\oplus(2,2)\oplus(2,0)\oplus(0,2),
\]

with dimensions `1+9+25+5+5=45`: scalar degrees zero, two and four and the two tensor degree-two orientations. It establishes no forced coexact-vector representation for this selected aligned state under its equivariance assumptions. Trace, scalar-derived trace-free and transverse-traceless tensor coefficients still need their actual source projections, kernels and exceptional-mode treatment in the wider-state computation. A more general first state requires its own vector-sector analysis. The scalar identities alone do not close the coupled Einstein–scalar system or replace the existing complete construction.

## 3. Differentiate the physical junction laws

Use one material shell, common induced metric `h`, scalar trace `psi`, and two outward normals. With `p_i=n_i phi_i`, the trace-reversed junction residuals are

\[
\mathcal J_{ab}=\sum_{i=\pm}K_{i,ab}-\frac13\sigma(\psi)h_{ab},
\qquad \mathcal J_\phi=\sum_{i=\pm}p_i+\sigma'(\psi),
\]

alongside induced metric/scalar matching. Outward-normal contributions are summed. In the identical-copy sector these reduce to `K_ab=sigma h_ab/6` and `n phi=-sigma'/2`. The order-two equations obey

\[
\boxed{D\mathcal J[U_2]=-D^2\mathcal J[U_1,U_1],}
\]

plus twice-differentiated matching and gauge-boundary equations. In a family-adapted material Gaussian normal collar, using the same derivative convention, the physical rows are

\[
\boxed{
\sum_iK_{i,2,ab}-\frac{\sigma_0}{3}h_{2,ab}
-\frac{\sigma'_0}{3}\psi_2h_{0,ab}
=\frac{2\sigma'_0}{3}\psi_1h_{1,ab}
+\frac{\sigma''_0}{3}\psi_1^2h_{0,ab},}
\]

\[
\boxed{\sum_i p_{i,2}+\sigma''_0\psi_2=-\sigma'''_0\psi_1^2.}
\]

Here `K_i,2,ab` is the lower-index tensor derivative. Any scalar trace uses the varying induced metric. Material matching reads `h_+,2=h_-,2=h2` and `phi_+,2=phi_-,2=psi2`; graph-evaluation terms have already entered the material coefficients.

For the pure-tension action, material shell stress `S_ab=-sigma(psi)h_ab` has

\[
S_{2,ab}=-\sigma_0h_{2,ab}-\sigma'_0\psi_2h_{0,ab}
-2\sigma'_0\psi_1h_{1,ab}-\sigma''_0\psi_1^2h_{0,ab}.
\]

These sources are fixed by the action and `U1`. There is no independent quadratic boundary forcing available to fit a desired waveform.

## 4. Quantitative implementation and hypothesis audit

1. **Gauge and material pullback.** Fix the first transformation and reconstruction, second gauge conditions and primal boundary completion. For a tensor pulled back by a collar diffeomorphism with first generator `xi1` and generator derivative `xi2`, the material second coefficient is `Q2+2 L_xi1 Q1+(L_xi2+L_xi1^2)Q0`. Apply the same convention on both copies and retain embedding, tilt, normal-normalization and gauge terms. Material Gaussian and bulk hyperbolic gauges require an explicit transition.
2. **Shell clock and stress.** Fix the proper-time/material measurement clock, its reconstruction through order two, and the compared tensor components. Include clock-map and pullback derivatives. A coordinate component at a changing coordinate time is not yet the observable. Use the common tension and induced geometry in the stress rows above.
3. **Initial constraints and data.** Instantiate the inherited second-order Cauchy correction for the wider state and verify all twice-differentiated constraints, `D C[U2]=-D^2 C[U1,U1]`, with complete gauge data. A localized inverse retains its stated background, support, kernel and derivative-loss hypotheses. Fix any free second-order homogeneous data explicitly; `U1` does not determine them.
4. **Boundary constraints and corner compatibility.** Check physical junctions, matches and gauge boundary rows together, retaining the hypotheses of the existing boundary/constraint propagation result. A scalar boundary-condition count does not supply the tensorial boundary system. Check the necessary initial corner jets.
5. **Full projection and support.** Instantiate all required source projections from the existing full finite-angular construction and compute their coupled coefficients with validated actual data. Prove that the localized Cauchy correction and gauge construction retain support and that the source/evolution remain in the controlled collar. A local quadratic source shares the first field's support, but a constraint correction need not. Keep the uncontrolled background region outside the domain of dependence by a justified propagation argument.
6. **Verification and nonlinear remainder.** Reconstruct the quantitative `U2`, substitute it into every bulk, constraint and boundary row, and bound numerical/analytic approximation errors in the stated norms and clock. Retain the existing conditional coefficient construction. To justify `O(epsilon^3)` uniformly, identify a suitable exact family, match its first two parameter derivatives to the chosen coefficients, and prove a quantitative third-derivative or nonlinear defect/stability estimate in the observation norm. Parameter smoothness at finite spatial regularity does not automatically supply that uniform constant or derivative-loss budget. All nonlinear gauge, support, inverse and lifespan hypotheses remain explicit.

The next concrete delivery is a validated wider-state source/projection table and quantitative realization of the already constructed conditional `U2`, including material shell stress and clock response, with an explicit remainder estimate conditional on the required exact-family/stability hypotheses. Linear waveform closure and the identity receipt below do not compute these quantities or remove those hypotheses.

## 5. Separate planning-check scope

`check_next_second_order_identities.py` performs six exact symbolic groups: ambient harmonicity/homogeneity, ambient decomposition, unit-S3 eigenvalues, tangential-gradient projection, the independent Laplace-product route, and normalized S3 moments. It rejects deliberately wrong harmonic constants, an omitted degree-two term, S2 eigenvalues, a wrong gradient sign and a wrong Taylor-factor convention.

`SECOND_ORDER_HARMONIC_IDENTITY_RECEIPT.json` records only these planning checks, exact residuals, checker/note hashes and runtime versions. The checker is planning code outside the frozen waveform capsule and is not counted in its replay roster. Passing it certifies angular algebra and convention controls; it does not verify quadratic bulk/junction sources, solve a second-order equation, or prove evolution or a nonlinear remainder.
