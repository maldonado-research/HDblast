# Independent linear-mode control reference

No new shooting calculation was run. Values and equations below were read from the completed Chat 9 and Chat 13 sources. The algebraic conversion from mass to time growth was independently recalculated with the standard homogeneous de Sitter scalar equation.

For unit-curvature de Sitter time T, a homogeneous harmonic obeys

\[
Y''+3Y'+\mu^2Y=0,
\qquad Y\propto e^{sT},\qquad s=-3/2+\sqrt{9/4-\mu^2}.
\]

At the static shell, T=H₀τ. This growth rate can therefore be compared directly with the growth of the scalar perturbation at the brane in the conformal-time nonlinear code while the perturbation remains linear.

| Detuning | Archived μ² | Expected s |
|---|---:|---:|
| .001 | −7.717871625176294 | 1.6571936312453648 |
| .1 | −7.528262476157398 | 1.6270213424531335 |
| .03 | −7.662148417066783 | 1.6483564628337088 |

The .001 numbers are in `untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/stability/independent_GN/RESULT.json`, which reports an uncertainty estimate 5×10⁻⁹ on μ² from resolution variation; that is a numerical uncertainty estimate, not a certified bound. Chat 9's separately reported interval enclosure is −7.71788<μ²<−7.71786. The .1 and .03 values are in `untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/runs/linear_predictions.json`.

## Longitudinal-gauge radial equation

Define H=ρ'/ρ, s_b=φ' (the background derivative, not the growth exponent), and G=φ''/φ'. The scalar master equation implemented in Chat 9 `stability/orchestrator_derivation/t_scan_orchestrator.py`, lines 13–31, is

\[
\psi''+2(H-G)\psi'+\left[-\frac43(\phi')^2-4HG+\frac{2+\mu^2}{\rho^2}\right]\psi=0.
\]

The regular cone data have ψ∝y^α with α=1/2+√(9/4−μ²). Reconstruct the scalar perturbation as

\[
\chi=-\frac3{\phi'}(\psi'+2H\psi).
\]

The scalar boundary condition is

\[
\chi'+2\phi'\psi+\frac{\sigma''}2\chi=0.
\]

This is a useful separate shooting reference, not the same gauge as the nonlinear conformal-coordinate variables. One must transform the mode before using it as constraint-satisfying nonlinear initial data.

## Independent Gaussian-normal formulation

Source: Chat 9 `stability/independent_GN/GN_DERIVATION.md`, sections 2–5. With Q=E', the bulk system can be integrated in (χ,χ',Q,ψ):

\[
\psi'=Q-\phi'\chi/3,\quad Q'=-4HQ-2\psi/\rho^2,
\]
\[
\chi''+4H\chi'+(4\psi'+\mu^2Q)\phi'+(\mu^2/\rho^2-U_{\phi\phi})\chi=0.
\]

The gauge-invariant shell mismatch is

\[
\mathcal M(\mu^2)=\chi'+\tfrac12\sigma''\chi+
 (\phi''+\tfrac12\sigma''\phi')\rho_b^2 Q\quad\text{at the shell}.
\]

These radial equations and junction identities were subsequently checked symbolically in Chat 11. The larger claims of completeness, decay, nonlinear stability, or physical reheating do not follow from matching this one mode. A nonlinear solver should also check both Einstein constraints, boundary residuals with the correct differentiated conditions, grid convergence, and causal separation from its truncated outer boundary.
