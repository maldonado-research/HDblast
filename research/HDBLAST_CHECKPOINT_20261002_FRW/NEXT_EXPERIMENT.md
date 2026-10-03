# Next experiment: a smooth curved control with an exact incoming vacuum

This is a proposed follow-up, not an executed result or completed preregistration. First resolve the finite curved-action matching, then register the numerical matrix before executing it.

Use the C-infinity step B defined in FRW_COUNTERTERMS.md, with a(eta)=exp[A B(eta/T)] and x(eta)=r[2B(eta/T)-1]^2. It is exactly static before eta=0 and after eta=T and crosses x=0. This removes both finite-order spline knots and finite-time approximate initial-state ambiguity. A useful initial control is r=4, T=1, A=0.2, with A=0 serving as the independent flat-limit check; these are prescribed-background test parameters, not fitted HDBLAST parameters.

Before evaluating a physical source:
1. Derive the complete relation between the positive-reference subtraction and a common covariant action. Explicitly state F(r)R and curvature-squared finite constants. Flat agreement alone does not fix them.
2. Use exact incoming static modes and a stable physical all-k mode solver; never substitute a possibly negative truncated WKB frequency as the physical low-k state.
3. Compute energy density, pressure and the paired mass-squared current from the same prescription. Test the local exchange identity including source work.
4. Register momentum cutoff, quadrature, time resolution and independent-solver refinements, with absolute as well as relative residuals. Separate coherent stress from occupation-only estimates.
5. Only after these tests pass consider a smooth background reconstructed under the original shell equations with controlled derivatives, initial prehistory and consistent backreaction.

The current source-free shell, and especially its eventual recollapse, cannot be declared changed by this proposed test. The later physical task is to demonstrate source consumption in coupled evolution, an energy ledger, decay and thermalization over a sustained expanding interval.
