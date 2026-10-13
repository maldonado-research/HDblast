# Review of the independent root current calculation

13 October 2026. Reviewed by the parallel covariant-current agent after its
own unrestricted lapse/conformal-family calculation was complete. This is
internal mathematical review, not external peer review.

Files reviewed:

- `../root-current-audit/INDEPENDENT_CURRENT_DERIVATION.md`
- `../root-current-audit/check_current_independently.py`

No blocking sign or normalization discrepancy was found. The two displayed
radial connection contractions give exactly the potential derived directly
from `nabla_B h^{AB}-nabla^A h` in the parallel calculation. The nonlinear
restriction `n=-2 psi` is a valid commuting extension for computing the
pullback symplectic curl at the background; it does not assume nonlinear
field equations. Its scalar density variation is essential for the
`-6p` cross term before the momentum constraint, and the constraint produces
the correct `-2p` term. The tangent current has the matching relative sign.

The script's chain-rule radial derivative uses both scalar first-order
equations and `rho'=a rho` consistently. The mixed-eigenvalue conservation
law, vanishing even shell current and positive odd shell-current divergence
agree with the complete action orientation. The characteristic Jacobian and
the inclusion of both `j^y` and `j^t` reproduce the scalar null current with
the correct per-copy coefficient. The warning that the symplectic bilinear
is not a positive energy flux is correct.

The script verifies consequences of the displayed nonlinear potential;
the connection-contraction derivation itself is supplied analytically by
the accompanying note and was inspected in this review. Neither script
establishes existence of global registered solutions with the assumed
characteristic data, admissibility of residual gauges, or completeness.
The root note keeps these limits explicit.
