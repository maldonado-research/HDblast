# Prospective finite-collar diagnostic

This selection was written before reading the instrumented rerun values. It is a new numerical diagnostic, separate from the archived junction-inferred result.

Use the two reruns of `main_ps0.5_G100_y1_bcut_dzf{1e-3,5e-4}`. Use the fixed computational collar Z in [-0.04,0], covered by the saved grid [-0.081,0], and the same shell clock window s=H0 tau in [1,4]. Clip its endpoint computational times by linear interpolation of the saved monotone (s,T) pairs. Do not optimize the collar or interval using residual size.

Use the scalar stress-energy current for xi=exp(-B) partial_T, with common spatial normalization exp(-3 A_b(0)) and the code's kappa5-squared energy multiplier. Set q=phi_T, r=phi_Z, and

    E = exp(3(A-A0)-B) [ (q²+r²)/2 + exp(2B) U(phi) ],
    F = -exp(3(A-A0)-B) q r,
    G = exp(3(A-A0)-B) [ (3A_T-B_T)r²/2 - (3A_T+B_T)q²/2
                         +(3A_T+B_T)exp(2B)U(phi) + B_Z q r ].

The exact continuum equation is partial_T E + partial_Z F = G. Therefore compare the endpoint collar energy change against the time integral of F(-0.04)-F(0)+integral G dZ. This includes the bulk potential and non-Killing geometry work. The energy need not be positive. It is not a global five-dimensional conserved energy.

Compute phi_Z and B_Z directly from the saved field values, never from a Neumann boundary value or an old derivative routine. Use separate degree-two and degree-four local polynomial stencils, called derivative orders two and four. At Z=0, all derivative stencil nodes must be strictly interior, Z<0. At the inner endpoint Z=-0.04, use nearby nodes on both sides from the saved extension. All other points use the closest allowed nodes. Interpolate state fields at a non-node inner endpoint with the same polynomial degree. Retain the actual wall state values, not an extrapolated wall value.

Use composite trapezoidal integration on the actual spatial nodes with the exact two collar endpoints added. Use all saved temporal nodes between the two clipped endpoints, plus separate linear interpolation of the already spatially integrated quantities at those endpoints. The words second/fourth order refer to the derivative stencils only; trapezoidal integration remains second order. Report every-other-time-node sensitivity and both stencil orders, without converting their agreement into a certified error bound.

Report the terms separately: Delta E_coll, integral inner flux, integral wall flux, integral metric work, their combined residual, the residual divided by the sum of absolute integrated terms, and the residual divided by max(abs(Delta E_coll),abs(predicted change)). The first denominator exposes large cancellations; it must not replace the second. The comparison is per one retained copy. At the shell, F(0)dT=-a³wv d tau; twice this wall flux has the same orientation and normalization as the total shell influx in the archive note.

Any failure of the complete balance or strong sensitivity must be retained. No independent bulk donor conclusion may be inferred merely from agreement of a boundary derivative with its imposed junction.
