# A finite-collar scalar-stress identity for the instrumented B1 run

This identity specifies the new diagnostic. It does not itself assert that the sampled numerical solution satisfies it.

For ds²=e^(2B)(-dT²+dZ²)+e^(2A)dX², the scalar equation is

    phi_TT - phi_ZZ = -3 A_T phi_T + 3 A_Z phi_Z - e^(2B) U_phi.

The scalar stress is T_AB=[phi_A phi_B-g_AB((partial phi)²/2+U)]/kappa5². Choose the explicit unit-observer time vector xi=e^(-B) partial_T, not a presumed Killing field. Its current is J^A=-T^A_B xi^B. Per coordinate three-volume, with the common constant normalization e^(-3 A_b(0)) and with energies multiplied by kappa5², the densitized components are exactly the E,F in `COLLAR_PROTOCOL.md`:

    E = e^(3(A-A0)-B) [(phi_T²+phi_Z²)/2+e^(2B)U],
    F = -e^(3(A-A0)-B) phi_T phi_Z.

Differentiate these two products, use phi_TZ=phi_ZT, and substitute the scalar equation. The remaining terms are

    partial_T E + partial_Z F
      = e^(3(A-A0)-B) [ (3A_T-B_T)phi_Z²/2 -(3A_T+B_T)phi_T²/2
                        +(3A_T+B_T)e^(2B)U+B_Z phi_T phi_Z ] = G.

This is the explicit form of div(J)=-T^(AB) nabla_(A xi_B). Retaining G is essential; the spacetime is time dependent and the observer field can accelerate. Integrating on the fixed collar [Zc,0] gives

    E_coll(T1)-E_coll(T0)
      = integral_T0^T1 [ F(Zc)-F(0)+integral_Zc^0 G dZ ] dT.

Positive F(0) is energy leaving the retained interior through the shell. At the wall, phi_T=e^B v and phi_Z=e^B w, so

    F(0) dT = -a³ wv d tau.

Thus twice the outward wall flux equals the total shell influx, with no extra lapse or factor of two. The inner flux uses the opposite outward normal of the collar, explaining its positive sign on the right-hand side.

For comparison, choosing xi=partial_T instead gives E_Tcoord=e^(3(A-A0))[(phi_T²+phi_Z²)/2+e^(2B)U], F_Tcoord=-e^(3(A-A0))phi_T phi_Z and G_Tcoord=e^(3(A-A0))[3A_T(phi_Z²-phi_T²)/2+(3A_T+2B_T)e^(2B)U]. That current has an extra shell lapse in its boundary integral. Mixing its collar energy with the unit-observer shell flux is incorrect.

Neither current supplies an observer-independent, positive, globally conserved scalar-plus-gravity reservoir. The selected diagnostic is a local stress balance with an explicitly declared observer field, inner boundary, metric-work term, and energy normalization.
