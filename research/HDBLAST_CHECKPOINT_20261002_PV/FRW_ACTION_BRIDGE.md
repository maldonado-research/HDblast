# Action-level bridge to an expanding shell

This is an independently checked analytical identity for the next calculation. No FRW quantum-source integration, new junction solver or bulk evolution has been implemented here.

With signature (-+++), R=6(Hdot+2H²), and the common action term
Gamma_F=integral sqrt(-g) F(phi)R,
variation before imposing the background gives
T_munu,F=-2[F G_munu+(g_munu box-nabla_mu nabla_nu)F],
J_phi,F=-F_phi R,
where the current convention is J_phi=-delta Gamma/(sqrt(-g) delta phi).

For spatially flat FRW:
rho_F=-6F H²-6H Fdot,
p_F=2F(2Hdot+3H²)+2Fddot+4H Fdot.
Direct differentiation verifies
rhodot_F+3H(rho_F+p_F)=-R Fdot=J_phi,F phidot.
At H=0, this reduces to rho_F=0, p_F=2Fddot, J_phi,F=0, precisely the pressure term tested in the flat benchmark.

A shell implementation must vary the same covariant brane action into both metric and scalar junctions. A tension-only replacement misses the induced Einstein tensor, derivatives of F, and the curvature-dependent scalar current. Local matching stress must not be counted twice between the brane action and a separate quantum source.

For constant alpha, the action alpha R² contributes
T_munu=-2alpha[2R R_munu-.5g_munu R²+2(g_munu box-nabla_mu nabla_nu)R].
It vanishes on flat and constant-curvature Einstein backgrounds, but generally contributes higher derivatives on time-dependent FRW. Vanishing Weyl curvature does not remove R². A constant-coefficient Euler term in four dimensions is topological subject to boundary terms; a phi-dependent coefficient is not.

The next physical step requires smooth background/state preparation, the remaining curvature matching conditions, a common absolute stress/current prescription, and only then consistent coupling through all shell junctions. The current finite matching choices are benchmark definitions, not measurements of HDBLAST parameters.
