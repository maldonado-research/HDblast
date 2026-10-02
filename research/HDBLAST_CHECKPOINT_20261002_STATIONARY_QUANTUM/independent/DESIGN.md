# Independent stationary closure review: prospective design

Prepared before new quantum source evaluations or new shooting roots. This review treats the proposed family as an explicit dimensionless mathematical model; it has no empirical mass or gravity normalization.

## Common action and units

Restore a length L only to keep dimensions visible: yhat=y/L, Rhat=R/L, xhat=L^2 m^2, zhat=L^2 H^2, rhohat=L^4 rho, Jhat=L^4 J_phi. The inherited potential and tension functions are dimensionless, with phi dimensionless. Freeze N=1 canonical free field and vary gamma=kappa_5^2/L^3; it multiplies rhohat and Jhat in the respective dimensionless junctions. Arbitrary positive gamma is a continuous mathematical amplitude, not a statement that a gravitational loop expansion is controlled. Gamma=0 is the classical source-free control, while the largest gamma values are formal stress cases with EFT control unestablished.

For f=1+b(phi-phi_b0)/2 and x=r f^2, hold r=2 H_0^2 fixed once the source-free branch has supplied H_0 and phi_b0, and require f>0. Then x_phi=b r f, x_phiphi=b^2 r/2. Evaluate phi-phi_b0 as eta-eta_ref to avoid subtracting nearly equal doubles. At b=±1 this is the original generic quadratic interaction with m0=0, Ghat=sqrt(r)/2, phi_chi=phi_b0-2/b; at b=0 it is constant mass m0^2=r, Ghat=0. The massless global endpoint is outside the positive-factor domain. Changing H at a trial shell changes x/z and r/z; it must not reset r to 2z. Recomputing a reference for every trial is a different action and invalidates the source derivatives already derived.

The complete fixed-reference Euclidean action is Gamma_E=8 pi^2 W/(3z^2). Differentiate at fixed x,r when varying the S4 radius. Under a constant rescaling of its metric, delta Gamma_E=-4 Vol rho delta(log H), giving rho=W-z W_z/2. The scalar source follows independently from Q=2W_x, J_phi=x_phi Q/2. An invariant state has T_mn=-rho h_mn and p=-rho; this gives the full static stress. W alone is not rho, and J_phi is not generally partial_phi rho.

For independent evaluation, use the resolvent finite part rather than differentiating the producer determinant numerically. With u=x/z, v=r/z, nu^2=9/4-u and Psi=psi(3/2+nu)+psi(3/2-nu), retaining the s=1 zeta pole gives

    FP Z(1;u)=[(u-2)Psi-u+4/3]/6,
    Q=z/(16 pi^2)*[(u-2)(Psi-ln v)-u+v+4/3].

Direct metric variation of the fully normalized spectral action (including all reference polynomials), rather than imposing its trace, gives

    rho=z^2/(64 pi^2)*[u(u-2)(Psi-ln v)-u^2+4u/3
                           +2uv-2v-v^2/2-D(u)],
    D(u)=u^2/2-2u+29/15.

The independent proper-time route will evaluate the metric integrand itself at selected converged endpoints, retaining both the spectral zeroth moment and first moment. The trace identity is a subsequent diagnostic, not the source definition. The useful common-action cross derivative rho_x=Q/2-z Q_z/4 is nontrivial even though the static Ward identity is 0=0.

## Independent radial equations

Integrate four second-order acceleration variables (R,P,eta,w), phi=1+eta:

    R'=P,  P'=-R*(w^2/4+U/6),
    eta'=w,  w'=U_phi-4 P w/R.

This does not impose P=sqrt[1+R^2(w^2/12-U/6)] at every step. Its Hamiltonian defect

    C=P^2-1-R^2(w^2/12-U/6)

has C'=0 analytically after the displayed equations are substituted. Thus an independent constraint monitor tests both the regular initial series and the acceleration evolution. Use the exact polynomial in eta to retain the displacement of order 10^-23 at the cone. Put a=-U_h/36, B=U_h^2/4320-U_phi,h^2/750, d=U_phi,h*(U_phiphi,h/280+U_h/630). The regular cone data are

    R=y+a y^3+B y^5, P=1+3a y^2+5B y^4,
    eta=eta_h+U_phi,h y^2/10+d y^4,
    w=U_phi,h y/5+4d y^3.

Solve the two actual residuals in coordinates ell=log10(abs eta_h), y_b, eta_h<0:

    E1=P_b/R_b-sigma(phi_b)/6-gamma*rho(z_b,x_b,r)/6,
    E2=w_b+sigma_phi(phi_b)/2+gamma*x_phi,b*Q(z_b,x_b,r)/4.

A regular bulk and both shell equations are necessary. The algebraic endpoint identity H_b^2=(sigma+gamma rho)^2/36-(sigma_phi+gamma J)^2/48+U/6 alone neither determines the regular scalar profile nor fixes the bulk Weyl sector.

## Prospective independent checks and gates

Root will freeze the complete primary grid and selected independent points before execution. No numerical source or root evaluation is authorized by this design alone. Unless the root registration specifies stricter gates, use independent DOP853 acceleration solves at (rtol=1e-11,y0=1e-3,max_step=.1) and (rtol=3e-14,y0=1e-4,max_step=.05), with absolute tolerances [1e-14,1e-14,1e-100,1e-100]. Record every attempted solve including failures. The selected points must be chosen before primary outputs are read; independent starts are inherited vacuum coordinates and continuation only within the independently computed registered sequence.

Require |E1|/delta and |E2|/delta <=2e-12, |C|/max(1,P^2) <=2e-10 throughout, and agreement between independent refinements within 2e-9 in phi and 2e-8 relative in positive H^2. Compare against the producer within those same observable gates, and within 2e-7 in each shooting coordinate. Independently compare source values against direct proper-time integrals at a preselected subset using source units z^2 for rho and z for Q; tolerance 1e-9 absolute + 1e-7 relative in these normalized source units. Report empirical refinement and bounded spectral/IR tails separately from uncertified quadrature and UV truncation.

Negative controls are evaluated at the correct solutions without tuning new branches: reverse J in the scalar junction; drop J for nonzero b; use W in place of rho; use p=+rho in the temporal Israel condition; permit r to track x under differentiation; and change the metric-variation curvature sign. A sign/drop control is meaningful only when its expected defect exceeds the predeclared absolute floor. The wrong r variation contributes W_r*x_phi to J, with W_r=-[(x-r)^2/2-2z(x-r)+29z^2/15]/(32 pi^2 r); this makes the control analytical and avoids redefining the model by accident. Report detection/nonapplicability separately at gamma=0 or b=0.

## Scope and closure dangers

1. De Sitter invariant vacuum polarization supplies p=-rho, not a radiation fluid. Static-patch thermality is not reheating or a particle population.
2. A finite-amplitude root resums source feedback in a semiclassical one-loop action. Success does not establish a small loop expansion, a physical EFT cutoff, or omitted graviton/higher-loop control.
3. Constant curvature-squared counterterms have zero metric variation on four-dimensional exact de Sitter, but field-dependent loop coefficients retain their scalar derivative. The retained exact determinant already includes them; adding or subtracting them again double counts the action.
4. De Sitter invariance and tracelessness force the projected Weyl tensor E_mn to vanish in this exact warped ansatz; an arbitrary dark-radiation contribution is incompatible with that symmetry. Integrating a regular radial interior still selects the scalar/bulk solution. Hamiltonian and shell checks do not prove a complete Lorentzian continuation, nonlinear stability, or a trajectory from the original unstable branch. Relaxing the symmetry requires bulk Weyl dynamics again.
5. Nonlocal effective-action variations away from this invariant family are needed for perturbations and time evolution. The static source formulas and root Jacobian are insufficient to infer quantum dynamical stability.
6. Small gamma*rho relative to sigma is not a sufficient perturbativity criterion near the detuning cancellation. Track Delta H^2/H_0^2, Delta phi/abs(phi_b0-1), the independently frozen linear prediction, and the size of its nonlinear remainder.
7. Benchmark grid failure must remain failure evidence. No post-result tuning of b, gamma, r, c, delta, stopping boundaries, or acceptance gates is allowed.
