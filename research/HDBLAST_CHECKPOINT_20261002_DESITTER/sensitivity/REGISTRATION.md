# Prospective stationary-branch sensitivity registration

Registered 2026-10-02 before this directory's first radial integration.

Scope: numerical differential response of the already selected corrected +1,
delta=0.001 classical branch. No shell mass, coupling, finite quantum matching,
quantum source, dynamical bulk, or new physical branch is selected here. Supplied
sources s=K rho and t=K J are independent infinitesimal inputs.

Frozen archived background: e_h=-1.3366923651588084e-23,
y_b=30.276680395885563, c=0.5975949350280132; ell=log10(abs(e_h)).
Use exact shifted polynomial U about phi=1 and sigma=2W+.001(1+c phi).
Pin archived BACKGROUNDS.json, gi_core.py, spectrum.py, and static bridge note
by SHA256 before calculation. Keep all runs and failures; no source/output edits.

## Derivation registered before numerical evaluation

Let f=partial_ell phi, z=partial_ell R/R, v=z'=partial_ell q.
Differentiate R'=sqrt(1+R²(w²/12-U/6)), phi'=w,
w'=U_phi-4qw and the analytic regular-cone series through R(y^5), phi(y^4).
The independent weighted identity is

    S=v+w f/3,   S'+4q S=-2z/R²,
    S_b=-2 R_b^(-4) integral_0^yb R² z dy.

The cone boundary term vanishes. Off-shell the exact derivative is
E1_ell=S_b-E2*f/3. This representation avoids subtracting nearly equal
metric/scalar terms and determines whether the archived numerical-floor entry
is generically nonzero. Independently derive this identity from the constraint
and acceleration equation before use. At a vacuum shell, E1_y=-H² and
E2_y=Bw. Also report exact off-shell derivatives to avoid hiding baseline error.

## Fixed numerical sequence and acceptance

Primary: double-precision DOP853 variational equations with integrated weighted
quadrature, max_step=.1, rtol in {1e-10,1e-12,3e-14}, y0 in
{1e-2,1e-3,1e-4,1e-5}; tiny scalar/sensitivity absolute tolerances 1e-100,
R absolute tolerance 1e-14. Refine max_step to .05 at the finest settings.
Independent methods: Radau rtol {1e-10,1e-12}, y0 {1e-3,1e-4}; analytic
complex-step shooting derivative at h={1e-15,1e-20}, DOP853 rtol=3e-14,
y0=1e-4. Use the same endpoint and check complex-step output derivatives
without importing the variational Jacobian.

Report the full run envelope plus conservative 10x observed cross-method/spread
estimates. Required absolute agreements: 2e-10 in phi_ell,
2e-14 in (H²)_ell, 2e-10 in each ordinary Jacobian entry;
for tiny E1_ell require direct variational and weighted-integral agreement
<=max(1e-19,1e-5*abs(E1_ell)), plus stable sign and start refinement.
For supplied-source observable susceptibility require <=1e-7 absolute agreement
in phi coefficients and <=1e-10 in H² coefficients. These are numerical
acceptance gates, not interval bounds or formal ODE truncation certification.

Negative controls: omit the scalar junction contribution to E1_ell; reverse
the scalar-current sign in the linear solve; omit the moving-endpoint terms.
Each must be detected by the analytic integral or differentiated constraint,
or by a nonzero registered observable discrepancy. Also compare the archived
finite-difference E1_ell to the newly resolved value without using it to fit.

Deliver phi_ell, (H²)_ell, residual Jacobian, d(ell,y_b)/d(s,t),
and d(phi_b,H²_b)/d(s,t), their refinement/error scales and the bound
|Delta phi| <= |C_phi,s s|+|C_phi,t t|,
|Delta H²|/H² <= (|C_H,s s|+|C_H,t t|)/H².
Do not claim the response is small without actual supplied source magnitudes.
