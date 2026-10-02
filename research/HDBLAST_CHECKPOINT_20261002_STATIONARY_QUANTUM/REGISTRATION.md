# Prospective stationary quantum shell continuation

This is a prospective, dimensionless mathematical model study. No source
evaluation, new radial integration or root solve is authorized until this file,
the exact producer and validator, EXPERIMENT.json, and their hash registration
have been committed and published by the coordinating researcher. The producer
requires an explicit public-freeze reference and verifies all registered hashes.
Preparation may read inherited public outputs, perform symbolic algebra and
compile source syntax. It may not evaluate new physical sources or integrate.

## Model fixed before calculation

Use the inherited corrected negative eta_h=phi_h-1 branch at delta=.001 and
c=.5975949350280132. The archived shooting seed is
eta_h=-1.3366923651588084e-23, y_b=30.276680395885563. The inherited reference
endpoint is eta0=-.00008405268308670538 and z0=H0²=.00005924014794328956.
Fix r=2z0=.00011848029588657912 once, including in every scalar and metric
variation. Let f=1+b(eta_b-eta0)/2 and x=r f², with b=-1,0,+1.
For b=±1 this is the original quadratic interaction family with m0=0,
Ghat=sqrt(r)/2 and phi_chi=phi0-2/b; b=0 is a constant positive mass squared r.
The fixed choices imply x_phi=b r f and x_phiphi=b² r/2. Require f>0 throughout
the local continuation domain. The reference x=r and x_phi=b r agree for every
integration method, without fitting the mass law to computed source values.

There is one free minimally coupled massive shell scalar, N=1, in its Euclidean
(Bunch-Davies) de Sitter invariant state. gamma=N kappa5²/L³=g5 labels the
declared mathematical gravitational model family. gamma=0 is a classical
control. gamma={0,.01,1,100,10000,1000000} is fixed before source values. These
are chosen research parameters, not empirical measurements. Large gamma probes
formal stationary closure; it does not establish a physical EFT hierarchy or
the validity of neglecting quantum gravity, higher loops or extra operators.
FIXED_MODEL.json records the shared exact decimal model inputs. For N=1 define
M5 L=gamma^(-1/3). Separately screen Ehat*gamma^(1/3)≤.1, with
Ehat=max(sqrt(z_b),sqrt(x_b),abs(q_b),max_sampled sqrt(abs(kT)),
max_sampled sqrt(abs(kN))), kT=U/6-w²/12, kN=U/6+w²/4. Samples include all
adaptive integrator nodes and 1001 equally spaced dense-output samples. Only
the endpoint q_b is used because cone slicing q diverges as a coordinate effect.
This is a sampled scale diagnostic and declared tenfold hierarchy screen, not
a certified cutoff. Screen failures do not stop or invalidate the mathematical
root tests; all high-gamma formal closure results are retained and labeled.

The finite common-action matching convention is inherited unchanged from the
public DESITTER checkpoint, evaluated with fixed positive r. There are no
additional local terms, retuned counterterms or adjustable subtractions. The
sources are rho=W-z W_z/2, Q=2 W_x, j=x_phi Q/2. Both s=gamma rho and t=gamma j
enter the shell junctions. rho and j are evaluated at each current endpoint.
They are not frozen to their reference values, and j is not rho_phi.

## Constrained bulk and nonlinear junction solve

Integrate the regular cone in eta=phi-1 and positive R' branch using the pinned
DESITTER sensitivity implementation. The integrated state is
(R,eta,w,R_ell,phi_ell,w_ell,I), ell=log10(-eta_h), I'=R R_ell. Its inherited
regular series includes R through y^5, eta through y^4 and I through y^7.
The constraint fixes R'=sqrt[1+R²(w²/12-U/6)] and the scalar equation is
w'=U_phi-4(R'/R)w. The shooting unknowns are (ell,y_b).

Solve C=B=0, with B=w+sigma_phi/2+t/2 and

    C=z+w²/12-U/6-(sigma+s)²/36.

Avoid subtracting two order .01 metric terms: W=1/3+eta²+eta³/3,
p=W_phi=2eta+eta², tau=delta(1+c phi)+s, and evaluate the algebraically exact

    C=z+(w-p)(w+p)/12-W tau/9-tau²/36.

Require positive q and sigma+s, so C=0 is equivalent to the original unsquared
metric junction E1=q-(sigma+s)/6. The validator reconstructs the unsquared
junction with high precision and checks the exact endpoint constraint.

The full Newton Jacobian contains all x_phi=b r f, x_phiphi=b² r/2 and source z derivatives at
fixed r. Source derivatives come from the analytic second derivatives of the
same convergent Hurwitz-series action, with rho_x=W_x-z W_xz/2,
rho_z=W_z/2-z W_zz/2, Q_x=2 W_xx, Q_z=2 W_xz. A separate high precision
Richardson action check validates these derivatives. The radial q_ell derivative
uses the cancellation-resistant identity q_ell=-2I/R^4-w phi_ell/3; the moving
endpoint has q_y=-z-w²/3. The exact off-shell product rule differentiates
C=E1[q+(sigma+s)/6]. Direct differentiation of the constraint is a check.

The explicit Newton cap is 12 iterations with 12 halvings, coordinate step caps
(.3,.2), |C|≤2e-18 and |B|≤2e-14. Seeds use only the archived source-response
matrix and, after publication, the reference sources and last successful branch
point. Source/residual evaluation always retains both junctions. The residual
scales are fixed z0 and delta*c/2. All rejected trial steps and failed cases are
retained. A failed point does not license changed tolerances or a new grid.

## Coverage, checks, and stopping boundaries

EXPERIMENT.json is the machine-readable specification. DOP853 primary and
refined resolutions solve all 18 registered points. Radau and an independent
cone-start refinement each solve gamma=0 and 1e6 for all three b values. Thus
48 intended successful solves are enumerated, including duplicate classical
controls. A separate reviewer will use a second-order bulk integration at the
six declared points b=±1 and gamma={0,100,1e6}; its source and integration
implementation are separate from this producer.

Hard endpoint domain: .75≤z/z0≤1.25 and |eta-eta0|≤.1|eta0|. Require
0<x/z<4.5, the positive-gap convergent Hurwitz domain. Shooting bounds are
-28<ell<-18 and 25<y_b<35. Leaving any boundary is reported as a bounded
failure; no positive-eta_h branch switch or endpoint search is allowed. The
domain bounds do not claim uniqueness, locate a fold, or establish the maximal
physical branch. A failure is not evidence of nonexistence.

Primary sources use 60 decimal digits and 128 terms; validation uses 80 digits
and 256 terms with inherited analytic truncation bounds. Every solved point
requires |C|≤2e-16, |B|≤2e-13, positive total tension and finite nonsingular
Jacobian. Source refinement uses 1e-22 absolute plus 1e-10 relative. Analytic
source derivative comparisons use 1e-20 absolute plus 1e-8 relative. The two
forms of constraint derivative agree to absolute 2e-17. Across integrations,
absolute differences must be ≤2e-14 in H², ≤2e-11 in eta and ≤2e-8 in shooting
coordinates. These are empirical floating-point checks, not rigorous interval
error bounds. They cannot resolve arbitrarily small differences by themselves.
For the action Hessian also bound the omitted series second derivative: with
a=3/2, theta=|q_s|/a², q_s=x/z-9/4, and
common=a³[1+a/(2N-2)]theta^(N+1)/(1-theta), take eB=common/(N+1),
eD=common/|q_s| and eDD=common/q_s² [N+theta/(1-theta)]. At q_s=0 the
omitted derivatives vanish for N≥2. Then eQx=2C0 eDD,
eQz=2C0(eD+|x/z|eDD),
eWzz=C0[2eB+2|x/z|eD+(x/z)²eDD], C0=1/(16pi²).
Triangle inequalities propagate these bounds to rho_x,rho_z,rho_eta,j_eta,j_z;
each recorded derivative tail must be below 1e-20. These are analytic series
tail bounds, while arithmetic and radial integration remain empirical checks.

Each resolution's shifts use its own solved gamma=0 control, cancelling the
shared archived seed mismatch. Report each shift against ten times its observed
cross-resolution spread plus a conservative floating-point floor (1e-18 for z,
1e-17 for eta). For gamma=.01 especially, unresolved shifts are reported as
unresolved. Do not call all registered shifts detected merely because roots
converged. Compare with the archived differential-response prediction. A
nonlinear correction is classified resolved only if it is at least 1e-3 of the
observed shift and larger than ten empirical refinement spreads plus the same
floor. This detection criterion is separate from numerical root acceptance;
failure yields a bounded linear-response-compatible result, not extra sampling.

Three wrong models are re-solved at b=+1,gamma=1e6 using the refined solver:
flip the scalar current sign, omit the metric source, and freeze both sources
at the archived reference endpoint. The correct residual at each wrong root
must exceed ten times at least one corresponding root-acceptance threshold.
The original wrong-model convergence and the correct residual are both saved.

## Provenance and permissible conclusions

Freeze hashes of this file, EXPERIMENT.json, FIXED_MODEL.json, both codes and requirements, and
all inherited action, source, sensitivity and bridge files before calculation.
Every execution uses a fresh output directory, logs its public freeze reference,
registration and producer hashes, observed checkout HEAD, library versions,
successful rows, rejected steps, failures and mutation results. Existing outputs
are never overwritten. The validator exits nonzero on any acceptance failure;
its checks do not rely on Python assert and work under python -O.

Permissible result: a bounded numerical stationary regular-cone solution family
for this declared one-field quantum-action benchmark, with empirical integration
and source checks and a separately measured nonlinear response. No uniqueness,
quantum stability, nonlinear attraction, time evolution, heating, particle
production, exit history, observational fit or physical parameter determination
follows from these stationary solves. A negative or unresolved outcome remains
a result and will be preserved with the same scope.
