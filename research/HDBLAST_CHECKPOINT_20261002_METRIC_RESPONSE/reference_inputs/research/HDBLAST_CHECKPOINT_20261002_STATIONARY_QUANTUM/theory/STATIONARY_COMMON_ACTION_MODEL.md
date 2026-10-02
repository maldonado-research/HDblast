# A constrained stationary shell model from the fixed-reference de Sitter action

Prepared 2 October 2026. This document defines a **mathematical model family**
and derives its stationary equations. It contains no new source evaluation,
radial integration, root search, or inferred physical parameter. The inherited
empty-shell solution and matched free-field action are inputs. A numerical
amplitude grid, stopping gates, methods, and files must be frozen separately
before numerical execution. The equations define a one-loop stationary
closure; solving them at finite amplitude does not calculate higher-loop
effects or establish dynamical stability.

## 1. Dimensions, normalization, and the model choice

Restore the inherited bulk length unit explicitly as a positive constant L.
Write dimensional quantities in terms of dimensionless quantities with hats:

    y = L yhat, R = L Rhat, H² = L^-2 z,
    U_phys = L^-2 Uhat, sigma_phys = L^-1 sigmahat,
    x_phys = L^-2 x, r_phys = L^-2 r,
    W_phys = L^-4 W, Q_phys = L^-2 Q,
    rho_phys = L^-4 rho, J_phi,phys = L^-4 j.

The remainder uses the dimensionless variables z,x,r and drops hats from the
bulk functions. The dimensionless geometric scalar phi has canonical bulk
normalization Phi = phi/kappa5, where kappa5 = sqrt(kappa5²) and
[kappa5²] = length³. Thus Phi has mass dimension 3/2. A canonical shell scalar
chi has mass dimension one. There are N identical, independent canonical
shell scalars; regulator multiplicities never count as physical N.

The doubled bulk and single-shell action is

    S = sum_two_copies {1/(2 kappa5²) integral sqrt(-g)
                          [R5-(nabla phi)²-2 L^-2 U(phi)]
                      +1/kappa5² integral_shell sqrt(-h) K}
        -1/(kappa5² L) integral_shell sqrt(-h) sigma(phi)
        + lambda N Gamma_1,ref[h, L^-2 x(phi); L^-2 r].

lambda is an optional formal loop-counting factor. Define the only combination
entering these stationary equations by

    gamma = lambda N kappa5²/L³.

For a literal N-species semiclassical interpretation set lambda=1 and specify
N, kappa5²/L³, and the physical EFT cutoff independently. Continuous gamma is
a mathematically legitimate source amplitude without any such physical
identification. A number assigned to gamma does not measure those separate
inputs. The common loop factor multiplies the entire matched determinant,
hence BOTH rho and j. The classical local shell coefficients are fixed while
gamma changes. Sources have dimensions

    [rho_phys]=[J_phi,phys]=mass^4, [Q_phys]=mass²,
    [kappa5² rho_phys]=[kappa5² J_phi,phys]=mass,

which match q_phys, w_phys, sigma_phys, and sigma_phys,phi in the two junctions.
The canonical scalar source is J_Phi = kappa5 J_phi. Changing normalization of
phi without also changing its mass-law derivative changes the model.

The final registered family specializes to **N=1, lambda=1**, and
gamma=kappa5²/L³. The gamma=0 point is its formal classical control. Large
positive gamma values are mathematical closure stress cases; they are not
automatically physically valid semiclassical gravitational parameters.

Use the inherited corrected regular branch at its frozen detuning and c, with

    B(phi)=1-phi+phi³/3,
    U=B_phi²/2-2B²/3,
    sigma=2B+delta(1+c phi).

Here B is the bulk superpotential, distinguished from the shell determinant W.
Let (phi_b0,z0) denote the unperturbed endpoint of the chosen regular branch.
The following mass law is a declared family, not a conclusion about nature:

    r = 2 z0 > 0,          a(phi)=1+b(phi-phi_b0)/2,
    x(phi)=r a(phi)²,      b in {-1,0,+1},
    x_phi=b r a,           x_phiphi=b² r/2.

In dimensional terms m²(phi)=L^-2 r a(phi)². For b=+/-1 this is the existing
matter-action interaction m0²+G²(phi-phi_chi)² with m0=0,
G=sqrt(r)/(2L), and phi_chi=phi_b0-2/b. For b=0 take G=0 and m0²=r/L².
The canonical-field coupling in G²(phi-phi_chi)² is G kappa5 multiplying
Phi-Phi_chi. These values are declared new couplings within the inherited
quadratic interaction family, not measured parameters. Positivity is ensured
on the registered local domain; the distant zeros at phi=phi_chi for b=+/-1
are excluded. During all endpoint and metric/scalar variations, r and
phi_b0 are fixed parameters. The equality r=2z0 is imposed once on the
reference branch; **r is not set equal to 2z or x at new endpoints**.

## 2. Common action and fixed finite matching

Use the regular Euclidean state on S4, continued to invariant Lorentzian de
Sitter, for positive x and z. The inherited finite Euclidean action density is

    C0 = 1/(16 pi²), h=x-r,
    A(t)=6t² sum_l [(l+1)(l+2)(2l+3)/6] exp[-l(l+3)t],
    c2 = h²/2-2zh+(29/15)z²,

    W(z,x,r) = -C0/2 integral_0^infinity ds/s³
       {exp(-xs) A(zs)
        -exp(-rs)[1+(2z-h)s+c2 s²]}.

This is a convergent complete action with the normalization fixed; a term
independent of x but dependent on z must not be dropped. Its physical sources
per species and its scalar source in the chosen normalization are

    Q=2 W_x,
    rho=W-z W_z/2,       p=-rho,
    j=x_phi Q/2=b r a Q/2.

The metric derivative includes Vol(S4)=8pi²/(3z²) and holds x,r fixed. The
scalar derivative holds the metric and r fixed. Only AFTER these independent
variations are taken are x=x(phi_b) and z=1/R_b² substituted.

The finite local functions in the inherited determinant include

    V=[x² ln(x/r)-3x²/2+2rx-r²/2]/(64pi²),
    F=[x ln(x/r)-x+r]/(192pi²),
    W_local=V-12zF+(29/15)z² ln(x/r)/(32pi²).

They obey V(r)=V_x(r)=V_xx(r)=F(r)=F_x(r)=0 and vanishing local geometric
logarithmic coefficients at x=r. These conditions specify a finite reference
model, not measured shell matching. The full curved-space action and its
sources generally remain nonzero at x=r. No additional potential, induced
Einstein coefficient, curvature-squared constant, extrinsic-curvature operator,
or other shell matter is added in this declared family. If such an operator
is added later its matched coefficient and both of its variations are new
model inputs; it must appear exactly once.

The exact action identities inherited from independent variations are

    rho_x = Q/2-z Q_z/4,
    rho_phi = j-z j_z/2,
    W_r=-C0 c2/(2r),
    -4rho=-xQ+C0 c2.

The second identity is particularly useful: j is generally NOT rho_phi at
fixed curvature. The trace equation is a check of the action-derived sources,
not a substitute for deriving the metric source.

## 3. Exact derivatives needed for source feedback

Put u=x/z, v=r/z, nu²=9/4-u, and

    Psi(u)=psi(3/2+nu)+psi(3/2-nu),
    A0=Psi(u)-ln(v),
    Psi_u=[psi_1(3/2-nu)-psi_1(3/2+nu)]/(2nu).

For imaginary nu take the real conjugate sum. At nu=0 the analytic limit is
Psi_u=-psi_2(3/2). The exact current and its partial derivatives at fixed r are

    Q=C0 z[(u-2)A0-u+v+4/3],
    Q_x=C0[A0+(u-2)Psi_u-1],
    Q_z=C0[-2A0-u(u-2)Psi_u+u-2/3].

For an action-first numerical implementation use

    rho_x=W_x-z W_xz/2,
    rho_z=W_z/2-z W_zz/2.

Independent exact formulas, suitable for cross-checks, are

    rho_x=Q/2-z Q_z/4,
    rho_z={x Q_z-C0[-2(x-r)+(58/15)z]}/4,
    rho_phi=x_phi rho_x,
    j_z=x_phi Q_z/2,
    j_phi=(x_phiphi Q+x_phi² Q_x)/2.

The x_phiphi Q/2 contribution in j_phi is essential for the quadratic mass
law. Differentiating only Q misses it. Here x_phiphi=b² r/2, so replacing this
term by its value for an exponential mass law changes finite-amplitude
feedback even though the two laws have identical reference x and x_phi.

The inherited convergent determinant implementation writes, with q=u-9/4,

    P(q,z)=B1(q)-B0(q) ln z,
    B0=zeta_H(-3,3/2)-zeta_H(-1,3/2)/4+q/8+q²/4,
    W=-C0 z² P
      +C0[rx-r²/4-2rz-(x²/2-2xz+29z²/15)ln r]/2.

Here B1=3 Z'(0;u) in the inherited series convention. Derivatives must include
P_q=B1_q-B0_q ln z and P_qq=B1_qq-(ln z)/2. Direct differentiation gives

    Q_x=-2C0 P_qq-C0 ln r,
    Q_z=-2C0[P_q-u P_qq-B0_q]+2C0 ln r,
    W_zz=-C0[2P-2uP_q-3B0+u²P_qq+2uB0_q]
         -C0(29/15)ln r.

Thus all feedback derivatives can be obtained from the same determinant
series through its second mass derivative. A tail estimate for W and Q alone
does not automatically control this second derivative. If a remainder bound
is built by term ratios, the extra (k-1) factor in the second derivative must
be included in its majorant. Special-function arithmetic errors remain
separate from analytic truncation bounds.

## 4. Bulk, two residuals, and a stable metric equation

For the unchanged classical bulk, y now denotes yhat, and

    R''=-R(w²/4+U/6),       phi'=w,
    w'=U_phi-4qw,           q=R'/R,
    q²=z+w²/12-U/6,         z=1/R²,
    q'=-z-w²/3,             z'=-2zq.

Use the regular cone R~y, w~0 and the outward normal +partial_y. The selected
branch has q>0 near the endpoint. Its shooting coordinates may be
u_coord=(ell,y_b), ell=log10|phi_h-1| with phi_h<1.

Define S=sigma+gamma rho and T=sigma_phi+gamma j. The TWO residuals are

    E1=q-S/6,
    E2=w+T/2.

Stress and current are both evaluated at the new endpoint. Even for b=0,
where j=0, rho retains curvature feedback. Simultaneously imposing the two
residuals and the radial constraint gives the exact endpoint identity

    z=S²/36-T²/48+U/6.

This constraint alone does not impose regularity at the cone or select the
radial profile. In particular it cannot replace the two-dimensional shooting
problem.

Direct subtraction in E1 can lose a residual below floating-point accuracy.
For the positive-q branch use the exactly equivalent residual

    D=q²-S²/36=z+w²/12-U/6-S²/36,
    E1=D/(q+S/6).

Do not evaluate D by subtracting the separate O(B²) terms. With eta=phi-1,

    B=1/3+eta²+eta³/3,     B_phi=2eta+eta²,
    tau=delta(1+c phi)+gamma rho,

the exact stable expression is

    D=z+(w-B_phi)(w+B_phi)/12-B tau/9-tau²/36.

The vacuum B²/9 terms have canceled algebraically. This retains the small
positive curvature against detuning-size terms, instead of subtracting two
extrinsic curvatures of order unity. Use eta internally so the tiny regular
cone displacement is represented rather than rounded out of phi=1. Still
record arithmetic and integration uncertainty; algebraic conditioning is not
an error certificate.

Squaring introduces another algebraic branch. Require q>0, S>0, and a
denominator q+S/6 bounded away from zero in the registered domain. A solution
of D=0 with q=-S/6 is not a solution of the selected unsquared junction.

## 5. Full two-coordinate Jacobian and endpoint variations

For a coordinate index a in {ell,y_b}, let phi_a,z_a,q_a,w_a be total endpoint
derivatives at fixed model parameters gamma,b,r,phi_b0. Then

    rho_a=rho_z z_a+rho_phi phi_a,
    j_a=j_z z_a+j_phi phi_a,

    M_1a=q_a-sigma_phi phi_a/6-gamma rho_a/6,
    M_2a=w_a+sigma_phiphi phi_a/2+gamma j_a/2.

Every source derivative is a partial derivative of the SAME fixed-reference
action. In particular one must neither freeze rho as phi changes nor replace
j_a by a density derivative.

The moving endpoint column is fixed by the bulk equations:

    phi_y=w, z_y=-2zq, q_y=-z-w²/3,
    w_y=U_phi-4qw.

Hence, writing E2vac=w+sigma_phi/2,

    M_1y=-z-w E2vac/3
          -gamma[-2zq rho_z+w rho_phi]/6,
    M_2y=U_phi-4qw+sigma_phiphi w/2
          +gamma[-2zq j_z+w j_phi]/2.

At a sourced solution, pairing simplifies the first expression to

    M_1y=-z+gamma z q rho_z/3+gamma z w j_z/12.

For the ell column let f=phi_ell, k=R_ell/R, v=q_ell, and define
I(y)=integral_0^y R R_ell dy. The unchanged bulk identity gives

    v+w f/3=-2I/R^4,
    (E1vac)_ell=-2I/R^4-E2vac f/3,
    M_1ell=-2I/R^4-E2vac f/3-gamma rho_ell/6,
    z_ell=-2zk.

This resolves the delicate classical metric derivative without subtracting
nearly equal endpoint sensitivities. Do not put E2vac=0 on a sourced root;
there E2vac=-gamma j/2.

If D is used as the first numerical residual, its Jacobian is

    D_a=(q+S/6) M_1a+E1(q_a+S_a/6),
    S_a=sigma_phi phi_a+gamma rho_a.

At E1=0 this reduces to D_a=2q M_1a. The off-shell term matters in Newton
iteration. Alternatively differentiate the factored D directly. The second
row remains M_2a. A condition number or determinant must say whether it
refers to (D,E2) or (E1,E2), since row scaling changes those diagnostics.

The endpoint constraint differential, including the varying quantum sources,
is

    [1-gamma(S rho_z/18-T j_z/24)] dz
      =[S(sigma_phi+gamma rho_phi)/18
         -T(sigma_phiphi+gamma j_phi)/24+U_phi/6] dphi
        +(S rho/18-T j/24) dgamma.

This is an additional on-shell check; it is not an independent closing
equation. Holding rho and j constant here discards finite-amplitude feedback.

## 6. Curvature response, existence, and what is not inferred

At gamma=0 the chosen mass ratio is exactly x=r=2z0, where the independently
derived sources simplify to

    Q0=z0/(12pi²),
    rho0=11z0²/(960pi²),
    j0=b z0²/(12pi²),       j0/rho0=80b/11.

These are exact symbolic reference values; no numerical source is evaluated
here. If the exact classical residual map is C1 and its exact Jacobian M0 is
nonsingular, the implicit function theorem gives a locally unique stationary
branch u_coord(gamma) in a neighborhood of zero, with

    du_coord/dgamma|0=M0^-1 (rho0/6,-j0/2)^T.

Thus the previously resolved supplied-source susceptibilities may be used
directly at first order. Write A_rho,A_j for the classical endpoint curvature
responses to independently supplied dimensionless (s,t), and B_rho,B_j for
the scalar responses. Then

    z(gamma)=z0+gamma z0²/pi² [11 A_rho/960+b A_j/12]+O(gamma²),
    phi_b(gamma)=phi_b0+gamma z0²/pi² [11 B_rho/960+b B_j/12]
                   +O(gamma²).

The action-derived vacuum curvature correction scales as z0² at the reference
point. The exact classical A_j is small but resolved in the inherited
sensitivity work; it is not identically zero. Finite-amplitude continuation
must use the full M above. Iterating the one-loop source self-consistently
resums its feedback; it does not supply all terms of a true two-loop theory.

The mathematical IFT conclusion is conditional until exact nonsingularity
and a neighborhood of validity are established. Floating-point roots and
convergence envelopes provide numerical evidence, not an existence theorem.
Certification would require validated cone initialization and ODE integration,
interval bounds for all source derivatives and radial sensitivities, and an
interval-Newton/Krawczyk or equivalent inclusion with an explicit domain.
An analytic small-delta argument would also need error bounds at the selected
finite delta. A finite determinant/condition number alone certifies neither
the exact residual map nor its neighborhood.

The quadratic x(phi) is an INPUT. Its flat local potential is exactly

    V=r²[a^4 ln(a²)-(3/2)a^4+2a²-1/2]/(64pi²),
    F=r[a² ln(a²)-a²+1]/(192pi²).

These functions and the geometric/nonlocal action contributions show why a
scalar-dependent mass cannot be treated as a constant supplied density after
the endpoint moves. They do not demonstrate that quantum corrections select
an exponential or quadratic mass law, relax vacuum energy, connect an unstable branch to
the corrected branch, produce particles, or heat radiation. Here p=-rho and
phi is constant on the shell. There is no time-dependent production,
energy-transfer, thermalization, causal-response, or attraction calculation.

## 7. Uniform domain and prospective amplitude protocol

A useful declared compact trial domain is

    3z0/4 <= z <= 5z0/4,
    |phi-phi_b0| <= dphi <= 1/2,
    b in {-1,0,+1}, gamma>=0,
    q>=q_min>0, S>=S_min>0.

It implies, uniformly,

    r(1-dphi/2)² <= x <= r(1+dphi/2)²,
    (8/5)(1-dphi/2)² <= x/z <= (8/3)(1+dphi/2)²,
    9/10 <= x/z <= 25/6 < 9/2.

The positive gap and fixed positive r ensure uniform smoothness of the action
and its differentiated sources on compact subsets. The last bound keeps the
inherited unshifted determinant series strictly inside its convergence disk
|x/z-9/4|<9/4. For uniform quantitative bounds, use the actual extremal ratio
to register theta_max<1. The allowed profile domain also needs regular R>0
away from the cone and a smooth positive-q choice; endpoint bounds alone do
not control all bulk profiles.

For b!=0, the excluded quadratic mass zero also gives a precise obstruction
to extending this Euclidean-vacuum family at fixed z>0. The inherited
homogeneous mode obeys Q~3z²/(8pi²x). With d=phi-phi_chi,

    x=r b² d²/4,    x_phi=r b² d/2,
    j=x_phi Q/2 ~ 3z²/[8pi²(phi-phi_chi)].

The leading current diverges with a coefficient independent of r and b.
The verifier checks this coefficient exactly. The zero of x_phi at phi_chi
does not cure the Euclidean zero-mode obstruction because Q diverges faster.
For finite positive gamma, fixed z>0 and finite classical junction data,
this singular source prevents smooth continuation through the mass zero
within the same stationary state prescription. This is a local mass-zero
obstruction; it makes no assertion about a simultaneous z->0 limit or a
different quantum state, and no heating process is inferred.

A prospective registration can use one geometric positive-amplitude grid
for every b, together with gamma=0 and predetermined precision/step refinements.
The declared grid is gamma in {0,10^-2,1,10²,10⁴,10⁶}; freeze it BEFORE
evaluation. These are mathematical amplitudes, not calibrated physical
couplings. No adaptive mass,
reference, b, or grid tuning to obtain a desired sign is permitted. Stop a
branch if a registered domain, positive-gap, denominator, convergence, or
response gate fails; retain the failure rather than silently extrapolating.

The final trial domain requires |phi_b-phi_b0|/|phi_b0-1|<=1/10 in addition to
3/4<=z/z0<=5/4 and the strict determinant disk 0<x/z<9/2. Register also the permissible
departure from the leading gamma response, and source/derivative precision
tolerances. The accuracy needed to resolve a correction must be set relative
to that correction, not only to the background residual. A grid may pass
the mathematical domain while failing a chosen perturbative-response gate.

Physical EFT validity is separate. If a cutoff Lambda_phys is declared, check
H/Lambda_phys and m/Lambda_phys, the appropriate bulk curvature and scalar
gradient scales, the fixed higher-operator truncation, and gravitational loop
counting in a specified N and kappa5² model. The crude combination
N kappa5² E³ is dimensionless but is not by itself a derived universal cutoff
or a full error bound. Without this additional matching there is no physical
EFT-validity claim, regardless of how well the stationary solver converges.

For the stated N=1 family, a conservative power-counting diagnostic can use
M5 L=gamma^(-1/3) and

    kT=U/6-w²/12,       kN=U/6+w²/4,
    Riemann5² L^4=24 kT²+16 kN²,
    Ehat=max[sqrt(z_b),sqrt(x_b),|q_b|,
             sup_bulk sqrt(|kT|),sup_bulk sqrt(|kN|)],
    epsilon_G=(Ehat/(M5 L))³=gamma Ehat³.

The frozen tenfold hierarchy diagnostic is Ehat/(M5 L)<=1/10, equivalently
epsilon_G<=1/1000, and is labeled a power-counting warning. The production
implementation's bulk maxima are sampled diagnostic maxima, not rigorous
supremum enclosures. This is not a proved quantum-gravity cutoff, species
bound, or error estimate. The bulk supremum refers to regular intrinsic
sectional curvatures; use the shell extrinsic curvature q_b, not
sup_bulk|q|. At the smooth cone q~1/y diverges purely because the de Sitter
slices contract, while kT and kN are finite. Treating that coordinate
extrinsic-curvature divergence as physical would reject every regular cone.

## 8. Why intrinsic higher derivatives do not add a radial equation here

For a fully covariant intrinsic shell functional in an invariant state,
symmetry forces its first metric variation on de Sitter to be proportional
to h_munu. Varying the S4 radius, including the volume, therefore determines
that single coefficient exactly. Uniform scalar variation determines j.
This justifies using W-zW_z/2 and x_phi W_x on the stationary family; it does
not identify W itself with a stress.

Intrinsic higher-curvature terms ordinarily give higher tangential derivatives
in general time-dependent shell equations. On exact de Sitter with constant
phi these derivatives vanish. Constant R4² has zero metric first variation
there, Weyl² vanishes, and constant Euler is topological, but scalar-dependent
coefficients still have scalar variations. For example a Euclidean term
f(x) z² gives zero rho and Q=2 f_x z²; dropping it from the action solely
because its metric source is zero would lose a real current. All such loop
terms are already included in the complete inherited determinant.

The determinant is intrinsic to the shell and introduces no normal bulk
derivatives or extra radial bulk degrees of freedom in this stationary
boundary-value problem. The classical Einstein-scalar bulk equations remain
second order, with the same regular-cone data. The induced metric and scalar
variation modifies their junction conditions. Varying y_b changes z and phi
and is consistently included by the total endpoint derivatives above; it
does not license varying r or differentiating the mass law during a pure
metric variation. The scalar and gravitational bulk constraints also enforce
the usual dependent normal-displacement condition on a glued shell.
Exact de Sitter invariance also forces the projected Weyl tensor, which is
invariant and traceless, to vanish in this ansatz. There is no independent
dark-radiation integration constant to vary inside this stationary family.

Choosing a different finite reference or adding matched local operators can
change the stationary solution; it is not a scheme-independent prediction
without a corresponding physical matching prescription. One can re-express
the same physical action in a different finite scheme only by compensating
the local couplings. Setting r equal to a new endpoint mass or curvature
during a variation instead changes the action and loses the reference
derivatives. The radius variation above must therefore keep r fixed even
when it happened to equal x or 2z at the original reference endpoint.

This argument depends on the declared intrinsic action and symmetric state.
An extrinsic-curvature operator, quantum bulk field, noninvariant state, or
off-de-Sitter perturbation would require additional equations or data. The
sphere-restricted W does not determine the full causal Hessian needed for
time-dependent stability; no such conclusion is claimed here.

## 9. Exact symbolic verification

`verify_stationary_model.py` is a portable SymPy verifier with an explicit
output path. It checks the dimensional bookkeeping, fixed-reference common
action derivatives, quadratic mass current, determinant Hessian, both
junctions, the stable metric residual, full feedback Jacobian, weighted bulk
identity, endpoint constraint differential, and local curvature variations.
Its deliberate wrong-formula controls include moving the subtraction
reference, removing curvature variation or current terms, losing the second
mass derivative, omitting endpoint motion, accepting the wrong squared root,
and double-counting local terms. All checks use explicit exceptions and
remain active under `python -O`. It performs no numerical source evaluation,
radial solve, or root search. Its execution report is evidence about the
displayed symbolic algebra, not about numerical or physical existence.

The final quadratic-model verifier passed **51 exact identities and detected
20 wrong-formula controls**, both normally and under `python -O`. The files
`FINAL_SYMBOLIC_CHECKS.json` and
`FINAL_SYMBOLIC_CHECKS_OPTIMIZED.json` record source hashes and runtime
versions. These are pre-registration symbolic algebra checks. During verifier
development, an initial check mistakenly substituted the differentiated
constraint for f' instead of w f'; the explicit exception caught that error.
The final check substitutes a polynomial dummy for w f', avoiding division
by w and remaining valid at w=0. No source calculation or radial run was
performed during that correction. That failure occurred during an earlier
unregistered exponential-mass model draft; it is a symbolic-verifier
development error, not a failed registered numerical experiment. The final
quadratic model and its exact verifier are the sole model specification here.

    python verify_stationary_model.py --output /tmp/stationary-symbolic.json
    python -O verify_stationary_model.py --output /tmp/stationary-symbolic-O.json

Choose unused output filenames on replay. Python and SymPy suffice; the
script does not import the repository's radial or quantum numerical code.
