# Smooth ambient reconstruction and actual local Einstein–scalar cone flux

13 October 2026. Independent continuation calculation. All files are confined
to this scratch directory. The prior physical-domain capsule is an unchanged
input. This is a sufficient local solution family, not a completeness theorem.

## 1. Result and hypotheses

There is an exact regularization of the scalar master equation that removes
its apparent cone singularity. In five regular Cartesian coordinates
`X^A=(mathsf T,X^1,...,X^4)`, with flat metric `eta_AB`, put

    s=eta_AB X^A X^B,     E=X^A partial_A,
    g_AB=Omega(s)^2 eta_AB,    phi=Phi(s),    Omega(0)=1.

Assume the background is an actual analytic Einstein–scalar solution on the
causal neighborhood considered, `Omega>0`, and `Phi_s!=0`. The registered
analytic regular-cone branch admits this local continuation: its regular
even/odd expansions give analytic `Omega(s),Phi(s)`; the background equations,
valid for positive s, continue analytically across s=0. A merely arbitrary
smooth continuation of expansion coefficients would not establish that it
solves the background equations. The local analytic background continuation
is an explicit hypothesis here; the reconstructed perturbations themselves
need only be smooth. A more general nonanalytic background would require
the corresponding direct ambient operator calculation.

Define

    L(s)=Omega_s/Omega + 2 Phi_ss/Phi_s,
    mathscr L F = Box_eta F - 2 L(s) E F
                  - [8 L(s)+(8/3)s Phi_s^2] F.                  (1)

This is a scalar normally hyperbolic operator with smooth lower-order
coefficients. For **every** smooth solution `mathscr L F=0`, set

    h_AB = 2s F eta_AB - 6 F X_A X_B,
    chi  = -3(E+4)F / [2 Omega^2 Phi_s].                       (2)

Then `(h_AB,chi)` is a smooth solution of the complete linearized bulk
Einstein–scalar equations, including at every point of the cone and at its
vertex. No spectral separation or inversion of `Box_gamma` or
`Box_gamma+4` is used. Smooth compact Cauchy data for (1) therefore produce
actual smooth finite-energy local coupled perturbations.

An explicit choice of such data below gives two perturbations that are
both non-gauge, with nonzero null presymplectic flux through a small
positive cone segment. This advances the prior conditional null-asymptotic
identity to an existence statement for this local class. Compact support
also supplies a finite-time shell-compatibility construction before the
support reaches the shell. It does not solve the unrestricted shell
initial-boundary problem or prove global stability.

## 2. Exact regularization in the exterior patch

On `s>0`, write `r=sqrt(s)` and `X^A=r n^A`, with `n^A` on unit dS4.
The conformal radius satisfies

    dy=Omega dr,    rho=r Omega,    p=phi_y=2r Phi_s/Omega.

For the registered scalar master set

    alpha=(rho^2 p^2)^-1,    w=(rho^4 p^2)^-1,
    q=2/(3rho^2),           P=rho^2 Psi.

The master equation is equivalent to `M[P]=0`, where

    M[P]=P_yy-(2a+2p_y/p)P_y
          +rho^-2(Box_gamma+4)P-(2/3)p^2 P,
    a=rho_y/rho.                                             (3)

Choose **exactly** `P=r^4 F`. This uses the conformal radius, not a claim
that `r=y` away from the cone. Since

    Box_eta = partial_r^2+(4/r)partial_r+r^-2 Box_gamma,

direct differentiation yields the exact identity

    M[r^4 F]=Omega^-2 r^4 mathscr L F.                         (4)

The two potentially singular algebraic terms cancel:
`(r^4)''-4(r^4)'/r+4r^4/r^2=0`. The first-derivative coefficient becomes
`4/r`, precisely the radial coefficient of the five-dimensional ambient
wave operator. The choice `P=r^3 F` fails this cancellation.

The longitudinal fields satisfy

    Psi=P/rho^2 = r^2 F/Omega^2,
    h=2Psi g-6Psi dy^2,
    chi=-3P_y/(rho^2 p).

Using `dy=Omega X_A dX^A/r` gives (2) exactly. The apparent division by
`p~2Phi_s(0)y` has disappeared. The resulting tensor and scalar are smooth
in the **full ambient topology**, including near-null approaches to the
vertex, because (2) uses polynomials in `X` and smooth functions of `s`.
No bounded-de-Sitter-time qualification is needed for this smoothness claim.

## 3. Full bulk residuals without a harmonic assumption

Take the field equations in trace-reversed form

    E_AB = R_AB - phi_A phi_B - (2/3) U(phi) g_AB = 0,
    S    = Box_g phi - U_phi = 0.                             (5)

The background equations are

    p_y=U_phi-4ap,     a_y=-rho^-2-p^2/3,
    a^2=rho^-2+p^2/12-U/6.

For arbitrary functions `Psi(y,x)` and `chi(y,x)`, choose
`h_yy=-4Psi`, `h_y mu=0`, `h_mu nu=2rho^2 Psi gamma_mu nu`. The exact
linearized Ricci components are

    delta R_yy = -4Psi_yy-16aPsi_y+2rho^-2 Box_gamma Psi,
    delta R_y mu = -3 D_mu(Psi_y+2aPsi),
    delta R_mu nu = gamma_mu nu {
      -Box_gamma Psi-rho^2[Psi_yy+10aPsi_y
                            +6(a_y+4a^2)Psi]}.                (6)

These are equations for arbitrary spacetime functions, not coefficients
identified only for a separated harmonic. For example, the tangential
trace-free Hessian terms from the lapse and conformal slice metric cancel
because `N=-2Psi` was imposed in the explicit reconstruction. No tensor
harmonic has been divided out.

The scalar equation variation is

    delta S = chi_yy+4a chi_y+rho^-2 Box_gamma chi-U_phiphi chi
               +6pPsi_y+4U_phi Psi.                          (7)

Now impose the reconstruction `chi=-3P_y/(rho^2p)`, which implies the
identity `Psi_y+2aPsi=-p chi/3`. Put

    C=p chi_y-p_y chi-3rho^-2(Box_gamma+4)Psi+2p^2Psi.

Substitution in (6), using only the displayed background equations, gives

    delta E_y mu=0,
    delta E_yy=-2C/3,
    delta E_mu nu=(rho^2/3)gamma_mu nu C,
    C=-3 M[P]/rho^2.                                         (8)

Direct differentiation in (7) gives the independent scalar identity

    p delta S = (partial_y+4a)C,
    delta S = -3/(rho^2p) (partial_y+2a)M[P].                  (9)

It is not necessary to rely only on a Bianchi implication for the scalar.
Combining (4), (8) and (9) gives the particularly transparent ambient
factorization

    delta E_AB = Omega^-2(3X_A X_B-s eta_AB) mathscr L F,
    delta S = -3/[2Omega^4 Phi_s] (E+6)mathscr L F.            (10)

For the registered analytic background, (10) is an equality of linear
differential operators with analytic coefficient functions. The coefficient
identities established on `s>0` continue analytically to the connected
ambient neighborhood. They therefore act identically on the jets of every
smooth F there, including the open region `s<0`. Continuity at the cone
alone would not justify an identity on that negative-s region. Equation (10)
proves sufficiency of (1) for all reconstructed fields. At the vertex the
Einstein prefactor vanishes, but this causes no gap: the reconstructed
solution has `mathscr L F=0` on the neighborhood, and the scalar identity
holds smoothly there as well.

This argument must not be reversed into completeness. Other perturbations
may not have this representative, and exceptional harmonics can have
additional solutions not described by a generic master reduction.
Constant and Hessian-degenerate subspaces were exceptions to deriving
necessary equations by harmonic division; they do not invalidate a direct
sufficient solution map verified by (10). The entire mass values zero
and minus four have not been removed.

## 4. A gauge-invariant physical diagnostic

Let

    K(s)=(nabla phi)^2=4s Phi_s^2/Omega^2.

Because `Phi_s!=0`, `K` is locally a smooth function of the background
scalar value. The first variation of the identically zero relational
background scalar `K_actual-K_background(phi_actual)` is

    Z=delta[(nabla phi)^2] - (K_s/Phi_s) chi.                 (11)

Under any smooth infinitesimal diffeomorphism the two variations are
`xi.K` and `(K_s/Phi_s)xi.Phi`, so they cancel. A nonzero Z is therefore
a sufficient obstruction to being pure gauge. For (2), direct contraction
gives

    Z=-6 Omega^-4 {
       [E-1-2sL](E+4)F-(8/3)s^2 Phi_s^2 F}.                 (12)

In null coordinates `u=mathsf R-mathsf T`, `v=mathsf R+mathsf T`, the
future cone is `u=0`, `v>0`. With `K0=Phi_s(0)`, its exact restriction is

    chi_cone=-3(v partial_v+4)F_cone/(2K0),
    Z_cone=4K0(v partial_v chi_cone-chi_cone).                (13)

At the vertex `chi(0)=-6F(0)/K0` is itself invariant under smooth gauge
transformations, since the background scalar gradient vanishes there.
Unlike a local claim that some fields have diffeomorphism form, (11)–(13)
prove physical nontriviality of the explicitly selected solutions below.

## 5. Smooth Cauchy data and a non-gauge flux pair

The standard smooth Cauchy theorem for a linear normally hyperbolic operator
with smooth lower-order coefficients applies to (1), for example Bär,
Ginoux and Pfäffle, *Wave Equations on Lorentzian Manifolds and Quantization*,
arXiv:0806.1036v1, Theorem 3.2.11 (printed pp. 85–86). On a small ambient
spacelike slice `mathsf T=0`, choose smooth compact data. If necessary,
extend the coefficients smoothly outside the physical causal neighborhood
to pose a usual Cauchy problem, then use finite propagation speed to
restrict to the region in which the coefficients are the actual background.
The characteristic cone is that of `eta`; all local data-to-field operations
in (2) preserve support. The reconstructed fields and their initial
derivatives are smooth and compact, have finite local Sobolev Cauchy energy
in a smooth inertial frame, and satisfy all linearized constraints by (10).
Finiteness here does not assert a positive conserved gravitational
Hamiltonian in a background with no chosen global timelike symmetry.

Here is an explicit two-parameter construction. Let a smooth radial outer
cutoff be identically one near the origin. Multiply the following data by
that cutoff, and write `L0=L(0)`:

    F1|T=0 = c=-K0/6,        partial_T F1|T=0=0,
    F2|T=0 = 0,              partial_T F2|T=0=d+e mathsf R^2,
    d=-2K0/15,               lambda=7e-5L0 d !=0.              (14)

Rotational symmetry and time reversal of the operator imply that F1 is
even in time and F2 is odd in time. The Taylor expansion of the wave
equation at the origin fixes

    F2 = d T + e T mathsf R^2
           + (4e-5L0 d)T^3/3 + O(|X|^5),
    F2_cone = d v/2 + (7e-5L0 d)v^3/24 + O(v^5).             (15)

Here `|X|` denotes an auxiliary Euclidean coordinate norm used for Taylor
estimates, not the Lorentzian radius. The coefficient follows by equating
the degree-one terms: `Box_eta(e T R^2+f T^3)=(8e-6f)T`, whereas the
lower-order terms on `dT` give `-10L0dT`. Smoothness supplies the needed
differentiated remainders. Equations (13)–(15) give

    chi1_cone = 1+O(v^2),
    chi2_cone = v/2 - 7lambda v^3/(16K0) + O(v^5),
    Z2_cone = -(7/2)lambda v^3 + O(v^5).                     (16)

The first perturbation is non-gauge because `chi1(0)=1`; equivalently
`Z1_cone=-4K0+O(v^2)` is nonzero on a sufficiently small positive segment.
The second is non-gauge on sufficiently small positive cone segments
because `lambda!=0` and (16). The cubic adjustment is intentional: the
leading linear cone scalar alone could agree with a translation-type
gauge variation and would not prove physicality of the second solution.

The prior action-derived Lee–Wald current applies directly to this
reconstruction: its pullback used the longitudinal field form and momentum
identity, both satisfied here independently of harmonic classification.
For one interior copy its exact cone restriction is

    omega^u_cone = (v^3/8) sqrt(gamma_S3)
       [chi1_cone partial_v chi2_cone
         -chi2_cone partial_v chi1_cone].                    (17)

The metric–scalar mixing term vanishes in this limit as established by
direct reconstruction; (2) supplies its uniform smooth asymptotics.
The Wronskian in (17) is `1/2+O(v^2)`. There is therefore an `epsilon0>0`
such that it is strictly positive on `0<v<epsilon0`. For any
`0<a<b<epsilon0`, the exact integrated +u flux is nonzero, with expansion

    integral_a^b omega^u_cone dv dOmega
       = Vol(S3) (b^4-a^4)/64 + O(b^6-a^6).                 (18)

The remainder has a uniform small-segment bound from the smooth Taylor
expansions. The sign reverses for outward flux from the exterior wedge.
This is a nonzero **antisymmetric symplectic pairing** of two actual,
non-gauge linearized solutions. It is not a positive-energy theorem, an
instability, or a new source of cosmological radiation.

Unlike the previous flat test-scalar comparator, these are coupled
Einstein–scalar perturbations of the specified local extended background.
Their small-cone flux resembles that comparator because the regular
normalization is Minkowskian at s=0; the exact equations and nonlinear
background coefficients have not been replaced by flat-space equations.

## 6. Finite-time compatibility with a distant shell

Suppose the registered background functions and nonvanishing `Phi_s` cover
the causal region from the data support up to a collar inside the shell.
The shell sits at fixed conformal radius `r_b`, hence

    mathsf R_shell(T)=sqrt(r_b^2+T^2).

Choose compact Cauchy data supported in `mathsf R<=R0<r_b`, and a time
interval `|T|<tau` with `R0+tau<r_b`. Finite propagation speed for (1)
implies support in `mathsf R<=R0+|T|`. Consequently F, h, chi and all
their derivatives vanish on a neighborhood of the shell throughout that
interval. The shell displacement can be set to zero. Induced-field
continuity and both linearized junctions are then identically satisfied.
Two copies may carry equal or opposite data without altering this argument.

This constructs a finite-time doubled-shell perturbation while its support
is causally separated from the shell, if the stated coefficient region is
available. It does not require a boundary condition for the master in the
auxiliary extended Cauchy problem; finite propagation isolates the arbitrary
auxiliary extension. It proves neither all-time shell reflection nor a
well-posed full shell initial-boundary problem after radiation arrives.

The small cone segment in (18) can be chosen with `b/2<tau` and within the
coefficient neighborhood, so the nonzero flux occurs during this finite
time window. Retaining the smooth region beyond the cone interprets it
as internal flux between subregions; treating the cone as a physical
boundary instead requires a separate choice of boundary data and charges.

## 7. Verification, inherited input and limits

`check_ambient_reconstruction.py` computes all 15 symmetric components of
the Christoffel-derived linearized Ricci tensor in a flat de Sitter chart
for arbitrary `Psi(y,t,x1,x2,x3)`, as well as the full scalar wave variation.
It then checks (8)–(10), the exact master and field transformations, the
relational invariant, Cauchy jets, and the nonzero cone Wronskian. The Ricci
formulas are covariant in the de Sitter indices, so the chart calculation
checks their coefficients rather than restricting the solution to a
homogeneous or separated mode. Deliberately wrong terms, including a
corrupted radial power and an omitted drift, are rejected.

Normal and optimized runs and injected-failure behavior are recorded in
the accompanying receipt. Symbolic checks verify local differential
identities and finite Taylor coefficients. The smooth normally hyperbolic
Cauchy theorem, finite propagation, analytic continuation of the background
and the stated scope arguments remain mathematical proofs in this note,
not mechanically certified PDE computations.

The unchanged input is
`HDblast/research/HDBLAST_CONE_PHASE_SPACE_20261013`, especially its
manuscript equations (4), (8)–(10), (20)–(24), the covariant-current note,
and the prior-checkpoint operator audit. Its radial positivity and threshold
claims are unchanged. We do not transfer the historical sixteen-field
normal form on the different near-phi=-1 branch; it is unnecessary for
the independent derivation above.

The new result is sufficient smooth local reconstruction and actual local
nonzero flux, including a finite-time support argument for an untouched
shell. It is not a full physical phase-space classification, a uniqueness
theorem for gauge representatives, a new global spectral threshold,
nonlinear stability, long-time radiation balance, or evidence for a
higher-dimensional origin of the Big Bang. Prior-art comparison of this
specific construction is still required before any novelty claim.
