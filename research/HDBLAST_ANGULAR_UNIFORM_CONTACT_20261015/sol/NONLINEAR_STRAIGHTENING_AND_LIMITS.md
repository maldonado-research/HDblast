# Nonlinear field-space straightening: exact algebra and remaining theorem

15 October 2026. Independent assessment. The field-space transformation
below is exact. Its principal and boundary-flux calculations are exact.
This note does not claim that these calculations alone prove a finite
Sobolev nonlinear existence theorem. The prior checkpoint already supplies
an exact localized initial constraint family at each fixed high finite
regularity; its use in the shell-contact evolution remains a separate step.

## 1. The nonlinear boundary vector field

Let q_N=(g_yy,h_ab,phi), q_D=g_ya. At the fixed shell the exact equations
with q_D=0 are partial_y q_N=F(t,Omega,q_N), where

    F_yy = (4 sigma(phi)/3) N^3
                    - 2 a_b N^4 tr_h(h0),
    F_h  = (sigma(phi)/3) N h,
    F_phi= -(sigma'(phi)/2) N,
    N=sqrt(g_yy).

This is a smooth vector field on the open set g_yy>0 with Lorentzian h.
The background lies in that open set. Define the perturbation U_N by
subtracting q0,N. In inward z coordinates its exact boundary law is

    -partial_z U_N + B(t,Omega,U_N)=0,
    B(t,Omega,v)=-rho_b [F(t,Omega,q0,N+v)-F(t,Omega,q0,N)].

In particular B(t,Omega,0)=0 and B_v(0)=rho_b B_y=M, with the exact
linear matrix in the companion note. This sign follows from
partial_z=-rho_b partial_y and is not an independently adjustable choice.

## 2. A local flow which fixes the background

For each fixed (t,Omega), let Phi_l(v) solve

    dPhi_l(v)/dl=B(t,Omega,Phi_l(v)),     Phi_0(v)=v.

On a sufficiently small field neighborhood and for a sufficiently small
fixed collar parameter, this smooth ODE flow and its inverse exist.
Because B(0)=0, Phi_l(0)=0. Take the same radial chi as in the linear
construction and put

    U_D=V_D,           U_N=Psi(t,Omega,z,V_N)
                       :=Phi_chi(z)(V_N).

Psi is a smooth local diffeomorphism in v, is the identity at z=0 and
outside the collar, preserves the background, and has Psi_z(0,v)=B(v).
At z=0, its Jacobian J=Psi_v is the identity. Consequently

    -U_N,z+B(U_N)=-V_N,z.

More generally the chain-rule calculation gives the boundary law
-J V_N,z=0 if the transformation is normalized differently at the
boundary; J remains invertible. Thus the exact nonlinear shell boundary
operator becomes homogeneous Neumann, with homogeneous Dirichlet on V_D.

All this is pointwise field algebra. There is no inversion of a tangential
differential operator and no angular-frequency denominator. Parameters
(t,Omega) enter through the known smooth background and invariant tensor
operations. Their derivatives therefore do not consume extra derivatives
of the unknown. Derivatives of the Jacobian and its inverse are governed
by ordinary finite-dimensional flow derivatives.

## 3. The principal symbol survives

Normalize the metric and scalar reduced equations so their principal
part is the common scalar expression

    g^{mu nu}(q) partial_mu partial_nu q

times the identity on all sixteen components. Under U=Psi(V), second
derivatives of U equal

    J partial_mu partial_nu V
      + Psi_vv(partial_mu V,partial_nu V)
      + known-parameter/first-derivative terms.

Multiplication by J^{-1} leaves exactly

    g^{mu nu}(q0+Psi(V)) partial_mu partial_nu V

as principal part. Every new term is at most quadratic in first derivatives
of V with smooth field-dependent coefficients, together with known smooth
parameter terms. No new component-dependent principal wave operator appears.

Because q_D=g_ya vanishes on the boundary, its inverse metric satisfies
g^{ya}=0 there. Hence g^{zt}=g^{zi}=0 at z=0 for the transformed system.
This is an identity for every admissible boundary field, not just the
background. All time/angular derivatives of these coefficients also vanish
at the boundary.

## 4. The variable-principal flux is favorable

For a smooth solution sufficiently close to the background, g^{tt}<0,
g^{zz}>0, and the four-dimensional spatial matrix g^{IJ} is positive.
Using the same positive bundle metric, the standard multiplier V_t gives
the energy

    (1/2) integral [ -g^{tt}|V_t|^2
                        + g^{IJ}<V_I,V_J> + |V|^2 ].

Its normal flux is g^{z nu}<V_nu,V_t>. The mixed coefficients vanish at
the boundary. The remaining g^{zz}<V_z,V_t> vanishes by V_D,t=0 and
V_N,z=0. Derivatives of g enter the interior estimate only through
bounded first derivatives of the coefficients. This gives an H1 a priori
estimate for smooth solutions of any appropriately frozen variable-
principal problem, uniform with respect to angular approximation.

Tangential/time differentiation preserves the static mixed boundary
conditions. The vanishing g^{za} traces also persist under these
derivatives. Higher normal derivatives can be recovered from the equation
because g^{zz} stays bounded away from zero. For large integer spatial
regularity, for example H^s x H^{s-1} with s>=6 in four spatial dimensions,
the relevant first derivatives and low second derivatives are bounded,
and the usual Sobolev product estimates give the expected tame a priori
commutator bounds. This is evidence that an ordinary quasilinear scheme
can close without an angular-frequency obstruction.

## 5. What still has to be proved before a nonlinear theorem

The positive background linear estimate is a completed construction. The
following parts are not supplied by that construction or the algebra above:

1. A full variable-principal linear mixed-boundary existence proof with
   controlled higher regularity, including a compatible smooth approximation
   at the time/boundary corner. Angular equivariance is lost for a general
   nonlinear iterate, so the background isotypic limiting construction
   cannot simply be reused unchanged.
2. A quasilinear iteration or regularization scheme whose fixed initial
   data and modified/frozen equations retain the needed high-order corner
   identities. Freezing coefficients and retaining a source naively can
   change higher initial time jets and break compatibility.
3. A precise contraction/continuous-dependence estimate, a stated finite
   Sobolev class, and a continuation criterion. The smooth a priori flux
   calculation is not a substitute for these steps.
4. Carry the already available exact physical initial constraint family,
   its exact initial harmonic gauge completion, and an explicit parameter-
   regularity budget into the mixed-boundary construction. The previous
   first and second Taylor coefficients alone are not this input.
5. A fixed-slab estimate sufficient to reach the specified contact time.
   Conditional small-data persistence near a smooth background is plausible;
   a theorem must state the smallness and exact-compatible-data hypotheses
   and actually establish that persistence.

Once a sufficiently regular exact reduced solution is constructed, the
exact Bianchi identity, Codazzi identity and scalar junction give the
homogeneous subsidiary system with its normal Dirichlet and tangential
Robin conditions. The same local bundle argument then proves subsidiary
uniqueness. Zero exact initial constraints/gauge data imply C=0. This last
step is conditional on exact initial data and the constructed solution,
rather than on a formal Taylor series.

The prior input is
publish-exact-precontact-20261014/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md.
I read that note for this assessment. Under its reflection, annulus,
background and weighted-inverse hypotheses, it gives for each fixed k>=12
a parameter-smooth exact constraint family, metric in H^{k+2}, metric
momentum in H^{k+1}, with prescribed scalar sectors, tangent to the selected
reflection-invariant seed. Its corrections are supported in a fixed closed
annulus, and the fields equal the background on open shell and vertex
neighborhoods. It also gives an exact sufficiently short precontact family
with C2 parameter dependence in a lower piecewise spacetime C4 norm.
These achievements need not be reproved. They do not yet give mixed-boundary
evolution through contact, a common spatial C-infinity family, or an
identification of its second derivative with any separately chosen frozen
second-order correction.

There is no identified boundary-algebra or high-frequency counterexample
to the proposed nonlinear formulation. There is also no fully completed
finite-H^s existence/continuation argument in this independent note.
The honest checkpoint is the all-angular linear theorem, with an exact
nonlinear straightening mechanism and a concrete next proof obligation.
