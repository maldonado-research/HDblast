# Quadrupole data and the ten quadratic symmetry tests

13 October 2026 (America/Los_Angeles). New analytical continuation; the
previous frozen shell-contact capsule remains unchanged. This note evaluates
the necessary functionals derived in the accompanying nonlinear-shell and
rotational-constraint proofs. Vanishing of these functionals is not a proof
of second-order solvability or nonlinear integrability.

## 1. Fixed data space

Fix the actual nonconstant regular background and an admissible sufficiently
thin annulus from the shell-contact theorem. Choose a nonzero real
`Q ∈ C_c^∞(I)`, and keep one fixed reflection parity between the two interiors.
The prior radial estimate gives `K_rad(Q,Q)>0`, where angular integration has
been removed from the action-derived Cauchy bilinear. It is the same radial
functional for every angular harmonic in this construction.

Let `n=X/r ∈ S³`, and let `A,B` be real symmetric traceless four-by-four
matrices. The degree-two scalar harmonics are exactly

    H_A(n)=nᵀ A n,      H_B(n)=nᵀ B n.

Consider initial master data

    q=Q(r) H_A(n),      p=ω Q(r) H_B(n),                 (1)

with a fixed nonzero real constant `ω` (carrying inverse-time units if the
coordinates are dimensional). This is an eighteen-dimensional initial-data
space at fixed `Q,ω`, inside the previously constructed smooth finite-angular
linear solution family. It is preserved by spatial rotations. Boosts leave
this particular subspace, but their pairings can be evaluated in the enlarged
finite-angular solution family. No separate Lorentz-invariant phase space is
being assigned to (1).

## 2. Sphere integrals and the rotation action

Write `V₃=Vol(S³)=2π²`. Rotational invariance, or direct monomial integration,
gives

    ∫ n_i n_j n_k n_l dΩ = V₃(δij δkl+δik δjl+δil δjk)/24.

Consequently

    ∫ H_A H_B dΩ = (V₃/12) tr(AB).                      (2)

For an antisymmetric matrix `M`, use the vector field
`J_M=(Mn)·∂_n`. Direct differentiation gives

    J_M H_A = H_[A,M].                                 (3)

The matrix `[A,M]=AM−MA` is symmetric and traceless. Equations (2)–(3) imply

    ∫ H_B J_M H_A dΩ
       = (V₃/12) tr(B[A,M])
       = (V₃/12) tr(M[B,A]).                           (4)

These facts are standard harmonic and matrix algebra; no novelty claim is
made for them.

## 3. Necessary conditions become a commuting-matrix criterion

The exact nonlinear momentum-plus-Israel functional `C_J` vanishes on every
admissible exact coupled solution. The background is a critical point of this
functional. For the compact initial perturbations here, the independently
derived Hessian is

    D²C_J[v,v] = −Ω_phys(v,Jv)
                = 4 K(Jq,p),                           (5)

where the factor four includes the two interiors and the Hessian factor two.
The sign uses the published future Cauchy convention, opposite to the
canonical `δπ∧δγ` convention. Only the required zero is orientation independent.

Combining (1), (4), and (5) yields

    D²C_JM[v,v]
      = (V₃/3) ω K_rad(Q,Q) tr(B[A,M]).                 (6)

For the basis `J_ij=X_i∂_j−X_j∂_i`, use `M_ij=-1, M_ji=1`.
Then `tr(B[A,M])=-2[A,B]_ij`, so

    D²C_Jij[v,v]
      = −(2V₃/3) ω K_rad(Q,Q) [A,B]_ij.               (7)

Because `[A,B]` is antisymmetric, its six entries with `i<j` determine it.
The prefactor is nonzero. Thus all six rotational quadratic conditions
vanish **if and only if**

    [A,B]=0.                                           (8)

All four boost candidates vanish for every pair in (1). Both initial data
have even spatial antipodal parity. At `T=0`, a boost acts on master data as

    (q,p) ↦ (X_i p, ∂_i q+X_i A_wave q).

The background wave-acceleration operator `A_wave` has radial coefficients
and preserves antipodal parity. Both transformed components are odd. The
radial derivatives, weights, and complete sphere measure in the physical
pairing preserve parity, so each angular integral is zero. Compact support
also makes the initial shell cut term zero. The necessary-condition
identification for all ten common Lorentz symmetries uses the accompanying
nonlinear Noether proof; the explicit rotational obstruction already follows
from the independent momentum/Israel derivation.

Therefore, within the restricted data space (1), **all ten Lorentz quadratic
tests vanish exactly when `A` and `B` commute**. Real commuting symmetric
matrices can be simultaneously diagonalized by an orthogonal rotation.
This is an algebraic classification of these necessary tests only. It does
not characterize the full nonlinear constraint range, exclude additional
adjoint-domain issues, or prove any commuting pair is integrable.

## 4. Explicit nonintegrable candidate and the original zero family

Take

    A=diag(1,−1,0,0),
    J=J_12,
    B=[A,M_12],     B_12=B_21=−2,

with all other entries of B zero. Then

    H_A=n₁²−n₂²,     H_B=JH_A=−4n₁n₂,
    J²H_A=−4H_A,
    ∫H_A²=V₃/6,     ∫(JH_A)²=2V₃/3,
    [A,B]_12=−4.

Equation (7) becomes

    D²C_J[v,v] = (8V₃/3) ω K_rad(Q,Q)
                = (16π²/3) ω K_rad(Q,Q) ≠ 0.           (9)

Subject to the classical-family assumptions of the independent obstruction
theorem, this excludes a parameter-C² exact coupled solution family with
first derivative v. A second-order correction cannot cancel (9), since
`DC_J[w]=0` for every allowed correction, including corrections nonzero on
the initial shell cut. This conclusion does not assume that the second-order
correction belongs to the scalar master image.

The initial data generate actual smooth non-gauge linear solutions by the
prior finite-angular theorem. The displacement component and velocity
component project onto orthogonal nonzero harmonics `H_A` and `JH_A`.
Their separately established shell activations cannot cancel as complete
angular traces. Thus the example can be chosen among the contacting linear
solutions, rather than merely among formal unconstrained initial tensors.

The degree-two specialization of the previous same-real-harmonic
displacement/velocity span has `A=aC, B=bC` for one fixed real symmetric traceless C, including pure
displacement (`b=0`) and pure velocity (`a=0`). It therefore passes these ten
necessary tests. Its nonlinear extension is still
unproved. The obstruction shows why solving the linear equations alone does
not establish nonlinear admissibility for an arbitrary superposition.

## 5. Verification boundary

The companion checker directly differentiates arbitrary traceless quadratic
polynomials, integrates their monomials over S³ using exact even moments,
and compares the results with the matrix identities and the explicit witness.
It includes deliberately incorrect normalization/sign/cancellation formulas.
Those algebra checks do not certify the PDE existence theorem, the nonlinear
Noether identity, or the differentiable-family hypotheses. Their proofs and
input provenance are separate.

This is a restriction within a specified classical gravitational model. It
is not a dynamical instability, a stability theorem, a cosmological energy
source, or evidence establishing a higher-dimensional origin of the Big Bang.
