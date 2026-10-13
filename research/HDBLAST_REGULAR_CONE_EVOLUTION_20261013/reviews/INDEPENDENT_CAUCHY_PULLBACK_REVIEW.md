# Independent review of the physical Cauchy pullback

13 October 2026. This is a bounded internal mathematical review by a separate
derivation process, not external peer review. Its target is section 6 of
`boundary-evolution/BOUNDARY_LAWS_AND_LOCAL_CAUCHY_PHASE_SPACE.md` and the
corresponding checker identities. Exact reviewed bodies are pinned in
`CAUCHY_PULLBACK_REVIEW_RECEIPT.json`.

**Verdict:** the one-copy physical pullback, its measure and conformal
factors, the integration-by-parts identity, and nondegeneracy on the stated
smooth compact-data image are correct. These arguments do not identify a
completed physical Hilbert space, a global gauge quotient, or the evolution
after the perturbation reaches the shell. The rest of the boundary note is
outside the scope of this review.

## Direct Cauchy calculation

Use `T=r sinh(t)`, `R=r cosh(t)`, `dy=Omega dr` and the prior action-derived
one-copy coordinate-density current

    j^t = -sqrt(gamma_S3) rho^2 cosh(t)^3
          [chi1 partial_t chi2-chi2 partial_t chi1
             +6(Psi1 partial_t Psi2-Psi2 partial_t Psi1)].

The absolute coordinate Jacobian is
`|partial(y,t)/partial(T,R)|=Omega/r`. At `T=0`,
`partial_y T=0`, `partial_t T=r`, so the ambient vector density is
`J^T=Omega j^t`. Also `partial_t=r partial_T` there. Choosing the
future-Cauchy pairing `sigma=-integral J^T dR dOmega` gives exactly

    sigma = integral Omega^3 r^3
       [chi1 chi2_T-chi2 chi1_T
         +6(Psi1 Psi2_T-Psi2 Psi1_T)] dr dOmega.                (1)

Reversing the global symplectic orientation changes the overall sign only;
it does not change the nondegeneracy result. No lapse or volume factor is
supplied by analogy with the auxiliary scalar action: (1) comes from the
Einstein–scalar current itself. The shell current vanishes on these
compactly supported initial data because the fields and embeddings vanish
on a shell neighborhood.

Let `q=F|T=0`, `p=F_T|T=0`, `D=r partial_r`, `f=Phi_s`, and
`b=-3/(2Omega^2 f)`. The exact reconstruction gives

    chi=b(D+4)q,       chi_T=b(D+5)p,
    Psi=r^2 q/Omega^2, Psi_T=r^2 p/Omega^2.                    (2)

The second formula is not obtained by simply changing q to p in the first:
`[partial_T,E]=partial_T` supplies its additional one. Derivatives of the
background coefficient b and of `r^2/Omega^2` vanish at T=0, since each is
a function of `s=R^2-T^2`.

Substitution in (1) produces the reviewed bilinear

    sigma((q1,p1),(q2,p2))=K(q1,p2)-K(q2,p1),
    K(q,p)=integral [m(D+4)q(D+5)p+6r^7 Omega^-1 qp] dr dOmega,
    m=9r^3/(4Omega f^2).                                      (3)

Both the scalar factor m and the gravitational coefficient agree with the
target note. Fixed-parity doubling multiplies (3) by two. Cross-parity
pairings cancel between the two copies.

## Positivity of the diagonal and the exact kernel argument

Since `L=Omega_s/Omega+2f_s/f`,

    (r m)'/m=4-2sL.

For real smooth compact data, the radial integration of `9qDq` yields

    K(q,q)=integral [m{(Dq)^2+(2+9sL)q^2}
                     +6r^7 Omega^-1 q^2] dr dOmega.            (4)

The outer boundary term is zero by compact support. At the vertex it is
`r m q^2=O(r^4)`, so it vanishes for smooth q. Angular dependence creates
no extra boundary term because the integration is over the closed S3 and
this calculation has no angular integration by parts.

Choose the support ball small enough that `2+9sL>=1`. The weights in (4)
are positive for `r>0`; hence `K(q,q)>0` for every nonzero smooth q in
the class. The bilinear K itself is not required to be symmetric. If a
pair `(q,p)` pairs to zero with every allowed compact test pair, use
`(0,q)` to get `K(q,q)=0`, hence q=0, then use `(p,0)` to get
`-K(p,p)=0`, hence p=0. These are admissible independent Cauchy data for
the normally hyperbolic master equation. This proves the stated
nondegeneracy on that image, rather than only positivity of a convenient
formal weight.

An ordinary smooth zero-charge diffeomorphism is a presymplectic degeneracy
when the complete Cauchy surface and shell contributions are included.
Thus no nonzero element of this particular reconstructed compact-data image
can be such a gauge direction. This inference assumes the complete
zero-charge gauge interpretation stated in the target. It does not classify
charged boundary symmetries, singular gauge generators, weaker completed
gauge groups, or physical solutions outside the image.

Finally, (4) is a diagonal auxiliary bilinear used to test the kernel of an
antisymmetric physical form. It is not a positive value of the symplectic
form on equal arguments and is not asserted to be the gravitational energy.
