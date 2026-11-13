# Evaluated all-mode weighted initial inverse

This continuation follows merged PR #138. It quantifies the selected
background constraint inverse on the entire reflection-invariant weighted
space, including arbitrarily high angular modes and both degenerate annulus
ends. The preceding estimates for the actual second source's finite
harmonic support did not provide this operator norm.

Keep the registered background, fixed canonical density `rho0`,
`ell=10^-6`, correction annulus `ell/2<z<5ell/2`, and

```
r=(z-ell/2)(5ell/2-z)/ell^2, kappa_w=ell/100=10^-8,
q=kappa_w^2 r^2, w=exp(-1/r), Gamma=(Z2)^4,
||u||_Hk(q,w)^2=sum_(j=0)^k integral w^2 q^(2j)|nabla_gamma0^j u|^2 dmu0.
```

The derivative multiplier is q, rather than its square root. Every original
and derivative tensor slot uses the fixed physical metric and full positive
Hilbert norm. The dot denotes compact-support closure in the displayed norm,
with no added ordinary endpoint trace on a potential.

| Selected inverse | Normalized source → potential | Evaluated upper bound |
| --- | --- | --- |
| Scalar `L_H^-1` | `dot H12(q,w) → dot H16(q,w)` | `10^5000` |
| Momentum `L_M^-1` | `dot H12(q,w) → dot H14(q,w)` | `10^2700` |
| Their Hilbert direct sum | Both normalized source blocks → both potentials | `10^5000` |

These are conditional, deliberately loose analytic ceilings in the fixed
model normalization. They are not measured effects or blast energies.
Existence, uniqueness and the smooth background retain the stated analytic
premises. The estimates quantify the same selected inverse, without
replacing it by a convenient right inverse or an archived York model.

## The complete operator and source convention

The actual normal operators are

```
L_H=w^-2 P_g w^2 q^4 P_g^*,
L_M=w^-2 P_m w^2 q^2 P_m^*.
```

Here `P_g` is the retained scalar metric-constraint variation. For momentum,
`P_m=-2 div_gamma0` and its adjoint is the ordinary Lie derivative of
`gamma0`, with the prescribed raised indices. It is a symmetric-gradient
operator. The canonical correction remains `delta pi=rho0 w^2 q^2 P_m^*Z`
with the original fixed density.

The source in the table is normalized. For raw scalar and canonical
momentum-density residuals, the pair is
`f=(w^-2 f_H,w^-2 f_M_density/rho0)`: the momentum source is first divided
by the fixed rho0. If `f_phys` already denotes that scalar/one-form pair,
then `f=w^-2 f_phys`. The exact radial unitary transform uses `ds=dy/q`,
`M=w sqrt(q rho^3)` and source `g=Mf`. Substituting the unnormalized
source into the weighted range norm would change the claim.

## How the all-mode estimate closes

[SCALAR_ALL_MODE.md](SCALAR_ALL_MODE.md) extends the scalar coercivity
argument uniformly to all allowed scalar degrees, then proves the H16/H12
estimate using `A=1-D_s^2+zeta^2`, `zeta=q sqrt(lambda)/rho`.
The unbounded angular frequency is part of the positive principal operator.
It is not bounded as a lower-order radial coefficient.

[MOMENTUM_ALL_MODE.md](MOMENTUM_ALL_MODE.md) includes all scalar/gradient
blocks and both coexact vector families. It proves a uniform weighted
energy gap and H1 Korn estimate. The six Killing directions have no
invariant part. [MOMENTUM_REGULARITY_ALL_MODE.md](MOMENTUM_REGULARITY_ALL_MODE.md)
then evaluates the full matrix graph estimate, commutator and spectral
absorption required for H14/H12. Mixed angular/radial principal terms and
the complete canonical pairing are retained.

[WEIGHTED_NORM_BRIDGE.md](WEIGHTED_NORM_BRIDGE.md) gives both directions of
the actual covariant weighted norm comparison through order sixteen,
including all derivative slots and the source. Its common explicit
constant is `13*8^16*16!<10^29`. Natural SO(4) Casimir identities handle
all angular modes without a harmonic-polynomial cutoff.

The low-order gap is squared in the normal inverse estimate by the energy
identity. A coercivity constant is not itself the inverse norm. All graph,
coefficient, connection, norm and domain constants used above are charged
in the four complete proofs.

## Verification and exact limits

[REPRODUCE.md](REPRODUCE.md) describes the authenticated exact replay and
failure controls. Algebraic oracles compare principal tensors, generic
row-energy operators and actual background polynomials with independently
expanded formulas; rational ledgers check the numerical allocations.
Normal and optimized Python reports must be identical. Finite checks
support the written analytic arguments and do not prove their continuous
existence or regularity premises.

[INTERNAL_REVIEW.md](INTERNAL_REVIEW.md) records the independent internal
review of all four proofs. [LITERATURE.md](LITERATURE.md) checks the exact
primary weighted-space theorem, source conventions and S3 spectrum, with
coverage and access limits. No unspecified constant from a primary theorem
is inserted into the numerical result. An apparent version-specific
formula mismatch in a cited Korn identity is documented with its corrected
substitution; the matching weighted isomorphism remains applicable under
its stated hypotheses.

This is a linear inverse at the selected background. It does not quantify
the nonlinear derivative defect on a positive state/parameter ball, an
admitted nonzero amplitude, a uniform nonlinear evolution tube or an
observable remainder. No physical-origin or observational confirmation
follows. [CONTINUATION.md](CONTINUATION.md) identifies the next input.
