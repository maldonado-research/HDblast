# Scalar domain theorem: result and verification scope

Dated continuation, 13 October 2026. The accompanying proof concerns an
abstract radial operator. It is not a proof of complete physical stability.

Fix one finite background interval `(0,T]`, with `rho>0`, nonvanishing
background momentum `p`, and regular-cone expansions
`rho~y`, `p~A_c y`, `A_c!=0`. Set

    alpha=1/(rho^2 p^2), q=2/(3rho^2), w=1/(rho^4 p^2).

Assume a positive reference `P0~C0 y^4` solves
`-(alpha P0')'+q P0=4w P0`, and assume `beta>0` with

    r0=alpha(T) P0'(T)/P0(T)-4/beta > 0.

The proof constructs the self-adjoint operator in `L2(w dy) plus C`, with
endpoint weight `1/beta`. The endpoint component equals the shell trace on
the form and operator domains. Its eigenvalue equation is

    -(alpha u')'+qu=Lambda w u,
    alpha(T)u'(T)=Lambda u(T)/beta,    mu^2=Lambda-4.

The cone is limit-point. The ground-state identity extends to the full closed
form domain. Define

    a_g=alpha P0^2,     m_g=w P0^2,
    alpha_*=inf a_g/y^4 > 0,     m_*=sup m_g/y^2 < infinity,
    M_g=integral m_g,   b_g=P0(T)^2/beta,   k_g=r0 P0(T)^2 > 0.

A proved weighted Hardy inequality yields

    spectrum(A-4I) is contained in [epsilon_delta,infinity),
    epsilon_delta=[4m_*/(9alpha_*)+(M_g+b_g)/k_g]^(-1) > 0.

This is an **explicit function of a fixed background**, not a numerically
evaluated or detuning-uniform lower bound. The result excludes all nonpositive
spectrum, including continuous spectrum, for this realization. Under the
stated differentiated cone expansions, its essential spectrum in `mu^2`
is exactly `[9/4,infinity)`. The proof also gives an ordinary-Robin companion
with no endpoint mass, conditional on the separately derived Robin sign.

The physical-state-space identification, exceptional harmonics and full
kinetic norm require the separate physical audit. No uniform dimensionless
gap, nonlinear stability or complete no-ghost theorem is claimed. Physical
masses scale with `H_b`, which tends to zero in the detuning limit.

## Portable files and exact checker

- `SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md`: complete proof, copied unchanged.
- `reviews/INDEPENDENT_DOMAIN_REVIEW.md`: independent internal review, copied
  unchanged and pinning the exact proof version.
- `check_liouville_and_gap_identities.py`: portable CLI adaptation with exact
  symbolic identities and deliberate wrong-formula controls.
- `provenance/`: unchanged original checker and receipt, plus separately
  labelled original-authoring pins. The original checker is historical
  evidence; use the CLI adaptation for portable replay.
- `PORTABLE_MANIFEST.json`: actual portable-file hashes and adaptation notes.

The checker requires Python and SymPy. From this directory, for example:

```sh
python check_liouville_and_gap_identities.py --output domain-normal.json
python -O check_liouville_and_gap_identities.py --output domain-optimized.json
```

The output argument may point outside the capsule. Identical environments
produce deterministic receipts, including the actual CLI script hash.
Correctness checks use explicit failures and remain active under `python -O`.
Deliberately incorrect algebra must be rejected for the run to pass.

The computations check algebraic identities and sensitivity to selected
wrong formulas. They are **not a formal verifier** of closure, self-adjointness,
Hardy coercivity, essential-spectrum arguments, or physical completeness.
Those claims are established and reviewed in the written proof.
