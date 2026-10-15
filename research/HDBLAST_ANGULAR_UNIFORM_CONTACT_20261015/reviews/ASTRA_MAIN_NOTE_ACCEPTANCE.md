# Astra review of the final frozen-linear checkpoint

**Accepted within its stated scope. No blocking mathematical finding.**

Reviewed main note SHA-256:
`cc559ed4a16c920310fca4c97c3af869508414b6d0ffc9fd8b4612e82989dff1`.

Reviewed exact nonlinear-flow companion SHA-256:
`da221f388e0397e4a3b06ce40be1feff1d170311f3b7eaf79f3644b667dd6882`.

The main proof uses the actual full boundary matrix. Its block-preserving
exponential removes the nonsymmetric Robin term, retains the common scalar
principal operator and gives zero flux in a positive split-orthogonal control
norm. Equation (7) is a sufficient high-Sobolev source norm: its initial jets
and time integrals bound the instantaneous source terms used in normal
recovery. Invariant full-bundle isotypic projections preserve all these
structures, so the estimates do not hide an angular-cutoff constant. The
fixed causal collar, extension of the compatible finite-mode construction,
corner recursion and compatible-forcing qualification are now explicit.

The nonlinear companion proves the exact state-flow reduction only. The
remaining mixed quasilinear existence, parameter-comparison and prescribed
contact-lifespan arguments are clearly marked as unproved. The prior exact
initial-data and short precontact families at fixed finite Sobolev regularity
are retained; no common spatial C-infinity parameter interval is inferred.

The separate expanded derivation also gives a precise sufficient full
Sobolev estimate and optional invariant Poisson boundary lift. A further
internal referee independently checked the source indices, projection
argument and half-order boundary-lift exponents.

The supporting checker passes twenty exact algebra tests under ordinary and
optimized Python. All six executions of three deliberately wrong controls
are rejected: wrong Robin sign, wrong lapse factor, and a false full-matrix
symmetry claim. These checks support the matrix, principal and flux algebra;
they are not a formal verification of the PDE theorem.

This is an internal mathematical review with an independent linear
derivation. The nonlinear companion was authored by this reviewing agent,
so its review here is not an independent authorship claim. External novelty,
exhaustive literature certification and external peer review are outside
scope. The accompanying JSON pins the reviewed files and enumerates limits.
