# Internal review of the formal reflection-invariant constraint theorem

14 October 2026. Reviewers are independent OpenAI Codex agent processes,
not external human referees. No novelty certification is supplied.

The mathematics-continuation reviewer independently checked the source
induction, nth-derivative factorial factors, absolute-density reflection
pairing, copy parity and finite-order scope. The literature reviewer
independently checked the same mechanism, the absence of invariant local
KIDs and the support argument. Both accepted the proposition after the
same support wording repair: apply the localized solver on each
predetermined `𝒜_n` itself, with `closure(𝒜_{n−1}) ⊂ 𝒜_n ⊂ 𝒜_∞ ⋐ 𝒜`.
The final text contains this repair.

The exact finite-group verifier completed in normal and optimized Python:
**266 positive checks and four wrong-claim controls in each mode**. The
controls reject central inversion as sufficient for eliminating rotations,
antipodal evenness as sufficient for rotational orthogonality, independent
reflection invariance for arbitrary quadrupoles, and confusion between a
second derivative and an `ε²` coefficient.

The central-inversion counterexample is an exact normalized rotational
sphere moment `−1/3` for an inversion-even noncommuting pair. All six
rotations remain invariant under central inversion alone; none is
invariant under the full sixteen-element reflection group. These checks
support the algebraic mechanism. The smooth localized solvability theorem
is a mathematical input, not a machine-verified theorem.

The reviewers agree that no convergence, exact nonlinear initial-data
family, shell-contact evolution or Big Bang inference follows. Standard
equivariance plus the localized solver produces formal coefficients
through every finite order for this symmetry sector.
