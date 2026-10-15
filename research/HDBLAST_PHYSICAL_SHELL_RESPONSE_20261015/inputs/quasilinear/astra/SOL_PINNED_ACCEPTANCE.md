# Independent acceptance of the Sol reduced-evolution proof

15 October 2026. This is an analytic crossreview, not a machine-checked PDE proof.

**Verdict:** accepted for the fixed m=14 prescribed-interval reduced mixed-boundary theorem, from sufficiently small H14 × H13 data with a common open shell collar and a common protected precontact seed interval. No analytic blocker was found in that scope.

## Exact reviewed sources

- `sol/INDEPENDENT_QUASILINEAR_MIXED_PROOF.md` — SHA-256 `76b83bf96b1f6cb82f63c8b2ea3acab9e9c9394e3f961db350a05d6ea68bb1f8`.
- `sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md` — SHA-256 `fcd7dc55b9459f0ebc56b1e6a77584be6a101878126efb692614378b4884af86`.

The final main-file change after the full read corrected three TeX backslashes only. The final linear-audit change clarified that the full geometric estimates are uniform in angular cutoff, while the intermediate fixed-angular radial construction may have constants depending on that cutoff. Both corrections were checked.

## Specific analytic checks

1. The angular projection commutes with the selected generators, time and normal derivatives, and preserves the fixed split. Coefficient multiplication is never commuted through it. The radial mass and stiffness retain their coercive L2 bounds.
2. The radial Galerkin time-energy correction cancels the differentiated weak spatial-form term. Its remaining time derivative contains at most the current time energy or lower ones. A positive weighted sum closes without assuming the radial H2 regularity that is recovered only after the weak limit.
3. A spacetime derivative of order r of u has spatial differential-polynomial order at most r−1 in Y=(u,u_t,D_xu). A derivative of g(u) of order q≥1 contributes order q−1, so the principal commutator total is at most m−1. At most a fixed finite number of recent iterates occurs. The protected precontact field is retained as a known field when the recursion reaches the cutoff seed; no cutoff-derivative ratio is bounded incorrectly.
4. The spatial product estimate follows from generalized Hölder and Gagliardo–Nirenberg with derivative orders summing to q, followed by weighted Young. It is linear in the high norms. The low A6 bound controls the needed undifferentiated Y and second derivatives, using the evolution identities.
5. Pointwise normal-recovery forcing costs one less derivative and is explicitly present. Purely spatial interpolation makes its preceding-iterate high coefficient arbitrarily small. Low energy and a sufficiently fast exponential envelope then close the whole prescribed interval; no nonzero-corner restart theorem is smuggled into the argument.
6. Difference energy closes in A1 because coefficients depend on u, not Du, and the second derivative of the preceding solution is bounded. The difference sequence converges factorially or in an exponentially weighted norm.
7. At order m−1, the difference forcing requires only the uniformly controlled A_m norm of the comparison solution. At order m, its unique extra principal derivative is D^(m+1)v multiplied by the difference. This gives the stated lower Lipschitz and top smoothing-comparison estimates.
8. Supported smoothing preserves a smaller fixed initial zero collar. Its lower error is o(N^−1), its high norm is O(N), and its top error tends to zero. Put the coarser smoother solution in the extra-high factor. The resulting Bona–Smith comparison gives strong top Sobolev continuity, not merely weak compactness.

The accepted theorem concerns reduced evolution. Exact harmonic completion, subsidiary constraints, both-copy matching, and parameter differentiability remain separate interfaces. No general arbitrary-corner continuation existence theorem, top-norm differentiability, common spatially smooth family, global stability, novelty, or Big Bang conclusion is accepted here. The source-specific verdict must be reconsidered if either pinned file changes.
