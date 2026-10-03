# Independent prospective tail-module review

Reviewed 2 October 2026, before the stress run's public freeze. Scope:
`stress_tail_bounds.py`, `STRESS_TAIL_DERIVATION.md`,
`pulse-norms/pulse_derivative_budgets.py`, their exact polynomial/root and
global rational proof sources, and recorded proof JSON. No deferred source,
interval, response, mode, or bound-point function was called in this review.
No checkout edit was made.

**No blocking defect found.** The following checks were independent static
algebra and enclosure reasoning, not numerical execution evidence.

- The q-prime and q-second bounds retain derivatives of the moving subtraction
  scale. The second contact factor is exactly
  `(1-v)(3v^4+3v^3+v²+v+1)=1+2v³-3v⁵`.
- The normalized anomaly coefficients are correct. In particular
  `a^4 d_ddot/epsilon=f''-5Lf'+4L²f`, so its first coefficient is
  `f''-12L²f`; the second is `30L²f+10Lf'`.
- Density normalization, pressure reconstruction, and the b=1,r=2 current
  `J=q+Q0 f/4` agree with the independently reviewed stress formulas.
  Combining the exact anomaly and density contact tails before taking the
  interval absolute value preserves a valid, potentially tighter enclosure.
- The stable factorizations of the three moment tails and positive reference
  tail are correct. `1-v=M²/[sqrt(K²+M²)(sqrt(K²+M²)+K)]` avoids cancellation.
  Interval dependency may widen these expressions, but cannot invalidate
  containment. Final binary64 upper endpoints are advanced outward.
- The pulse module encloses every primitive extremum, in increasing u order,
  including the parity-dependent zero root. Exact rational comparisons order
  the observation endpoint against each root bracket; an unresolved comparison
  raises an exception. Interval sums of absolute successive differences
  enclose partial total variation, including the final endpoint contribution.
- The separate polynomial verifier checks derivative recurrence, complete
  root count, simple roots, parity, disjoint brackets, and all registered
  endpoint orderings. Thus no additional unsampled extremum is assumed away.
- The rational global proof uses interval polynomial arithmetic, enclosing
  square roots, and reciprocal exponential Taylor bounds with a valid
  geometric remainder. For any endpoint-zero primitive g,
  `|g(t)|+TV(past)<=TV(full)` follows from `TV(future)>=|g(t)|`.
  Its resulting ceilings therefore bound every partial derivative budget.
- The global stress certificate consistently uses pi>3, L<=2/3,
  `1-v<=M²/(2K²)`, the correct positive-polynomial maxima, and the common
  source envelopes. Its rational caps concern omitted momentum tails only.

Recorded evidence reports PASS in ordinary and optimized Python for the
pulse proof: 18 identities, six complete critical-root certificates, and
seven mutations. Both stress-tail algebra records report PASS: 25 identities
and six mutations. Both rational pulse-envelope records report PASS with
global `(N0,N1,N2)` ceilings `(97,2926,161866)` for B and
`(96,2760,154078)` for uB. `GLOBAL_STRESS_CEILINGS.json` reports the exact
rational upper enclosures below the declared K=256 analytic ceilings:

    normalized density tail < 3/100000,
    normalized pressure tail < 1/1000,
    normalized current tail < 3/1000000.

These recorded proof results were inspected, not rerun in this review.
They establish analytic feasibility and omitted-band envelopes. Actual
deferred interval values, finite-band quadrature accuracy, independent mode
agreement, and numerical Ward ledgers remain for the registered post-freeze
run. The module and its document state this separation accurately.
