# Pre-registration preparation history

No real bump source callback or saved quantum-array decode occurred during
these preparations. Repository source and frozen packet bytes were read only.

The first fabricated CLI attempt used system `python`, which lacks python-flint,
and failed at import. The pinned existing interpreter
`/workspace/hdblast-research-work/environment/bin/python` was then used.
Its first attempt failed because Python's Fraction constructor does not accept
python-flint's `fmpq` object; conversion now uses its exact rational string.
Both failures occurred before any real callback, and the latter during
fabricated output serialization only.

The first sparse-polynomial affine test used `actual.contains(expected_arb)`.
Two rounded balls can contain the same exact rational without one containing
the other's full rounding interval. The test now checks the exact rational
against exact exported endpoints. This was a test correction; the polynomial
recentring implementation did not change.

Benchmark 001 retained the first successful pre-review route. Source guard
payload keys were subsequently aligned with the root wrapper; exact source
coefficient endpoints were added to output rows; phase/source/cap contributions
were combined into one explicit final inflation with a proved exact export
radius relation. Benchmarks 002/003 cover that version. The source callback's
common payload is now exactly source, center, and the prescribed scope.
Fabricated callback-shape tests stop before source coefficient construction.

The common outer schema represents the shared inflation factor as a scalar,
so the earlier per-U/W pair was simplified without changing any scientific
numbers. Benchmarks 004/005 cover the final reviewed-candidate implementation,
with identical data SHA256 and matching 238 normal/optimized fabricated
checks and twelve rejected controls. The separate outer validator passed the
final fabricated primary output with the actual `1e-20` width gate enabled.

`READY_FOR_REVIEW_V2.json` pins the final code, method and two complete
benchmarks. Review readiness is not a public freeze or authorization to call
the real source. The parent owns independent review, registration, public byte
readback and the authenticated bounded physical runner.
