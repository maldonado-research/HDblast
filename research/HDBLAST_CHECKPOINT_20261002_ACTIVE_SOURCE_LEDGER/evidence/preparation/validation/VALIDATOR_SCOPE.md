# Authoritative comparison validator

`validate_active.py` is a standard-library-only result reader. It authenticates
its registered entry bytes and integrity helper before importing the helper,
then verifies the complete frozen checkpoint before reading route results.
The independent wrapper's complete frozen-verification receipt must match the
validator's own receipt. The primary and independent producer hashes, manifest,
registration and public freeze are bound separately.

The field inventory is generated from the validator's literal contract tuples
in `FIELD_INVENTORY.json`. Unknown, omitted and additional scientific fields,
controls, contexts or cases are rejected. Integer metadata rejects Boolean and
floating-point substitutes. Scientific values are finite decimal strings read
as exact fractions; rational serialization receipts remain exact fractions.

All twelve cases, primary GL24/32 controls, independent GL16/16, GL16/24 and
GL24/24 controls and both MP80/100 contexts are retained. Every primary numeric
row and every independent MP scalar/control field enters the precision
comparison. Independent native profiles are represented once; this does not
claim an MP reevaluation of their source, phase or trajectory arithmetic.
All six cross-route control pairs are compared at both contexts. Independent
ledger, forcing and joint control comparisons cover their complete reported
scalar and profile universes.

Exact definitional closures, signed decompositions, reported maxima, precision
gaps, serialization and fixed stored-input agreement use `1e-12`. Scientific
profile, control, flow, operator and reconstruction gates use `2e-7`. A complete
scientific consistency failure remains a retained negative result with exit
zero. Integrity, schema, arithmetic, serialization and authoritative resource
failures raise an exception. Strict attribution requires every consistency
gate to pass and the same canonical case to satisfy the registered strict
criteria in both routes.

Saved density differences are independently recomputed from endpoint profiles.
The original native global-prefix Simpson ledger, reset ledger and fine
double-step ledger are checked between the two independently implemented
readers. The result JSON contains only interval profiles, so this validator
does not claim a third replay of the complete native global-prefix sums.

Raw continuous finite-K versus discrete-momentum contact/profile differences
remain explicit full-profile reports. Own discrete contacts and physical versus
scaled baseline normalization are audited at initial, midpoint and endpoint
only. Those three audits do not prove full-profile target matching.

`test_validate_active.py` manufactures complete result records and resource
receipts. Its attribution fixture checks classification policy; it is not a
physical solution, numerical trajectory or scientific evidence. It verifies
retained negative outcomes, fatal mutations and bootstrap rejection before a
deliberately marked helper can import. These tests import no numerical
producer, array decoder or source model and evaluate no saved study values.
The production driver separately supplies actual nonroot execution, resource
and process-containment evidence. Existing metric status remains `FAIL`.
