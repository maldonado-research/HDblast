# Scoped manufactured operation rehearsal edit

Verdict: **PASS** for the explicitly requested pre-freeze fabricated-entry
change. No production source evaluation or registered-source tail certificate
was performed by this rehearsal.

Old wrapper SHA256:
`334a1cdee4992493231c4646ea2dcf242196ca67542560054a3ddd45424b456f`

New wrapper SHA256:
`b32e1a515cfa3f62cd23f520fd201c502a856a82dc81f405b0c0daf6a292e6a6`

Engine unchanged:
`de995d2edeca3944400fd90b892521f05c941b6d8a1f75c2a820740b82d995c0`

## Exact edit boundary

The only changed existing function is `run_fabricated`. The only added
function is `manufactured_arb_operation_rehearsal`. Comparing all other
top-level module nodes by Python AST confirms that imports, production
`run_primary`, production `run`, the arithmetic/source budgets, the fabricated
polynomial provider, and the CLI remain unchanged.

The new helper performs 128 unrelated manufactured Arb jet computations.
Their centers are positive exact rationals 2+(2j+1)/128, j=0,...,127; they are
outside the registered negative physical source interval. Each manufactured
source is exp(t/(2+t^2)), t=c+H*x. The helper explicitly uses the inverse of
2+t^2 and an exponential, and the existing pure forcing algebra computes the
inverse defining L, both eta derivatives, and the Lg coefficients through
degree 24. Each iteration requires exactly 25 finite g coefficients and 25
finite Lg coefficients. These are finite-operation checks, not analytic
remainder tests for the registered bump source.

The helper temporarily sets 256-bit Arb precision and series cap 27. Its
context manager and try/finally restore both caller precision and caller cap.
It never calls either registered source callback. `run_fabricated` invokes
the rehearsal first, then its unchanged manufactured polynomial fixture, and
adds only constant metadata to that fixture's budget:

```
scope: MANUFACTURED_ARB_OPERATIONS_ONLY
registered_source_calls: 0
manufactured_operations: 128
g_coefficients_checked_per_operation: 25
Lg_coefficients_checked_per_operation: 25
source_tail_certification: false
```

## Verification

Both of these checks pass; the check script uses explicit exceptions rather
than assertions, so its verification remains active under Python -O:

```
PYTHONDONTWRITEBYTECODE=1 /workspace/hdblast-cloud-setup/venv-arb/bin/python \
  cauchy-proof/manufactured-rehearsal-check/check_rehearsal_entry.py --mode normal
PYTHONDONTWRITEBYTECODE=1 /workspace/hdblast-cloud-setup/venv-arb/bin/python -O \
  cauchy-proof/manufactured-rehearsal-check/check_rehearsal_entry.py --mode optimized
```

Instrumented pure forcing calls confirm exactly 128 manufactured operations,
cap 27, precision 256, degree 24, and positive unrelated centers. Both physical
callbacks are replaced with rejecting traps; their invocation count is zero.
The checks start with caller precision 177 and cap 13 and confirm exact
restoration of both after `run_fabricated` returns.

Normal and optimized fixture payloads are byte-identical to one another and
unchanged from the saved pre-edit payload. Existing budget fields are unchanged
from the pre-edit budget after removing the one newly requested constant
metadata field. The full new normal and optimized budgets are byte-identical;
no timing/resource value was added to the stable evidence.

Payload SHA256, before and after, normal and optimized:
`186055875b3504d8513217fdb6662d68d2c4e03cbb87f0ed2b810bd02f37b22c`

Full new budget SHA256, normal and optimized:
`fc371db5ffd04275edc350c1e26cd2a7e931f558b647f8c354c9906cb465e262`

Evidence directory:
`cauchy-proof/manufactured-rehearsal-check/`

It contains the pre-edit source snapshot, pre-edit payload and budget, normal
and optimized post-edit payloads and budgets, verification receipts, the
explicit-exception check script, and `PRIMARY_ROUTE_SCOPED_DIFF.patch`.

This change leaves the prior production mathematics review applicable. Root
still owns the authenticated physical bootstrap, public freeze, production
event counting, protocol, and any scientific run.
