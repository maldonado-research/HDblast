# Offline independent recovery controls

This is a control plan, not an execution or authorization receipt. The helper
candidate is not yet available. All proposed tests use disposable manufactured
evidence, asset bytes, dummy credentials and subprocess/transport doubles. No
live endpoint, real publication state, repository, package, or scientific input
is changed or executed.

The acceptance boundary is one PUT to the same initialized destination using a
validated immutable asset capture. Invalid evidence must cause zero PUTs. Once
the persistent attempt latch exists, every later failure and every subsequent
invocation must preserve it and permit zero additional PUTs. No test authorizes
publish, metadata mutation, deletion, another initialization, or fallback upload.

| Area | Independent failures and boundary checks |
| --- | --- |
| Evidence | Change each original pending/type/draft/publish/source/input pin separately; remove fields; use wrong types, duplicate JSON keys, stale bytes, or swapped evidence. |
| Initialization | Change file, version, bucket or inherited IDs separately; mix different initialization receipts; inject path/query destinations or coercible values. |
| Unavailable gate | Accept only the exact locked HTTP 400 body; reject other statuses, near-match bodies, appended bytes, truncation, malformed encoding and body-only status text. |
| Capture | Reject symlink/nonregular assets; alter size, SHA256, MD5 or required stat identity; replace or modify source during capture; exercise partial/interrupted reads. |
| Sealed upload | After capture replace source; confirm upload uses exact captured bytes from offset zero; attempts to write or resize sealed memfd must fail. |
| Exclusivity | Existing latch file, directory or symlink blocks; concurrent reservations permit at most one attempt; cwd/path changes cannot select a different latch. |
| Persistence | Reserve before any PUT can leave; keep latch after timeout, partial transmission, lost connection, malformed response, collector failure and success. |
| Token custody | Distinctive dummy token absent from argv, logs, artifacts and exceptions; reject header injection; check reflected response errors do not leak it. |
| Attempts | Connection errors, partial send, 401, 429 and 5xx do not cause replay, credential refresh, alternate transport or fallback upload. |
| Destination/TLS | Fixed HTTPS destination and normal certificate/hostname verification; no redirects, input overrides or environment weakening. |
| Collection | Header/body/stderr overflow, endless stream, timeout, split boundaries, informational responses, conflicting final status, stale collector and truncated capture fail closed within fixed limits. |
| Integration | Valid fixture makes exactly one PUT; invalid fixtures make none; original evidence stays unchanged and no publish/delete/metadata operation occurs. |

Transport assertions will bind to the stable helper design. In particular,
curl retry suppression alone must not be treated as proof against curl's
internal protocol replay. Any approved curl command must suppress Expect and
redirect/authentication replay routes and preserve the persistent unknown-outcome
latch. This plan makes no server-side compare-and-swap or delivery guarantee.
