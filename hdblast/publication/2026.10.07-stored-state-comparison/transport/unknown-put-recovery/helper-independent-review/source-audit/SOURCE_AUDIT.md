# Recovery helper final source audit

Accepted with no remaining blocker in the assigned source-audit scope. Final helper SHA256: `39ec701ae3133b699bbad1f74031c549a0aedbc2116109cf5be199d3c7b0916d`, 19,883 bytes. The final candidate snapshot is pinned alongside this report. This audit read source only: no helper import/main, live request, real asset/state read, mutation, or repository edit occurred. Parent owns unit controls, manufactured loopback curl controls, and any authorized live attempt.

Acceptance covers the specifically authorized identical PUT to the existing initialized ZIP key in draft 23228395. The helper retains the unknown original intent and requires independent controller GET hash reconciliation afterward. It cannot authorize a repeat, clear pending, commit, or publish.

## Production URL and curl transport

Fixed single HTTPS URL targets draft 23228395 and the unchanged initialized ZIP content key. Absolute /usr/bin/curl receives -q first, --http1.1, one URL, explicit PUT/upload-file, empty Expect, retry 0, max-redirs 0, HTTPS protocol restrictions, no redirect following/auth negotiation/insecure TLS/debug options. Timeout bounds and response-file bound are set. Required proxy/trust environment is preserved; ZENODO_ACCESS_TOKEN and SSLKEYLOGFILE are removed from child environment. Bearer is passed only through an ASCII-validated, quote/backslash-escaped stdin config. Actual curl wire behavior is validated by parent controls, not this source-only review.

## Source and local evidence bindings

External helper SHA256 and explicit attempt flag required. Controller aa7e6cfe, inventory f6b7db33, metadata 0d089617 and initialized-file response 7b4e819d are complete-byte pinned. Original state/journal external hashes are checked under EvidenceLock and before/after upload. State requires exact draft, unpublished false, source/input hashes, create attempted true, absent/false publish latch, no ZIP content verification, and exact pending PUT shape/URL/body hashes/size/context. Controller constructor retains original family/prior/baseline/receipt binding validation.

## Original controller writes suppressed

View.save is a no-op, View.event writes only SIDE/GUARD_JOURNAL.jsonl, View.capture writes only exclusive sidecar files, and View.write always raises. Main invokes current, owner_list, draft and a GET to content, not reconcile/write/prepare/publish. Original STATE and JOURNAL are hash-checked unchanged. EvidenceLock protects the original evidence directory. No original pending clearing, commit or publish call is present.

## Immutable byte custody before fresh guards

Complete ZIP is captured before current/owner/draft/file/missing-content guards. O_NOFOLLOW and single-link regular-file checks, full byte count/MD5/SHA256, unchanged open-descriptor stat, unchanged source-path stat, and exact asset pin are required. Captured anonymous memfd has all four immutable seals, verified by GET_SEALS, and rewinds to offset zero. Curl receives only that sealed fd via pass_fds and /proc/self/fd. Linux-only gate and UAPI fallback ADD1033/GET1034/seal bits15 match installed Linux headers.

## Fresh remote guards and identity

After capture, public prior/latest checks require no completed publication, owner listing is exhausted and matches the known draft, draft family/owner/state/core metadata are authenticated by the pinned controller, exact 27 inherited plus one pending ZIP membership is required, inherited pins and initialized immutable identities/bucket are unchanged, new key/file/version/bucket and content/commit links are exact. Content GET requires bounded complete HTTP integer400 and exact integer400/string-message JSON unavailable response. Source and original evidence hashes are rechecked after guards and before latch.

## Durable exclusive single attempt guard

Existing fixed SIDE refuses another invocation and exclusive mkdir prevents replacement; parent HERE is fsynced immediately after mkdir. EvidenceLock wraps guards and upload. ONE_PUT_ATTEMPT_LATCH uses exclusive no-follow 0600 creation, flush/file fsync and SIDE-directory fsync before the sole run_curl call. A guard/transport failure leaves SIDE/latch for manual review; no automatic subsequent invocation, reinitialization, commit, publication or pending clear is attempted. This is a client-side attempt guard, not a proof about proxy/server internal delivery.

## Bounded, redacted response custody

JSON string values and keys are redacted before serialization; both literal and normal JSON-escaped bearer representations are removed from raw buffers. Three separate streams are each bounded at 8MiB; overflow or 190-second wall deadline kills curl, then drains at most5seconds before stopping. Curl max-time180/max-filesize8MiB further bound transfer. Metrics must appear exactly once at stderr end with strict numeric syntax and status000 or100..599. Missing metrics/overflow/deadline/exit are retained as outcomes; generic result explicitly requires reconciliation, never declares upload verification success. Original evidence is rechecked, redacted raw outputs have hashes, and receipt/STOPPED are exclusive fsynced sidecar writes. Parent retains responsibility for proxy-auth environment scope and complete GET receipt validation.

## Resolved findings

- Post-serialization raw token replacement misses permitted quote/backslash-bearing tokens: Recursive pre-serialization key/value and direct/JSON-escaped raw buffer redaction.
- New SIDE entry parent not fsynced before attempt: Parent directory fsync immediately after exclusive SIDE mkdir; durable latch file+SIDE fsync preserved.
- Full asset capture occurred after remote freshness checks: Capture now precedes the complete fresh current/owner/draft/file/content guard chain.
- Child environment inherited unnecessary bearer variable and curl TLS key-log hook: Explicit child env removes only ZENODO_ACCESS_TOKEN and SSLKEYLOGFILE while preserving proxy/trust.

## Limits

This review proves the source guards are present and ordered as described. It does not establish the real server cause of the earlier unknown PUT, validate a real upload, certify every libcurl or proxy internal resend path, or establish server-side compare-and-swap publication semantics. Empty Expect, fresh single-URL HTTP/1.1 and disabled redirect/retry/auth-negotiation controls are appropriate for the root-approved route; parent transport controls supply observed behavior. Proxy/trust settings are deliberately retained. Raw diagnostics are redacted for the known bearer; this audit does not claim comprehensive identification of unrelated environment secrets or arbitrary encodings. No result status bypasses the original pending latch or complete content readback requirement.
