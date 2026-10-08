# One identical archive upload recovery

The original archive PUT has an unknown outcome. Its saved intent is 240,135,519
bytes, SHA256 `ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189`,
MD5 `e424502f12dbaacd36abaaabdd3fb74e`, to the already initialized archive in
owned draft 23228395. The original controller, STATE and JOURNAL remain intact.

GET evidence shows the archive entry still pending with its original file and
version IDs and no size/checksum. Both modern content and exact-version native
bucket content return HTTP400: the file is unavailable. This establishes current
readback state, not the HTTP response to the first PUT. The first transport hid
its underlying exception, so a proxy upload limit is not yet established.

Use the exact vendor media type for the modern draft record GET. Generic JSON
and legacy serialization of this pending record returned HTTP500; the exact
vendor response returned HTTP200 with the expected owner and family. The helper
uses the accepted controller's original GET functions and media type.

The candidate is independently reviewed before any root invocation. It captures
the sealed archive into an immutable Linux memfd, then performs fresh prior,
latest, ownership, draft, inherited-file identity and exact unavailable-content
guards while holding the original evidence lock. Its subclass directs evidence
to a separate sidecar and disables original state writes. Any failed guard stops
before a PUT. The helper authenticates the accepted controller and final input
pins, original STATE/JOURNAL pins, and initialized file-identity snapshot.

The fixed recovery directory is durably created and permits no rerun. A durable
exclusive one-attempt latch is saved before starting curl. Curl uses one fresh
HTTP/1.1 operation, verified TLS, the configured proxy and trust environment,
no redirects or retries, and an explicitly empty Expect header. `-q` suppresses
implicit curl configuration. The Bearer header goes through in-memory stdin
configuration; it is absent from argv, saved configuration and child environment.
SSLKEYLOGFILE is removed from the child. Header, response and error streams are
independently bounded and redacted before saving.

The earlier proposed Expect:100-continue route was held because libcurl can
resend following HTTP417. The empty-header route requires independent installed
curl controls to establish that its manufactured417 case does not repeat.

Root must authenticate the reviewed helper bytes before first execution. An
external pinned captured-byte bootstrap may read the file once, compare SHA256,
compile only those captured bytes with the intended `__file__`, and call
`main(['--helper-sha256', REVIEWED_SHA, '--attempt-identical-put'])`. The helper's
own pin check is an additional guard. The preparation agent never invokes this
live path.

Every upload result, including HTTP2xx or curl exit0, still requires complete
authenticated readback. The helper never clears the original pending action,
commits a file, saves metadata, publishes, deletes, reinitializes, switches
routes or retries. Root runs the original accepted controller's GET-only
reconciliation after the attempt. Only its complete byte count, MD5 and SHA256
check may clear the original PUT pending action. A failed readback preserves
both original and recovery histories and requires a separately reviewed next
step.

Multipart links advertised by a native bucket do not establish that deployed
Zenodo supports multipart writes. Current Zenodo source overrides the upstream
multipart resource; a fragment PUT to that native route can replace the whole
file. No such endpoint is used here. Archive, inventory, metadata and scientific
package bytes are unchanged.
