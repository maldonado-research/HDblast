# Independent one-shot recovery review

Accepted for the manufactured control boundary: helper SHA256
`39ec701ae3133b699bbad1f74031c549a0aedbc2116109cf5be199d3c7b0916d`.
The final independent suite passes 80 controls in normal Python and the same
80 controls under `-O`, with zero failures. This is an offline review receipt,
not authorization, upload evidence, or a claim that any live PUT occurred.

The suite never calls helper `main`, `load_controller`, `controller_view`, the
real curl executable, a live endpoint, or any scientific/source callback. The
module import contains only its stdlib definitions. Socket creation and
connection are blocked. Inputs, assets, tokens, response objects, environments
and latch directories are manufactured; response-collector controls use actual
anonymous pipes fed by an in-process subprocess double. Temporary fixtures are
removed. No original publication state/journal, repository, package or input
asset is read or edited by these controls.

The controls cover original pending intent, draft/published/create/publish stage,
controller/inventory/metadata bindings, initialized and inherited file/version/
bucket identities, unavailable response status/body precision, pinned read
leaves, exact SHA256/MD5/stat asset capture, immutable memfd write/resize refusal,
partial writes, changes during capture, preservation after capture, exclusive
and concurrent persistent latch reservation, token header-injection rejection,
recursive escaped redaction, fixed curl options and destination, metrics
ambiguity/numeric/trailing checks, three-channel 8MiB collection bounds, one
subprocess launch, retained proxy/CA settings with token/keylog variables removed,
and wall/held-pipe drain termination.

The original `4db55b72...` snapshot and its original test harness and receipts
are retained. That first run showed 57 PASS and 12 FAIL in each mode. Two
expectations were subsequently corrected to the documented contract: the
original state may omit its publish-attempt flag, and curl status `000` denotes
an unknown no-HTTP outcome. Three failures exposed missing Linux seal constants
in the actual Python 3.12.14 runtime. The remaining failures covered permissive
publish flag values, changed initialization identity, numeric HTTP status,
ambiguous metrics and malformed duration. The corrected helper supplies reviewed
Linux UAPI seal constants, strict guards, unique anchored metrics, recursive
redaction and bounded post-kill draining. Added controls then verify those
boundaries explicitly.

The curl command uses `-q` first, HTTP/1.1, one fixed URL, a PUT from the sealed
descriptor, empty `Expect:`, retries zero, no redirect following, HTTPS-only
protocols and ordinary certificate verification. These argument and single
process controls do not by themselves prove one wire request for the actual
curl binary. The separate transport reviewer owns its manufactured socket-level
417/empty-Expect check. This receipt makes no server-side compare-and-swap or
successful-delivery claim. Every actual outcome still requires fresh read-only
reconciliation; original pending remains uncleared and commit/publish counts
remain zero in the helper's design.

No remaining blocker was found within this independent offline control scope.
