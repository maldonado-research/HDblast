# Focused independent local adoption review

Source SHA256 `f3cd239eba834806961a6ed246ab183ee63c149221f1fca2a9dba7a2cb5db4da`
passes 17 focused controls in normal Python and the same 17 controls in `-O`,
with zero failures or remaining blockers found in this scope. Source and tests
are stable; no helper changes were made by this reviewer.

The suite uses manufactured old/new STATE, journals, creation and cleanup proof
bytes in disposable directories. It copies the accepted `aa7e6cfe...` controller
source bytes solely to satisfy the immutable byte-digest binding and never imports
that controller. The adoption helper import contains stdlib definitions only.
No helper main, network, protocol, source callback or actual state/journal read or
mutation occurs. Temporary recovery evidence is not an authorization receipt.

The valid application changes only NEW `draft_id` from null to `23228395`; all
other fields remain equal. It preserves the NEW journal's byte prefix, appends
the intent and completion events, retains the durable exclusive latch and leaves
protected original bytes unchanged. Dry-run writes nothing. Repeated application
is refused. An injected atomic state-replacement failure retains the latch and
original state bytes, preserves the journal prefix and forbids automatic repeat.
Other cases reject missing roles, stale source pins, nonfresh state, binding/input
mismatch, wrong creation/cleanup proof, incomplete journal boundary, aliased input,
symlink state and recovery paths inside execution directories.

Main authenticates a hardcoded complete INPUT_PINS manifest hash and its nineteen
fixed input roles. The temporary author-pinned API fixtures exercise the local
implementation; they do not create a live authority path. This review therefore
does not expand into generic parser or protocol tests. Source review confirms the
helper has no controller import or transport, takes both execution locks, checks
all pins again after its durable intent, atomically replaces NEW state and checks
the original state/journal after replacement. The unchanged controller's fresh
remote guards remain a separate mandatory next step under root ownership.

Root reported that its actual local adoption and separate readback had completed
while these isolated controls were being prepared. This receipt records only the
independent manufactured tests; it neither performed nor independently read that
actual adoption. The original pending ZIP history remains outside the mutation
scope.
