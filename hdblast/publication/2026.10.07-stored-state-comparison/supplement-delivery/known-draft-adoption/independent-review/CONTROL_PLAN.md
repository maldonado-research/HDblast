# Focused local draft adoption control plan

Preparation only; root owns actual invocation. No helper main, protocol behavior,
network, protected original state/journal, source callback, package or real new
execution state is invoked or mutated. Forthcoming controls use disposable local
files and explicitly manufactured evidence, not authorization receipts.

The only permitted NEW state change is `draft_id: null -> "23228395"`. Every
other field/value remains identical. NEW pending stays null and published stays
false; the old genuine creation/unknown-PUT evidence stays unchanged. Adoption
adds one auditable NEW journal event after a durable exclusive persistent latch.

Focused cases will cover valid adoption and sole-field comparison; exact journal
prefix preservation; original-byte preservation; source/input/state/journal/
creation/cleanup-proof pin mismatches; binding mismatch; wrong draft identity;
nonnull pending/published/previous draft; wrong cleanup absence proof; exclusive
latch refusal and repeat refusal; aliases into protected originals; and partial
local-write failure retaining the latch. Network/protocol implementation is out
of scope; root runs the accepted remote guards after the local step.

The final receipt must bind a stable source SHA256 and executed temporary-file
controls. This plan does not establish actual adoption or authorize repetition.
