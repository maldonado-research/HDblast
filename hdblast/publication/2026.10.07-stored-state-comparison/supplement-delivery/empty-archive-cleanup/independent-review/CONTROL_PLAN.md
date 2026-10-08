# Offline exact empty-entry cleanup controls

Preparation only. The root owns any live DELETE; these controls never invoke
helper main, a live transport, original state/journal, repositories, package,
real ZIP inputs or scientific/source callbacks. All subsequent requests and
evidence are manufactured. No archive size-cap or causal diagnosis is inferred
from the two authenticated integration 401 outcomes.

The target is solely the original initialized unavailable ZIP file entry under
draft 23228395, with its exact pinned file/version/bucket identities and strictly
None size/checksum. The only permitted mutation is one DELETE to its fixed modern
file-record URL, after a persistent durable exclusive latch. Native-version URLs
are GET-only. Inherited/user files, drafts and publication endpoints are outside
the mutation scope. Original pending STATE and JOURNAL remain intact.

Planned independent controls:

- Change target draft/file/version/bucket IDs, status, size/checksum, destination
  links, owner/family/latest or any inherited pin/identity independently: reject
  before DELETE, including completed or zero-byte materialized entries.
- Require fresh modern and exact-version native unavailable GET responses; reject
  missing/partial/truncated/error/near-match/wrong-version or present content.
- Authenticate both prior sources and each prior failure/latch/body/header/stderr
  byte pin; reject stale/missing/substituted evidence and changed operation scope.
- Check fixed DELETE URL, no request body or alternate/native/delete endpoints,
  durable O_EXCL latch before transport, one attempt and persistent refusal after
  lost response, malformed success, timeout or failure.
- Exercise existing latch files/directories/symlinks and concurrent reservations.
- Unknown reconciliation performs GETs only. Absence requires complete fresh
  draft/list evidence preserving all inherited entries and identities. Presence,
  changed identity or ambiguous listing stays unresolved and never repeats DELETE.
- Valid and failed flows preserve original pending state and never publish,
  create/delete another draft, delete inherited files, or change metadata.

This plan is not an acceptance, upload, deletion or authorization receipt. Final
review must bind the stable candidate source SHA256 and executed offline controls.

Production-shape clarification received before controls were bound: the pinned
initialization and fresh modern file metadata legitimately omit size/checksum
for this empty pending entry. The earlier "strictly None" key-presence suggestion
is historical and is not the accepted boundary. The final controls must accept
either omission or JSON null, reject every nonnull value including zero, and
also require both modern/native exact unavailable 400 responses and the final
fresh identity guard before DELETE. No real metadata was read by this reviewer.
