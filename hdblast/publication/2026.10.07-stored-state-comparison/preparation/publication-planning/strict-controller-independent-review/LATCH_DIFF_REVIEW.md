# Controller source diff and publication latch review

The current controller is accepted with no blocker in this static source diff and latch review. The review does not execute the controller, import it, call the network, edit repositories, or change publication state. Behavioral validation of the revised comparisons belongs to the parent independent review.

Archived source: `27fd844cc903a445688a177339b7af633c0a3581f71eeedf7e20bb015bee36a7` (46,109 bytes). Current source: `aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b` (49,240 bytes). Both were independently calculated from complete source bytes.

The complete `Controller` source segment is byte-identical, SHA256 `5c461b56768594c53defed553f1ec1afbfe2c3ebadb171cd6b6697cfc7e66307`. `Transport`, `main`, `EvidenceLock`, `revision_if_match`, `assert_pins`, `immutable_file_identity`, and every other top-level definition outside the comparison changes are also byte-identical. All 34 module nodes outside function/class definitions are byte-identical, covering constants, imports, and entry-point dispatch.

Only `HtmlEvents`, `normalized`, `value_equal`, `metadata_equal`, and `file_map` changed; `vocabulary_decorations` and `json_structure_equal` were added. Within `HtmlEvents`, `handle_starttag` and `handle_data` changed; `handle_decl`, `handle_pi`, `unknown_decl`, `parse_endtag`, `close`, and `VOID_TAGS` were added. All other HTML methods are unchanged. These edits strengthen rich-text event preservation, complete JSON type/shape checking, required metadata/access presence, and remote file-list name consistency. No mutation, recovery, transport, stream, journal, latch, source-binding, ETag, or command-line function changed.

## Persistent create latch and explicit safe 406 retry

create_post_attempted is saved before create POST; a repeat without retry is stopped. Retry requires original attempted latch, no pending/draft/publication, no prior retry, exactly original empty-body WRITE_INTENT plus definitive integer HTTP 406 WRITE_RESPONSE, and saves create_406_retry_attempted before retry. write persists pending intent before transport request and keeps ambiguous outcome pending. Evidence: `Controller.create`, `Controller.write`, `Controller.reconcile`.

## Persistent publish latch

publish rejects published state, absent/tampered verified receipt, pending/reconciled publication, or publish_post_attempted; sets and saves publish_post_attempted before write. Published state and post-attempt latch prevent repeated POST. Evidence: `Controller.publish`, `Controller.write`.

## Read-only recovery and known owner/family fallback

Recovery uses GET evidence, owner listing must be exhausted, exact owner/family checks and at most one owned draft apply, and latest-only fallback requires the already journaled numeric record and exact published/draft booleans. Pending mutations are cleared only after the corresponding evidence passes. Evidence: `Controller.reconcile`, `Controller.owner_list`, `Controller.current`, `Controller.draft`.

## 29 fresh full streams before and after publication

FINAL_COUNT remains 27+2. verify traverses self.old+self.add, validates exact membership/size/MD5/completed state, record-bound content links, immutable identities, and complete byte count/MD5/SHA256 via stream_check for every file. All 29 receipt bases must be FRESH_COMPLETE_CONTENT_STREAM. publish calls verify immediately before latch and verify(published=True) after publication polling; public streams use public=True. Evidence: `Controller.verify`, `Controller.stream_check`, `Controller.publish`, `assert_pins`.

## Fresh immutable identity map before publication

publish compares fresh file_id/version_id mapping to the completed stream receipt immediately before publication; no prior-release identity shortcut replaces the new edition full streams. Evidence: `Controller.publish`, `immutable_file_identity`.

## Source and input/state binding

Continuation binding contains prior/family IDs, actual controller SHA256, baseline, published baseline, and prior receipt canonical hashes; mismatched binding and orphan journal are refused. Input SHA256 state pins and sealed inventory/completed scientific replay marker, inherited 27 pins, exactly 2 valid local addition pins, and final totals remain required. main still authenticates inventory/metadata by explicit SHA256 and requires explicit mutation/publish flags. Evidence: `Controller.__init__`, `Controller.validate_inputs`, `frozen_json`, `main`.

## Raw ETag equality precedes decimal revision conversion

publish requires the fresh raw ETag equal verified[etag] before revision_if_match; conversion requires canonical quoted nonnegative decimal and exact nonboolean integer revision_id, returning its decimal digits for If-Match. Metadata save uses the same unchanged converter. These client checks do not establish server-side publication CAS semantics. Evidence: `Controller.publish`, `revision_if_match`, `Controller.prepare`.

## Public publication readback

After publication, bounded anonymous public polling waits for published status; verify(published=True) checks new record/family, expected files link, complete metadata and all 29 public content streams; current then rechecks the public prior and current latest edition. Evidence: `Controller.publish`, `Controller.verify`, `Controller.current`.

## Limit of acceptance

This proves the publication safeguards above have not changed in source and remain present in the revised controller. It does not establish remote Zenodo behavior, prove server publication compare-and-swap, or independently certify every rich-text/JSON equivalence case. There are no blockers in the assigned source diff and latch scope.
