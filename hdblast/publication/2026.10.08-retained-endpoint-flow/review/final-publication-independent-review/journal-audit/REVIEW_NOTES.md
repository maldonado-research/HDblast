# Independent saved publication journal review

The final saved journal passes the independent data-only audit. It records 197
events and exactly 93 successful file streams: two complete 31-file draft rounds
before publication and one complete 31-file public round after publication.
Each round binds the same 453,646,943 bytes, with exact inventory filenames,
byte counts, MD5 and SHA-256. No duplicate, unknown, incomplete or unscoped stream
appears in the terminal journal.

The exact nine recorded writes comprise one new-version creation; initialization,
content upload and commit for each of the two additions; one complete metadata
save; and one publish. Their target URLs, body sizes, body hashes, ordering and
successful responses were independently checked. No mutation followed the last
complete draft round except its one publish. No unknown or failed write outcome,
truncated write response, reconciliation or retry appears in this journal.

The publish intent is dated 2026-10-08T17:49:16.409885+00:00; the HTTP 202 response
is dated 17:49:19.454259. The final public content receipt is dated
17:52:21.331322 and the terminal journal event 17:52:22.725316. Saved modern record
creation is 17:49:16.754647 and its metadata publication date is 2026-10-08.
These are publication dates, separate from prior scientific execution dates.

Draft, public and latest saved identities match new record 23244754, DOI family
17088132 and owner 1386319. Version and concept DOI identifiers, latest flags,
semantic version, completed file listings, embedded file identities, and all 31
draft/public file-and-version identity pairs were checked. The final draft
receipt binds revision 9 and the persisted controller state binds that receipt.

External terminal pins supplied by the root publisher were authenticated:

- `JOURNAL.jsonl`: `20386fdd5734494bbb2dd92516af3099ab4c8eb8ade8d3a0da7d616e36cf6d1f`
- `VERIFIED_PUBLIC.json`: `71804ff951b4a0aa6643b939fc133afc3352f9545162a141443c1ab26d267e0b`

The accepted controller source, sealed inventory and sealed metadata were
authenticated against their independent external pins before audit. The final
journal remained identical during capture. The earlier observed partial journal
is an exact prefix of the final journal and remains preserved. A later partial
snapshot records the ongoing public stream round; neither partial observation
was treated as final acceptance.

The independent checker rejected 16 manufactured corruptions and accepted its
partial baseline in both normal and optimized Python. Those 34 control outcomes
and their preserved artifacts are bound by the control evidence manifest.
No actual publication failure was observed in the accepted terminal journal.

## Evidence limits

This review performs no network calls, attachment streams, remote mutations,
publisher execution/imports, or scientific source execution. It verifies saved
evidence. Parent review supplies fresh metadata/listing observations and complete
metadata/historical preservation checks.

Stream events omit the content URL and HTTP headers. Their full-200, exact URL,
unencoded/nonpartial and complete EOF/checksum semantics follow from the pinned
controller, the scoped file-list sequence, and the saved receipt/listing
identities. This is not an independent second attachment fetch or remote
cryptographic attestation.

The accepted source flushes and fsyncs journal rows; the journal is not hash
chained and does not independently detect general forged or missing history.
Earlier GET captures and verification receipts reuse filenames, retaining the
latest JSON value; journal rows retain earlier response hashes and times.
Per-write response JSON objects use distinct names. JSON captures are parsed,
redacted and reserialized, so their saved-byte hashes cannot generally recreate
the raw HTTP response hashes recorded in the journal. No claim of account-wide
single-publisher exclusion or absence of unrecorded external actions is made.
