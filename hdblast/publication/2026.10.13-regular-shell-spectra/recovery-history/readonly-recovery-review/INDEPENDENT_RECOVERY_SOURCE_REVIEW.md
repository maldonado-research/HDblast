# Independent read-only recovery source review

**Source assessment: PASS for the bytes pinned below.** Live recovery and the publisher child's offline controls were not completed at the time of this source assessment. A separate final review receipt will record the controls gate. No full 35-file public-readback success is claimed here.

Reviewer: `scalar_domain_20261013`, independent of the wrapper author. Date: 13 October 2026. Internal AI-assisted review is not external peer review.

| Source | SHA256 |
|---|---|
| `../readonly_public_recovery.py` (9683 bytes) | `55c0a5fc0c527caceb19d990f743da4f6a99a7ffe5185a35d6a0486793b2a458` |
| `../regular_shell_edition_controller.py` | `2c53e32a24fdc922d07f711cb67b8451da0be3325f740658d91d1edc0d97ca34` |

The complete wrapper and the relevant original controller paths were read: transport requests/streams, URL and metadata checks, original `stream_check`, `verify`, `current`, persisted-state and continuation bindings, publication/write methods, and evidence locking. No substantive unresolved source defect was found within this bounded task.

The wrapper hashes the base controller before loading its exact source bytes and requires its own externally supplied source pin. It does not edit the controller or create a base bytecode cache. Before network activity, it requires the existing fixed evidence directory, record 23335485 already marked published, concept 17088132, no pending write, and the fixed controller/inventory/metadata pins. The original controller constructor then verifies its full continuation binding and frozen inputs. The existing evidence lock prevents simultaneous use by the controller.

Remote writes are blocked at three relevant layers. The controller's write, create, prepare, publish, reconcile, owner-list, draft and draft-adoption entrypoints are replaced with immediate failures. The transport requires GET, no request body, `_public is True`, and no Authorization header. Its actual HTTP opener checks GET, no body and no Authorization again before dispatch and retains the base HTTPS destination checks. The base no-redirect opener, proxy handling and default certificate verification remain unchanged. The entrypoint still reads the configured token for the base transport's constructor and redaction; **this is not a zero-token-read runner**. The checked path does not transmit that token as Zenodo Authorization, and this review did not read its value.

The original `Controller.verify(published=True)` executes unchanged. It rechecks the prior/current public identities, latest record, metadata, complete 35-file inventory, each expected public content URL, and remote file/version identities. It requires 35 fresh completed body checks before writing the original public-readback receipt. It does not accept the older incomplete rounds or reuse their bytes.

Only `stream_check` adds retry behavior. For each file it performs at most three calls to the original `Controller.stream_check`, each with fresh counters, MD5 and SHA256 state. It retries only the exact reasons `CONTENT_STREAM_FAILED` and `CONTENT_COMPLETE_MD5_SHA256_DIFFERS`; other failures, including size overflow and identity/state/URL guards, stop. All three unsuccessful attempts stop the run. A one-second delay separates retryable failures. The original success event is emitted only after complete count/MD5/SHA256 equality, so a rejected or partial attempt cannot substitute for an accepted body.

Attempts are recorded separately in an append-only sidecar with run identifier, timestamp, wrapper/base/input hashes, record identifier, filename, attempt number, failure reason, retryability, observed byte/hash tuple, EOF indicator, and a bounded selection of HTTP headers when available. Failed attempt data is diagnostic evidence, not accepted content. The sidecar records overall completion with the exact original public-receipt hash, or terminal failure. Logging failures themselves cannot manufacture success.

Successful original stream events remain in the original journal, while retry observations stay in the sidecar. This lets the existing packet builder select one actual anonymous file-list event followed by 35 actual original success observations. Delivery must explicitly include the wrapper identity and failure/recovery history: the existing saved-witness checker does not itself authenticate the extra retry sidecar. A completed recovery requires the original 35-file receipt **and** the wrapper's matching completion evidence. Root/operator delivery should not silently present this as an uninterrupted first attempt.

`PRE_RECOVERY_READONLY_SNAPSHOT.json` independently records the local pre-run state: unchanged base controller, published record 23335485, no pending write, 206 original journal events, nine write intents, exactly one publish intent and one publish response, and no completed public-readback receipt. It is a baseline for the later evidence comparison, not a new network observation.

The review does not identify the cause of the earlier generic failures. A later byte-exact success would establish the accepted response's integrity; it would not establish that all earlier mismatches or transport failures were harmless, nor prove storage corruption. This is bounded read-only recovery, with no scientific, repository, approved APS, or proof-script changes.
