# Portable stored-state comparison replay

This helper belongs to delivery preparation after the separately registered scientific source freeze. It does not change the candidate/source registration, construct a public GO, download inputs, decode arrays or execute a physical source during preparation. The root has now completed the real prospective public-source GO and started actual normal execution. Actual saved normal/optimized outputs, their complete package seal and the first real portable replay remain pending root completion and review; this helper preparation does not read active actual outputs.

The final archive is planned as `HDBLAST_STORED_STATE_COMPARISON_PACKAGE_20261007`. Preserve the repository paths below inside that extracted directory:

```
replay_checkpoint.py
PAYLOAD_MANIFEST.json
REPLAY_PINS.json
PUBLIC_GO.json                       # outside the scientific checkpoint
research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON/...
research/...                        # all four original + four compact input files
research/...                        # producer, old input manifest and static audit
 evidence/actual/normal/...          # use final paths declared in REPLAY_PINS
 evidence/actual/optimized/...
```

The root package includes all eight required input files (212,975,098 uncompressed input bytes), the immutable registered checkpoint, original explanatory mathematics/literature, actual results and execution proofs, and the original standalone replay helper. Those inputs remain at their pinned repository paths, so the worker's `--repository-root` points to the extracted package. The producer and original capsule files are authenticated as opaque complete bytes; this helper never imports that producer or regenerates native arrays. This is independent of the numerical source reconstruction explicitly authorized by the registered public GO.

`PAYLOAD_MANIFEST.json` has exact schema `{"schema_version":1,"files":{"relative/path":{"bytes":123,"sha256":"..."}}}`. Include every package file except the three root control files `replay_checkpoint.py`, `PAYLOAD_MANIFEST.json` and `REPLAY_PINS.json`. Include the external public GO, all candidate files including `FULL_REGISTRATION.json`, all eight inputs, producer/old manifest/static audit, original reports, saved normal/optimized outputs and complete entry/custodian receipts/logs. Do not include private PDFs, third-party paper bodies, token/credential material, bytecode, native libraries or setup environments.

The three control files are authenticated through independent external SHA256 arguments. The manifest excludes them to avoid circular self hashes. Every other file must match its manifest size and SHA256. Extra files, empty directories, missing files, unsafe paths, symlinks, hardlinks, FIFOs/devices/sockets, bytecode/native code, unstable file identities and directory changes fail before study imports. Complete inputs are streamed for hashing in bounded chunks. The fixed limits are 256MiB per package member, 1GiB total manifest payload, 8,192 entries and depth24; these are package-authentication limits, distinct from the worker's resource limits.

Obtain all three external pins from the independently sealed release/immutable delivery record; do not derive expected trust pins from the extracted package. Start the package helper itself in isolated Python mode and suppress bytecode. From the extracted package, after setting `PAYLOAD_SHA256`, `HELPER_SHA256` and `REPLAY_PINS_SHA256` to those externally supplied nonsecret values:

```sh
python -I -B replay_checkpoint.py --package "$PWD" \
  --payload-manifest-sha256 "$PAYLOAD_SHA256" \
  --helper-sha256 "$HELPER_SHA256" \
  --replay-pins-sha256 "$REPLAY_PINS_SHA256" --verify-only
```

The exact required runtime is Linux, Python3.12.14, python-flint0.9.0, sympy1.14.0 and mpmath1.3.0, as pinned by the scientific registration. The helper is stdlib-only through complete package authentication; registered saved-output verification subsequently requires these three numerical packages. The recorded preparation executable/native package hashes are runtime provenance, not a cross-host byte equality requirement. Fresh execution also needs `libseccomp.so.2` and a nonroot real/effective user. Install dependencies separately using their exact versions and verified package/TLS sources; no installation or network request occurs in this helper.

`--verify-only` calls the frozen registration guard's actual `authenticate(candidate, registration_sha256, public_go, public_go_sha256)`. Before importing the saved reader, it independently checks each saved custodian PASS status, exact mode, strict integer zero exits/wait4 outcome, absence of stop/monitor errors, fixed limits, nonroot/measured resources, ENTRY_RECEIPT hash linkage and held child-log hashes. It then passes the authenticated verified context to `readback.validate_output` for both saved modes, compares their exact stable science files, and reauthenticates afterward. No retained array is opened for numerical decoding and no real source/target callback is invoked. The external public GO must be a genuine independently verified, complete immutable public-source receipt. A local integrity seal or manufactured fixture is not authorization, and there is no fabricated public-GO bypass in this helper.

For fresh replay, choose a nonexistent output directory outside the package, with a real existing parent:

```sh
python -I -B replay_checkpoint.py --package "$PWD" \
  --payload-manifest-sha256 "$PAYLOAD_SHA256" \
  --helper-sha256 "$HELPER_SHA256" \
  --replay-pins-sha256 "$REPLAY_PINS_SHA256" \
  --output-root /absolute/external/fresh-replay-directory
```

After saved-output authentication, the helper runs normal Python and optimized Python as two fresh custodian processes. The custodian is isolated and stdlib-only before its own authentication; it is not invoked in the helper process that previously loaded numerical saved-output verification libraries. Credential/environment bindings are not forwarded; the child environment contains only fixed PATH and locale settings. The registered launcher enforces 900 wall/CPU seconds, 512MiB address space/RSS, 128MiB file/aggregate output, NPROC0/core0, one thread, nonroot and fail-closed syscall filtering. Aggregate output is a monitored/final acceptance limit, not a filesystem quota. The helper does not expose resource/gate/node/degree overrides.

Fresh output directories are exclusive. Worker failures retain their authoritative `EXECUTION.json`, entry receipt and logs; the helper reports failure and does not overwrite a previous attempt. The supervisor's wait4 exit/signal/resource outcome remains authoritative. On helper interruption, SIGINT is sent to the custodian so its cleanup can terminate and reap the worker.

The exact scientific byte identities are `DATA.json`, `SCIENCE_SUMMARY.json`, `SOURCE_CERTIFICATE.json`, `DIAGNOSTIC_CROSSCHECK.json`, `NODE_TARGETS.jsonl.gz`, `SOURCE_ATTEMPTS.jsonl` and `DECODE_ATTEMPTS.jsonl`. The deterministic `OUTPUT_READBACK.json` is also compared. Each fresh mode must match the saved normal output and the other fresh mode. `RESOURCE.json`, runtime/path observations and `EXECUTION.json` timing/RSS/wait4 fields remain fresh execution evidence and are not byte identities. The entry receipt's complete output pins and the authenticated readback still validate each run's own resource file.

The final helper receipt is `PORTABLE_REPLAY_RECEIPT.json` in the fresh external directory. Publication or a passing helper test is not evidence that this real replay ran. Root must complete actual saved-output verification, complete payload/control seals, independent review and real normal/optimized extraction replay before assigning any release PASS status.

Preparation controls use small manufactured nonexecutable byte-closure packages, fake saved-export text and an explicitly invalid `NOT_AN_AUTHORIZATION_RECEIPT` public-GO file. They do not manufacture a usable GO or invoke real arrays/sources. All24 controls pass in normal and optimized Python. Independent review found that the first helper revision merely hashed saved custody receipts without enforcing their PASS/mode fields. That revision and its two failing reproductions are preserved; the correction validates both saved and fresh custody records before accepting their execution claims. Initial fixture bytecode contamination was caught by the helper, retained as a development failure log, and fixed by suppressing bytecode before fixture module imports. The production invocation already requires `-I -B`.
