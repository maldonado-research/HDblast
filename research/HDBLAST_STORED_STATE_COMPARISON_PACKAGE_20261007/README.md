# HDBLAST stored-state comparison: portable package

This package preserves the registered scientific source, its exact input bytes, the initial actual normal/optimized outputs and the evidence needed to verify or repeat the finite stored-state comparison. The registered normal and optimized executions passed, as recorded in the included entry and custodian evidence. Read the included results report and figures for the actual comparison values and their interpretation.

The later portable extraction replay is a separate delivery check. Final package control hashes, publication identity and fresh replay evidence belong to the external delivery record. This README supplies neither a package seal nor a public authorization. At preparation of this text, final delivery sealing and the first portable extraction replay remain pending root evidence.

## Package contents and source identity

The registered source freeze is Git commit [`bedcca7e86995da1230c12a20fae3755b31f4a94`](https://github.com/maldonado-research/HDblast/tree/bedcca7e86995da1230c12a20fae3755b31f4a94/research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON). The complete scientific checkpoint lives at `research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON/`. It retains all 117 registered frozen files and `FULL_REGISTRATION.json`; post-execution result and evidence additions accompany those unchanged registered bytes.

The extracted package is intended to have this layout:

```text
HDBLAST_STORED_STATE_COMPARISON_PACKAGE_20261007/
  README.md
  replay_checkpoint.py
  PAYLOAD_MANIFEST.json
  REPLAY_PINS.json
  PUBLIC_GO.json
  research/
    HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON/
      FULL_REGISTRATION.json
      REGISTRATION_CONTRACT.json
      provenance/INPUT_SPEC.json
      execution/...
      decoder/...
      source/...
      theory/...
      transport/...
      RESULTS.md
      ... root-authored figures and result documentation
      ... registered proofs, documentation, reviews and actual evidence
    ... exact original repository paths for retained inputs and provenance
  evidence/actual/normal/...
  evidence/actual/optimized/...
  preparation/helper/...
  preparation/root-independent-review/...
```

The package includes all four original archives and all four compact incoming-state capsules. Their eight repository-relative paths remain exactly those declared in `provenance/INPUT_SPEC.json`. The original producer, input manifest and prior static audit also remain at their declared paths. The producer and input/archive bytes are authenticated as complete opaque files; the helper does not import the producer or regenerate arrays. Use the extracted package root as the worker's repository root.

The top-level `PUBLIC_GO.json` remains outside the scientific checkpoint. It binds the complete registration and reviewed source to the immutable public freeze. The package also carries original mathematical proofs, methods and literature summaries, internal reviews, actual saved outputs, and complete entry/custodian receipts and logs. Private source caches and third-party paper bodies are not required for replay.

`preparation/helper/` preserves the helper's preparation controls, documentation and independent review, including the earlier finding and corrected acceptance evidence. These manufactured controls test delivery safeguards; they are not actual study runs or public GO authorizations. `preparation/root-independent-review/` preserves the root's independent output and interpretation review. The root-authored `RESULTS.md`, figures and associated result documentation present the actual comparison and its limits; they accompany the frozen source rather than change its registered inputs or algorithms.

## Mathematical scope and normalization

The registered universe contains four separately identified source/resolution capsules, 49,152 capsule-node occurrences, 12 finite positional prefix cases and 20 selected exact-reader arrays. Capsule identities remain separate even when their momenta coincide. The exact reader preserves the represented binary80 values, including the retained momentum nodes and quadrature weights; it does not replace them with binary64 values, regenerated quadrature or continuum weights.

Saved first-order response components are `u_1` and `w_1`. Their normalized mathematical coordinates are

```text
U_saved = u_1 / epsilon_represented
W_saved = w_1 / epsilon_represented
```

Both divisions use the same exact represented nonzero epsilon specified by the frozen InputSpec. Exact represented momentum `k`, quadrature weight `dk` and represented Pi also come from that contract. Use those rational representations throughout; decimal approximations are not substitutes. The comparison convention is saved minus the certified prescribed Bunch–Davies prehistory target. No stored state is projected or renormalized.

The direct density and pressure quantities `R` and `P`, including their saved-minus-target differences, use the inherited `a0(eta)^4 / epsilon_represented` scaling of the first-order physical stress response. The background scale factor is `a0(eta)=L(eta)=-1/eta`; it is distinct from the incoming time anchor `a=-9/2`. A finite prefix sums the direct per-mode operators with the exact positive weights `mu_j=dk_j*k_j^2/(2*Pi_represented^2)`. Reported `R/P` bounds therefore have these normalized stress units and this retained finite-weight measure.

The target enclosure gate applies to every exported U rectangle and every exported W rectangle. Its L1 radius is

```text
[(Re_hi - Re_lo) + (Im_hi - Im_lo)] / 2 <= 1e-18
```

That gate concerns uncertainty in the certified target enclosure. The actual stored-minus-target errors are distinct quantities reported in the results and node exports. The fixed source, target construction, degree, precision, node universe and gates are defined by the registered contract and are not replay options.

The phase-uniform finite-weight transport bounds concern each selected discrete prefix throughout `eta in [-4.5,-3.5]`, with `L` ranging from `2/9` to `2/7`. They bound the incoming saved-minus-target difference under identical subsequent forcing and complete contacts, whose inhomogeneous differences cancel. They do not supply a continuum interpolation theorem, time or momentum quadrature remainder, later retained numerical trajectory enclosure, full contact/source bound or ultraviolet completion. The inherited full twelve-case pressure/contact continuum certificate remains **UNRESOLVED**; metric calibration remains **FAIL**. Higher-dimensional Big Bang causation remains **NOT_ESTABLISHED** and external novelty **NOT_ASSESSED**. Internal AI-assisted review is not external peer review or proof-assistant formalization.

## Runtime and externally supplied trust pins

Use Linux with Python 3.12.14, python-flint 0.9.0, sympy 1.14.0 and mpmath 1.3.0. Fresh execution additionally requires `libseccomp.so.2` and nonroot real/effective user IDs. Install the runtime separately using verified package sources. The helper performs no installation, download or network request. Recorded interpreter/distribution byte hashes are host provenance; runtime versions and the registered native guard requirements remain the portable prerequisites.

Set `PAYLOAD_SHA256`, `HELPER_SHA256` and `REPLAY_PINS_SHA256` to the exact nonsecret SHA256 strings supplied by the independently sealed external delivery record. Those external values authenticate `PAYLOAD_MANIFEST.json`, `replay_checkpoint.py` and `REPLAY_PINS.json`, respectively. Final control hashes are kept outside this manifested README to avoid circular seals.

Before any study import, the helper authenticates the complete package byte closure. It rejects missing or extra files, unregistered empty directories, unsafe paths, symlinks, hardlinks, special files, cached/native code and changing file identities. The package must remain unchanged during verification or replay. The frozen registration guard then authenticates the genuine external public GO and registered checkpoint. A local byte-integrity seal or a manufactured fixture does not authorize real source construction or retained-array decoding.

## Verify the included actual outputs

From the extracted package root, with the three externally supplied variables set:

The bootstrap below captures all three control files with the standard library, checks their externally supplied SHA256 values before executing any helper code, then executes the captured authenticated helper bytes. It does not load the helper from disk a second time. The helper subsequently authenticates the complete payload before importing study modules.

```sh
python -I -B - "$PWD" "$PAYLOAD_SHA256" "$HELPER_SHA256" "$REPLAY_PINS_SHA256" --verify-only <<'PY'
import hashlib, os, pathlib, re, stat, sys
root = pathlib.Path(sys.argv[1]).absolute()
for parent in (root, *root.parents):
    if parent.is_symlink():
        raise SystemExit("Symlink in package path")
controls = dict(zip(("PAYLOAD_MANIFEST.json", "replay_checkpoint.py", "REPLAY_PINS.json"), sys.argv[2:5]))
captured = {}
for name, expected in controls.items():
    if re.fullmatch(r"[0-9a-f]{64}", expected) is None:
        raise SystemExit("Invalid external SHA256")
    fd = os.open(root / name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or before.st_size > 8 * 1024 * 1024:
            raise SystemExit("Unsafe control file")
        raw = stream.read(8 * 1024 * 1024 + 1)
        after = os.fstat(stream.fileno())
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise SystemExit("Control file changed")
    if len(raw) != before.st_size or hashlib.sha256(raw).hexdigest() != expected:
        raise SystemExit("External control SHA256 mismatch")
    captured[name] = raw
entry = str(root / "replay_checkpoint.py")
sys.argv = [entry, "--package", str(root), "--payload-manifest-sha256", controls["PAYLOAD_MANIFEST.json"],
            "--helper-sha256", controls["replay_checkpoint.py"], "--replay-pins-sha256", controls["REPLAY_PINS.json"], *sys.argv[5:]]
exec(compile(captured["replay_checkpoint.py"], entry, "exec"), {"__name__": "__main__", "__file__": entry})
PY
```

This authenticates both saved custody receipts, including successful actual status, exact normal/optimized mode, authoritative outcome, fixed limits, entry/log links and measured resource bounds. It passes the guard's verified context to saved-output readback, checks both modes, compares their stable scientific bytes and reauthenticates the package afterward. It performs no new retained-array decode or physical source/target construction.

## Run a fresh normal and optimized replay

Choose a nonexistent output directory outside the package, with a real existing parent. The following path is illustrative; set it for your machine:

Use the same captured-byte bootstrap above, replacing its trailing `--verify-only` argument before the here-document with `--output-root /absolute/external/new-replay-directory`. Keep the package path and all three externally supplied pins unchanged. Perform the bootstrap checks before every invocation, including a fresh replay.

After verifying the saved package, the helper launches normal and optimized runs through separate isolated custodian processes. Credentials are not forwarded. The registered worker limits remain 900 seconds for wall/CPU time, 512MiB address space/RSS, 128MiB per-file and aggregate accepted output, one thread, no child processes and fail-closed syscall restrictions. Aggregate output is monitored and checked at acceptance; it is not a filesystem quota. No CLI option changes the numerical or resource contract.

Each fresh run must pass custody and saved-output checks. The exact scientific byte identities are `DATA.json`, `SCIENCE_SUMMARY.json`, `SOURCE_CERTIFICATE.json`, `DIAGNOSTIC_CROSSCHECK.json`, `NODE_TARGETS.jsonl.gz`, `SOURCE_ATTEMPTS.jsonl` and `DECODE_ATTEMPTS.jsonl`; `OUTPUT_READBACK.json` must also agree. Each fresh mode is compared with the saved normal output and the fresh modes are compared with one another. Timing, RSS, runtime/path observations and custodian execution provenance remain fresh evidence rather than scientific byte identities.

Fresh directories are exclusive. Failure evidence and logs are retained, and an earlier attempt is not overwritten. A successful fresh replay records its receipt in the chosen external output directory. Keep that receipt with the external delivery evidence; its result is not embedded in this package README. A passing preparation control or published archive alone does not establish that this fresh portable replay ran.
