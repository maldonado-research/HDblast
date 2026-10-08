# Retained endpoints and causal bulk response: completed research

The [results](../../../research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW/RESULTS.md)
enclose all 98,304 saved later-endpoint comparisons from the same exact
represented incoming state, with 24 finite weighted cases. Maximum Cartesian L1
discrepancies satisfy U<2.520e-16 and W<3.158e-15. All 48 signed finite R/P
intervals exclude zero. The 1e-18 gate concerns target-enclosure precision;
complete target radii satisfy U<7.574e-29 and W<1.263e-28.

Original normal/optimized runs and one independently authenticated fresh portable
pair reproduce all nine scientific files exactly. Fresh normal execution took
195.074 seconds and 71,864 KiB peak RSS; optimized took 193.158 seconds and
70,308 KiB. Both used the unchanged 900-second, 512-MiB and 128-MiB worker caps.
See [independent portable acceptance](portable-replay/FINAL_PORTABLE_REPLAY_REVIEW_RECEIPT.json)
and [the complete compact evidence manifest](portable-replay/FINAL_PORTABLE_REPLAY_REVIEW_MANIFEST.json).
The source replay package remains byte-for-byte sealed; these later witnesses
are deliberately outside its authenticated closure.

All 103 package files, the scientific manuscript/presentation and both bulk-model
packets passed anonymous complete public streams at immutable commit
`52126b98d3e4b24680017fa2537b8449bd75d94f`: 155 files / 75,623,575 bytes.
The numerical CI at that commit also passed. The first bulk workflow had a YAML
format error; its [preserved failure and successful correction](GITHUB_WORKFLOW_FAILURE_001.json)
do not alter the research bytes. Final delivery checks are tracked separately.

The [conditional scalar model](../../../research/HDBLAST_CHECKPOINT_20261008_CAUSAL_BULK_RESPONSE/README.md)
derives a stable half-space response and an exact four-dimensional continuum
replica of its reduced Gaussian measurements. Its reviewed incident-packet
extension reflects the incoming energy: pure continuum packets, with both
bound-mode projections absent, leave no permanent local brane energy in this
static linear model. A mechanism for lasting energy transfer still needs a
specified dynamical or gravitational action and a complete energy ledger.

Continuous numerical paths, momentum/contact/UV completion remain UNRESOLVED;
historical metric calibration FAIL; higher-dimensional Big Bang origin
NOT_ESTABLISHED; external novelty NOT_ASSESSED. Internal AI-assisted checks are
not external peer review. The original source freeze, input bytes, failed
attempts and older DOI records retain their historical identities.

## Authenticate and read back the saved actual package

Use Linux/nonroot Python3.12.14 with python-flint0.9.0, SymPy1.14.0,
mpmath1.3.0 and libseccomp.so.2. The eleven original input/provenance leaves,
213,114,777 bytes, originate at commit
`13a30ef46f5c90c1b01ae83b609db65fd0f8a709`; a later checkout is accepted only if
every consumed leaf has exactly the declared bytes. No original NPZ is copied
into the new source bundle. The package preparation README's statement that
fresh replay was not yet executed remains a dated assembly fact; the independent
evidence above establishes its later completion.

These external SHA256 values were checked against the independent review and
immutable public bytes:

- Payload manifest: `976e92a7a14ecb13ed04cbb9f7ed15fff3a285bd40175a28929dff53677d0752`.
- Replay helper: `5de25511aea5579b7e0ebe94208ee1584d0869ca248fa8939cdbb29ed25ff6fd`.
- Replay pins: `1407dc4bcf0a312d1bc23a29e93e3c1149071f4ce2705ccbbeef430e82446d68`.
- Captured bootstrap: `438a0478fa185f46ca06372036caa7aee1f9e9659a2778b9d1a9e98651122f0a`.

Capture/authenticate the bootstrap before executing it. With the registered
environment Python, the following performs saved-output readback only:

```sh
python -I -B - <<'PY'
from pathlib import Path
import hashlib, sys
repo = Path('/workspace/HDblast')
package = repo / 'research/HDBLAST_RETAINED_ENDPOINT_FLOW_PACKAGE_20261008'
entry = package / 'bootstrap_replay.py'
raw = entry.read_bytes()
if hashlib.sha256(raw).hexdigest() != '438a0478fa185f46ca06372036caa7aee1f9e9659a2778b9d1a9e98651122f0a':
    raise SystemExit('External bootstrap hash differs')
sys.argv = [str(entry), '--package', str(package),
    '--helper-sha256', '5de25511aea5579b7e0ebe94208ee1584d0869ca248fa8939cdbb29ed25ff6fd',
    '--payload-manifest-sha256', '976e92a7a14ecb13ed04cbb9f7ed15fff3a285bd40175a28929dff53677d0752',
    '--replay-pins-sha256', '1407dc4bcf0a312d1bc23a29e93e3c1149071f4ce2705ccbbeef430e82446d68',
    '--', '--repository-root', str(repo), '--verify-only']
exec(compile(raw, str(entry), 'exec'), {'__name__': '__main__', '__file__': str(entry)})
PY
```

Expected status is PASS_VERIFIED_SAVED_RETAINED_ENDPOINT_OUTPUTS
with zero new physical callbacks or retained-array decodes. Fresh repetition of
the unchanged study is documented by the package helper and its reviewed
bootstrap; routine environment setup performs saved verification only.
