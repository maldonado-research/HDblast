# Replay the conditional boundary-flow and metric patch

Python 3.12 or later suffices. The standard-library replay independently
evaluates exact rational formulas from nine immutable public input
snapshots. It imports no private producer proof or precomputed report.
The report contains 220 flow checks and 119 metric checks, including source
hash/size checks and retained formula markers. These finite checks support
the analytic arguments in the proof notes; they do not evolve the family.

Obtain the complete capsule at the reviewed commit and the external
`MANIFEST.json` SHA256 pin from
`.github/workflows/boundary-flow-patch-20261106.yml` at that same commit.
Set `PIN` to that literal 64-character value. A pin recomputed from
untrusted inputs is not an authenticity anchor. From the repository root,
choose fresh output paths outside the closed capsule:

```sh
python -I -S -B research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/boundary-patch-normal.json
python -O -I -S -B research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/boundary-patch-optimized.json
cmp /tmp/boundary-patch-normal.json /tmp/boundary-patch-optimized.json
python -I -S -B research/HDBLAST_BOUNDARY_FLOW_PATCH_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/boundary-patch-controls-fresh
```

Both modes must produce byte-identical canonical JSON. The checker
authenticates the manifest against the externally supplied pin, rejects
unlisted or missing files and all symlinks, and authenticates every listed
file before loading the evaluator. Output paths must be fresh and outside
the capsule. `MANIFEST.json` lists every capsule file except itself.

The controls directory must not exist beforehand. Its receipt records a
passing baseline and six intended failures: changed source bytes, an extra
file, a missing file, a symlink, a stale external manifest pin, and an
internally rehashed false mathematical input. The last changes the actual
inner-image radius from `101/100000` to `1/100`, rehashes the capsule and
passes the new pin. It must reach and fail the computed reconstructed
physical-metric bound `<1/40`, specifically the diagnostic
`Arithmetic failed: metric:literal_inner_state_embeds_full_metric_patch`.
It is not a metadata-only failure.

The flow is the literal background-subtracted proper-y flow with `B(0)=0`.
Its mixed normal-coordinate/state jets are auxiliary chart derivatives.
Metric energy and causal-speed margins are conditional on the pointwise
physical patch and exact frozen-iterate boundary traces. No evolved norm,
nonlinear amplitude interval, observer clock, remainder bound, common
observation cylinder or physical confirmation follows from this replay.
