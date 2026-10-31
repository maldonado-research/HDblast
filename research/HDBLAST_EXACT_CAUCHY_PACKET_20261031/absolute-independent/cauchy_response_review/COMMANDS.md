# Portable independent Cauchy response audit

The checker and control driver use only the Python standard library. Every
root is supplied explicitly. `--help` performs no computation or writes.
Each output path or control output directory must be absent before a run.

Define the path variables for the relocated artifact layout, then run:

```sh
python check_cauchy_response.py \
  --source-root "$TASK_SOURCE_ROOT" \
  --same-response-root "$TASK_SAME_RESPONSE_ROOT" \
  --absolute-root "$TASK_ABSOLUTE_ROOT" \
  --velocity-root "$TASK_VELOCITY_ROOT" \
  --independent-velocity-root "$TASK_INDEPENDENT_VELOCITY_ROOT" \
  --producer-source "$TASK_CAUCHY_ROOT/checks/response_cauchy.py" \
  --producer-freeze "$TASK_CAUCHY_ROOT/checks/FREEZE.json" \
  --producer-report "$TASK_CAUCHY_ROOT/reference/RESPONSE.json" \
  --producer-report-sha256 22554f77c72f0f489bbd1daa5c40274411d0cd71faff71fe50c38cbbd2c814c2 \
  --output "$TASK_REVIEW_OUTPUT"
```

The optimized invocation is identical with `python -O`. The control driver
accepts all the same root, source, freeze and report arguments, replacing
`--output` with `--output-dir`:

```sh
python run_controls.py \
  --source-root "$TASK_SOURCE_ROOT" \
  --same-response-root "$TASK_SAME_RESPONSE_ROOT" \
  --absolute-root "$TASK_ABSOLUTE_ROOT" \
  --velocity-root "$TASK_VELOCITY_ROOT" \
  --independent-velocity-root "$TASK_INDEPENDENT_VELOCITY_ROOT" \
  --producer-source "$TASK_CAUCHY_ROOT/checks/response_cauchy.py" \
  --producer-freeze "$TASK_CAUCHY_ROOT/checks/FREEZE.json" \
  --producer-report "$TASK_CAUCHY_ROOT/reference/RESPONSE.json" \
  --producer-report-sha256 22554f77c72f0f489bbd1daa5c40274411d0cd71faff71fe50c38cbbd2c814c2 \
  --output-dir "$TASK_CONTROL_OUTPUT"
```

The executed control driver performs two normal/optimized successful
checks, two successful checks after sparse relocation of all selected
scientific and proof files, and 34 rejected subprocesses. The rejections
cover 13 separately rehashed report corruptions in both modes, three disk
source/table corruptions in both modes, and a report symlink in both modes.
Every rejected run has a nonzero exit and leaves no success artifact.

The canonical independent report is
`executed/INDEPENDENT_RESPONSE_NORMAL.json`; the optimized and both sparse
relocated reports are byte identical at SHA256
`c96bb97129c23d21852609f6f0dad88796e386a5f7778217be091a288c99038a`.
The checker source SHA256 is
`c388404391f1aafe804dfd8e708de8f9dcecc8186af0163fc5502a33188c4714`;
the control driver source SHA256 is
`ca794821321b01f0e3582c83ae8aa2eb74d3afb6346d54d42b56b57af3cdc746`.
The executed `CONTROLS.json` SHA256 is
`93b3c06eca0fc1aa6a2aa7985fd3413416b9bad98ab3709762cf6dc6cb6000f6`.

The report independently reproduces all 35 public response ledger nodes,
every operator and source constant, every advertised cap/error/role/flag,
all comparison ratios, the selected-source maps, and all nine formula or
stale-velocity witnesses. It authenticates 1489 selected files before and
after the calculation and imports no scientific module. Both independently
reaggregated local velocity witnesses use all 23 positive-Gram components.

The native velocity caps give rr upper approximately `3.699886190370885e24`.
Replacing only the velocity inputs with the separate independent velocity
enclosures gives rr upper approximately `3.699886213379224e24`. Both exact
rational ledgers pass each second-error upper at most half its own PR109
upper and the C1 no-increase gate. These are global upper-certificate
comparisons. The separately fresh local velocity witnesses provide their
own local upper-versus-lower comparisons. Retained conditional source,
coefficient and corner hypotheses remain necessary; no nonlinear shell
contact result is established.
