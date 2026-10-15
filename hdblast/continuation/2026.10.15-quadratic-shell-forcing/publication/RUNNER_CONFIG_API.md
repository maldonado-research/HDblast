# Build-time runner configuration

`build_replay_runner.py --config publication/RUNNER_CONFIG.json --source-root .
--output-dir publication/replay-tools` reads existing scientific files and complete
goldens without executing them. Use `--replace` explicitly for an accepted
prefreeze revision. It writes `reproduce.py`, `REPLAY_ROSTER.json`, and a separate
`RUNNER_BUILD_RECEIPT.json`. The last file belongs in publication evidence unless
root deliberately includes it in source. The configuration is embedded as literal
Python `CHECKERS` and `CONTROLS`; runtime never reads a configuration file.

Required top-level fields:

- `schema`: `hdblast-build-time-replay-config-v1`.
- `status`: `APPROVED_FOR_FROZEN_REPLAY` (requires actual reviewer acceptance).
- `capsule`: `HDBLAST_QUADRATIC_SHELL_FORCING_20261015`.
- `python_version`: the tested exact version, currently `3.12.14`.
- `dependency_versions`: exact installed distributions, e.g. `{"sympy":"1.14.0"}`.
- `checkers`: nonempty list below; total enumerated groups must be positive.
- `controls`: nonempty list of explicit external formula/input controls below.

Each checker requires `name`, `script`, `golden` and exact positive `status`.
All source paths are canonical relative paths under `--source-root`.

Optional checker fields:

- `dependencies`: copied source/input paths, excluding the checker itself.
- `arguments`: CLI string list; default `[]`.
- `output`: fresh result JSON path inside the copied workspace; `null` means
  complete stdout JSON. It cannot overlap any copied input or reserved logfile.
- `optimized_golden`: separate golden only when execution metadata differs.
- `status_path`: JSON key/index path; default `["status"]`.
- `byte_equal`: require generated JSON byte equality as well as complete strict
  JSON equality. Use true for deterministic producers.
- `additional_outputs`: list of `[fresh_relative_path, source_golden_path]` pairs;
  every complete output is compared byte for byte, including full CSV files.
- `checks` / `embedded_checks`: `null` for no enumerated inventory, otherwise
  `{"path":["check_names"],"kind":"names"}` for unique string lists or
  `{"path":["checks"],"kind":"true_map"}` for dictionaries of named `true`
  results. Actual group/control names and counts are derived from the complete
  normal golden and must match the optimized golden exactly. Invented/manual
  group totals are not accepted.
- `check_count_path` / `embedded_count_path`: optional JSON path to a declared
  integer total, independently matched to the actual named inventory.
- `mode_fields`: explicit execution metadata only, e.g.
  `[{"path":["runtime","optimization"],"normal":0,"optimized":1}]`.
  Each value is type-checked in its mode; only these fields are removed for the
  normal/-O scientific comparison. They are never removed for golden comparison.
- `required_fields`: additional exact typed gates, e.g.
  `[{"path":["cell_data","rows"],"equals":4353}]`.

Each control requires `checker`, `name`, exact `reason`, and `mechanism`.
The mechanism is copied identically into the frozen roster and process receipt.
It also requires `failure_format`: `fail_line`, `traceback`, or `json`:

- `fail_line` requires the entire diagnostic to be `FAIL: REASON\n`.
- `traceback` requires Python traceback prefix and the exact final
  `EXCEPTION: REASON` line including its final newline. Specify `exception`.
  Optional `diagnostic_prefix` permits only those exact configured bytes before
  the traceback header; default is empty. Extra leading/trailing diagnostics,
  missing prefixes, wrong reasons and trailing blank lines are rejected.
- `json` requires complete stdout JSON with exact failure status (`FAIL` by
  default, optional `status` and `status_path`) and exact reason at the required
  `failure_path`. It must still exit1 and leave any expected result file absent.

A CLI formula injection supplies distinct `arguments`. A copied-input corruption
supplies `target`, which must be a declared dependency rather than the checker,
and exactly one of:

- `append`: nonempty literal text.
- `before` plus `after`: exactly one occurrence of changed formula text.
- `json_path` plus `value`: exact JSON field replacement.
- `json_path` plus rational-string `factor`: scales one scalar or list of rational
  strings, emitting rational strings. It must change the input bytes.

No corruption may be counted from an unrelated exception, a mere nonzero exit,
or a stale passing result. Both modes execute each positive/control identity in
fresh independently authenticated workspaces. Their log/result records retain
`path`, `bytes`, `sha256`; extra outputs preserve the existing evidence schema.

The final builder also requires process fields `finite_positive_check_names`
and `embedded_control_names` to match the frozen roster names exactly. Named
inventories, not just matching totals, are therefore retained end to end.
