# Final replay roster contract

Do not create an approved roster until the scientific scripts, output schemas,
meaningful corruption controls and copied data dependencies are stable.
`REPLAY_ROSTER.json` belongs inside the assembled scientific source and is
therefore authenticated by its external manifest pin. The evidence builder
reads that committed source file; a separate unpinned count setting is rejected.

Required fields:

```json
{
  "schema": "hdblast-replay-roster-v1",
  "status": "APPROVED_FOR_FROZEN_REPLAY",
  "capsule": "HDBLAST_COORDINATE_SHELL_SOURCES_20261015",
  "modes": ["normal", "optimized"],
  "counts": {
    "successful_checker_runs": 0,
    "external_corruption_processes_rejected": 0,
    "embedded_wrong_formula_controls_rejected": 0,
    "finite_positive_check_groups_each_mode": 0
  },
  "positive_processes": [],
  "external_controls": []
}
```

The zeros above are documentation placeholders, not an approved roster. The
builder rejects zero positive processes or finite positive groups. Replace all
counts with actual verified counts when the final runner is accepted. The required `positive_processes` list names each checker and script with its
`finite_positive_groups_each_mode` and `embedded_controls_each_mode`, plus exact
`finite_positive_check_names` and `embedded_control_names` lists. These names
are derived from complete source goldens at build time and verified at runtime.
The required
`external_controls` list gives each `checker`, `name`, `expected_reason` and
`mechanism`. The runner compares both lists exactly with its execution plan;
the evidence builder binds every result to those approved names, scripts, group
counts and failure reasons. Additional explanatory top-level fields are allowed.

The fresh receipt retains the previous pipeline's `steps` structure: each
process has `kind` (`positive` or `external_corruption`), `checker`, `mode`,
`exit_code`, `log` and `output`. Result/log records contain `path`, `bytes` and
`sha256`. Extra generated files use `additional_outputs: [{"output": record}]`.
Every positive process has `scientific_output_matches=true`,
`finite_positive_check_groups` and `embedded_wrong_formula_controls`, plus
`finite_positive_check_names` and `embedded_control_names` matching the frozen
roster exactly. Every
external control has a unique `control`, an `expected_reason`,
`expected_reason_detected=true`, and `expected_failure_exit_code=1`.
Its `control_mechanism` must equal the frozen declaration's `mechanism` exactly;
mechanism aliases from older checkpoints are not accepted.

The builder independently tallies every process and both execution modes,
requires the same process/control identities in both modes, rejects duplicates
and incorrect exits, and verifies the replay's manifest, runner and verifier
pins against the published immutable source. Numerical checks, conditional
analytic conclusions and unavailable nonlinear/observational quantities remain
separate. Future generated full JSON/CSV output must be compared completely.
