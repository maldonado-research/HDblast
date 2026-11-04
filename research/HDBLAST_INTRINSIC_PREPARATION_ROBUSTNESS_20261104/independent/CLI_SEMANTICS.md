# Portable independent replay interface

This release preserves the accepted private rational calculation. It adds
portable paths, explicit source guards, complete report comparison and a
single isolated replay command. It creates no new direction or numeric
campaign and imports no producer code.

Every executable requires Python `-I -S -B`; `-O` is supported. Each requires
exactly one `--repo-root`, `--candidate` and `--output`. Paths must be absolute,
normalized and free of symlink indirection. The repository root and candidate
must exist. The output's parent must exist; the output must be absent.
Abbreviations, duplicated arguments, unknown arguments and `--help` are
rejected. `recompute.py` alone also accepts `--controls`.

The candidate has a specific role in each command:

| Program | Candidate |
| --- | --- |
| `recompute.py` | Preserved old `HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json`, checked against its fixed byte pin |
| `compare_candidate.py` | Portable bounds producer report, compared with the full frozen scientific schema and exact declared portable provenance |
| `witness.py` | Preserved old `HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/RAW_WITNESS.json`, checked against its fixed byte pin |
| `validate.py` | The frozen portable bounds producer report, SHA-256 `16efc39b531ea3aea14aa6157f3f804df5d335a34294aa7cae7fb3fb9a6d00d5` |

Success writes one canonical JSON report to a fresh output and emits one
complete one-line JSON receipt on stdout, with exit code 0 and empty stderr. Rejection
emits one complete JSON object with status `REJECTED`, exits 2, leaves stderr
empty and creates no new output artifact. An already existing file or symlink
is preserved. No result is inferred from the last stdout line: the validator
requires exactly one stdout line, parses its entire payload as one JSON object,
and checks each of eight rejection objects against its complete frozen
status-and-reason witness.

`compare_candidate.py` uses `importlib` to load its local `recompute.py` only
after checking the helper's source digest, both old source JSON pins and both
local JSON input pins. It checks the same thirty independent grouped
scientific consequences. It also requires equality of the entire frozen
scientific projection, including all gates and unknown-key rejection. Only
the four declared top-level provenance fields may differ from the private
report, and their portable values are themselves exact and byte-pinned.

Run `validate.py` once. It retains all six canonical artifacts in a sibling
directory named `<output-stem>.artifacts`, with names `BOUND_NORMAL.json`,
`BOUND_OPTIMIZED.json`, `COMPARISON_NORMAL.json`, `COMPARISON_OPTIMIZED.json`,
`WITNESS_NORMAL.json` and `WITNESS_OPTIMIZED.json`. It records explicit
`-I -S -B` flags, adds `-O` for optimized children, and checks byte equality
within each pair. It then executes exactly eight real CLI rejection controls.
Its control receipt does not require a second parent-wrapper mode run.

The six parsed-source and nine candidate-report in-process mutant cases are
separate from the eight actual CLI rejection processes. Across the canonical
mode pairs those case suites execute twelve and eighteen times respectively.
Canonical witness mode pairs replay one existing target analytically twice;
these are not new bump queries. All new background, bump, evolved-field
oracle and PDE-solve counts are zero. Root auditing of twenty-eight preserved
closed capsules checks byte closure and does not replay those old proofs.
