# Independent preparation certificate review

Both checkers use Python's standard library, exact `Fraction` arithmetic,
and explicit checks that remain active under `python -O`. They do not import
the numeric or theory producers, load scientific modules, or call oracles.

The cell checker rederives all 242 profiles from frozen raw input boxes,
checks the five initial entries and retained background bindings, and matches
41 exact scientific fields. Its fourteen in-process negative controls test
coverage, derivatives, adapters, Gram/error/history factors and the affine
witness. The coarse checker independently establishes the continuous bump
bounds and reconstructs 23 exact fields with eight in-process controls.

Run with fresh report paths outside the source repository:

```
python -I -S -B independent/reaggregate.py \
  --repo-root REPO --numeric-root CAPSULE/numeric \
  --raw NUMERIC_RUN/RAW_WITNESS.json \
  --producer-report NUMERIC_RUN/REPORT.json \
  --comparison-raw NUMERIC_OTHER_MODE/RAW_WITNESS.json \
  --comparison-report NUMERIC_OTHER_MODE/REPORT.json \
  --output FRESH_CELL_REPORT.json

python -I -S -B independent/check_coarse.py \
  --repo-root REPO --coarse-report COARSE_RUN.json \
  --output FRESH_COARSE_REPORT.json
```

Repeat with `-O`; each result must match its corresponding frozen reference
byte for byte. Successful stdout is one JSON line with status and field/control
counts. Failures produce one JSON line, empty stderr and exit 2. The output
must be fresh and outside the source repository; the cell output must also be
outside the numeric source component. Canonical reports contain no workspace
paths or elapsed timings.

The component manifest binds this entire directory. The root checkpoint must
authenticate it and the retained dependency closure before and after running
the checkers. The cell checker independently authenticates the frozen campaign,
witness rule and three numerical source files; the coarse checker authenticates
its two retained source files. A direct invocation does not independently
reprove or authenticate the entire inherited analytic/PDE framework.

Both methods certify the same nonzero slopes and fixed lambda=-10^21 under
the retained conditional framework. The background/operator/development and
gauge premises remain shared. The raw route also assumes the registered bump
jet enclosures; the elementary route supplies independent continuous bump
bounds. Neither route is a new PDE observation, external peer review, a
proof-assistant formalization or an effective nonlinear amplitude result.
