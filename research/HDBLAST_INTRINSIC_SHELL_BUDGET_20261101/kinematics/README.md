# Actual shell normal-clock intrinsic budget

This component certifies a positive second parameter derivative of the
volume-normalized intrinsic three-curvature defect and a positive Buchert
kinematic coefficient on the full scaled-time interval `[31/25,63/50]`,
including the selected event `(x,r)=(5/4,0)`. The result is conditional on the
retained exact-family and actual-first-field premises. It uses the actual
shell scalar source, its complete inherited errors, and the whole history
from zero. It uses no second-response field.

Run from any directory, with a source repository containing the two closed
capsules named in `SOURCE_PINS.json`:

```sh
python -I -S -B check_shell_budget.py --repo-root /path/to/repo --output-dir /path/to/fresh-output
python -I -S -B -O check_shell_budget.py --repo-root /path/to/repo --output-dir /path/to/another-fresh-output
```

The output directory must be fresh and outside both this code directory and
the source repository. Each ordinary production run makes exactly 73 actual
oracle calls: 64 history cells, eight continuous target cells, and one selected
point. It writes `ACTUAL_SHELL_BUDGET.json`. Normal and optimized reports are
byte-identical and contain no absolute source, code, or output paths.

The recorded reports are in `results/`. Their shared SHA256 is
`42949cc9f4dc1c2377ac5b6497d79712c0c0d8550a524ad04d83cb799f7df475`.
The exact scientific projection has SHA256
`99d48f7340e347f0846405d9e0afc2025e223bc77c169ba3e906f01b150f85ad`.
This is exactly the frozen private trace campaign's mathematical projection;
portability changes only metadata and file organization.

Coarse outward enclosures throughout the continuous interval are:

| Actual second parameter derivative | Enclosure |
| --- | --- |
| Intrinsic curvature defect `C2` | `[3.73e-24, 1.024e-22]` |
| Expansion variance `Var2` | `[2.55e-9, 2.953e-8]` |
| Buchert coefficient `Q_B2` | `[1.70e-9, 1.969e-8]` |
| Mean squared shear `S2` | `[0, 1.734e-33]` |

The exact rational endpoints and independently checked outward decimal
endpoints are in the reports. The shear enclosure does not certify a strict
positive mean squared shear on the whole interval.

`METHOD_AND_SCOPE.md` gives the units, history bounds, formulas, limitations,
and protocol history. The preserved protocol includes the original global
master method that contained zero, the coarse trace pilot that contained
zero, and the successful 64-cell trace assessment. None is presented as the
final fresh-source certificate.

Run all strict subprocess controls once; the driver tests both interpreter
modes and relocates the code and both source capsules:

```sh
python -I -S -B run_controls.py --repo-root /path/to/repo --raw-report results/ACTUAL_SHELL_BUDGET_NORMAL.json --output /path/to/fresh-controls.json
```

The driver makes two fresh positive campaigns, 73 actual oracle calls each,
and twenty strict rejection subprocesses with no new oracle calls. The
`--raw-report` and `--inject-error` producer options support this explicit
review/control path; ordinary production uses the first two commands above.
