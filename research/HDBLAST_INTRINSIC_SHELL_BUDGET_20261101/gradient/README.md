# Portable conditional actual shell-gradient campaign

This component re-evaluates the frozen nine-call campaign on the physical shell
`r=0`: eight whole cells covering `x in [31/25,63/50]`, followed by the fixed
target `x=5/4`, with `ell=10^-6`. No point, interval, or refinement changed.
Every complete metric, scalar, and background jet and every derived enclosure
is equal to the original frozen scientific result.

The intrinsic scalar formula and its angular factor were already derived in
`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/shell-invariant/SHELL_INTRINSIC_SCALAR.md`.
This component provides the frozen current actual-error enclosures and their
join to the selected exact-family and intrinsic-clock premises; it does not
claim a new scalar formula. The analytic exact-family and actual-field
existence/uniqueness/error premises remain conditional inputs.

At the fixed target, outward integers give

    -803004132630938124 <= mean Q2 <= -659922099878760930 < 0.

Across the entire closed interval, they give

    -1203216564043213024 <= mean Q2 <= -408804027527669570 < 0.

Here `Q2` is the second parameter derivative of the Lorentzian intrinsic
shell scalar `(D phi)^2`. It is only one contribution to a shell Gauss scalar
curvature budget. It does not determine full shell curvature, expansion,
radiation, force, or an effective finite-amplitude parameter radius. A
negative angular mean does not imply negativity at every angle.

## Run

Use Python **3.12.14**, with no third-party packages. In these commands
`COMPONENT` is this directory and `REPO` is a source root containing the two
public capsules named in `SOURCE_PINS.json`. No Git checkout metadata, original
private artifact, private manifest, or absolute private path is required.
Requested output directories must be fresh and outside both code and source.

    python3.12 -I -S -B COMPONENT/producer.py --repo-root REPO --output-dir NEW_RESULT
    python3.12 -I -S -B COMPONENT/audit_aggregation.py --repo-root REPO --report NEW_RESULT/REPORT.json --output-dir NEW_AUDIT
    python3.12 -I -S -B COMPONENT/test_controls.py --repo-root REPO --output-dir NEW_TESTS

Repeat with `-O` added to the interpreter options. Normal and optimized
producer reports are byte-identical, as are the independent aggregation
reports and the control-suite reports. Relocating both component and minimal
source closures also preserves those bytes.

The producer authenticates the closed public capsule inventories and every
member before scientific imports and after the nine calls. It compares every
scientific output to the pinned `FROZEN_REFERENCE.json`. The independent
aggregation auditor imports neither scientific nor producer code: it uses a
separate exact-Fraction interval implementation, reauthenticates all public
source members, verifies complete raw jets and coverage, and checks all
temporal/spatial squares, factors, denominators, ratios, and hulls. The final
report is also bound to its canonical byte hash.

`SOURCE_PINS.json`, `FREEZE.json`, and output report hashes changed solely to
remove private metadata and describe the portable producer. The original
freeze/report hashes remain provenance strings, never runtime dependencies.
The old campaign directory remains untouched. Upstream origin descriptions
are normalized to relative dependency names, lengths, and SHA-256 values in
the portable report.

## Retained evidence

`REPORT_NORMAL.json` and `REPORT_OPTIMIZED.json` are identical campaign outputs.
`AGGREGATION_NORMAL.json` and `AGGREGATION_OPTIMIZED.json` are identical audit
outputs. `TESTS_NORMAL.json` and `TESTS_OPTIMIZED.json` are identical strict
control reports. `COPY_PLAN.json` identifies the essential files, optional
evidence, and the public source closures. `MANIFEST.json` pins this component.

Each control suite runs four positive producer subprocesses and four positive
auditor subprocesses: normal and optimized, each in original and relocated
layouts. These are separate from nine negative cases in both modes:

| Controlled mistake | Mutation | Required final witness |
| --- | --- | --- |
| Second-derivative factor | Halve the computed K2, S2, Q2 entries | `SECOND_DERIVATIVE_FACTOR` |
| Angular eigenvalue 8 | Replace spatial contribution by 7/8 of its value | `ANGULAR_GRADIENT_FACTOR` |
| Physical shell radius | Replace exact shell label by bulk r=1/8 | `PHYSICAL_SHELL` |
| Ell scale | Replace background ell by 10^-3 | `ELL_UNITS` |
| Background induced metric | Replace retained rho by one | `BACKGROUND_JET_REFERENCE` |
| Scalar/metric sector | Substitute metric L,Lt for scalar s,st | `SCALAR_JET_REFERENCE` |
| Raw data integrity | Corrupt a retained metric A interval | `METRIC_JET_REFERENCE` |
| Report integrity | Add an unregistered field to the canonical report | `REPORT_INTEGRITY` |
| Source capsule integrity | Alter an oracle file in disposable source copy | `SOURCE_MEMBER_BYTES` |

Every rejection checks the exact final JSON object, exit code 2, empty stderr,
and absence of an output directory. A substring found elsewhere in a traceback
or earlier output cannot satisfy a control.
