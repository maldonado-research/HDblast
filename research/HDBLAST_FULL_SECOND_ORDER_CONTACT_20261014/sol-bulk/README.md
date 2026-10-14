# Independent scalar bulk continuation

The reviewable derivation is [FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md](FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md).
It gives the complete factorial-convention quadratic l=4 source, the general
two-base scalar harmonic linear operator, compulsory scalar gravitational
completion, and exact vanishing of coexact-vector forcing for the aligned
quadrupole. The first longitudinal gauge is fixed throughout.

The inherited single conformal master cannot solve the full second equation:
its linear angular scalar-TF residual is zero while the quadratic source is
generally nonzero. A trace repair alone also fails the temporal mixed row.
The note supplies the necessary extra amplitudes and their explicit sourced
constraint reconstruction, rather than claiming a proved scalar contact
theorem.

Run the standalone direct metric checker:

```sh
python checker/direct_metric_checker.py --include-controls --output checker/normal.json
python -O checker/direct_metric_checker.py --include-controls --output checker/optimized.json
```

Run the global exact angular decomposition and vector-selection checker:

```sh
python verify_angular_selection.py --output ANGULAR_NORMAL.json
python -O verify_angular_selection.py --output ANGULAR_OPTIMIZED.json
```

The direct checker requires only Python's standard library. The angular
checker requires SymPy; the frozen environment is Python 3.12.14 and SymPy
1.14.0. Both use exceptions rather than assertions and accept `--output`.
The direct checker's embedded controls must all reject; the standalone
angular `--wrong-formula` choices must each return exit code 1, also with
`python -O`.

The receipts distinguish exact local sample checks from globally symbolic
identities. They establish formula verification, not a shell boundary
estimate, constraint propagation theorem, two-master reduction, or full
nonlinear continuation. No inherited source file or Git state was modified.
