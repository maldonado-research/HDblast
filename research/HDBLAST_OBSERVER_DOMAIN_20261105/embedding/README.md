# Conditional global observer-domain criterion

5 November 2026, UTC. Research note.

[OBSERVER_DOMAIN.md](OBSERVER_DOMAIN.md) gives explicit sufficient conditions
for global noncrossing and spacelike cuts, conditional on an existing regular
observer flow on the entire proper-clock cylinder. It proves embedding across
different clocks, not merely on each cut. It supplies the compatible spectrum
`(1/2)I < G_e < 2I` under symbolic input inequalities.

Flow existence, tube retention, full physical shell fields (including the
moving embedding), and uniform input budgets remain assumptions. No PDE
theorem, evaluated nonlinear constant, or effective numerical amplitude is
supplied. The global representation uses `S^3 subset R^4`, not a global chart.

The continuous proof is primary. The standard-library Fraction checker
supports its finite algebra; it does not certify the analytic assumptions.
Unaltered public source copies in `inputs/` are pinned by `SOURCE_PINS.json`
with immutable source URLs, repository paths and byte identities. They comprise
the retained README, observer-flow budget, area-density bound and source
dependency note. `MANIFEST.json` closes the payload inventory.

Run with fresh absolute output filenames outside this capsule:

```
python3 /ABS/capsule/check_observer_domain.py --output /ABS/reports/normal.json
python3 -O /ABS/capsule/check_observer_domain.py --output /ABS/reports/optimized.json
cmp /ABS/reports/normal.json /ABS/reports/optimized.json
```

The output parent must exist. Existing files, relative or symlink paths,
and destinations inside the closed capsule are rejected. Success stdout is
exactly the canonical JSON report written to the output file, including its
final newline. The checker uses no assertions, floats, third-party packages,
network, subprocesses, legacy scientific programs, or oracle output.
