# Exact replay inventory

The three current checkers evaluate the following exact named inventory in
each of normal Python and optimized Python. A tensor identity counts once
when its report groups its components; the counts do not inflate those rows.

| Checker | Positive algebra checks | Wrong-formula controls | Reviewed result |
|---|---:|---:|---|
| Nonlinear shell Noether algebra | 6 | 8 | [noether.json](results/noether.json) |
| Independent rotational constraints | 46 | 9 | [rotational.json](results/rotational.json) |
| Quadrupole admissibility | 29 | 4 | [quadrupole.json](results/quadrupole.json) |
| Total per mode | **81** | **21** | [exact contracts](CHECK_CONTRACTS.json) |

Each positive row passes its exact identity or registered nonzero witness
criterion. Each wrong-formula row detects the deliberately incorrect
expression. The Noether and quadrupole CLI failures alter a reviewed corner
sign and rotation identity respectively. The rotational CLI failure exercises
an explicit nonzero guard; it is separate from that checker's nine scientific
wrong-formula controls and is not counted in this table.

In the fixed radial-profile, time-scale, and reflection-parity quadrupole
family, the four boost tests vanish by antipodal parity, and the six rotational
tests vanish exactly when the two harmonic matrices commute. This classifies
all ten tests in that restricted family.

The replay runner requires byte-identical complete scientific receipts in
both modes and rejects a failed check even if a checker writes an explicitly
failed JSON report. It tests copied-proof corruption, external-manifest-pin
corruption, and each fixed prior-manifest pin after resealing only the negative
test's outer manifest. All test outputs go outside the capsule. A successful
complete replay writes `REPRODUCTION_RECEIPT.json` there; the frozen input
capsule and its relocated authentic copy must remain unchanged.

These computations support the local Noether algebra and the second-order
restriction argument. The analytic assumptions and proof steps remain
reviewed mathematical statements. Passing the rotational tests does not
establish nonlinear continuation, sufficiency, physical energy positivity,
or a cosmological prediction.
