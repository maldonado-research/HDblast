# Exact-family domain and seed audit

The conditional gate passes. Read [DOMAIN_SEED_GATE.md](DOMAIN_SEED_GATE.md)
for the precise analytic qualifications, coordinate conversion, complete
support, and matching displacement/velocity/scalar seed.

The candidate is `ell/2<z<5ell/2`, `T_t=3ell/2`, `Z=6ell`, `c_star=2`,
with a fixed existential boundary-flow collar `0<h<min(flowradius,ell/4)`.
The strict causal margin is `ell/2`. In ambient ADM variables the exact
current profiles are `G=alpha q_ell/s^3`, `K=G/b`, with
`partial_t=R partial_ambient_T` on the common initial slice.

The source pin manifest SHA-256 is
`ce3d0ea33c758e2aafd54cdc6d8486054e198369f29e87c277a6eb60fe4c0b80`.
All 15 source paths are relative to an explicitly supplied repository root.
After authenticating that manifest, reproduce from this audit directory,
substituting the source checkout path:

```sh
python -I -B check_domain_seed.py --repo-root /path/to/checkout --pins-sha256 ce3d0ea33c758e2aafd54cdc6d8486054e198369f29e87c277a6eb60fe4c0b80
python -I -O -B check_domain_seed.py --repo-root /path/to/checkout --pins-sha256 ce3d0ea33c758e2aafd54cdc6d8486054e198369f29e87c277a6eb60fe4c0b80
python -I -B run_controls.py --repo-root /path/to/checkout --pins-sha256 ce3d0ea33c758e2aafd54cdc6d8486054e198369f29e87c277a6eb60fe4c0b80 --output-dir /tmp/hdblast-domain-seed-fresh
```

Python 3 and SymPy 1.14.0 are required. The control runner requires a new or
empty output directory outside the checkout, and writes `NORMAL.json`,
`OPTIMIZED.json`, copied mutants and `CONTROL_REPORT.json` there. Normal and
optimized receipts have identical scientific bytes. The control report records ten copied-source
mutation rejections and ten strict-inequality rejections across the two
modes. Each nominal run also verifies nine symbolic seed identities and
rejects three wrong clock, density, or scalar-velocity formulas.

The recovered exact-shell capsule's externally pinned manifest and all 51
members are reauthenticated before algebra. The analytic weighted inverse,
mixed evolution, finite propagation, and uniqueness remain written
conditional mathematical theorems. No finite-amplitude sign, effective
epsilon radius, transported observation, or change to a public artifact
is certified.
