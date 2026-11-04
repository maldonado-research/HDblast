# Frozen preparation-response numerical component

This component replays the registered one-direction `P=1` experiment at
`T=5/4`, `ell=10^-6`. It makes no background or evolved-first-field oracle
calls. It regenerates only the fresh bump/profile boxes and aggregates the
full homogeneous five-mean response bound.

```
python numeric/reproduce.py --repo-root REPO --output-dir FRESH_OUTPUT
python -O numeric/reproduce.py --repo-root REPO --output-dir FRESH_OUTPUT_O
python numeric/run_controls.py --repo-root REPO --output FRESH_CONTROLS_JSON
```

The output directory must be fresh and outside the repository and this
component. `REPO` is a source tree containing the retained `research/`
capsules. No Git metadata, HEAD value, environment-specific location, or
private worktree is used. The enclosing checkpoint authenticates the full
retained closure; this component authenticates every direct input listed in
`SOURCE_PINS.json` and every local scientific member in `LOCAL_PINS.json`
before scientific imports, then authenticates them again after calculation.

The original `CAMPAIGN_FREEZE.json` and `POSTPROCESS_FREEZE.json` are byte
unchanged. The portable refactor separates generation from aggregation without
changing the direction, boxes, exact endpoints, or registered sign tests.
`references/ORIGINAL_REPORT.json`, `ORIGINAL_RECEIPT.json`, and
`RAW_WITNESS.json` preserve the original artifacts exactly.
`SCIENTIFIC_REFERENCE.json` is the complete original report with only
`source_authentication` excluded. Every one of its fields must match exactly.
Portable reports explicitly version the replacement authentication/provenance
and retain the original artifact hashes. Their raw witness bytes remain
identical to the original raw witness. Normal and optimized outputs have the
same canonical report, raw witness, and receipt bytes.

The real negative controls alter scalar preparation, initial norm, resolvent
feedback, first forcing cancellation, target error, kernel comparison, proper
time factors, and lambda sign in the arithmetic. Each runs the actual CLI as a
subprocess, reaggregating the authenticated raw witness without regenerating
profile boxes. Each must exit 2 with a single exact JSON witness, empty stderr,
and no output artifact. These aggregate-only controls are distinguished from
the positive fresh box generation.

The raw protocol records shared `rho_b,p_y` and shell intervals, and for each
of 241 continuous cells and one target records `rho,p,a`, `Q` and three
dimensionless derivatives, `h,h_r,V_L,s,s_r`, and all five adapter vectors.
An independent checker can reconstruct formulas and all aggregate endpoints
using only standard-library rational arithmetic; it need not import this
producer or call an oracle.

See `NUMERICAL_REPORT.md` for the result and conditional scope. The finite
`lambda=-10^21` witness reverses the signs of the formal volume and expansion
second coefficients. No effective nonlinear epsilon radius is provided.
