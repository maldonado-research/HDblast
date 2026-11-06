# Replay the conditional quantitative evolution inputs

Python 3.12 or later suffices. The replay uses only the standard library,
including exact `Fraction` arithmetic and formal sparse rational
polynomials. It recomputes all finite calculations from fifteen immutable
public snapshots. It imports no private producer report or saved output.
The retained first-source Taylor engine is itself an original public
snapshot whose source bytes are authenticated before it is loaded.

The report has 447 named checks: 45 source hash/size/Git-blob checks,
248 coefficient checks, 65 formal gauge polynomial identities,
37 gauge norm checks and 52 energy checks. These finite checks support
the conditional analytic arguments in the proof notes. The reduced metric
and moving point-map tests use exact rational test jets; the general
analytic reduction is stated and derived in `COEFFICIENTS.md`.

Obtain the complete capsule at the reviewed commit and the external
`MANIFEST.json` SHA256 pin from
`.github/workflows/evolution-bounds-20261106.yml` at the same commit.
Set `PIN` to that literal 64-character value. A pin recomputed from
untrusted inputs is not an authenticity anchor. From the repository root,
choose fresh output paths outside the closed capsule:

```sh
python -I -S -B research/HDBLAST_EVOLUTION_BOUNDS_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/evolution-bounds-normal.json
python -O -I -S -B research/HDBLAST_EVOLUTION_BOUNDS_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/evolution-bounds-optimized.json
cmp /tmp/evolution-bounds-normal.json /tmp/evolution-bounds-optimized.json
python -I -S -B research/HDBLAST_EVOLUTION_BOUNDS_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/evolution-bounds-controls-fresh
```

Both modes must produce byte-identical canonical JSON. The checker
binds the manifest to the external pin, rejects unlisted or missing files
and every symlink, authenticates every listed file, and then loads the
pinned evaluator. JSON duplicate keys and nonfinite constants are rejected.
Output paths must be fresh and outside the capsule. `MANIFEST.json`
lists every capsule file except itself.

The controls directory must not exist beforehand. Its receipt records a
passing baseline and six intended failures: changed source bytes, an extra
file, a missing file, a symlink, a stale external manifest pin, and an
internally rehashed false mathematical input. The last changes the actual
square-root energy forcing multiplier in `DATA.json` from `95` to `94`,
rehashes the capsule, and passes the new pin. It must reach and fail the
computed inequality `multiplier^2 > 1/(2*lambda)`, with `lambda=1/18000`,
specifically `Arithmetic failed: energy:95_square_root_force_multiplier`.
Since `94^2=8836<9000`, this is a false mathematical claim that passes
source and manifest authentication before arithmetic rejects it.

Coefficient jets concern fixed background and pointwise state variables;
they do not bound evolved spacetime first jets. The initial Cauchy ceiling
is in the explicitly defined positive unstraightened tensor norm. The
energy propagator still requires complete transformed lower-order and
first-jet budgets. Finite Sobolev products, commutators, normal recovery,
actual-family bootstrap, a nonzero amplitude interval, and the common
proper-time/actual-area observation domain remain unclosed.
