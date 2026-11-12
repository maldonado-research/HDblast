# Reproduce the forced second-wave checkpoint

Python 3.12 or later with only the standard library suffices. Obtain the
literal external manifest pin from `.github/workflows/second-wave-evolution-20261112.yml`
at the same reviewed commit and set PIN to that SHA256. A pin recomputed from
untrusted inputs is not an authenticity anchor. From the repository root,
use fresh output paths outside the closed capsule:

```sh
python -I -S -B research/HDBLAST_SECOND_WAVE_EVOLUTION_20261112/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/second-wave-normal.json
python -O -I -S -B research/HDBLAST_SECOND_WAVE_EVOLUTION_20261112/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/second-wave-optimized.json
cmp /tmp/second-wave-normal.json /tmp/second-wave-optimized.json
python -I -S -B research/HDBLAST_SECOND_WAVE_EVOLUTION_20261112/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/second-wave-controls-fresh
```

The 438 named checks comprise 24 source size/SHA256/Git-blob checks, 218
formal identities and 196 finite constant/allocation checks. The evaluator
exhausts all 74,256 mixed source derivative allocations. Both modes must
produce byte-identical canonical JSON. Changed/missing/extra files, symlinks,
stale manifest pins, duplicate JSON keys and nonfinite constants are rejected
before code imports. Eight immutable public snapshots are the only source
inputs; no private producer report or saved numerical result is imported.

The formal engine uses a two-fiber time/space Taylor system with independent
nonzero source jets to check twelve forced initial updates, plus the opposite
source sign in normal recovery. It checks mixed bivariate quadratic products,
including u*u_tt and u*u_zz, and independently extracts the complete second
amplitude equation with all factors two and the declared energy sign. These
finite identities support the general sixteen-field arguments; they are not
a numerical PDE simulation or independent proof of existence.

The seven corruption controls include five source/inventory/pin failures and
two internally rehashed false mathematical inputs. Authentication passes for
the mathematical mutants before the formal arithmetic rejects them:

- Omitting initial forcing must fail at
  `Arithmetic failed: formal:forced_initial_recurrence:0:0`.
- Omitting a mixed-product binomial term must fail at
  `Arithmetic failed: formal:mixed_Leibniz:value_value:0:2:0`.

Huge exponentials are kept as symbolic expressions with exact rational and
integer comparisons. Sobolev multiplication, geometric/natural-action
identities, energy arguments, and retained analytic premises remain written
proofs. The test count does not establish external novelty, physical origin,
a finite nonlinear amplitude interval or an actual-area observable remainder.
