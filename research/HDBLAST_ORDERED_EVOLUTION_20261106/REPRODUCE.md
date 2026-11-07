# Replay the ordered-background and forcing checkpoint

Python 3.12 or later and only the standard library are required. Obtain the
external manifest pin from `.github/workflows/ordered-evolution-20261106.yml`
at the same reviewed commit and set PIN to its literal SHA256. A pin computed
from untrusted input files is not an authenticity anchor. From the repository
root, use fresh output locations outside the closed capsule:

```sh
python -I -S -B research/HDBLAST_ORDERED_EVOLUTION_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/ordered-evolution-normal.json
python -O -I -S -B research/HDBLAST_ORDERED_EVOLUTION_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/ordered-evolution-optimized.json
cmp /tmp/ordered-evolution-normal.json /tmp/ordered-evolution-optimized.json
python -I -S -B research/HDBLAST_ORDERED_EVOLUTION_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/ordered-evolution-controls-fresh
```

The canonical reports must be byte-identical. There are 173 named checks:
24 source size/SHA256/Git-blob checks, 50 formal polynomial checks, 44 ordered
background bound checks and 55 forcing-map arithmetic/allocation checks.
Eight immutable public snapshots are authenticated before any evaluator code
is loaded. The checker rejects changed, missing, extra or symlink files,
stale pins, duplicate JSON keys and nonfinite constants. No private producer
report or saved numerical result is an input.

The formal engine uses 85 field/jet indeterminates, two fibers and three base
directions with noncommuting connection matrices, nonzero connection
curvature and state-dependent coefficient maps. It compares the literal
coordinate equation with the fully substituted bundle equation, checks
background maps, Hessian ordering and both signs of the complete second
variation. These are exact finite symbolic identities supporting the general
sixteen-field derivation, not a numerical simulation or an independent
existence theorem. The enormous exponential bound is kept as a symbolic
expression; all comparisons are rational/integer, without floating-point
overflow or a claim to have evaluated that exponential numerically.

The seven corruption controls include five authentication failures and two
false mathematical inputs. The latter are internally rehashed with a new
valid manifest pin, so they reach actual formal arithmetic:

- Omitting the coordinate linear substitution must fail at
  `Arithmetic failed: symbolic:coordinate_to_bundle_row:0`.
- Replacing the principal-variation factor two by one must fail at
  `Arithmetic failed: symbolic:second_variation_row:0`.

The Sobolev embedding, tensor connection formulas, derivative-allocation
argument and conditional energy theorem remain written mathematical proofs
and retained analytic premises. A passing report does not supply an evolved
N14 bound, the complete source-time/corner hierarchy, an actual-family
amplitude interval or a physical observable remainder.
