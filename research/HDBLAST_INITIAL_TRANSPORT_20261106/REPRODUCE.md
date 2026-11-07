# Replay the initial transport checkpoint

Python 3.12 or later suffices; only the standard library is used. Run from the
repository root, obtaining the literal external manifest SHA256 from
`.github/workflows/initial-transport-20261106.yml` at the same reviewed commit.
Set `PIN` to that value. Recomputing a pin from untrusted inputs does not supply
an authenticity anchor. Select fresh output locations outside the capsule.

```sh
python -I -S -B research/HDBLAST_INITIAL_TRANSPORT_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/initial-transport-normal.json
python -O -I -S -B research/HDBLAST_INITIAL_TRANSPORT_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/initial-transport-optimized.json
cmp /tmp/initial-transport-normal.json /tmp/initial-transport-optimized.json
python -I -S -B research/HDBLAST_INITIAL_TRANSPORT_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/initial-transport-controls-fresh
```

The canonical reports must be byte-identical. There are 317 named checks:
27 source size/SHA256/Git-blob checks, 34 formal polynomial checks and 256
finite bound checks. The checker authenticates the closed file inventory
before loading capsule code or the retained source engine. Missing, extra,
changed or symlink files, stale manifest pins, duplicate JSON keys and
nonfinite JSON constants are rejected. Reports must be fresh and outside
the capsule. The manifest lists every capsule file except itself.

The global first-packet majorants are recomputed from the original pinned
rational engine, without importing a saved producer report. The formal
polynomial check composes a nonlinear root with a non-diagonal three-state
quadratic map, with indeterminate Cauchy data and nonzero explicit time
coefficients. It verifies both the general time-dependent chain rule and the
initial even-frame specialization. This is a symbolic check of finite jets;
it does not numerically integrate or certify the actual continuous flow.
The general sixteen-field chain identities and tensor estimates are proved
in INITIAL_TRANSPORT.md.

The seven rejection controls include five authentication/inventory failures
(changed source, extra member, missing member, symlink, stale pin) and two
false mathematical inputs. For the latter, the control rehashes the modified
capsule and supplies its new external pin, so authentication passes:

- Root second derivative `-1/4` changed to `0` must fail at
  `Arithmetic failed: symbolic:root_second_derivative`.
- Second-velocity cross factor `2` changed to `1` must fail at
  `Arithmetic failed: symbolic:root_second_velocity`.

These checks reach actual formal arithmetic, not a metadata-only comparison.
The unchanged old evolution-input capsule can be replayed separately with
its own original workflow pin; this checkpoint does not replace that report.

Finite checks support the conditional proof. Retained analytic existence,
weighted construction, flat extension and all-point background bounds remain
source premises. No evolved norm, admitted amplitude, observer domain or
physical-origin conclusion follows from a passing report.
