# Reproduce the closed subsidiary-closure capsule

Use Python 3.12.14 with the exact binary-wheel hashes in
`requirements-checks.txt` (SymPy 1.14.0 and mpmath 1.3.0).

```
python3.12 -m venv /tmp/hdblast-subsidiary-env
/tmp/hdblast-subsidiary-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/requirements-checks.txt
```

Obtain the full manifest SHA256 from
`.github/workflows/subsidiary-closure-20261116.yml` at the trusted tested
commit. Do not replace that external expected digest by a digest computed
from an untrusted local capsule. Assign it to PIN, then use fresh outputs:

```
/tmp/hdblast-subsidiary-env/bin/python -I -B research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-subsidiary-normal.json
/tmp/hdblast-subsidiary-env/bin/python -O -I -B research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-subsidiary-optimized.json
cmp /tmp/hdblast-subsidiary-normal.json /tmp/hdblast-subsidiary-optimized.json
/tmp/hdblast-subsidiary-env/bin/python -I -B research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-subsidiary-auth-fresh
/tmp/hdblast-subsidiary-env/bin/python -O -I -B research/HDBLAST_SUBSIDIARY_CLOSURE_20261116/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-subsidiary-math-fresh
```

The main checker authenticates its external manifest, exact closed file
roster, every member's size/hash, safe paths, regular files and strict JSON
before importing the pinned local evaluator. All 41 source snapshots are
checked against size, SHA256 and Git blob identity, with 92 qualified
ancestor-manifest member joins. Distinct manifest namespaces are retained;
an identical basename is not treated as a unique ancestor member.
The P1--P9 proof-source receipt and retained parameter ledger are joined
to those exact inputs. No antecedent source producer is imported/executed.

Five protocol faults test changed source bytes, an unlisted member, a
missing member, a symlink and a stale external manifest digest. They
must reject before evaluator import. The separate 48 mathematical controls
change an oracle claim in a fresh complete capsule, rebuild its manifest
and pass its new external digest to a fresh full checker. Authentication
must succeed and rejection must reach the intended named mathematical
or boundary-contract check. A syntax/parser failure is not substituted
for the target. The harness preserves optimized Python in its children.

The local campaign runs main reports and the mathematical fault suite in
both modes and compares outcomes. Normal/optimized main reports are
byte-identical; no assertion controls validation. CI repeats both main
reports, protocol faults and optimized fresh mathematical controls with
the same pinned wheels and literal manifest digest.

The 32 replay checks include schema/engine guards and finite algebra
checks supporting their explicitly limited identities.
Auxiliary constant-metric/frozen-symbol/polynomial-conjugator/small-chart
witnesses do not prove full variable-coefficient Einstein evolution,
the exponential derivative formula or all actual 16-field coefficients.
The four written proofs and their internal review supply the conditional
continuum join. Source authentication establishes bytes and provenance;
it does not validate the retained background or other analytic premises.

Outputs must be fresh and outside this capsule. A completed evidence
index is immutable; start a new private evidence root for later work.
No physical energy calibration, observational fit, parameter derivative,
new amplitude interval or empirical origin is computed by this replay.
