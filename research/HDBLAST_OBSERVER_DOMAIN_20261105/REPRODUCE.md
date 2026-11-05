# Reproduce the finite checks

Use Python 3.12 and its standard library. Obtain the trusted outer manifest
SHA-256 from the pinned `observer-domain-20261105.yml` workflow at the commit
being checked, and substitute it for `TRUSTED_MANIFEST_SHA256` below. Choose
fresh absolute output directories outside this complete checkpoint:

```
python -I -S -B research/HDBLAST_OBSERVER_DOMAIN_20261105/checks/run_checkpoint.py --manifest-sha256 TRUSTED_MANIFEST_SHA256 --output-dir /tmp/hdblast-observer-domain-replay-new
python -I -S -B research/HDBLAST_OBSERVER_DOMAIN_20261105/checks/run_controls.py --manifest-sha256 TRUSTED_MANIFEST_SHA256 --output-dir /tmp/hdblast-observer-domain-controls-new
```

The outer replay authenticates its entire closed inventory and runs three
complementary components in normal and optimized Python modes. All six actual
report files must match their pinned references byte-for-byte, and stdout must
equal the corresponding report. Finite guard counts per mode are embedding
**50**, normal-clock **38**, and coefficient-bootstrap **383**. These are finite
algebra checks accompanying continuous proofs, not machine verification of the
global topological or analytic theorems.

The embedding component has its own nine-file closed capsule and four public
source copies. Its two canonical reports are deliberately OUTSIDE that inner
capsule in `references/`. Adding reports or extra files under `embedding/`
would break its exact inventory. The outer manifest includes the inner
manifest as an ordinary authenticated member; it excludes only itself.

The separate controls run two authentication-only positives and six actual
early source/CLI rejections: wrong manifest pin in each mode, changed public
input, changed reference, changed producer, and occupied output. These test
source transport, not the scientific hypotheses or mathematical proof mutants.
Rejected controls run no scientific producer.

No prior campaign, nonlinear PDE solver, field oracle or physical numerical
estimate is executed. Independent internal review is documented separately
from these runtime checks. The common coordinate-field neighborhood,
model-specific constants and a numerical admissible amplitude remain inputs.
