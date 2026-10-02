# Reproduce the registered scientific failure and complete saved data

Use normal Python3.12 and the exact five packages in requirements-replay.txt. All inputs are contained in this checkpoint; no parent checkout or credentials are needed. Keep source/package files unchanged and use fresh external output paths.

For a published split package, run code/reassemble_package.py with --parts-manifest PACKAGE_PARTS.json, --output a fresh external ZIP path and --expected-sha256 the published full ZIP SHA. It authenticates every part, full ZIP, payload hash, CRC and freeze binding. Extract the verified ZIP into another fresh external directory.

From its checkpoint directory, execute:

```python
import json, subprocess, sys
from pathlib import Path
root = Path.cwd()
receipt = json.loads((root / "FREEZE_RECEIPT.json").read_text())
subprocess.run([sys.executable, str(root / "code/reproduce_metric_expected_scientific_failure.py"),
    "--checkpoint", str(root), "--output", "/tmp/hdblast-metric-unused-control",
    "--public-freeze-commit", receipt["public_freeze_commit"],
    "--registration-sha256", receipt["registration_sha256"],
    "--independent-manifest-sha256", receipt["independent_manifest_sha256"]], check=True)
```

The wrapper invokes the unchanged32-command replay with all frozen gates. It requires completed primary12rows, all four mode evolutions and complete631-member raw archives within900seconds/256MiB, then the anticipated internal Ward/scientific failure. A successful control prints PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL while the underlying experiment remains FAIL. Unsupported failures are fatal and preserved. No cross-route validator PASS or nine accepted-calibration figures are promised.

`python code/replay_metric.py --output /tmp/hdblast-metric-unused-plan --plan-only` authenticates prospective and final package membership and lists32commands without physical evaluation. Direct full replay deliberately exits nonzero on failed science. Saved-only audits may later reconstruct operators and compare routes while retaining FAIL.

FULL_REGISTRATION protects prospective bytes; the post-freeze FREEZE_RECEIPT binds the exact public commit. Final package/replay receipts are kept outside the checkpoint to avoid self-reference. Prior original and high-precision failures remain separate immutable evidence. This special-point conformal channel does not complete actual-root/lapse/bulk/state response or establish stability, heating or cosmological support.
