# A certified auxiliary response after shell contact

This reproducible checkpoint bounds one conditional auxiliary response mean at
`(a,b)=(3/2,1/4)`, including the shell history and the complete feedback into its
five closed mean components. The outward interval is

`55,843,882.796267 <= mean(W_yy) <= 737,593,144,520.766258`.

The result concerns an auxiliary forced increment. A full evaluated numerical
waveform, physical/invariant interpretation, nonlinear integrability and a
cosmological origin are not established. The retained total `U` already has a
conditional analytic construction and conditional physical completion. See
[the manuscript guide](manuscript/README.md),
the [exact proof](analytic-review/MEAN_BLOCK_POINT_CERTIFICATE_PROOF_FINAL.md)
and [independent internal numerical review](numeric-review/NUMERICAL_REVIEW.md).
Internal AI-assisted review is distinct from external specialist review.

## Reproduction

Use Python 3.11 or newer; the numerical replay uses only the standard library.
No installation, network request or third-party Python dependency is needed.
Keep the extracted capsule files read-only (`0444`) and directories read-only
(`0555`). Obtain the expected manifest SHA-256 from the separately published
`SOURCE_FREEZE.json` and run:

```text
python3 -I -S -B reproduce.py --manifest-sha256 EXPECTED_MANIFEST_SHA256 --output-dir /tmp/hdblast-boundary-fresh
```

Replace `EXPECTED_MANIFEST_SHA256` with that 64-character hash and choose a new
absolute output directory outside the capsule. Allow several minutes and
roughly 30 MB of output/scratch space. The resulting
`REPRODUCTION_RECEIPT.json` records the actual child processes, raw and typed
golden comparisons, targeted controls and preservation checks. The PDF is
included with editable TeX; PDF typesetting is not part of this replay.

Before any retained numerical import, the runner authenticates the full new
source inventory, the retained 49-file forced-response archive and its nested
169-file geometric archive, including every manifest member. It extracts the
two retained capsules into read-only scratch trees. Their older scientific
campaigns are not run again. The precise retained premise and source-certificate
scope remain conditional.

The unmodified `compute/certify_boundary_response_v2.py` is executed through
`portable_actual.py`. That adapter changes only its three input-root globals;
its source hash and numerical algorithm are preserved. The fresh normal and
optimized source reports are joined, without changing their bytes, in one
read-only scratch folder for the five-mean-block calculation. The independent
reviews use the freshly generated reports.

`compute/BOUNDARY_N16_NORMAL.json` and its optimized counterpart intentionally
retain the earlier coarse full-state allowance, which contains zero. The
separate five-mean-block certificate proves the positive interval using the
same source integral and a sharper complete coupled feedback bound. These are
different bounds on the same auxiliary quantity, not different source runs.

The seven targeted numerical corruptions, eleven synthetic analytic controls
and six independent formula distinctions are reported separately. They are
in-process checks and do not establish rejection of every possible incorrect
certificate or replace external mathematical assessment.

Private editorial history, credentials and publication operations are absent
from this capsule. No new journal submission or Zenodo version is made here.
