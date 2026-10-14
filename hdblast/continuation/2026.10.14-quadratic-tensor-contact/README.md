# Quadratic tensor checkpoint: publication evidence and next target

The [research packet](../../../research/HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014/README.md)
contains the complete fixed-longitudinal quadratic tensor source and the
specified even-master forced tensor coefficient through controlled shell
contact. The coupled scalar/vector/gauge contact system remains open.

- [Frozen source and ZIP pins](FROZEN_SOURCE_RECEIPT.json)
- [Publication receipt: merged PRs, successful deployments and preserved work](PUBLICATION_RECEIPT.json)
- [All three successful checks on the reviewed research PR head](PRE_MERGE_CHECKS.json)
- [Fresh normal/optimized replay: 12 successful runs and four actual corruption rejections](fresh-replay/REPRODUCTION_RECEIPT.json)
- [Live readback: every frozen research file and four website pages](FINAL_LIVE_READBACK.json)
- [Next calculation: forced even l=4 scalar equations and shell boundary matrix](NEXT_COUPLED_CONTACT_TARGET.md)

The readback verifier uses only Python's standard library. With the reviewed
public and shared-site checkouts available, run:

```sh
python verify_tensor_published.py --hdblast-root /path/to/HDblast --shared-root /path/to/maldonado-research.github.io --output /tmp/tensor-live-readback.json
```

It compares 49 current public targets with those local bytes. The saved
readback is a dated witness; later intentional website edits may change its
hashes. Source reproduction is described in the research packet separately.
Historical process outputs retain original workspace paths, with distributed
evidence bytes listed in PUBLICATION_EVIDENCE_MANIFEST.json.

These checks authenticate published files and replay exact algebra. They do
not certify the analytic proof from check counts, establish a cosmological
origin or constitute external peer review. No new Zenodo version or APS
submission was made for this GitHub checkpoint.
