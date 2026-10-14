# Verified final publication

The [final live readback](FINAL_LIVE_READBACK.json) verifies **21 targets**:
the public English/Spanish pages, shared homepage and HDBLAST page, both
earlier frozen research capsules, the exact-family proof packet including
its refreshed manifest/ZIP, completed replay evidence, and fresh referee
reports. Every target returned HTTP 200 and matched the reviewed local bytes.
It follows the documentation clarifications in merged research PR50.

[Publication receipt](PUBLICATION_RECEIPT.json) · [Successful CI and deployments](CI_AND_DEPLOYMENT_RECEIPT.json) · [Fresh referee reviews](FRESH_REFEREE_REVIEWS.md) · [Theorem and limits](../../../../research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/README.md).

To repeat the readback using current checkouts of `HDblast` and the shared
`maldonado-research.github.io` repository, run the standard-library verifier:

```sh
python hdblast/continuation/2026.10.14-exact-localized-precontact/publication/verify_published.py \
  --shared-root /path/to/maldonado-research.github.io \
  --output /tmp/hdblast-live-readback.json
```

The default HDBLAST root is the repository containing this script; use
`--hdblast-root` to override it. Both checkouts must contain the reviewed
files. A future source or website change can legitimately make a comparison
fail; the saved receipt records the observed publication at its timestamp.

The result is conditional exact localized initial data and sufficiently
short precontact evolution at high finite regularity. Through-contact
evolution remains open. These are internal AI-assisted reviews, with no
external peer review or novelty certification. No new Zenodo version or
APS submission was created.
