# Conditional nonzero physical shell response

This checkpoint proves qualitative nonvanishing of a physical scalar response
for one explicitly specified aligned quadrupole in the registered HDBLAST
Einstein–scalar model. Both displacement and velocity are nonzero. The proof
uses the complete coupled reconstruction and physical shell laws, followed by
the previously constructed conditional exact nonlinear family.

The [main theorem](PHYSICAL_SHELL_RESPONSE.md) specifies its background,
localized exact-data, support, gauge and finite-time hypotheses. A smooth
annular bump is registered in the conjugated master variable and converted
back to actual constraint-compatible ambient data. A near-front Duhamel
argument proves that the physical first scalar trace cannot remain zero in
any interval extending beyond the linear arrival infimum. It does not assume
that every nonzero bulk perturbation produces this observable.

For a proper-time clock launched normally from a marked initial shell section,
define the squared scalar difference on each physical spatial section and
integrate it over a fixed finite clock window. The exact response has

\[
\mathcal R_\varepsilon(\rho_bT_1)
=\varepsilon^2\mathcal C_{T_1}+o(\varepsilon^2),\qquad
\mathcal C_{T_1}=\frac{\pi^2\rho_b^4}{3}
\int_0^{T_1}\cosh^3t\,S(t)^2dt>0.
\]

Consequently the response is positive for sufficiently small nonzero
amplitudes. Its leading coefficient depends on the selected first tangent,
not on the second exact initial correction. No numerical waveform,
coefficient value or effective amplitude threshold is provided. The exact
nonlinear first-arrival time remains open.

- [Main result, full hypotheses and proper-time observable](PHYSICAL_SHELL_RESPONSE.md)
- [Independent aligned positive-velocity derivation](sol/ALIGNED_POSITIVE_VELOCITY_SHELL_RESPONSE.md)
- [Independent pure-displacement comparison, a different tangent](astra/PHYSICAL_SCALAR_ACTIVATION.md)
- [Critical analytic acceptance](referee/FINAL_PHYSICAL_RESPONSE_ACCEPTANCE.md)
- [Astra joint review](astra/PRIMARY_AND_SOL_ACCEPTANCE.md) and [Sol cross-review](sol/INDEPENDENT_CROSS_REVIEW.md)
- [Retained analytic inputs and preservation](RETAINED_INPUT_PROVENANCE.json)
- [Accessible archive/source map](archive/SHELL_RESPONSE_SOURCE_MAP.md)
- [Bounded primary-literature comparison](literature/SHELL_RESPONSE_PRIMARY_LITERATURE.md)

## Reproduce the finite checks

Use Python 3.12 and the pinned dependency in `requirements.txt`. The external
[frozen-source receipt](../../hdblast/continuation/2026.10.15-physical-shell-response/FROZEN_SOURCE_RECEIPT.json)
records the manifest, verifier, runner and ZIP hashes. Compare those hashes
before executing downloaded Python. The scoped GitHub workflow performs this
authentication before installing dependencies or running scientific code.

After extracting the matching ZIP and obtaining the manifest hash from that
receipt, run from this directory:

```sh
python -I -B verify_manifest.py --expected MANIFEST_SHA256_FROM_FROZEN_RECEIPT
python -m pip install -r requirements.txt
python -I -B reproduce.py --manifest-sha256 MANIFEST_SHA256_FROM_FROZEN_RECEIPT --output-dir /tmp/hdblast-physical-response-replay
```

Choose a new output directory for each execution. The runner verifies every
source byte and runs each positive/control process from a fresh checker copy,
normally and with Python optimization. Its expected totals are six positive
processes, thirty-two rejected wrong-formula processes and forty-one finite
positive check groups per mode. There are no embedded control processes in
this checkpoint. The root mode flag is the only positive JSON field normalized
for comparison, and it is separately checked against the actual process mode.

The controls target omitted embedding displacement, gauge signs, the Hessian
factor two, angular normalization, proper-time measure, second-jet contamination,
double counting the one common shell, bulk conjugation, shell-normal signs,
normal lapse terms and the ambient/global velocity conversion. Each control
must fail for its specified reason. Retained prior scripts and receipts are
unaltered analytic inputs and are not counted as newly executed checks.

The symbolic identities and corruption controls check specified algebra and
normalizations. They do not computer-certify the analytic existence,
activation or nonlinear remainder proof. The published literature ledger and
dated provenance receipt describe inspected sources; raw primary-paper bodies
remain internal and are not part of this executable replay.

This is conditional mathematics with internal AI-assisted derivations and
review. External peer review, numerical magnitudes, global stability,
observational agreement and a higher-dimensional origin of the Big Bang are
not established. This checkpoint makes no new Zenodo version or APS submission.
