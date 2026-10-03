# Complete active-source ledger checkpoint package

The two payload parts reconstruct the unchanged repaired checkpoint ZIP: 32,195,260 bytes, SHA256 `0cd3dd2da604444541c06e98167f537f769c0ce430108bf02354e5090f328c4f`. Each part is at most 24 MiB. `PACKAGE.json` records their ordered lengths and hashes, the embedded payload manifest, repaired public freeze commit and registration hash.

Use Python 3.12 or newer and fresh output/receipt filenames in an existing directory:

```bash
python reassemble.py \
  --parts-manifest PACKAGE.json \
  --expected-sha256 0cd3dd2da604444541c06e98167f537f769c0ce430108bf02354e5090f328c4f \
  --output /tmp/HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER.zip \
  --receipt /tmp/HDBLAST_ACTIVE_SOURCE_REASSEMBLY_RECEIPT.json
```

The standalone helper verifies every part, complete archive membership, regular-file modes, the embedded manifest, and every payload member's SHA256 and CRC before atomically creating the output. It performs no extraction, physical array decoding or scientific calculation. `transport_evidence/` records its manufactured normal/optimized guards and byte-exact reassembly of all 238 canonical archive members.

The ZIP preserves the registered inputs, methods, outputs, initial execution failure, and separately disclosed repair chronology. Fresh replay evidence placed beside these distribution files is outside the immutable ZIP. Transport verification does not assign a scientific classification or establish the higher-dimensional-blast hypothesis.
