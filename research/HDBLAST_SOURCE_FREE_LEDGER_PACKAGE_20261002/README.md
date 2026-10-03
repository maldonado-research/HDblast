# Complete source-free ledger package

The two parts reconstruct a **48,470,467-byte** standalone ZIP containing all 172 payload files plus its manifest: frozen inputs, code, registration, original calculated results, reviews, and figures. They are transport chunks, not independent scientific packages.

```bash
python reassemble.py --parts-manifest PACKAGE.json --expected-sha256 b98291b7d8a2a40d5ddf79d877bb0177c03b748becc726d63a94456adc4830a0 --output /a/fresh/ledger.zip
```

Extract the verified ZIP into a fresh directory and follow its `REPRODUCTION.md`. Python 3.12.14 and the pinned dependencies are required for the numerical replay. Both complete routes include all twelve cases at 80 and 100 decimal digits; each route has one 900-second/256-MiB budget.

The measured classification is **LEDGER_ERROR_DEMONSTRATED** on the fixed source-free interval. The earlier metric calibration remains **FAIL**. The fresh ZIP replay receipt and hosted publication evidence are external to the ZIP so they do not alter the payload they verify. No Zenodo DOI is implied by this GitHub package.
