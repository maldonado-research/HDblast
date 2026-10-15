# Authenticated fresh reproduction

Use Python **3.12.14** and SymPy **1.14.0**. The scientific arithmetic is
integer, exact rational or symbolic; plotting additionally uses Matplotlib
3.10.8 and is presentation only. No old scientific program is executed.

Download the source ZIP and the separate `SOURCE_FREEZE_RECEIPT.json` from the
dated delivery. Verify the ZIP hash against that receipt before extraction.
Read its external `manifest.sha256`, check the extracted `MANIFEST.json` hash
with a trusted SHA256 tool, and then run:

```sh
python -I -S -B HDBLAST_QUADRATIC_SHELL_FORCING_20261015/verify_manifest.py --expected EXTERNAL_MANIFEST_SHA256
python -m pip install -r HDBLAST_QUADRATIC_SHELL_FORCING_20261015/requirements.txt
python -I -B HDBLAST_QUADRATIC_SHELL_FORCING_20261015/reproduce.py --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir fresh-quadratic-replay
```

Replace the uppercase placeholder with the receipt's exact hash. Keep the
output directory outside the source. It must be absent or empty. The runner
authenticates the entire frozen inventory before importing or executing any
scientific program. It checks the separate verifier, exact Python/dependency
versions and the full named replay roster. Each process receives a distinct
copy containing only its declared dependencies. Outputs begin absent.

Both normal and optimized Python modes run every declared positive calculation
and corruption control. All JSON scientific fields are compared with the frozen
goldens; only individually type-checked execution-mode metadata may differ
between modes. The producer regenerates the complete report and CSV byte for
byte, and independent checkers reconstruct every continuous coefficient cell.
The source roster names each group, embedded control, external mutation,
expected diagnostic and mechanism. A nonzero exit from an unrelated exception
does not qualify as a successful control. The unchanged source is authenticated
again after replay. The receipt retains process logs and full output hashes.

This reproduces finite calculations and checks conditional enclosures. It does
not prove analytic existence/regularity hypotheses or calculate the complete
second-order coordinate source, full `U2`, total stress or nonlinear remainder.
