# Independent exact binary80 codec review

Review date: 2026-10-08. Scope: the declared stdlib-only exact binary80 NPY/NPZ
storage codec and adversarial fabricated fixtures in this directory. This is
an independent subagent review of the candidate implementation, not root
approval, study registration, publication, or authorization to decode retained
scientific payloads.

## Reviewed artifacts

| Artifact | SHA256 |
| --- | --- |
| `exact_binary80.py` | `6556d52228ce3aa2f8604d29e3dcc9d3a81ad947dc4f126300128162ce15bd87` |
| `review_codec_independent.py` | `18dd8725ea8b4101b15ddb8f4e40dac0b82343696c99a35ecc00f4ad8edb87b6` |
| `test_exact_binary80.py` | `625b5e1dfaee5c873605ba7428f80ac16fb55d5358784de0fc3e80e20fd71b07` |

The verdict applies to these bytes. A changed implementation requires renewed
review of the affected behavior.

## Findings and resolution

1. **Resolved strict-spec defect.** NPY version tuples containing Booleans or
   floats compared equal to `(1, 0)`. The implementation now requires exactly
   two exact integer components. Independent regressions reject all fabricated
   Boolean/float variants.
2. **Resolved blocking complete-member verification defect.** `ZipExtFile`
   truncated decoded member output to the declared uncompressed size. A crafted
   stored or deflated member with an extra output tail passed when local and
   central sizes/CRC described only the prefix. The codec now reads the exact
   compressed region itself, bounds actual output, and verifies exact output
   length/hash/CRC. Independent regressions reject both compression methods
   with zero scalar decode calls.
3. **Resolved related raw-deflate defects.** Extra compressed bytes and a
   truncated stream whose output/CRC happened to be complete were accepted by
   the former `ZipExtFile` path. The raw reader now requires the deflate
   end-of-stream marker and rejects unused or unconsumed trailing compressed
   data. Both independently fabricated cases are rejected before scalar decode.
4. **Resolved central-directory allocation concern.** Metadata parsing formerly
   preceded the explicit roster cap. The bounded manual preflight now checks
   the exact central record roster/region before creating `ZipFile`. A crafted
   trailer declaring one member with 64 central records is rejected while a
   patched `ZipFile` constructor has zero calls.

## Independent validation

Commands run from this directory:

```text
python -m unittest -v test_exact_binary80.py review_codec_independent.py
python -O -m unittest -v test_exact_binary80.py review_codec_independent.py
```

Both final-candidate executions completed successfully: **34 test methods
passed** (15 author-suite methods and 19 independent-review methods), with
multiple boundary and malformed-input subcases. The normal run reported
0.064 seconds and the optimized run 0.047 seconds; these timings are observations,
not performance guarantees. Fixture construction intentionally emits Python's
duplicate-name warning in the duplicate-member rejection case.

The checks cover exact minimum/maximum subnormals, minimum normal, maximum
finite, one significand ULP, negative values, both signed zeros, arbitrary
padding, malformed/nonfinite binary80 encodings, complex component order,
checkpoints and exact integer ranges, strict bounded literal headers, version
pins, file/header/member hashes, all-member failures before any decode,
stored/deflated/local-ZIP64 containers, CRC/roster/trailer/resource failures,
opaque metadata decode flags, and isolation from pathname changes after
verification. Canonical boundary expectations are explicit exact Fractions;
no host float conversion occurs.

Static re-review confirms that `verify_inputs` verifies an opaque whole-file
snapshot, checks/copies every complete archive member, and validates every
header before returning the object through which selected numeric iterators
are exposed. Raw inflation output is bounded on every call. All binary80
numeric decoding remains exact integer/Fraction arithmetic under the explicit
x87 little-endian padded-16 layout. Canonicality is checked lazily during the
consumer pass, as documented.

## Verdict and scope limits

**No outstanding blocking findings in the reviewed candidate bytes for the
declared storage contract.** The discovered archive counterexamples are retained
as independent regressions, and pass in normal and optimized Python.

No historical/retained NPY or NPZ payload was opened or decoded. No producer or
scientific source callback was imported or invoked. No Git or network operation
was performed. The independent harness creates disposable fixtures in this
decoder directory; the author harness uses disposable system temporary files.
Both harnesses clean up their fixtures.

This review does not authenticate a supplied registration or layout-evidence
digest, prove the producer ABI, certify scientific arithmetic/error enclosures,
or approve a live comparison. The caller must authenticate the registered pins
and evidence independently. Root review approval and the required public
registration/freeze remain separate gates before any retained payload decoding.
