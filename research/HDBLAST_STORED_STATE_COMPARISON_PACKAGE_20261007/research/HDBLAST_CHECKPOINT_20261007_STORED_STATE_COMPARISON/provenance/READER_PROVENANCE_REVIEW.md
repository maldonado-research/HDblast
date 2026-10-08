# Independent storage/provenance review

The candidate reader at SHA256
`6556d52228ce3aa2f8604d29e3dcc9d3a81ad947dc4f126300128162ce15bd87`
supports the explicit inherited ABI and the four candidate InputSpec records.
This conclusion is conditional on authenticating the registered InputSpec and
layout-evidence bytes in the real caller. It authorizes no retained decode.

The binary80 conversion uses all 64 significand bits and the explicit integer
bit. Normal values are signed M*2^(E-16446); subnormals are signed
M/2^16445. Signed zero is retained separately from the rational zero. NaN,
infinity, pseudo-special, pseudo-denormal and unnormal encodings fail closed.
Complex elements are real sixteen-byte storage followed by imaginary
sixteen-byte storage. Padding is ignored numerically and remains in all
complete-file/member hashes.

The whole source file is copied into an owned bounded snapshot and its exact
size/SHA256 verified before ZIP parsing. The parser preflights the complete
ordered central-directory roster, rejects duplicates and unsupported names
and modes, and validates contiguous local member regions. It reads complete
raw stored/deflated regions rather than relying on ZipExtFile's declared-size
truncation. Actual member size, full SHA256 and CRC are verified; deflate must
end exactly at its region boundary with no hidden output or trailing
compressed bytes. All nineteen members pass these checks before any header
parse or selected value iterator is exposed.

The NPY parser uses a bounded ASCII AST and accepts only the three distinct
literal dictionary keys, an exact supported descriptor, literal False order,
and bounded tuple of exact nonnegative integer dimensions. It never evaluates
the AST and accepts no expressions, object/structured layout, pickle or
arbitrary-code path. Header hashes and dimensions/lengths are explicitly
pinned. The other fourteen members remain authenticated opaque metadata; the
candidate spec selects exactly the five incoming-state members.

`INPUT_SPEC_HYDRATION_REVIEW.json` confirms metadata-only construction and
validation of all four reader dataclasses with explicit capsule-compatible
limits. The hydration helper opens only the candidate metadata JSON, and
never calls verify_inputs or iter_values on a retained file.

`NATIVE_BINARY80_FIXTURE_REVIEW.json` independently validates the historical
ABI premise using newly fabricated native NumPy NPY serializations on the
current Linux x86_64, little-endian, sixteen-byte long-double, nmant=63,
iexp=15 platform. Ten real values cover signed zero, a binary80 ULP, finite
normal extremes, and subnormal boundaries. Four complex values confirm
component ordering and signed zero. Every decoded component agrees with the
native scalar's exact integer ratio. These fixtures contain no retained
input and are tied to the reader hash above. The actual historical value
encodings remain to be checked after GO.

The public producer/runtime/protocol evidence establishes the explicit x87
premise; `<f16` alone would not. Independent fabricated reader tests and
regressions owned by the decoder reviewer are separate supporting evidence.
The whole-file capsule digest also authenticates its expected timestamps,
attributes and extra-field bytes; the candidate full inventory preserves
these exact metadata, even where the generic decoder checks container
semantics rather than exposing individual timestamp pins.

No retained payload has been decoded in this review. No physical producer
callback executed, no producer module was imported, and no prior checkpoint
was edited. Actual input authorization, finite/canonical full-array checks,
positive/increasing full momentum checks, exact prefix-band equivalence and
all numerical comparisons remain behind the reviewed public freeze and GO.
The full twelve-case certificate remains **UNRESOLVED** and the proposed
higher-dimensional cause remains **NOT_ESTABLISHED**.
