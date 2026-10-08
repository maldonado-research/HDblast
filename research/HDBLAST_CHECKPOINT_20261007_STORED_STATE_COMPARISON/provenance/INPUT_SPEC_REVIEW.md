# Candidate incoming-state input specification

`INPUT_SPEC.json` covers all four inherited capsules and all twelve fixed
source/grid/cutoff cases at eta = -9/2. It remains a candidate specification:
there has been no GO for decoding the retained state in this audit.

The current opaque audit passes the published manifest SHA256, both identical
producer text files, all four whole original archive hashes, all four whole
compact capsule hashes, and all nineteen complete member streams per capsule
against both the capsule and its original archive. Header-prefix hashes,
member lengths and CRC32 also pass. `OPAQUE_AUTHENTICATION_RECEIPT.json`
records this evidence. No producer or saved-array reader was imported, no
retained value was decoded, no source callback executed, and no Git command
ran. The used prior inputs retain their before/after byte identities.

Each compact capsule contains nineteen members. Exactly five members are
selected for a future incoming-state decode: `k.npy`,
`momentum_weights.npy`, `observation_eta.npy`, `u_1.npy`, and `w_1.npy`.
The other fourteen remain metadata-only. The full compact inventory and all
its member/hash/header pins live in each `reader_inputspec`; the original
631-member inventories live in `ORIGINAL_ZIP_INVENTORIES.json`. There are
49,152 nodes across the four separate capsules and 86,016 node appearances
across the twelve nested prefixes.

| Source | Grid | K | Prefix count | Full count |
| --- | --- | ---: | ---: | ---: |
| positive_B | coarse | 64 | 2048 | 8192 |
| positive_B | coarse | 128 | 4096 | 8192 |
| positive_B | coarse | 256 | 8192 | 8192 |
| positive_B | fine | 64 | 4096 | 16384 |
| positive_B | fine | 128 | 8192 | 16384 |
| positive_B | fine | 256 | 16384 | 16384 |
| signed_uB | coarse | 64 | 2048 | 8192 |
| signed_uB | coarse | 128 | 4096 | 8192 |
| signed_uB | coarse | 256 | 8192 | 8192 |
| signed_uB | fine | 64 | 4096 | 16384 |
| signed_uB | fine | 128 | 8192 | 16384 |
| signed_uB | fine | 256 | 16384 | 16384 |

The producer uses the positional prefix `int(LD(K)/width)*16`, rather than a
literal `k<K` mask. This follows from authenticated producer text lines
340–346. Its momentum rule at lines 201–208 maps promoted `leggauss(16)`
nodes within full panels of width 1/2 or 1/4. That source recipe does not
establish actual retained node values or source-to-source equality. The
authoritative values remain the stored represented values.

After GO, the reader must verify full finite positive increasing `k`, full
finite positive aligned weights, and exact equivalence of each positional
prefix to the strict band `k<K`, including rejection of cutoff equality.
These are integrity requirements derived here, rather than assertions
present in the old producer. In particular, an exact equality
`sum(momentum_weights)==K` is not an inherited invariant and must not be
invented. Prefix sums and any other totals are computed deterministically
from the exact decoded ratios only after GO.

`momentum_weights` is the dk quadrature weight. The physical node weight is
`mu_j=momentum_weights_j*k_j^2/(2*Pi_represented^2)`, as established by producer
lines 299–304 and 340–344 and active-source protocol lines 94–95. Epsilon is
exactly `3777893186295716171/37778931862957161709568`; Pi is exactly
`14488038916154245685/4611686018427387904`. Source-text constants,
mathematical pi, exact decimal epsilon, or fresh native evaluation cannot
replace those represented ratios. Stored `u_1,w_1` already include epsilon;
the normalized state is `U=u_1/epsilon`, `W=w_1/epsilon`, without projection.

Producer lines 27–28 and 441–453 establish the six-entry observation map
`(-11/2,-9/2,-4,-7/2,-5/2,-3/2)`. Thus suffix `_1` identifies eta=-9/2,
`_2` identifies -4, and `_3` identifies -7/2. The incoming anchor is history
index 192 coarse or 384 fine; the prehistory source-cell count is 64 coarse
or 128 fine. These remain source-derived metadata until values are read.

The `<f16` and `<c32` descriptors identify storage sizes and byte order; they
alone do not distinguish binary80 from another sixteen-byte format.
`LAYOUT_PROVENANCE_BASIS.json` binds the public original-results runtime:
Linux x86_64, Python 3.12.14, NumPy 2.2.6, 63 fraction bits, and epsilon
2^-63. The active-source public protocol additionally specifies x87 binary80
ABI at line 334. The candidate reader layout is little-endian x87 binary80
in sixteen-byte slots, with meaningful bytes 0–9 and six padding bytes.
Complex storage consists of the real sixteen-byte slot followed by the
imaginary sixteen-byte slot. Padding remains part of the byte hash even
when ignored numerically. Fabricated native serialization and independent
reader review must validate this interpretation before a real decode.

Safe access requires the whole capsule hash before ZIP interpretation, exact
ordered membership and ZIP metadata, bounded streams, then all full member
and header hashes before any selected payload decode. There is no filesystem
extraction, pickle, eval, producer import, or native ndarray load. No future
cell coordinate is assumed to be a binary64 value; no stored source-cell
coordinate has been decoded here.

The full twelve-case pressure/contact certificate remains **UNRESOLVED**;
metric calibration remains **FAIL_UNCHANGED**; the higher-dimensional
Big Bang cause remains **NOT_ESTABLISHED**. This audit authenticates inputs
and defines a conditional reader boundary; it does not certify state error
or upgrade scientific conclusions.
