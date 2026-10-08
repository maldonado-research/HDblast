# Continuation input provenance, 2026-10-08

The candidate `INPUT_SPEC.json` selects nine numeric streams from each of the
four existing 19-member compact capsules: `k.npy`, `momentum_weights.npy`,
`observation_eta.npy`, and `u_1.npy`/`w_1.npy`, `u_2.npy`/`w_2.npy`,
`u_3.npy`/`w_3.npy`. The complete archive roster, order, byte lengths, file
hashes, member hashes, header hashes, CRCs, descriptors, shapes, and storage
layout pins are preserved from the frozen stored-state specification. The
other ten members remain opaque and are excluded from numeric decoding.

The candidate enables a later reader selection. It supplies no retained-decode
GO, public registration, or scientific error certificate. In this work no
retained NPY member stream or header was opened, no retained scalar was decoded,
and no physical source or producer function was imported or called.

| Source | Grid | Nodes | Capsule bytes | Capsule SHA256 |
| --- | --- | ---: | ---: | --- |
| positive_B | coarse | 8192 | 2339164 | `adb5728a77da58b99fb4009d1cb88c5298ff9c4a3bb59d5ad0cf9f4a2c81eb4e` |
| positive_B | fine | 16384 | 4671836 | `7f36f2be2b9c345b556668aaf1574bfa07ad01c394d5b1acd1c4b793fa0c2eaa` |
| signed_uB | coarse | 8192 | 2339164 | `bb5b7ade887db7a9b86c29cbd5c88ebad8c3b0217a9db716ecb57ca05265db84` |
| signed_uB | fine | 16384 | 4671836 | `16f271b8cb668adb151f48cf44788947339a1532e702ff41d7f0ad8456056d59` |

All four compact and all four original archive hashes and sizes were checked
by opaque streaming, and their ZIP central metadata matched the frozen
inventories exactly. The compact selected-member pins also matched the
active-source `INPUT_MANIFEST.json`. This pass did not recompute complete NPY
member/header hashes or reauthenticate Git publication objects. Those pins and
the claim of published-blob identity are inherited from the frozen metadata.
`STATIC_PROVENANCE.json` distinguishes these scopes and records zero member
streams opened. The linked artifacts in the candidate's inherited
`layout_evidence.path` and `original_zip_inventories.path` remain relative to
the upstream stored-state checkpoint; their complete upstream repository
paths and hashes are recorded in `STATIC_PROVENANCE.json`.

The producer observation map at
`research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py:27`
is `(-11/2, -9/2, -4, -7/2, -5/2, -3/2)`. Its mode-write branch at lines
429–443 uses observation indices, so the continuation selections mean:

| Index | Eta | Selected state | Role |
| ---: | --- | --- | --- |
| 1 | -9/2 | u_1, w_1 | incoming anchor |
| 2 | -4 | u_2, w_2 | interior snapshot |
| 3 | -7/2 | u_3, w_3 | terminal snapshot |

This is source-derived metadata; the stored `observation_eta` values have not
been inspected in this work. A future retained reader must check the exact
six-entry map and all selected state vector alignments after its authorized
decode boundary.

The full inherited per-step mode trajectory is unavailable in these archives.
Each of the four original 631-member inventories contains only `u_0` through
`u_5` and `w_0` through `w_5`, at six observation times. The compact capsules
contain only indices 1, 2, and 3. Producer lines 384–395 overwrite the current
`u,w` at each step; lines 403 and 429–442 write those modes only when the step
is an observation. There is no array of intermediate per-node `u,w` states in
these inventories.

The `history_*` arrays carry aggregate cutoff observables, the ledger,
baseline contact, source and forcing jets, and geometry. Specifically,
`history_values` has shape `(577,3,9)` or `(1153,3,9)`: nine integrated
quantities for three cutoffs. This data does not retain the modes for 8192 or
16384 individual momentum nodes. It cannot supply inherited mode comparisons
at each continuation step. No inference about other unrelated archives or
unrecorded transient memory is made.

| Grid | Step | Full history count | Anchor index | Interior index | Terminal index | Continuation cells |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| coarse | 1/128 | 577 | 192 | 256 | 320 | 128 |
| fine | 1/256 | 1153 | 384 | 512 | 640 | 256 |

These counts follow from producer lines 27–31 and 356–376. The 129/257
continuation grid points include both endpoints; only the three stated
observation modes are retained. Source-cell coordinates are not asserted
from decoded data, and native regeneration would not authenticate an input.

The inherited twelve strict-band cases retain their positional prefix rule
`int(LD(K)/width)*16` from producer lines 340–346. Coarse prefix counts are
2048, 4096, 8192; fine counts are 4096, 8192, 16384 for K = 64, 128, 256.
There are 49,152 separate capsule nodes and 86,016 appearances across those
nested prefixes. Separate capsules remain separate even if source recipes
agree. A future reader must check positive increasing k, positive aligned
weights, and exact equivalence of every positional prefix to `k<K`, rejecting
cutoff equality. No exact equality of weight sums or cross-capsule values is
assumed.

`momentum_weights` is the inherited dk quadrature weight. The physical weight
is `mu_j = momentum_weights_j*k_j^2/(2*Pi_represented^2)`, from producer lines
299–304 and 340–344. Every selected mode retains epsilon already; normalization
is `U_i=u_i/epsilon_represented`, `W_i=w_i/epsilon_represented` for i = 1, 2, 3,
without projection. The pinned represented constants remain
`epsilon=3777893186295716171/37778931862957161709568` and
`Pi=14488038916154245685/4611686018427387904`.

The copied `exact_binary80.py` is byte identical to the frozen upstream file,
SHA256 `6556d52228ce3aa2f8604d29e3dcc9d3a81ad947dc4f126300128162ce15bd87`.
`UPSTREAM_HASH_RECEIPT.json` pins it and the exact copied fixture/test files.
Its explicitly pinned layout is little-endian x87 binary80 in 16-byte slots,
with ten meaningful bytes and six numerically ignored but byte-hashed padding
bytes. Complex slots store real then imaginary. Finite canonical values are
returned as exact `Fraction` values; signed zero is preserved separately.
All complete members and headers are verified before numeric iterators are
exposed; each selected slot is checked for canonicality when iterated. This
storage contract does not itself validate scientific arithmetic.

Validation used only manufactured bytes. From this directory run:

```sh
python -B input_spec_adapter.py
python -B -m unittest -v test_exact_binary80.py review_codec_independent.py test_continuation_selection.py
python -B -O -m unittest -v test_exact_binary80.py review_codec_independent.py test_continuation_selection.py
```

The metadata command passed all four hydrations. Each test command passed
36 test methods: the 34 frozen codec tests plus two continuation-selection
integration checks. Those checks exercise a manufactured 19-member archive
with all nine selected streams, preserved signed-zero complex components,
arbitrary padding, and unselected history streams that fail closed on decode
requests. They also hydrate the actual candidate JSON with snapshot and scalar
decode functions patched to fail if called. The duplicate-ZIP-name warning is
an intentional malformed-fixture control. `MANUFACTURED_VALIDATION.json`
records commands, outcomes, and artifact hashes.

Metric calibration remains **FAIL_UNCHANGED**, the full twelve-case
pressure/contact certificate remains **UNRESOLVED**, and a higher-dimensional
Big Bang cause remains **NOT_ESTABLISHED**.
