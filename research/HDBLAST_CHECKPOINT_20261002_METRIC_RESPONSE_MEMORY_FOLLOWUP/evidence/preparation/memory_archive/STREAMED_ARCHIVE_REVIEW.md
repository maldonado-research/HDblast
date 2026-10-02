# Scientific invariance review of the proposed archive-only followup

The proposed independent producer changes archive allocation, serialization
and hashing only. This review found no scientific source, grid, arithmetic,
gate or resource-budget change. It imported or executed neither producer.
The old producer SHA256 is
`a9f6afdb4c1d484299e5bc09b99597560f266a4769133a09c8a388c2e7d8a147`;
the proposed producer is
`4fa625f2b534d741e2a48d6f0bfe84a89ee17bf306805906f962a8552e005906`.

Eighteen function bodies are byte-identical, including the source/forcing,
transport, direct observables, cutoff sums, Simpson history, combination,
resource checks and main execution. Top-level constants, grids, gates and
entrypoint are exactly the same AST. The provenance gate only requires the
new `stream_npz.py` helper in the frozen independent manifest. The SHA helper
only replaces a whole-file byte read with bounded one-MiB chunks.

The evolution loop and all retained raw-field expressions are exactly the
same AST after removing a narrowly enumerated set of archive operations:
the old retained dictionary construction, duplicate `archived_*` references,
deferred bulk save and output-path placement; the new synchronous initial
member writes, output-path placement, explicit close, and release of already
written observation dictionaries. One archive context wrapper ensures closure
on normal return or exceptions. Numerical calculations, gates, observations,
mode copies, contacts and histories are retained. The independent review
passes both normally and with optimization, and rejects six fabricated
mutations to amplitude, transport sign, gate comparison, raw history name,
saved complex mode and closure.

The helper writes each actual ndarray synchronously with NumPy's standard
NPY formatter into ZIP_DEFLATED/ZIP64, preserving its dtype and shape. It
retains member names and ZIP metadata but no ndarray references; object
arrays, duplicate names and path names are rejected. Completed members remain
readable if an exception interrupts a later calculation. Extended real and
complex precision witnesses contain a 2^-60 increment near one and survive
without binary64 conversion. The author's synthetic full-schema test checks
all 631 NPY members' exact bytes, dtype, shape and values using fabricated
arrays with the prior fine archive's read-only header schema. Its measured
I/O-only peak is 34804 KiB; no saved physical array was loaded or recomputed.

These checks establish archive format and scientific invariance, not a future
physical memory budget or scientific PASS. The earlier recorded resource
failure remains FAIL. In particular, the earlier saved fine positive_B Ward
declaration exceeds the unchanged endpoint gate; synchronous serialization
cannot improve the numerical response. A separately registered completion
attempt must retain every gate and may transparently qualify its expected
scientific failure after all four runs are saved within budget. Such a
control's PASS would mean reproducible scientific FAIL with complete data,
not a validated metric calibration.

The review receipts are `STREAMED_INVARIANCE_NORMAL.json` and
`STREAMED_INVARIANCE_OPTIMIZED.json`; the source-only verifier is
`verify_streamed_archive_invariance.py`. Actual shifted-root, lapse/gauge,
bulk, state and coupled-evolution prerequisites remain unchanged.
