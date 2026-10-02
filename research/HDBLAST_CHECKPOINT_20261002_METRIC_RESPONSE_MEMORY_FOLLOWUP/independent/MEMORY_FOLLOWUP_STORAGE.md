# Archive-only memory followup

This is a new prospective experiment. It does not replace the failed original
or the failed first followup. The first followup used public freeze
`5ff571cdb7729398834ff3c025b50191c4ee29f2`; its independent producer exceeded
the unchanged256MiB resource gate at306088KiB after153.80774094seconds. Both
completed positive-B raw archives and that failed result remain preserved.

Read-only ZIP/NPY-header inspection found631members in the fine archive,
168327096 uncompressed NPY bytes, and65892627 compressed bytes. The old
producer retained all six snapshots until final serialization, and copied
their references into failure-state bookkeeping. Its SHA function then read
the entire compressed archive while `evolve` still held those snapshots.
The internal resource check preceded that return-dictionary hash; the next
main-level check recorded the excess peak.

`stream_npz.py` writes each observation array synchronously as the exact
standard `.npy` member used by `numpy.savez_compressed`, inside a
ZIP_DEFLATED/ZIP64 NPZ. It retains member names and ZIP metadata, without
retaining ndarray references. Existing numeric/complex precision, shape,
Fortran-order metadata, member names and ordering are preserved. Completed
members remain readable after a context-managed interruption. Source/metric
jets, modes, every bare/subtraction/contact piece, baselines, and the complete
direct scalar histories retain the same631-member schema per complete run.

The producer deletes only the observation `contact_archive` and `integrands`
dictionaries after synchronous writing. Failure bookkeeping retains current
modes and direct scalar histories; completed observations already reside in
the partial raw NPZ rather than duplicate `archived_*` aliases. Each archive
closes on exception. SHA256 uses1MiB chunks rather than reading the full file.

`ARCHIVE_ONLY_INVARIANCE.json` proves18 producer functions byte-identical,
all model/grid/gate/budget declarations identical, and all remaining numerical
and retained-field expressions in `evolve` identical after explicitly
removing only the storage/lifetime operations. The separately audited WKB
and stable-baseline modules remain byte-identical. `PRIOR_FAILED_PRODUCER.py.txt`
retains the exact prior source for this read-only comparison.

`STREAM_FORMAT_CHECKS.json` records27 synthetic controls, including extended
precision increments of2^-60 near1 for real/complex longdouble, exact NPY-byte
agreement, normal NumPy loading, array-reference release, bounded hashing,
and completed-member preservation on synthetic interruption.
`FULL_SCHEMA_STREAM_CHECKS.json` records all631 fine-sized synthetic members
checked for exact NPY bytes and dtype/shape/value loading; its archive-only
peak was34804KiB. These are synthetic storage tests, not physical evaluations
or a certified bound on the next physical producer's peak.

Physical formulas, evolution, source/grid/cutoffs, gates,900-second wall and
256MiB memory budgets remain unchanged. In particular the repair does not
change Ward, refinement or cross-route acceptance. A successful memory check
would not imply those scientific gates pass. No new physical source values,
modes, stresses, quadratures or response evaluations occurred during this
preparation. Physical execution requires the new public registration first.
