# Provenance of `basis/`

These files were produced on 13–14 October 2026 by a companion Claude Code
session on this Mac (workspace `~/HDblast-research-work/second-order-20261013`)
and copied here unchanged; hashes are in `../MANIFEST.json`.

Corrections to `SECOND_ORDER_LOCALIZED_CONSTRAINTS.md` (left unchanged as a
historical record):

* §0 item 3 says the exact verification of `P_H[w] = S_H` was done. At copy
  time `verify_cc.py` had not finished; it later terminated without output.
  The identity is instead established by `new/fast_verify.py` (step 1 and the
  momentum identity, exact), `new/exact_random_point.py` (exact rational
  identity tests) and `independent/indep_check.py` (independent numerics).
* §0 item 4 ("must deform the monopole metric") holds only for the packet in
  the reconstruction gauge; see Theorem 3 and Corollary 3.1 of the technical
  note, as corrected after the referee report.
