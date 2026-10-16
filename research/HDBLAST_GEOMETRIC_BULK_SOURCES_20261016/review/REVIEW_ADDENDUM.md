# Isolated replay portability addendum

This addendum preserves REVIEW_RECEIPT.json and its original scientific replay
reports. It accepts only the subsequent import-path fixes and records current
portable replay bindings in FINAL_REVIEW_RECEIPT.json.

Running a script with Python -I removes the script directory from the implicit
module search path. The geometry producer's checker therefore now explicitly
adds its own directory before importing geometric_sources, and its nested
control subprocesses use the same -I -S -B flags. The indexed curvature core,
frame implementation, complete coefficient maps and numerical source formulas
are unchanged. The producer's updated freeze binds its isolated replay.

Two independent review entry points also imported independent_coordinate_oracle
from their own directory. Each now explicitly inserts that directory into
sys.path. Removing exactly that one added line reproduces its previously
accepted file hash. All other files bound by the original review receipt are
byte-identical to their accepted versions. Import-only smoke checks with
-I -S -B establish that the current independent entry points load correctly.

PORTABLE_FRAME_EXPECTED.json and PORTABLE_CORE_EXPECTED.json are **expected
replay outputs rebound only for the checker hash**, not new mathematical run
receipts. Their scientific payloads are byte-for-byte the same JSON values as
the earlier successful actual runs. check_portable_review.py checks that exact
relationship. The portable capsule's fresh normal/optimized execution compares
the current checker output with these complete expected outputs. The other
scientific replay goldens are unchanged.

The independent full-domain numeric review separately binds the dyadic interval
implementation. This addendum does not extend the geometric review to a newly
implemented arithmetic ring. The existing independent exact Fraction endpoint
evaluation of the retained interval source maps remains valid.

REPLAY_ROSTER.json gives the portable contracts. PUBLIC_FILES.json is the
review-directory allowlist, while the original receipt's accepted_files map is
a cross-directory scientific acceptance manifest. The archived seal_review.py
is retained for provenance; its old-file consistency check is intentionally
superseded by check_portable_review.py after these explicit path-only changes.

Acceptance remains internal independent mathematical review with no blocking
finding in the stated scope. This addendum makes no new mathematical claim.
