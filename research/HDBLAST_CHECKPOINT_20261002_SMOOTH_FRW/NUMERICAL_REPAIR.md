# Preserved first attempt and mechanical repair

The first registered execution used source commit 983c472dfe70d2a830b8e9b538995732dc38626f and the unchanged REGISTRATION.md hash. The A=0 primary physical modes completed and were saved before subtraction. NumPy then treated ndarray-left addition with the Taylor-jet object as an object ndarray. Calling sqrt on that ndarray raised AttributeError. No renormalized source, cutoff gate or matrix verdict was produced by this attempt.

The original failure log, failure metadata, provenance and incomplete summary are in outputs/first_attempt. The completed numeric modes remain in the cloud execution artifacts, with size and SHA256 recorded. The original source is retained in code/history and the public registration commit.

The repair adds only Jet.__array_priority__=1000 and its explanatory comment, so NumPy dispatches array-left operations to the reflected jet operators. An independent reviewer checked the exact diff. Physical evolution, counterterms, background, state, action convention, cutoffs, quadrature, time grids, tolerances and acceptance gates are unchanged. REPAIRS.json records both source hashes. The protocol has not been rewritten.

The added no-mode preflight checks array dispatch, shapes, finite counterterms, complete subtraction exchange, exact-static subtractions and the explicit flat formulas. Its executed PASS report is retained. It does not stand in for the physical matrix. A fresh full matrix is executed only after this repair is committed; its source hash and outputs are separately recorded.
