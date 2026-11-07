# Quantitative initial-data transport through the retained boundary chart

6 November 2026. Conditional mathematical checkpoint, continued from PR #133.

The previously bounded raw harmonic Cauchy data are now transported through
both the proper-y radial component square root and the actual inverse boundary
flow. In the fixed product bundle and H14/H13 norm defined in
[INITIAL_TRANSPORT.md](INITIAL_TRANSPORT.md), the selected first and second
amplitude derivatives at zero satisfy

```
T1 = ||V1||H14 + ||V1,t||H13 < 10^5200,
T2 = ||V2||H14 + ||V2,t||H13 < 10^10000.
```

The deliberately large rounded bounds cover one copy and the sum-of-integrals
norm on two identical reflected copies. They concern the specified gauge and
weighted canonical curve, conditional on the retained background and initial
construction. They do not bound a finite-amplitude evolution.

The proof charges the physical-to-product connection and volume change,
constant shell normalization, symmetric tensor norm, root Hessian, inverse-flow
jets and the factor-two term in the second time velocity. It uses the known
smooth first packet as a multiplier; it assumes no new high supremum norm for
arbitrary evolved solutions. Flat support lies inside the established physical
collar, so no outside cutoff derivative is used.

[REPRODUCE.md](REPRODUCE.md) describes source authentication, exact arithmetic,
normal/optimized replays and corruption controls. Nine immutable public inputs
are pinned to main commit `66ce4c07cd9b16cb67f5e9f147402d16781c676e` in
[SOURCE_PINS.json](SOURCE_PINS.json). The proofs, source snapshots and checker
form a new capsule; previous checkpoints are unchanged.

The finite checks support the written conditional argument. They are not an
independent proof of its analytic premises, external peer review or evidence
for a higher-dimensional origin of the Big Bang. These constants have no SI
energy interpretation. [CONTINUATION.md](CONTINUATION.md) states the next gaps.
