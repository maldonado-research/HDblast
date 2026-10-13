# Cone-current checkpoint: verification

This GitHub-only continuation supplies the reviewed technical note and
[immutable manifest inventory](../../../research/HDBLAST_CONE_PHASE_SPACE_20261013/MANIFEST.json).
It does not update Zenodo 23335485, which remains the preceding spectral
edition, or submit anything to APS. The concept DOI family is unchanged.

The relocated replay passed 79 positive algebra/scaling checks and 23
wrong-formula controls in each of normal and optimized Python. Receipts
matched byte for byte; deliberately corrupted inputs, an altered pinned
manifest and injected checker failure were rejected. The source and
relocated capsules remained unchanged. These counts do not certify global
existence, phase-space completeness, novelty or a cosmological origin.

The public CI workflow authenticates the frozen manifest and replay source
before installing symbolic dependencies. Its manifest SHA256 is
`6021915456815936a80c864d5e503b05725325a2cb15cbd882051fade0890194`.

The capsule contains 64 files / 1,095,538 bytes including its manifest.
Its technical PDF has ten pages; internal mathematical, bibliography and
render-layout reviews are under `reviews/` in the capsule. Separate agent
review is not external peer review. The note remains subject to author review.

The supplied figure is optional for algebra replay. To regenerate it in a
scratch copy with NumPy and Matplotlib installed, run
`python figures/make_cone_figure.py` from that copied capsule; regenerated PDF
metadata can differ. Do not edit the frozen capsule in place.

CI run and complete remote-source/site readback evidence is retained with
the corresponding private research handoff after deployment.
