# Incoming-pulse energy test

For this static linear single-channel model, a specified asymptotic pure
continuum packet has exact unit reflection and its local brane energy vanishes
at late times. Both position and velocity bound-mode projections must initially
vanish; merely setting q and qdot to zero does not imply that condition. This
result distinguishes transient energy transfer from permanent capture or
reheating. It does not exclude nonlinear, interacting or cosmological HDBLAST
models.

The preserved theorem, contract and controls are under `model/`; independent
proofs and controls are under `review/`. The portable manifest maps all exact
accepted bytes. Original manifests retain historical absolute provenance paths.
The producer checked 25 exact identities and eight negative cases; the reviewer
checked 35 independent exact cases. Both modes passed all 68 controls. No actual
blast state, observational data, retained array or numerical source was used.

From this directory, with the pinned research Python, SymPy and mpmath:

```sh
python -I -B model/verify_incident_controls.py --output /tmp/incident-normal.json
python -I -B -O model/verify_incident_controls.py --output /tmp/incident-optimized.json
```

Metric calibration FAIL, full continuous certificate UNRESOLVED, physical origin
NOT_ESTABLISHED and novelty NOT_ASSESSED remain unchanged.
