# Causal bulk response: conditional scalar testbed

This checkpoint derives a stable flat 4+1-dimensional half-space scalar model's
retarded boundary response, positive continuum spectral density, dissipative
energy flux and filtered Gaussian noise. It also constructs an exact continuum
of four-dimensional fields with the same reduced response and measurement laws.
These observables therefore do not identify spacetime dimensionality or the
origin of the Big Bang. This model is separate from the de Sitter pulse numerics.

The declared domain is M>0, m>0, c>0 and 0<g²<m²(c+M). The model's coupling g
is unrelated to the endpoint calculation's time-dependent forcing g(eta).
Ground or thermal equilibrium is a state premise; thermalization and a blast
initial condition are not derived. Raw equal-time continuum force noise has an
ultraviolet divergence, so measurement smearing is required. Stable bound modes
and their occupation must be retained. The finite closed quadratic rival theorem
concerns an ideal exact function; no finite noisy discrimination bound is given.

Accepted files under `model/` and `review/` are copied byte for byte. Their
original manifests retain absolute historical paths as provenance, rather than
pretending those paths are portable. `PORTABLE_CONTENT_MANIFEST.json` records
the corresponding paths in this repository. The inherited numerical checkpoint
and all frozen historical files remain unchanged.

From this checkpoint, with Python 3.12.14, SymPy 1.14.0 and mpmath 1.3.0:

```sh
python -I -B model/verify_half_space_controls.py --output /tmp/half-space-normal.json
python -I -B -O model/verify_half_space_controls.py --output /tmp/half-space-optimized.json
```

The preserved independent script writes its receipt beside its source. To keep
the checkpoint unchanged, copy it into a fresh temporary directory and run it
there:

```sh
bulk_review_tmp=$(mktemp -d /tmp/hdblast-bulk-review.XXXXXX)
cp review/independent_controls.py "$bulk_review_tmp/independent_controls.py"
python -I -B "$bulk_review_tmp/independent_controls.py"
```

The producer checks 26 exact identities and four negative cases; five numerical
integrals are explicitly uncertified diagnostics. Independent review checks 14
exact identities, two negative cases and six uncertified diagnostics. Its
acceptance covers only the declared conditional analytic model.

Metric calibration remains **FAIL**; the full continuous numerical certificate
remains **UNRESOLVED**; higher-dimensional Big Bang origin remains
**NOT_ESTABLISHED**; external novelty remains **NOT_ASSESSED**. Gravity,
cosmological expansion, backreaction, reheating, a specified blast energy ledger,
and an observational likelihood still need a coupled physical model.
