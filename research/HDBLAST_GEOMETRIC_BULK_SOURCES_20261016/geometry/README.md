# Geometric Ricci and reduced Einstein source

`geometric_sources.py` completes the parent source inventory with the full
Ricci Hessian and natural/reduced Einstein source, using ordinary first-field
amplitude two-jets in the registered frame. `indexed_geometry.py` is the
reusable covariant-jet contraction core. Read
`GEOMETRIC_SOURCE_DERIVATION.md` for conventions, every index/sign conversion,
the harmonic assumption, and unsupported claims.

Import this folder on `sys.path`, then call
`evaluate_sources(m,s,b)`. The sibling `projection/angular_generic.py` provides
the default exact polynomial backend; an API-compatible backend can be
passed explicitly. Output tensor indices are 0=y, 1=time, 2..5=ambient
tangent sphere. Each component is an angular polynomial dictionary. The
sibling projector exports and reconstructs the 23 source coefficients.

From this folder, with new output paths:

```sh
python -I -S -B check_geometric_sources.py --run-controls --output /tmp/geometry-normal.json
python -I -S -O -B check_geometric_sources.py --run-controls --output /tmp/geometry-optimized.json
```

The two outputs must be byte identical. The checker and its nested mutation
subprocesses explicitly support isolated standard-library execution. Six direct mutations are available
through `--mutation`; see `REPLAY_ROSTER.json` for their required diagnostics.
The checker never uses Python assertions and rejects an existing output path.
The four retained parent inputs are authenticated by `INPUT_PINS.json`.
They fix conventions and provide the inherited matter/scalar/gauge adapter;
they are not silently substituted for the new Ricci calculation.

This is a mathematical source evaluator and exact coefficient result. It
does not claim an evaluated second waveform or exact nonlinear family.
