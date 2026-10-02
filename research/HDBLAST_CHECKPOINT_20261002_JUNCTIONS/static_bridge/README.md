# Static de Sitter feasibility bridge

This result derives conditional analytic response and reanalyzes the archived classical junction Jacobian. It does not compute quantum expectation values, select physical parameters, run new ODE/BVP experiments, or solve a coupled quantum shell.

Files:

- `STATIC_DESITTER_SOURCE_BRIDGE.md`: source conventions, regular-cone response, exact local static limits, an explicit candidate Euclidean state, and the proposed next registration.
- `verify_static_bridge.py`: 33 algebra/archive assertions and five deliberate wrong-formula controls.
- `SOURCE_PINS.json`: expected hashes for 12 public input files at `e17a01b428bb8049e919c42376ab0359e41d613c`.
- `STATIC_BRIDGE_CHECKS.json`: executed check results and matrix inversion of the archived Jacobian, including the numerical-floor warning for its tiny current-to-radius component.
- `PORTABILITY_CHECKS.json`: an extracted source tree without Git metadata passes; an altered public matching input is rejected by its expected source hash.

Use Python with NumPy and SymPy already supplied by the project workflow. No new environment configuration is required for this bounded verification.

```sh
python verify_static_bridge.py --repo /path/to/HDblast --output /path/to/replay-output
```

The verifier reads a checkout or extracted repository source tree and requires all 12 inputs to match the frozen hashes. It detects the source tree from current/script ancestors when `--repo` is omitted. It records the scientific source commit separately from any observed checkout HEAD; it works without Git metadata after source extraction. All outputs are confined to `--output`, which defaults to the script directory.

The next physical calculation needs a prospective registration fixing mass/coupling/gravity inputs, finite matching, Euclidean-state renormalization, independent source derivatives, regulator limits, and response/error gates. The note is a concrete study proposal, not that registration.
