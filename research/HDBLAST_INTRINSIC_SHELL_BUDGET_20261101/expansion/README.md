# Portable shell expansion geometry component

The main proof is `DERIVATION.md`. It applies established Gaussian
foliation, averaged Raychaudhuri and normalized scalar-curvature variation
machinery to the selected exact family's first shell field. It derives

```
(<R3>−6/Rvolume²)_2 = 10(u−H∫₀^τu)²/(3A²)
```

and the exact forced first shell trace equation used for scalar-source
recovery. It also isolates the second scalar mean still needed for a full
volume response. It supplies no new nonlinear existence theorem or actual
four-dimensional fluid interpretation.

The essential component files are:

- `DERIVATION.md`: clock, curvature defect, averaged balances and trace ODE.
- `GAUSS_IDENTITY_STANDARD_APPLICATION.md`: self-contained standard
  five-to-four-dimensional Gauss application; no private note is needed.
- `check_shell_expansion.py`: the unchanged 21 scientific identities and
  five wrong-formula controls, with portable source/output handling. Each
  control must end with its exact mutation-specific ValueError witness,
  exit nonzero, emit no success stdout and leave its fresh output absent.
  The canonical reports record these checks as structured control rows.
- `SOURCE_PINS.json`: repository-relative source inputs plus local proof
  hashes. The script authenticates these before importing SymPy.
- `REFERENCES.md`: claim-to-source mapping and scope of standard machinery.
- `TRACE_RECOVERY_AGGREGATION_AUDIT.md`: independent mathematical check of
  the numerical Volterra extraction, separate from the 21 identities.

From this component directory, with SymPy 1.14.0 available:

```
python -I -B check_shell_expansion.py --repo-root CHECKOUT --run-controls --output NORMAL.json
python -I -B -O check_shell_expansion.py --repo-root CHECKOUT --run-controls --output-dir optimized
```

The second command writes `optimized/SHELL_EXPANSION.json`. These reports
must agree byte for byte. Omitting both output options writes the same
canonical JSON to stdout. Source authentication uses file hashes rather
than requiring the checkout's HEAD to equal the recorded source commit;
the latter records provenance and permits an identical source tree in a
later commit. A containing capsule can authenticate its full closed source
set before invoking this component.

No absolute private workspace path is a dependency. Generated reports are
finite algebra certificates, not PDE solves or specialist proof review.
