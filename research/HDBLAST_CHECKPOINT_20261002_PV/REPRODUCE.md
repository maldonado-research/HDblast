# Reproduce

This benchmark needs no private archive, cloud-drive files or local computer access. The mass profile is specified analytically and every input parameter is in REGISTRATION.md plus the pre-execution PROTOCOL_ADDENDUM.md.

Use Python 3.12 with NumPy 2.2.6, SciPy 1.15.3, mpmath 1.3.0 and Matplotlib 3.10.1. Limit BLAS threads to one for the many small mode-vector operations.

```sh
python -m pip install numpy==2.2.6 scipy==1.15.3 mpmath==1.3.0 matplotlib==3.10.1
python code/verify_protocol.py
python code/verify_matching.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python code/pv_crossing.py --registration REGISTRATION.md --registration-sha256 fc72bf086bb7747af0a51848ef997bf11d9d4451ee894aef899ea0ea8ca2888d --source-ref local-replay
python analysis/summarize_results.py
python analysis/plot_results.py
```

The fixed matrix runs eight dynamic cases and one algebraic constant-mass null. The first runtime and control results are reported in AUDIT_AND_EXECUTION.md. Outputs are written into outputs/pv_actual/. Full per-case NPZ curves and spectra are regenerated; they use allow_pickle=False. The curated first-run JSON and primary curves are kept separately in outputs/ for evidence preservation.

The runner exports scalar energy, pressure, mass source, work, trace correction, direct physical-mode comparison, potential-only pressure control and normalization reconstruction bounds. Comparing errors scaled by m_infinity^4 or m_infinity² is not the same as a relative error on a small quantum source; raw magnitudes are included.

Download HDBLAST_CHECKPOINT_20261002_PV.zip from the checkpoint directory. From its extracted directory, verify the downloaded file using:

```sh
python code/verify_package.py ../HDBLAST_CHECKPOINT_20261002_PV.zip
```

The verifier checks unique safe paths, ZIP CRCs, inventory and every SHA-256 digest. The ZIP contains all necessary code, results, selected curves, plots and documents; no private correspondence or unrelated project files are included. GitHub Actions artifacts provide regenerated arrays but have finite retention.

A successful finite regulator/cutoff comparison is numerical evidence, not a mathematical proof of the limit or validation of an HDBLAST cosmology.
