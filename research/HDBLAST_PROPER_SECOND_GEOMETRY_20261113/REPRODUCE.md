# Reproduce the authenticated finite checks

Use Python 3.12.14 and the exact SymPy 1.14.0 / mpmath 1.3.0 wheels pinned
with SHA256 in `requirements-checks.txt`. Install in a separate environment:

```
python -m venv /tmp/hdblast-proper-check-env
/tmp/hdblast-proper-check-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/requirements-checks.txt
```

Run the authenticated protocol from the repository root. Replace `PIN`
with the independent literal manifest SHA256 in the checked workflow
`.github/workflows/proper-second-geometry-20261113.yml`; do not obtain the
trust pin by hashing an unauthenticated replacement capsule.

```
/tmp/hdblast-proper-check-env/bin/python -I -B research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/checks/check_checkpoint.py --manifest-sha256 PIN --output /tmp/hdblast-proper-normal.json
/tmp/hdblast-proper-check-env/bin/python -O -I -B research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/checks/check_checkpoint.py --manifest-sha256 PIN --output /tmp/hdblast-proper-optimized.json
cmp /tmp/hdblast-proper-normal.json /tmp/hdblast-proper-optimized.json
/tmp/hdblast-proper-check-env/bin/python -I -B research/HDBLAST_PROPER_SECOND_GEOMETRY_20261113/checks/run_controls.py --manifest-sha256 PIN --output-dir /tmp/hdblast-proper-controls-fresh
```

Outputs and the controls directory must be fresh and outside the closed
capsule. The runner rejects symlinks, a changed source, unlisted/missing
files, a stale manifest pin, duplicate JSON keys and nonfinite JSON values
before importing its pinned evaluator. Dependency wheels are authenticated
by the separate hash-required install; runtime versions are checked too.

The 119 named checks comprise seventeen source snapshots times three
size/SHA256/Git-blob checks, nine formal algebra comparisons covering 119
exact components, and 59 exact rational accounting checks. The component
count is distinct from the number of named comparisons. Direct oracles
differentiate scalar/path/area polynomials and a curved four-dimensional
Lorentzian metric, while production formulas use covariant connection and
Ricci variations. A nonconstant scalar background and an uncorrected
coordinate velocity provide explicit premise-failure witnesses. The
separate constant-conformal family is an algebra control, not the actual
registered initial family.

Thirteen fresh child-process rejection cases cover changed/extra/missing
or symlinked source, stale manifest bytes, seven changed mathematical
factors and an opposite connection-factor mutation. Mathematical mutants
have valid refreshed manifests and must fail a recomputed exact oracle,
so they are not merely authentication failures.

The arithmetic ledger never exponentiates the enormous B1/B2 ceilings.
It verifies their formal exponent comparisons and the rational trace,
observer, curvature and area coefficient accounting. These checks do not
prove the continuous Sobolev, PDE, geometric or existence arguments. Read
the three complete proof notes and their stated analytic premises.
