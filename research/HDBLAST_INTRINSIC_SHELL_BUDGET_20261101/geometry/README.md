# Portable independent shell geometry component

This component checks the first-only geometric identities used by the shell
trace numerical result. It carries independent proofs, exact symbolic
checkers and source pins with repository-relative paths. It has no private
workspace dependency and runs no numerical oracle.

The primary new application is the intrinsic curvature defect

    C=<R3>−6/Rvolume², C''=10E²/(3A²),

and the Buchert coefficient

    Q_B''=2E'^2+(32/3)E'J'+(16/3)J'^2,
    E=u−HP, P'=u, J'=P/A².

These follow from standard round-metric criticality, spatial diffeomorphism
invariance and the hypersurface averaging identities. The independent
proof is `CURVATURE_BUCHERT_AUDIT.md`.

The scalar-gradient Q'' formula is inherited prior work; `DERIVATION.md`
records that distinction, exact-family/clock joins, source-domain limits,
necessary background-anchor qualifications and qualitative finite-family
angular fraction consequences. `ANGULAR_FRACTION.md` proves the angular
measure formula and exact rational tail bound. No source supplies a
numerical finite-parameter radius, shell energy, reheating or new generic
mathematical theorem.

Run the authenticated wrapper from any working directory:

```sh
python -I -B run_geometry.py --repo-root /path/to/source/repository --output /tmp/geometry-normal.json
python -I -O -B run_geometry.py --repo-root /path/to/source/repository --output /tmp/geometry-optimized.json
```

The wrapper uses only the standard library to authenticate the component
and pinned repository-relative source files **before** importing SymPy or
executing a checker. It authenticates them again afterward. Each checker
uses explicit runtime guards, including under Python optimization. Normal
and optimized canonical reports match byte for byte. `requirements.txt`
records the exact symbolic dependency used here. The enclosing integration
wrapper must authenticate this component's complete manifest before
executing its wrapper, following the same source-before-import rule.

The five source pins establish the exact-family/common-domain, physical
shell, old intrinsic observable and intrinsic proper-clock conventions.
Their analytic theorem and actual-first-field premises remain conditional;
this component checks the finite geometric reductions and does not
computer-certify those analytic theorems.
