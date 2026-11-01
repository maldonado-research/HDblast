# Weighted-family monopole initial data

Continuation of the pinned exact-family TT transfer, 1 November 2026.

The [derivation](WEIGHTED_MEAN_INITIAL_DATA.md) supplies the same selected
family's explicit radial second-constraint sources, fourth-order weighted
metric equation, second-order weighted momentum equation, weak formulations
and physical Cauchy conversion. Its scalar field and coordinate-time scalar
velocity both have zero second parameter derivative initially; its physical
normal scalar velocity has a nonzero prescribed conversion term.

The published selection does not instantiate the radial defining function
or weight parameters. For every fixed admissible choice, the retained theorem
selects a unique weighted correction, but the source does not determine one
numerical mean profile. Explicit compact radial homogeneous corrections and
the weight-variation formula describe the remaining completion dependence.

A result independent of those missing choices is available for the actual
current packet. With `gamma2_mean=diag(A,B,B,B)` and first scalar amplitude
`v`, the exact areal-radius adjoint identity and retained coefficient boxes
give

```
integral rho^2(A+2B)dy
 <= -(49/150) integral rho^4 v^2dy < 0.
```

This is a gauge-fixed initial monopole moment. It is not a shell-volume
response, a radiation-energy measurement, a numerical scalar history or a
new model. The moment identity has an earlier source in the repository;
the strict estimate here applies it to the current displacement and velocity
packet, without choosing unrecorded weights.

Replay requires Python with SymPy 1.14.0:

```sh
python -I -B checks/check_weighted_mean.py --repo-root /path/to/HDblast --output /tmp/weighted-mean-normal.json
python -I -O -B checks/check_weighted_mean.py --repo-root /path/to/HDblast --output /tmp/weighted-mean-optimized.json
```

Both output paths must be new and outside the component and repository source
roots. Omit `--output` to print the complete canonical
JSON report to stdout. The source tree must contain the six unchanged public
repository files listed in [SOURCE_PINS.json](SOURCE_PINS.json); it need not
have the pinned commit as its current HEAD. No private source copies or
network access are required. Authentication occurs before loading SymPy and
again after all checks, including a check that the pin file stayed unchanged.

The script authenticates six repository source files and performs 34 exact
identity/bound checks and eight wrong-formula controls per mode. It includes
an independent full curvature Hessian calculation, canonical conversion,
formal adjoints, integral identities and rational sign bounds. The normal
and optimized reports are equal and retained in `results/`. Seven additional
entrypoint controls verify authentication before the SymPy import,
authentication after evaluation, pin-file stability, preservation of existing
reports, rejection of output paths inside code/source roots and rejection
of a deliberately incorrect radial operator:

```sh
python -I -B checks/check_portable_io.py --repo-root /path/to/HDblast --output /tmp/weighted-mean-io-normal.json
python -I -O -B checks/check_portable_io.py --repo-root /path/to/HDblast --output /tmp/weighted-mean-io-optimized.json
```

These controls use temporary copies of public inputs and leave the supplied
repository unchanged. Their two reports are equal. They are separate from
the 34 scientific checks and eight wrong-formula controls. Every rejection
requires exit code 1, an exact match to the intended final exception line,
empty stdout and an absent or byte-preserved output artifact. Each report
retains these witnesses in seven structured rows; temporary path prefixes
are normalized only after the exact raw comparison passes.

No weighted boundary-value solve or evolution was run. Analytic inverse and exact-family existence remain
conditional retained inputs. Internal AI-assisted checks are not external
specialist review.
