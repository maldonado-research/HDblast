# Frozen actual initial mean: portable replay

This component regenerates the actual initial mean for the explicit bump
preparation. Its scientific reference files preserve every rational interval,
profile row, tail, target and cap from the frozen calculation. It does not
evaluate an evolved actual field or a shell observable.

The campaign remains `N=1024`, with 994 continuous packet cells, flat endpoint
tails of width `1/64`, and point targets `r=3/2,5/4`. No portability run changes
the grid, targets, bump, background uncertainty or source profile. The second
parameter derivative convention is retained throughout.

From any location, with a repository containing the pinned input capsules and
the public gauge API, write fresh outputs outside the source and code trees:

```
python initial/reproduce.py --repo-root /path/to/repository --output-dir /tmp/prepared-initial-normal
python -O initial/reproduce.py --repo-root /path/to/repository --output-dir /tmp/prepared-initial-optimized
```

The default gauge API directory is
`research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/gauge` relative to
`--repo-root`. A relocated API directory may be supplied explicitly with
`--gauge-api /path/to/gauge`. Both `mean_initial.py` and
`conditioned_mean_initial.py` must have their frozen bytes.

Each run writes `RAW_BOXES.json`, `CONDITIONED_RAW_BOXES.json` and a
`RECEIPT.json`. Outputs are created only after source custody and complete
scientific equality checks succeed. An existing output directory is rejected.
There is no Git checkout, Git HEAD, private path, private manifest or network
requirement. The code uses Python's standard library.

Before importing scientific code, the wrapper authenticates the essential
local code, campaign freeze, scientific reference files, both gauge adapters
and all seven files of the initial-input source capsule. It checks those bytes
again after calculation. The scientific outputs are path-free and identical
in normal and optimized Python; mode and process paths are not scientific
fields. The receipts identify the input hashes and exact source call counts.

The two normative scientific files are
`references/RAW_BOXES.json` and `references/CONDITIONED_RAW_BOXES.json`.
They include all 994 rows; they are not sample summaries. Comparison is exact
Python data equality across every scientific field, including all exact
endpoints, display-only decimals, target values, global caps and protocol
choices. Runtime timing and private source-custody records are represented
by the new portable receipts rather than copied into the scientific files.

## Initial source and arithmetic

The authenticated source module is
`research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/interval_core.py`.
Its complete transitive data-loading closure comprises `INPUT_PINS.json`
and the five files it names. `SOURCE_PINS.json` pins all seven files by
repository-relative path, byte count and SHA256. The source loader authenticates
the same declared input bytes itself. The stored origin strings in those old
source files are provenance only and are never used as fallback paths.

Each primary generation makes **1990 actual initial-background calls**:
994 integral cells, two point profiles, and 994 full-cell profiles. The
conditioned adapter makes zero additional background calls. Both stages make
**zero evolved actual-field oracle calls**. This count is checked at runtime;
an authentication control performs no scientific imports or oracle calls.

The exact outward arithmetic grid is `2^112`. The bump is
`Q=exp(4-1/((r-1)(2-r)))` on `1<r<2`, zero otherwise, with `ell=10^-6`.
The source's validated degree-28 exponential enclosure, derivative recurrence,
and monotone flat-tail bounds are retained. The initial background uses its
actual retained shell box and differentiated ODE recurrences, with its existing
uncertainty. The full-cell integrals use `dy=-rho dz`:

```
I = ell integral_1^2 rho^5 QH dr,
J = ell integral_1^2 rho^3 Q dr,
f = (2I/(9J)) Q.
```

With `KH=integral_yL^y rho^4 QH ds` and `Kb=integral_yL^y rho^2 Q ds`,
the exact, stable primitive is

```
xi = (3/2)(f_y-a f)/p^2 - [KH-(I/J)Kb]/(rho^4 p^2),
A=f+2xi_y, B=f+2a xi, tau_angular=-QM/(6a).
```

This avoids quadraturing the much larger `f''` term. Every primitive cell
retains its width. The exact analytic moment identity supplies compact support
and exterior zero; the code does not infer exterior zero by subtracting
independent interval boxes for `I` and `(I/J)J`.

The conditioned representation is a checked second harmonic gauge of the
same frozen initial boxes. Its displacement is
`(f,0,-u^2/3,f,-p xi)`. The additional global `p_y` and retained `p_yy`
come from the authenticated old actual-background coefficient certificate,
explicitly identified in the conditioned output. No new oracle evaluates them.
The full nonlinear gauge naturality and existence assertions remain retained
conditional inputs.

The integral `I` is strictly negative. Outward decimal summaries are
`I in [-1.945e22,-1.500e22]`, `J in [0.8395,0.8439]`, and
`2I/(9J) in [-5.146e21,-3.952e21]`. At `r=5/4`, conditioned
`w0=-p xi in [1.131e30,1.977e30]`. On `[0,7/4]`, the Gram C1 sum
`sup(|U0|+|V0|+|U0_r|)` is at most `1.855e31`. The exact raw endpoints
are authoritative. Retained background uncertainty limits precision.

These are mathematical second parameter derivatives. They provide no SI
conversion, effective amplitude, energy, physical generator or evolved shell
sign. In the original gauge, the large angular term `2a xi` cancels against
radial displacement in the anchored observable; its positive sign alone is
not a shell-volume conclusion.

## Bounded controls and a fresh relocated positive run

```
python initial/checks/check_portability_controls.py --repo-root /path/to/repository --gauge-api /path/to/gauge --output /tmp/prepared-initial-controls.json --relocated-positive-output-dir /tmp/prepared-initial-relocated
```

The controls run six bounded negative cases in both Python modes. They require
child flags `-I -S -B` (plus `-O` for optimized children), exactly one JSON line
on stdout, empty stderr, the exact full structured `FINAL` rejection witness,
exit code 2, and absence of any
new scientific output. The cases cover source code, gauge code, essential local
code, the independent reference layer, an existing output and output inside a
source tree. The source/gauge mutations contain deliberate import failures,
so their authentication witnesses also demonstrate rejection before import.
The optional positive run copies the entire seven-file source capsule, code,
references and both gauge APIs to fresh locations and performs the full frozen
generation. Its count is reported separately from the negative controls.

This portable revision changes only the control protocol and its records.
The earlier portable revision is preserved. Scientific code, source and gauge
bindings, references, grid, targets and the 1990-call generation protocol are
unchanged. The inherited primary and relocation receipts remain the scientific
replay evidence; this revision reruns the twelve bounded negatives under both
normal and optimized control drivers, for 24 executed rejected children and
zero initial-background or evolved-field calls in those negative processes.
A separate fresh relocated positive generation verifies `-I -S -B` execution
and makes the unchanged 1990 initial-background calls, with zero evolved-field
calls. The canonical control roster remains twelve negatives plus one relocated
positive; the second control-driver run is recorded only as extra packaging
validation. The old weaker stdout and
child-flag protocol is superseded, rather than described as strict evidence.
