# Independent postcontact numerical audit

At the fixed point `(x,r)=(5/4,1/8)`, the independently computed actual-target
enclosure is strictly negative:

```
h in [
 -2099303583227573148515637989828051/5192296858534827628530496329220096,
 -75816962575598594402102064535893/5192296858534827628530496329220096
]
display only: [-0.4043111633297405, -0.014601815851683174].
```

Normal and optimized Python reports are byte identical, SHA256
`e670411b8fcb4f8a8ec335c4759fd0694007ffb5cc3b56ab8e90f842097327b9`.
This remains a conditional formal response enclosure under the inherited
actual first-solution error and response existence/regularity premises.
It establishes no finite-amplitude shell-contact solution or gauge-invariant
observable, and makes no comparison between finite W1 and W2 accuracy.

## Independent computation and common inputs

`independent_numeric.py` imports only the authenticated old
`G/numeric/actual_first_jets.py` as its scientific implementation input.
That old module supplies actual common-gauge ordinary physical two-jets,
including its complete inherited C3 errors, physical one-sided shell traces,
background enclosures, and fixed unit `ell=10^-6`. This audit does **not**
independently reprove those inherited numerical inputs.

The native postcontact TT producer is neither imported nor executed by the independent numerical checker. The independent checker
owns integer outward interval arithmetic on grid `2^112`, literal Laurent
evaluation of all 78 terms in the original `metric_forcing/ab/5` data map,
the physical shell forcing formula, cone clipping, characteristic endpoint
integration, coefficient and feedback evaluation, and the compact quadratic
screen formula. Both original TT source orientation rows are checked equal.
No original source evaluator, shell evaluator or screen evaluator is imported.

All five source capsules in `SOURCE_PINS.json` are authenticated as closed
inventories before the old scientific import and again after the computation.
This includes the restored exterior capsule, manifest
`0bc16f65598ae4977a7b8ab70a66c13ec85d844c871d6fee7739c4ec472d1c70`,
with 45 listed files plus the manifest, 46 total members. The selected
background field-equation/Hamilton-constraint note is separately pinned.
The checker checks its own bytes again at completion. The source tree is
read-only to this work; all outputs are outside it.

## Cone, boundary history, and complete feedback

The target cone is
`0<=s<=T, max(0,s-(T-R))<=v<=T+R-s`, with `T=5/4,R=1/8`.
The fixed 8 by 8 grid over `[0,T] x [0,T+R]` was chosen before evaluation.
For an original cell `[a,b] x [c,d]`, its projected intersection is enclosed by

```
x in [a,min(b,T+R-c,T-R+d)],
r in [max(c,0,a-(T-R)),min(d,T+R-a)].
```

Empty intersections are omitted. These are boxes enclosing the intersection,
so they may overcover the cone locally. Their union covers it completely.
The resulting 43 boxes require no evaluation outside the old actual-jet
domain. No adaptive refinement or point-budget tuning occurred.

The derivative uses two endpoint characteristics with **positive** endpoint
weights. The uniform 32 time subdivisions plus the exact reflected turn
`s=T-R=9/8` give 33 segments and 66 source boxes. Whole-cell source intervals
times exact rational segment lengths are integration enclosures, not point
quadrature estimates. The inhomogeneous Neumann contribution is retained:

```
G = -(4/3) ell rho (a A - 2p s/3) M,
(L G)_x(P) = -G(9/8),
V_x(P) = -G(9/8) + (1/2) integral [F(s,T+R-s)+F(s,abs(s-(T-R)))] ds.
```

Here the scalar field amplitude `s` inside the shell formula is distinct from
the integration-time symbol. Eight whole-history shell intervals and one
endpoint interval are evaluated. Their intervals include zero; this is not
a proof that the generally nonzero physical source vanishes. A rational live
implementation control detects the lift sign even though changing the sign
of the symmetric numerical endpoint interval would be invisible.

The coefficient calculation independently retains
`A=-3ell tanh(t)`, `B=-3ell rho a` and the full original `D` expression.
The original field equation `a_y=-rho^-2-p^2/3` and Hamilton constraint
`V=p^2/2-6a^2+6/rho^2` reduce `D` to `-8ell^2/cosh(t)^2`. The reduced interval
is recorded as a conditional cross-check; the primary error bound uses the
full unsimplified interval. Thus no sharpened correlated background boxes
are inferred from interval precision.

With `X=max(2 sup|u|/T,sup|u_x|,sup|u_r|)`, the cone is backward closed under
both endpoint families and `X(K_N f)<=T sup|f|`, `X(LG)<=2 sup|G|`.
The complete scalar TT feedback contraction is

```
sigma = sup|A| + sup|B| + T sup|D|/2,
q = T sigma < 1,
X(U) <= (2g+Tf)/(1-q),
|U_x(P)-V_x(P)| <= q (2g+Tf)/(1-q).
```

This bounds every omitted iterate. It is not a one-iterate error estimate.
Approximate displays of the exact rational bounds are `f=0.0826459738028`,
`g=1.32339556112e-21`, `q=5.42651760331e-5`, and feedback error
`5.60630212302e-6`. The complete `U_Tx` enclosure is
`[-0.00454993145347068,-0.00016432353726657007]` (display only).

## Fresh point quadratic and normalization

The fixed quadratic is reevaluated from fresh actual first jets at this point:

```
Rphys=rho c,
Fquad=((D+4M-L)Mt-2M Lt)/(4rho)-E(D+L)/(2Rphys)
      +(C/Rphys-My/2)(Rphys(Ct-theta C)/(2rho)-Rphys Ey/2+B/2),
h=3 U_Tx/(2ell rho^2 c)+4Fquad/Rphys.
```

These derivative suffixes are ordinary physical derivatives. The second
parameter derivative has no additional factor `1/2`. Gauge-dependent first
field amplitudes and all terms of the fixed quadratic are retained. The
fresh point quadratic contribution is positive,
`4Fquad/Rphys in [3.48148648036e-18,1.42037292928e-17]` (display only).
The opposite magnetic orientation is exactly `-h`; equal TT metric source
rows do not imply magnetic cancellation.

There are 119 old actual-jet calls: 43 source-cap boxes, 66 characteristic
boxes, eight shell history boxes, one shell endpoint, and one actual apex.
The maximum remained the frozen 140-call budget. The independent interval
and the producer interval differ by tiny outward-rounding amounts because
of different operation ordering and the clipped cover. Their displayed
endpoints agree; exact equality is neither assumed nor claimed.

## Separate report audit and real controls

`report-algebra/audit_report.py` performs 202 exact report-only checks against
the producer report, with no scientific producer import. It checks every
endpoint box, source and shell cap, integral/seed aggregation, tail, point
quadratic, scale and final orientation. Its normal and optimized reports
match. Its five altered-report controls are separately scoped to report
consistency, not independent numerical evaluation.

`run_controls.py` creates real checker-source mutations and runs each in both
Python modes. All twelve subprocesses reject without success artifacts:

| Alteration | Exact rejected witness |
|---|---|
| Boundary lift sign | `got=7/3; expected=-7/3` |
| Reflected endpoint sign | `got=0; expected=5/4` |
| Omit fixed gauge quadratic term | `got=0; expected=-1/2` |
| Omit B and D scalar feedback from sigma | `got=1/7; expected=2969/5544` |
| Change a closed source member | `CLOSED_MEMBER:.../README.md` |
| Replace that member with a symlink | `CLOSED_MEMBER:.../README.md` |

The first four controls exercise live implementation functions at exact
rational fixtures before scientific imports. The last two corrupt isolated
source copies and fail closed authentication before any scientific import.
They are altered-implementation checks, not purported alternate physical
solutions or fresh completed numerical campaigns with changed equations.
`CONTROLS.json` retains their full exact failed witnesses, exit codes, and
mutant checker hashes. `REPLAY_RECEIPT.json` binds the main artifacts and
separately records the comparison with producer report data.

A further code and mathematics review checked the independent checker and
its aggregation proof. `report-algebra/REVIEW_NUMERIC_CODE.md` proves the
cone projection and backward closure and checks every normalization.
`NUMERIC_CODE_CHECKS.json` records 4,016 focused exact arithmetic and cover
cases, including 210 negative reciprocal intervals and even powers. No
blocking issue was found. The included `check_numeric_code.py` harness
replays all 4,016 checks, including the inverse and even-power tests.


## Portable files and replay

The independent audit uses Python 3.9 or newer and the standard library.
Keep `SOURCE_PINS.json` beside `independent_numeric.py` and retain the
`report-algebra` subdirectory structure: the focused code-check harness
loads `../independent_numeric.py` relative to its own file. No installed
packages or network requests are needed. The source repository must contain
the five complete pinned capsules and selected background note at the
`research/...` paths recorded in `SOURCE_PINS.json`.

The public file set is exactly:

```
independent_numeric.py
SOURCE_PINS.json
NUMERIC_NORMAL.json
NUMERIC_OPTIMIZED.json
run_controls.py
CONTROLS.json
PUBLIC_NUMERIC_REVIEW.md
REPLAY_RECEIPT.json
report-algebra/audit_report.py
report-algebra/ALGEBRA_NORMAL.json
report-algebra/ALGEBRA_OPTIMIZED.json
report-algebra/check_numeric_code.py
report-algebra/NUMERIC_CODE_CHECKS.json
report-algebra/REPORT_ALGEBRA_REVIEW.md
report-algebra/REVIEW_NUMERIC_CODE.md
```

`REPO` below denotes the source checkout root, and `PRODUCER_REPORTS` the
directory containing the native producer's `POSTCONTACT_NORMAL.json` and
`POSTCONTACT_OPTIMIZED.json`. Run from the independent audit directory.
Every output name must be fresh. The numeric, report-algebra, and code-check
outputs need an existing parent directory; the controls output directory
must not exist and is created by its runner.

```sh
mkdir replay
python independent_numeric.py --source-root "$REPO" --output replay/NUMERIC_NORMAL.json
python -O independent_numeric.py --source-root "$REPO" --output replay/NUMERIC_OPTIMIZED.json
python run_controls.py --source-root "$REPO" --output-dir replay/controls-run
python report-algebra/audit_report.py --producer-dir "$PRODUCER_REPORTS" --output replay/ALGEBRA_NORMAL.json
python -O report-algebra/audit_report.py --producer-dir "$PRODUCER_REPORTS" --output replay/ALGEBRA_OPTIMIZED.json
python report-algebra/check_numeric_code.py --output replay/NUMERIC_CODE_NORMAL.json
python -O report-algebra/check_numeric_code.py --output replay/NUMERIC_CODE_OPTIMIZED.json
```

The numerical modes produce identical canonical reports. The report audit's
`negative_report_data_controls` dictionary contains five distinct altered
report cases, each with `rejected` and `first_failed_witness` fields. Each
mode evaluates those same five cases in process. These are five distinct
report-consistency controls, or ten mode evaluations; they are separate
from the twelve real implementation/authentication subprocess controls in
`CONTROLS.json`. The 4,016 focused arithmetic and cover checks likewise
belong to their included replay harness and are not counted as physical
source mutations or new physical samples.
