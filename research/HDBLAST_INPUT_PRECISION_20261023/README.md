# Regular-shell and initial-data precision

23 October 2026. New bounded continuation of the retained five-mean
Einstein–scalar shell calculation. Earlier capsules and manuscripts are frozen.

The same regular branch at detuning `delta=1/1000` has the certified necessary
shell enclosure

```
-0.00008405268308665527 <= beta <= -0.00008405268308663626
radius = 1901/200000000000000000000 = 9.505e-18
```

This sharpens the old analytic detuning interval of radius `3.196e-6` by more
than `3.36e11`. The proof inherits existence and uniqueness in the specified
tube from the authenticated retained branch theorem. It uses an exact
degree-four polynomial graph, its complete defect, a positive damping bound,
a finite-cone memory term, and one necessary-enclosure contraction. It does
not assume convergence of a formal series or solve a singular shooting IVP.

The companion initial calculation forms actual-minus-nominal jets before
taking interval ranges. It separately charges shell/Taylor-coefficient
differences, rational packet/cutoff/exponential errors, and the validated
fifth-derivative ODE remainder with outward arithmetic. The old candidate and
its failed upper-bound certificate stay frozen. Any new nominal coefficients
are separately registered and reconstruct a single rational C1 initial/free
wave rather than unrelated per-cell values.

The completed 96-cell calculation gives
`e_free < 8.518093e27 < 1e28` for one explicitly registered recentered
candidate. The original candidate, assessed with the sharper branch interval,
still has an insufficient upper cap of `1.758e30`. This does not prove a lower
bound on its actual error. The new candidate changes only the fixed nominal
degree-four background coefficients; its packet/cutoff, physical formulas
and boundary reflections are unchanged. Candidate coefficient SHA256:
`739240283c2bee53afa559c4ce71c1e3356cf2d0ad4e4b0f18b9621c077e20d3`.

The suggested component allocations were uncomputed sufficient designs. The
packet radial bound exceeds its `4e27` design row and the total radial bound
exceeds `8e27`. Acceptance uses the exact combined free-wave expression below,
which passes; it does not assert that every illustrative row passed.

| Contribution | Value H cap | Radial H cap | Scaled time H cap |
|---|---:|---:|---:|
| Shell/initial coefficients | <1.590e19 | <1.085e21 | <2.424e11 |
| Packet and representation | <1.932e24 | <8.519e27 | <4.639e16 |
| ODE Taylor remainder | <2.594e14 | <8.806e15 | <6.228e6 |
| Total | <1.932e24 | <8.519e27 | <4.639e16 |

Every displayed upper cap is rounded outward from the exact report.

The norm and error target are those of the retained response preflight:

```
X(V) = max(sup|V|_H/3, sup|V_x|_H, sup|V_r|_H)
e_free = max((delta0_H + (3/2)*deltav_H)/3, deltar_H + deltav_H)
provisional initial allocation = 1e28
H Gram weights = (1,2,1,3,1), means = (K0,L0,N0,J0,w0)
```

Certification here means a bound proved using exact rational or outward
interval arithmetic under explicit model and upstream theorem premises.
Internal independent AI-assisted review checks the proof and implementation;
it is not external peer review, a journal decision, or a proof-assistant
formalization. Displayed digits are mathematical enclosures, not physical
measurements.

The evolved numerical total response, an intrinsic observable value/sign and
an exact nonlinear physical family remain unestablished. This checkpoint is
not evidence that a higher-dimensional blast caused the Big Bang. No new
Zenodo version or journal submission accompanies it.

The methods comparison checks relevant 2026 primary literature alongside
established validated-numerics methods. It is a bounded search and establishes
no novelty or priority claim. No third-party PDF or private editorial history
is included.

## Files and replay

- `manuscript/INPUT_PRECISION.pdf` and the two TeX files form the four-page
  working note; compile editable source only in a writable copy.
- `branch/ANALYTIC_BRANCH_COMPARISON.md`, `GRAPH_POLYNOMIALS.json` and
  `BRANCH_COMPARISON_REPORT.json` contain the exact proof/data.
- `branch/portable` recomputes the complete finite polynomial/barrier/root
  arithmetic without project imports, while explicitly inheriting branch
  existence. Its controls reject wrong coefficients, damping, inward rounding
  and altered numerical references.
- `initial` contains the 96-cell native evaluator, explicit candidate and
  three-channel error ledger, fixed-candidate comparison, proof and controls.
  `initial/NUMERICAL_REPLAY.md` documents all five full upstream capsules.
- `review` contains separate internal mathematical audits. The authentication
  addendum updates only the producer identity; the numerical gate is unchanged.
- `literature` records a bounded primary-source methods comparison including
  relevant 2026 work. `SOURCE_PROVENANCE.json` binds selected native artifacts.

Use an externally retained manifest SHA256 from the publication evidence:

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_HASH --output-dir FRESH_OUTSIDE_DIRECTORY
```

This authenticates the exact capsule before and after replay, runs finite
branch arithmetic, initial ledgers, model identities and independent algebra
controls in normal/optimized Python, and tests relocated branch rejection
controls. It does not reexecute upstream PDE proofs. Add
`--native-dependency-map /absolute/deps.json` to execute both normal/optimized
96-cell interval evaluations against all complete P/N/G/C/T closures, compare
all report fields except the two declared location strings, and test rejection
of unlisted caches and symlink paths. The new CI workflow uses this native mode.

The retained native branch scripts are audit originals with recorded workspace
paths and optional historical archive lookup. Their portable counterpart is
the location-independent finite arithmetic replay. The optional independent
native branch audit similarly requires the retained complete historical
source tree; no private numerical dependency is required by the default replay.

No script should write into a frozen source capsule. Build logs/render images
are audit artifacts; rendered images are retained with the native work rather
than copied into this source archive.
