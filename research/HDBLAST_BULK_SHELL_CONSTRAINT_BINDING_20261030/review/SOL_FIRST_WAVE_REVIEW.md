# Independent first-wave normalization and initial-jet review

30 October 2026. Independent internal review. No frozen source is modified.
No producer result is imported by the independent checker.

The current independent result is **PASS**. No sign, raw-amplitude,
normalization, radial-drift, or initial-acceleration defect was found in the
registered first reconstruction and original first-gauge operator.

## Independent construction

`independent_wave_identity.py` uses its own standard-library `Fraction`
Laurent-polynomial engine and formal x,r differentiation. It isolates only
`actual_first_jets.frame_amplitudes` from authenticated source AST. It does
not import that numeric module, the producer checker, or any saved result.
The coordinate Lie derivative followed by parallel-frame conversion gives
all seven returned metric amplitudes and the plus scalar shift exactly.
The complete `R vartheta W` and `Bcoef vartheta(W_r/ell+kW)+pf`
reconstruction is checked separately with independent formal coefficient
symbols; `Bcoef` is distinct from the mixed metric amplitude Bmix. The
operator derivation then treats u as arbitrary, so it needs no master-wave
equation or scalar-reconstruction identity.
The actual `physical_derivatives` AST is separately checked in all six
0/y/t/yy/yt/tt rows, including `F_yy=L²F_rr+aLF_r`; variable rho is
differentiated by the independent engine rather than held constant.

The retained full-frame DC rows are independently transcribed from
`research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/theory/CONSTRAINT_FORMULA_NOTE.md`.
The first reconstruction uses an untraced Hessian, so the full23 notation is
explicitly converted by `J=L−8M/3`, `Htf=M`, and `beta=−10/3`.
This is essential in all three covector rows.

The resulting frame covector is compared to the raw operator. The original
operator's A,B,D matrices and force are extracted as isolated list-expression
ASTs from
`research/HDBLAST_FIRST_GAUGE_WAVE_20261016/operator/check_gauge_operator.py`.
All three raw residual rows and the source agree with the independent
transcription of the retained physical-half-line system.

Let `Res=U_xx−U_rr−P U_x−Q U_r−Rmat U`, with raw `U=(f,b,d)`, and
`s_f=−6 ell² rho² (u_y+4au)`. Exact collected coefficients vanish for

```
ell² rho² DC_y = s_f − Res_f,
ell² rho DC_0 = Res_b,
(ell² rho/c) DC_a = −Res_d D_aY.
```

The timelike normalization follows because a contravariant t operator
component becomes the frame covector `−rho (Pxi)^t`. The angular component
becomes `rho c (Pxi)^angular`; this accounts for its different sign and c.
The raw f radial drift is −3 ell rho a, whereas the b,d drifts are
−5 ell rho a. These differences follow from differentiating variable rho.

The formal reduction uses `rho_y=a rho`, `a_y=ay`,
`H_t=c^−2`, `(c^−1)_t=−H c^−1`, and `c^−2=1−H²`.
No background ODE simplification is needed for these identities.

## Same-operator initial specialization

The initial specialization sets the raw gauge values and velocities to zero
as radial functions, hence also their spatial/mixed derivatives. Its only
nonzero gauge second-time derivative is fixed by that same PDE:

```
f_xx=s_f,
f_tt=−6rho²(u_y+4au),  b_tt=d_tt=0.
```

The reconstructed metric then has `A=−4u,D=−2u,L=2u`, their corresponding
time velocities, and `Bmix_t=−6rho(u_y+4au)`; other mixed/angular initial
values and velocities vanish. The same DC rows vanish exactly.

All 21 metric value/y-derivative/t-derivative entries are also compared
directly against the completed initial-QC input table, extracted from the
`m`, `bt`, and `m.update` source ASTs at
`research/HDBLAST_REGISTERED_INITIAL_CONSTRAINTS_20261030/checks/check_registered_adm_binding.py:268`.
Every entry agrees. In particular setting f_tt=0 would change Bmix_t and
fails this join.

## Authentication and reproduction

Every directly read scientific source is compared to fixed SHA-256 pins
before evaluation and hashed again afterward. The reviewed source base is
Git commit `d98c641aabd658ff1bf456f29ae8b69886d52e54`, recorded as provenance;
Git is not required at runtime. The complete first-gauge
capsule manifest matches its copied manifest byte for byte; all 66 members
are authenticated before evaluation and read back after it. The tree is
closed to unlisted files and symlink members. Only isolated authenticated
function/list/table ASTs are evaluated; no scientific module is imported.

Normal and optimized final runs produce byte-identical JSON reports:

```
python -I -S -B independent_wave_identity.py --repo-root ROOT --output INDEPENDENT_NORMAL_FINAL.json --run-controls
python -I -S -B -O independent_wave_identity.py --repo-root ROOT --output INDEPENDENT_OPTIMIZED_FINAL.json --run-controls
```

Here ROOT is the repository/export root containing the retained research
capsules. Output paths must be new. The optional controls launch ten fresh
corrupted-formula subprocesses per mode, require the precise expected failed
identity and absence of a result file, and retain the exact nonzero rational
remainder in the report. They cover scalar sign, radial rho conversion,
raw/normalized f, angular clock term, radial drift, source sign, timelike
sign, angular c, variable-rho differentiation and zero initial acceleration.
The report additionally retains seven short independent formal witnesses.

## Producer review and fresh replay

The producer's separate checker is also reviewed. Its full five-dimensional
parallel-frame connection contraction is independent of this review's
retained-DC-row path. Its quadratic angular third-derivative rules and trace
elimination are consistent, the Lorentz contraction uses the negative time
inverse, and the same differential operator supplies the initial
acceleration. Its scalar coefficient Bcoef remains separate from Bmix.

Fresh reviewer executions of the producer checker pass all 77 exact checks
and all 15 specific corrupted subprocess controls in each Python mode, with
byte-identical reports. They include the requested raw-f normalization and
separate ell²/rho² source-scale traps. No code defect is found. The direct
21-entry completed-initial-QC-table join is supplied independently by this
review, rather than attributed to the producer's mathematical initial table.

## Dependency boundary

This check binds exact coefficient algebra to registered code. It does not
solve the gauge PDE or infer its exact validity from intervals. The actual
solution's neighbourhood PDE identity, zero Cauchy data, physical-half-line
existence and retained regularity remain analytic premises. The implication
`Res=s => DC1=0` is established by the exact identities above. It is stronger
than an interval containing zero, and it has a different role.

No shell/Codazzi source binding, subsidiary uniqueness, nonlinear existence,
cosmological-origin, finite-amplitude or physical-accuracy claim is made by
this independent first-wave review.
