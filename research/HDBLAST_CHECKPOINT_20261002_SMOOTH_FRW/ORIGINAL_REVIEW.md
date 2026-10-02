# Independent internal review: smooth prescribed FRW control

The registered matrix has overall status **FAIL**. The flat control A=0 passes
the registered gates after the K192 extension. The curved A=0.2 control fails
the final pressure cutoff gate. Independent mode evolution and artifact integrity
checks pass. These are distinct conclusions; the latter do not override the
registered pressure failure.

## Scope, chronology and provenance

The review prepared its own Radau implementation, analytic compact-step
derivatives, local-action variation checks, and formal WKB grading/finite-cutoff
trace checks before physical evolution. It imported no primary numerical source
into the Radau solver. Exact symbolic checks executed before registration are
algebraic derivations, not prospective empirical outcomes.

The numerical registration was committed at
983c472dfe70d2a830b8e9b538995732dc38626f before mode runs; its SHA256 is
37fda5b30332d0079ec44db708e266788fdf5859303707bb709750fdcf59d898.
The first integrated run completed physical mode evolution but failed during
counterterm evaluation because NumPy-left arithmetic returned an ndarray instead
of a Jet. The repair, committed at 57a97a2, adds array dispatch priority and leaves
the scientific expressions, settings, and gates unchanged. The first failure,
physical modes, and source provenance are retained.

Independent Radau results use the unchanged registered source. Two comparison
outputs retain both original and repaired primary-source hashes; their numerical
differences are identical. The repaired matrix source SHA256 is
32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a.

This review is internal independent implementation/review, not external peer
review. It does not establish priority or physical validity of HDBLAST.

## Analytic review

Independently verified:

- The full fourth-order graded WKB expansion, second-order Q, paired source sign,
  and finite-cutoff trace coefficient functions c5 through c13.
- The local covariant action variations for V, F R and alpha R^2, including direct
  pressure formulas regular at H=0. The alpha trace shift is -12 alpha box R.
- The scalar heat-kernel anomaly polynomial, its flat reference-mass shift,
  exact finite-cutoff moment primitives, and their continuum limits.
- Direct pressure logarithmic matching, the C^2/Euler decomposition, and the
  displayed background/w''''/u2'' derivative formulas.

The common-action bridge specifies a finite convention; it cannot eliminate the
freedom to add finite local terms. The FRW calculation cannot independently fix
the Weyl/Euler extension. Heavy-field decoupling is the stated standard asymptotic
argument for bounded smooth, gapped, exactly static-prepared backgrounds; an
explicit certified uniform remainder constant or full operator proof is not
supplied. The document preserves this qualification.

An isolated H=0 event is not flat space when derivatives of H are nonzero; the
flat-trace wording was corrected before execution. Also x=r(2B-1)^2 touches zero
without changing sign. A negative canonical frequency squared can instead arise
from the curvature term -a''/a; the physical ODE remains regular.

## Independent physical modes

All registered selected-mode Radau controls pass. Maximum internal raw
Wronskian error is 1.80e-14. Maximum DOP853/Radau absolute differences are
1.09e-13 in u, 1.66e-12 in u', 1.59e-11 in bare rho, 5.41e-12 in bare p and
9.90e-14 in bare Q. These are below the registered gates.

Exact incoming-vacuum subtraction, static-future coherent Bogoliubov
decomposition, and deliberately omitted current/pressure checks pass. Bare
exchange residuals from direct ODE derivatives check algebra/normalization and
must not be described as independent continuum-conservation evidence.

## Matrix result and pressure limitation

| Control | Final pressure cutoff change, K96 to K192 | Allowed change | Result |
| --- | ---: | ---: | --- |
| A=0 | 0.0008291985 | 0.0030376492 | PASS |
| A=0.2 | 0.0226242834 | 0.0112327284 | FAIL |

For A=0.2 the final rho and Q cutoff differences are 7.40e-5 and 6.20e-5,
respectively, and both pass their stated gates. Solver, quadrature, time, raw
Wronskian, finite-cutoff energy exchange, pressure trace, static null, and future
occupation-energy checks pass. The pressure cutoff failure is the outstanding
registered acceptance failure.

A read-only examination locates the maximum curved pressure shell difference
at eta=0.059375, early in the transition. The pressure shell maxima for
K24→48, K48→96 and K96→192 are 0.6176240, 0.1825564 and 0.02262428. Their last
ratio is 0.12393, which does not justify assuming a resolved K^-2 asymptotic tail.
The finite-K trace decomposition shows this shell difference is dominated by
the derivative-dependent Q response, not by its small energy or explicit
anomaly-shell terms. This is descriptive analysis of existing output, not an
additional registered experiment or a certified tail estimate.

The selected curved finite-K trace residual is 8.46e-13 using physical-equation
Q derivatives and 1.49e-5 using independently sampled Q derivatives. The
renormalized conformal-energy ledger residual is 1.46e-7. Their small size does
not cure the pressure bandwidth limitation: exact identities hold at every
fixed comoving cutoff.

Future pressure and Q retain substantial coherent contributions. Maximum
departures from occupation-only estimates are 0.05055 and 0.02348 in the curved
case, and 0.14742 and 0.01563 in the flat case. An occupation-energy plot cannot
replace these observables. These finite-K quantities are not a thermal-fluid
model or a demonstrated radiation era.

The static vacuum null residual maxima are 7.92e-9 in rho, 2.74e-9 in p and
1.17e-12 in Q. The refined matrix has much smaller solver/quadrature changes
than its pressure cutoff discrepancy. Increasing precision in the final
subtraction alone cannot recover information lost by a physical mode solver;
the registered tight and half-step checks appropriately examine that risk.

## Read-only output audit

The independent archive audit loaded all 13 completed cases with
allow_pickle=False, checked finite numeric arrays and matching physical/observable
archive grids, reintegrated saved mode quantities, recomputed raw Wronskians,
and independently recomputed future occupation energy. This integrity audit
passes while retaining matrix status FAIL.

Its first attempt assumed u,v were in arrays.npz; the committed implementation
stores them in physical_modes.npz. The first KeyError log is retained. The
reader was corrected to combine both files and check overlapping arrays agree;
no research data or numerical evolution was altered. This post-execution audit
is explicitly separate from registered numerical experiments.

No unregistered cutoff extension was executed. The bounded outcome supports a
verified implementation and a passing flat benchmark, with curved pressure
requiring a separately specified follow-up. It does not support claiming a
fully resolved curved continuum source, changing the archived shell, or asserting
self-consistent backreaction, decay, thermalization or cosmological success.
