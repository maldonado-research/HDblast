# Independent reduced-source and covector review

31 October 2026. Private finite-binding review. No old source was edited and
no publication was attempted. The registered source-role identities and the
eight retained covector wave rows pass a separate exact implementation. No
mathematical defect was found in the inspected reduced-Bianchi signs or
conditional specialization premises.

`check_independent_source_roles.py` uses only the Python standard library. It
implements its own sparse Laurent arithmetic, exact `S3` quotient reduction,
intrinsic covariant differentiation, normalized angular bases, full physical
frame connection and all index contractions. It imports no scientific
producer, historical arithmetic helper, sphere helper or Bianchi checker.
The source map is read as authenticated data, rather than reconstructed from
the identities to be tested. The new operator is obtained by directly
contracting two full covariant derivatives, including the derivative slot.

The executed positive checks are:

* All eight modes of `QC_actual-QC_background=2 q.C1`, with `C1` independently
  obtained from the linear connection contraction. The complete
  `QC_actual_d=-2 q^{ab}P_dab` contraction is independently checked too.
* All twenty retained tensor modes of `PC1=2 P.C1` and
  `gauge_gradient=sym nabla QC_actual`, including the two TT modes. Temporal
  index raising has the Lorentzian sign.
* All twenty coefficient rows of each registered role
  `QEred_general=QE-gauge_gradient+PC1` and `QEred=QE-gauge_gradient`.
* The three scalar rows of
  `Qphi_general-Qphi=2[(q^{cd}C1_d)(dphi0)_c-C1^c(deta)_c]`.
* Full angular reconstruction of the above quantities. Independent diagonal
  Gram matrices have nonzero entries in the three scalar, two gradient and
  seven tensor bases, so equality of full components binds every retained
  modal coefficient. Both PC1 TT coefficients have 67 collected terms;
  deleting a TT source is tested at the all-twenty role level.
* `B0(sym nabla Z)=(Box0+Ric0)Z/2` for eight arbitrary lower-covector modes
  with all 48 free ordinary jets through order two, with direct rough-wave
  contraction, the scalar-Einstein radial `p^2` term, physical rho scaling
  and the exact Einstein and de Sitter quotient relations.

The first metric source roles concern the retained seven-profile `l=2`
ansatz, with its 42 ordinary jets through order two free. They do not claim
that seven profiles span an arbitrary spacetime metric. The arbitrary
eight-mode operator calculation is separate from that first-field ansatz.

The directly contracted wave is independently projected and reconstructed.
Every one of its eight exact coefficient polynomials agrees with
`/workspace/hdblast-covector-wave-31a/WAVE_NORMAL.json`, after renaming the
independent profile symbols. The producer report is read as data only; the
producer wave helper is not imported. The compared report and helper bytes
are SHA256 pinned before and after comparison.

A coherent pointwise-lowering counterexample is also executed. For
`q=A(y)Y dy^2`, all other first metric profiles identically zero,
`A=0, Ay=1, Ayy=Ayyy=0` and vanishing time jets give equal actual/background
QC values at the point. At the consistent background point
`rho=p=1,a=1/2,ay=-4/3,V=5,Vp=0,c=1,theta=0`, direct contraction gives
`L0(QC_actual-QC_background)_y=10Y^2`, which is nonzero. This uses commuting
ordinary profile jets rather than invented independent covariant C1 jets.

Normal and optimized isolated subprocesses produce byte-identical canonical
reports. Twenty distinct corruptions per mode, forty actual subprocesses in
total, all fail with a collected nonzero witness and no successful output.
They alter lowering/PC1 factors and temporal raising, use the background QC
prematurely, remove symmetrization normalization, reverse the gauge sign,
drop PC1, delete TT source modes, select natural/principal/harmonic sources
at the wrong stage, omit radial `p^2`, change the temporal Ricci sign, omit
the operator half factor, remove derivative-slot/warp/FRW connections, use
positive time contraction or omit physical rho scaling. Authentication-only
failures do not count as these mathematical rejections.

The source map, its directly named provenance sources, all 29 files of the
retained natural capsule and the independent checker are authenticated before
and after computation. The sole nonstandard execution unit is this checker;
its AST has no assertions or dynamic-code entry calls and permits only its
explicit standard-library import set. The runner authenticates its pinned
checker and itself around the fresh replay batch; every child authenticates
its full selected source closure before and after calculation. The receipt records actual
commands, return codes, exact nonzero witnesses, absent rejected outputs and
the independent eight-row wave comparison. Runtime flags and receipt data
remain outside the canonical positive reports.

Reproduce with a fresh directory:

```sh
python -I -S -B /workspace/hdblast-sol-reduced-31a/run_independent_review.py \
  --root /workspace/hdblast-reduced-publish-31a \
  --wave-dir /workspace/hdblast-covector-wave-31a \
  --output-dir /tmp/hdblast-sol-reduced-fresh
```

Canonical reports:

| Artifact | SHA256 |
| --- | --- |
| `check_independent_source_roles.py` | `6c1384d380127c89149e5752da839dd280ad507f3cf9550de06d686f0cbe7cef` |
| `final-verified/SOURCE_ROLES_NORMAL.json` | `e28da8b19da1785bfbb19248168d1173cea145acb0a593f877ef16abe40e1188` |
| `final-verified/SOURCE_ROLES_OPTIMIZED.json` | `e28da8b19da1785bfbb19248168d1173cea145acb0a593f877ef16abe40e1188` |
| `final-verified/RUNTIME_RECEIPT.json` | `1d5781749fc2d3330ecbdca8a8e0927a2377b2c8c81c0a9745d2f015046ad95b` |
| `ANALYTIC_SIGN_REVIEW.md` | `56c297c1446e7ed7302a88940af2ab869bd738fc513294f5a0451919cdd7c4df` |
| Compared producer `WAVE_NORMAL.json` | `226117368a98a7a01255547147fcbf56d9019e1e2d29880136bd5592f8a3f549` |

The separate source producer's general residual check was inspected at
`/workspace/hdblast-reduced-bianchi-31a/checks/check_reduced_source_binding.py`
SHA256 `8863101fa6733836f50ea06d9067dc64f4ee639f5629c0928a4456df32fec017`.
Its `E1`, `nabla E1`, `S1` and covector wave are independently constructed,
and `QEred_general` comes from the registered map. Thus its target is not
made zero by defining one of those operands through the target identity.
The lower P/C1 contraction, scalar-general correction, full angular role
reconstruction, B0 map and final signs match the independent derivation.
Its wave helper is now included in its complete pre/post import/AST closure.

The general left side is
`B0(QEred_general)+p Qphi_general e^y+2 deltaB_q(E1)+2 S1 d eta+(1/2)L0 QC_actual-B0(PC1)`.
The conditional homogeneous companion needs `E1=0,nabla E1=0,S1=0` from
the first-natural neighbourhood identity and separately
`C1=0,nabla C1=0` from harmonicity. It retains QC_actual. Exchanging to
QC_background inside L0 additionally requires second C1 prolongations as a
sufficient premise. A pointwise equation value supplies no such unproved
prolongation. `ANALYTIC_SIGN_REVIEW.md` gives the independent derivation and
the precise sufficient-premise table.

This review supplies finite source/operator evidence for an existing analytic
subsidiary argument. It is not external peer review, a new PDE existence or
uniqueness result, a new nonlinear solution family, a numerical accuracy
improvement or a cosmological origin conclusion. The actual numerical third
jet adapter is a separately implemented/reviewed task; this checker neither
imports it nor claims to validate its numerical error certificates.
