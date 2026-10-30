# Conditional natural-equation residual from the W1q error certificate

This is a separate error-transfer route for the unchanged `W1_quadratic`
candidate. It does not reuse the filtered `Z_epsilon` candidate, filtered
source, or filtered candidate jet bounds. It retains the original target,
full original quadratic source, physical shell law and coefficient domain.
The resulting bound is conditional on an additional formal-target
subsidiary identity. That identity is named, not proved by an interval box.

## Exact natural equations and additional premises

Use the natural trace-reversed Einstein-scalar and scalar equations

```
E_ab=Ric_ab-dphi_a*dphi_b-(2/3)V(phi)g_ab,
S=Box_g phi-Vprime(phi).
```

The metric equation is not the Einstein tensor `G_ab`. With
`g=g0+lambda*q+(lambda^2/2)*W` and the analogous scalar expansion,
`E2[W]=DE[W]+QE` and `S2[W]=DS[W]+Qphi` denote second parameter
derivatives. Their `lambda^2` coefficients are half these quantities.

Let `La` be the actual registered unit-principal lower-order operator and

```
R[V]=V_xx-V_rr-La V-F.
```

The complete original source normalization is
`F_metric=-2ell^2rho^2 QEred`, `F_scalar=ell^2rho^2 Qphi`. The primitive
scalar assembly already includes the natural `-DC[V].dphi0` coupling.
The exact previously verified identities therefore give

```
E2[W]=R_metric[W]/(2ell^2rho^2)+sym(nabla DC[W])+(QE-QEred),
S2[W]=-R_scalar[W]/(ell^2rho^2).
```

Here `sym(nabla C)_ab=(nabla_a C_b+nabla_b C_a)/2`. The actual first
harmonic identity `C1=0` must hold on a neighbourhood; it is not an
independent pointwise selection from a source interval box. Under that
retained premise, the complete 20-row source-map identity is

```
QE-QEred=sym(nabla QC).
```

The further formal-target premise is

```
DC[U]+QC=0 as a C1 covector identity on the physical domain.
```

It permits differentiating the identity. The retained target and first
field have the required physical-side regularity. For `e=W1q-U`, shared
quadratic source cancellation then gives

```
E2[W1q]=R_metric[W1q]/(2ell^2rho^2)+sym(nabla DC[e]),
S2[W1q]=-R_scalar[W1q]/(ell^2rho^2).
```

If another report uses `e=U-W1q`, the gauge-gradient sign reverses and its
positive norm bound is unchanged. The preceding finite-candidate reduced
defect is recorded with the opposite sign `r1=-R`; its norm cap `F0`
therefore applies. No source-averaging discrepancy is present. There is
no `2F_H` charge and no filtered-candidate difference term in this route.
The equality itself does not assert that the candidate residual vanishes.

## Converting the coordinate error into a frame covariant derivative

The previously authenticated exact harmonic map is

```
P=ell*rho*DC[e]=Rx e_x+Rr e_r+Zmat e,
Zmat=alpha Za+eta Zeta+zeta Zzeta,
alpha=ell*rho*a, eta=ell*tanh(ell*x), zeta=ell*sech(ell*x).
```

The input is the full23 positive probability-normalized S3 H norm. The
covector output Gram has no mixed-tensor factor two. The exact derivative
map gives the rational joint bound `qder=15/8>sqrt(7/2)` when both derivative
inputs have the same cap. The inherited weighted Frobenius squares of
`Za,Zeta,Zzeta` are `245/2,131/2,68`, so the inclusive operator caps
`12,9,9` suffice. Write the resulting operator caps as `z0,zx,zr` for
`Zmat,Zmat_x,Zmat_r`.

For the retained C1 error norm
`X(e)=max(sup||e||_H/3,sup||e_x||_H,sup||e_r||_H)<E1`, use

```
V0=3E1, Vx=Vr=E1,
Vxx=Exx, Vxr=Exr, Vrr=Err.
```

This calculation deliberately uses the immediate `3X` value bound. It
does not introduce the stronger zero-wedge `3X/4` value estimate or need
an additional support-path argument. The new C2 error report supplies all
three second-error upper caps; an absolute candidate jet must not be
substituted for an error jet, or vice versa.

Consequently

```
P0=qder*E1+z0*V0,
Px=qder*max(Exx,Exr)+z0*E1+zx*V0,
Pr=qder*max(Exr,Err)+z0*E1+zr*V0.
```

The actual background identities `rho_r=-alpha*rho`, `rho_x=0` imply

```
e_y DC[e]=-(P_r+alpha*P)/(ell^2rho^2),
e_0 DC[e]=P_x/(ell^2rho^2).
```

The radial product term has a plus sign inside the numerator. Converting
component derivatives to covariant derivatives also retains the background
connection and angular derivatives. The retained rough angular eigenvalue
cap is24 and the positive covector connection Gram bound is
`sum_i |Omega_i C|^2<=3(a^2+h^2)|C|^2`. The resulting conservative charge is

```
connection=2alpha_abs+5ell+(7/4)(alpha_abs+eta_abs),
Gerr=(Px+Pr+connection*P0)/(ell^2rho_lower^2).
```

The norm of symmetrization is at most one in the positive Frobenius norm.
The complete original reduced-defect H cap bounds each metric/scalar
projection. Thus the conditionally natural residual majorants are

```
cap(E2[W1q])=F0/(2ell^2rho_lower^2)+Gerr,
cap(S2[W1q])=F0/(ell^2rho_lower^2).
```

These formulas are an explicit bridge from a C2 *error* bound to natural
equation coefficient residuals in the declared background frame. They are
not an independent proof of the full constraints. The scalar bound can
already be obtained directly from the original reduced defect; a tighter
metric C2 error affects `Gerr`, not that retained scalar cap.

## Numerical bound for the supplied unchanged-candidate error report

Using the new combined initial-data/C2 report, with SHA256
`85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923`,
and all retained C1, reduced-defect, actual-background and harmonic-map
inputs, the exact rational majorants give

```
sup_D+ ||E2[W1q]||_positive,background-frame,L2(S3) <1.476e35,
sup_D+ ||S2[W1q]||_background-frame,L2(S3)          <2.118e30.
```

The unrounded upper caps are approximately `1.4756669711e35` and
`2.1173167933e30`. They are conditional large absolute upper bounds.
They do not show that the natural equations are satisfied accurately.
No reduction factor is asserted against the older filtered-candidate
example, which uses a different candidate and a different direct route.

The new normal and optimized reports are byte identical and contain 65
explicit checks. All 20 metric source differences and 23 forcing
normalizations are retained. The scope-control executor rejects 18 actual
corrupted executions in each mode: two missing conditional premises,
eleven incompatible candidate/source/target/domain/norm/error inputs,
and five wrong-formula injections. Explicit premise flags record the
conditional theorem input; they are not a numerical proof of the premises.

## Source locations and what is inherited

Paths below are relative to the HDblast repository root:

- `research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/theory/NATURAL_RESIDUAL_NOTE.md`
  supplies the direct natural/reduced identities, angular/connection bounds
  and frame conversion; its numerical example is for a different candidate.
- Its `checks/check_natural_residual.py` verifies all20 source-map
  differences and the metric/scalar forcing factors. The old CLI expects
  filtered absolute jets and adds the filtered-source discrepancy. Passing
  W1q/error jets through that old schema would mix candidates or quantities.
- `research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/theory/CONSTRAINT_TRANSFER_NOTE.md`
  already separates direct constraint values from the premise-based
  `DC[W]+QC=DC[W-U]` route. Its source maps and positive norms are retained.
- `research/HDBLAST_COMMON_GAUGE_FORCING_20261016/evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`,
  section7, gives the conditional subsidiary argument: actual first
  equations/harmonic identity, trace-reversed Bianchi identity, physical
  initial Hamiltonian and four momentum constraints, initial gauge
  completion, Codazzi/junction compatibility and the homogeneous
  mixed-boundary subsidiary uniqueness argument yield `DC[U]+QC=0`.
- Its `inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md` and
  `inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md` provide the inherited
  analytic physical/initial gauge identities. Hash authentication verifies
  that those inputs have not changed. The new numerical error or residual
  calculation does not independently rederive every one of those identities.

The new checker independently verifies the full polynomial source
differences and forcing factors, authenticates the old matrix/background
reports, recomputes their charges and performs exact rational error
transfer. Explicit CLI premise flags state the two additional identities.
Their presence is not evidence that the physical constraints were proved
by that execution. Missing flags, changed candidate/source/target identities,
incomplete norms and corrupted formula variants must fail explicitly under
normal and optimized Python.

## Scope, normalization and the next physical task

The output is a uniform positive Frobenius norm of background-orthonormal
frame components, with probability-normalized angular L2. It is not a
Lorentz-invariant scalar, an observed curvature signal or a force. It is a
second-parameter coefficient residual, not a residual of an exact finite
amplitude spacetime. Large absolute caps do not establish useful physical
accuracy or relative smallness.

The proposed small ratio of global coordinate second-derivative seminorms
uses a different denominator, dominated by an initial endpoint. It cannot
be called an Einstein-equation, curvature or observational relative error.
It can mask errors at points where the target is small. No arbitrary
coordinate derivative scale is imported as a physical residual scale here.

The next concrete physical verification should separately audit the
actual initial Hamiltonian/momentum and shell/Codazzi/subsidiary premises
which enable this cancellation, and choose a declared natural-equation
reference scale or a fully anchored observable with its own error budget.
The present formula does not remove that work.

Projecting this conditionally bounded `E2[W1q]` onto normal constraint rows
would provide another conditional bound for the candidate. It would not be
an independent verification of the actual target's physical Cauchy
constraints: those constraints already enter the formal completion premise.
A noncircular audit must evaluate the physical Cauchy-constraint identities
directly, with the correct trace reversal, Cauchy variables and initial
metric/lapse/shift conversions.

An anchored proper-time shell-volume functional is already specified in
`research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/observable/SHELL_PROPER_TIME_OBSERVABLE.md`.
Its second-order mean uses `J0` directly and an integral of `N0` through
the clock; all first-field clock/quadratic terms must remain. For the same
first field those fixed terms cancel in a candidate-minus-target error,
but a uniform C2 bound alone does not evaluate either history or prove an
observable sign. The current sharper C2 estimate is most directly useful
for the natural-residual bridge above; old first-jet/point tables and clock
histories must not be silently replaced or inherited from a new candidate.
