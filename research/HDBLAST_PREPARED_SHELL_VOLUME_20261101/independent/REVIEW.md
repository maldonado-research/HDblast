# Independent review of the explicit mean and conditioning gauge

1 November 2026. Portable independent mathematical review. **The new five-mean initial
jet, its closed projected evolution, and the proposed second-only gauge
conditioning are consistent, conditional on the retained constraint,
regularity, and evolution theorems. The final independent reconstruction also
accepts the positive proper scalar mean and negative proper volume/expansion
second coefficients at the single frozen target `x=5/4`.**

This portable review uses repository-relative source inputs and requires no Git metadata. The selected source files read
for this review are hashed in `SOURCE_PINS.json`; that is not a fresh authentication
of every dependency closure in the public repository. The review read the
new producer `mean_initial.py` and `MEAN_INITIAL_GAUGE.md` and checked their
formulas independently, without importing them. It made zero oracle calls.

## Initial physical and harmonic jet

Write `A=f+2xi_y`, `B=f+2a xi`, and retain mean order `(K,L,N,J,w)` with
positive Gram diagonal `(1,2,1,3,1)`. All suffix-2 quantities are actual
second parameter derivatives, not coefficients of `epsilon^2`.

Direct inversion of the canonical density map gives

```
<K2_yy>=QM/(3a)+6u ut/rho,
<K2_angular>=-u ut/rho,
<Pi2>=-u j/3.
```

The initial displacement and scaled coordinate velocity are

```
U0=(A,0,-u^2/3,B,0),
V0=ell*(
  2rho QM/(3a)+(32/3)u ut,
  rho[A_y/2-3B_y/2+4aA-3aB-(13/3)u u_y-a u^2/3],
  -2rho QM/(3a)+(4/3)u ut,
  -(4/3)u ut,
  0).
```

The radial displacement derivative is `-ell rho partial_y U0`. There is
no extra proper-time factor in `V0=partial_x U0`.

An independent direct contracted-Christoffel calculation verifies the
time and radial gauge components for general radial second profiles.
The exact lapse prescription yields
`<dot N2>=rho^2[QM/(3a)-u ut/rho]`. Multiplying the shift velocity by the
perturbed spatial metric is essential: the mean `dot g_ty` is not the mean
`dot beta^y`. Omitting this product gives the wrong quadratic coefficients
`-37/3,-97/3` instead of `-13/3,-1/3`. The check explicitly detects this
omission. The lapse-square contribution `N=-u^2/3` is retained. The scalar
coordinate second velocity vanishes pointwise by the density/lapse
cancellation, even though the physical normal scalar velocity does not.

The explicit metric primitive and `D=-QM/(6a)` satisfy the two mean
canonical constraint equations with their full Hessians. Direct integration
by parts also checks the numeric producer's equivalent primitive:

```
xi=(3/2)(f_y-a f)/p^2-[KH-(I/J)Kb]/(rho^4p^2),
K_z=-rho^5 QH+(I/J)rho^3 b.
```

Positive `ell*Delta r` suffix sums from the target to `r=2` have the correct
orientation for increasing outward y. Interval dependency loss in the
zero-total-integral cancellation can widen an enclosure; it does not change
the underlying smooth compactly supported profile.

The exact initial harmonic prescription imposes `C=0`; the physical
constraints and reduced accelerations give its initial normal derivative
zero. Those accelerations are consequences of the reduced equation, not
extra data. This implication uses the trace-reversed normal constraint
projections and the retained analytic premises. Finite algebra checks alone
do not prove it.

## Complete mean projection and source independence

The Haar projector must act on natural tensor sections. Its invariant
metric range consists of the four slots `yy,y0,00,angular trace`; the scalar
supplies the fifth. There is no invariant angular one-form or trace-free
angular tensor. Radial background bulk, gauge, and physical boundary
operators commute with this projector. At second parameter order all
unknown second fields occur linearly. Thus every nonmean second component
decouples from the mean, while all quadratic first-field products must be
projected into the mean forcing. The same first tangent gives exactly the
same reduced and physical shell Hessian sources.

`algebra_checks.py` separately constructs the connection and curvature
actions on full five-dimensional symmetric tensors, then checks preservation
of the invariant subspace. `INDEPENDENT_MEAN_MATRICES.json` records the
resulting complete `A,B,D,M` matrices. The scalar derivative coupling is
retained; explicitly, its underlying contraction is

```
k:Hess(phi0)+p dC[k]_y
 = (p_y+4ap)K-3ph L
   +p[(K_y+N_y-3J_y)/2-L_t/rho].
```

The Robin frame rows are

```
K_y=-4aK-2aN+6aJ+(4sigma'/3)w,
N_y=-aK-sigma'w/3, J_y=aK+sigma'w/3,
w_y=pK/2-sigma''w/2,
```

with their retained forcing, and `L=0` is the single Dirichlet row. There is
no coordinate-frame `2a` tangential self-coupling after frame conversion.

The exact family's projected second derivative and the new five-mean
Volterra problem share these operators, sources, and initial jet. Retained
linear uniqueness therefore identifies them. No angular master inversion,
especially no exceptional `ell=0` substitution into an `ell=2` formula, is
used. The mean has the initial support `[1,2]`; its zero wedge follows from
that support and the unchanged first-field source support even if other
second angular data and higher nonlinear orders have larger support.

## Second-only gauge conditioning

Let the second pullback generator be
`Z=X partial_y+b partial_t`, with initial data

```
X=-xi, b=0, X_t=0, b_t=a xi.
```

Independently contracting the vector wave operator plus Ricci gives the
full monopole gauge equations, where `Ht=tanh(t)`:

```
X_tt=rho^2(X_yy+4aX_y)-3Ht X_t
     -2a rho^2(b_t+3Ht b)-rho^2(4a_y+8a^2)X,
b_tt=rho^2(b_yy+6ab_y)-3Ht b_t+(3+3Ht^2)b-2aX_t.
```

They imply initially

```
X_tt=rho^2[-xi_yy-4a xi_y+(4a_y+6a^2)xi], b_tt=0.
```

The shell conditions are `X=0` and coordinate `b_y=0`. They preserve the
fixed-shell representation and the vanishing mixed metric row. The initial
generator and all its jets vanish on the open shell collar. These are the
monopole vector equations themselves; no absent angular-gradient amplitude
is introduced by a limiting procedure.

For `k_new=k+Lie_Z g0`, `w_new=w+Z(phi0)`, the new initial pair is

```
U0_new=(f,0,-u^2/3,f,-p xi),
V0_new,L=ell rho[-f_y+a f-(13/3)u u_y-a u^2/3],
V0_new,K/N/J=V0_K/N/J, V0_new,w=0.
```

In particular `-p xi` is mandatory. Removing the radial Lie metric term
while retaining zero scalar data would change the physical family.
The lapse second value remains zero: the radial variation of `g_tt` and
the term `2g_tt b_t` cancel. The chosen vector wave keeps `dC[Lie_Zg0]=0`,
so the same homogeneous reduced operator applies. At second order its
addition leaves every quadratic source unchanged. It is a second-only
change of representation of the same physical jets; a smooth diffeomorphism
realization preserves the exact physical family and anchors. Exact
finite-parameter relative-harmonic gauge beyond second order is a separate
coordinate completion and is not needed to identify this second response.

At the shell, `X=0`, so the second scalar trace itself is unchanged. The
second tangential time reparametrization changes the total metric means by
`delta N=-2b_t`, `delta J=2Ht b`; its graph-clock change is
`delta Bclock=-rho_b b`. Hence
`(3/2)delta J+3H delta Bclock=0`: the anchored volume coefficient is also
unchanged. Conditioning can improve a norm bound, but its actual numerical
benefit still needs evaluation.

## Proper scalar mean and volume clock

Use the retained harmonic shell first amplitudes `X,Z,V,d`,
`Cclock'=X`, and `E=Z-H Cclock`, with prime denoting proper time. The
first scalar is `psi Y`; the background shell scalar is constant. The graph
clock displacement is `-Cclock Y`, and the first volume-density coefficient
is `(3E-8d)Y`. Expanding the normalized geometric average gives

```
M2=w0-Cclock psi'/3+(E-8d/3)psi.
```

The final measure term must not be omitted. A second measure term cancels
between numerator and denominator against the constant background scalar.
All these first-order corrections are unchanged by the second-only gauge.

The complete volume formula in the retained proper-time specification is
consistent with the determinant Hessian and graph clock. In particular
its direct second metric term is `(3/2)J0`, not `J0/2`, and its clock needs
the full history of `N0`. With the current intrinsic trace identity it can
also be recovered from

```
v2''+4H v2'-2v2/A^2
 =(F'_b/2)M2+F''_b psi^2/12-psi'^2/6+4psi^2/(3A^2)
   -5E^2/(3A^2)+Q_D,2/2.
```

The proposed equivalent kinematic term is correct:
`Q_D,2/2=E'^2-(16/3)E'D'+(8/3)D'^2` when `D'=-P/A^2`.
All history integrals must lie inside the proved shell domain. In particular
the older pilot's shell history ends at `x=5/4`; it does not by itself cover
the entire newer target interval ending at `63/50`.

## Error-budget implications and limits

The old homogeneous coefficient bounds, forcing bounds, and contraction
identities can be reused on a contained domain because the operator and first
sources are unchanged. The old initial error and total `U_cap` cannot be
reused for this family or its new gauge. On the cone
`0<=x<=T=5/4`, `0<=r<=T-x`, the zero wedge gives `|U|<=delta X(U)`
with `delta=1/4`; the same holds at the shell. The resulting bound is

```
q=T(2*10^-4+delta*10^-6)+delta/2000
 =1201/3200000 < 1.
```

A two-term candidate still requires a fresh complete ledger
`q^2 U_cap+(1+q)e0+e1`, including any omitted forcing exactly once and
all matrix/profile approximation errors. A full residual uses instead
`residual/(1-q)`. A loose upper bound does not prove a lower bound on actual
error, and visible cancellations do not supply a physical response sign.

The checks above support exact algebra, projection, gauge equivalence,
and observable normalization. The numerical follow-on below uses these
identities and retained continuous input boxes. Neither part supplies a
finite-amplitude radius, radiation or heating, an empirical signal, or a
cosmological-origin conclusion. This is internal AI-assisted review, not
external specialist review.

## Frozen numerical response: independent reconstruction accepted

The final campaign was frozen before execution. Current portable initial
and response reference hashes and byte counts are in `SOURCE_PINS.json`.
Their scientific fields match the previously reviewed frozen results.

`reconstruct_response.py` imports no scientific producer and makes no
oracle calls. It reads the retained original raw initial profiles and
independently derives the conditioned displacement, radial slope, and
velocity in exact rational arithmetic. It verifies 14,910 component
inclusions against the separate conditioned boxes, all partition joins,
the frozen target footprint, and the pinned necessary public source files.
The actual first shell metric/scalar/background raw oracle boxes retained
in the final report are reused as input; they are not reevaluated.

The independent code reconstructs the fresh Gram initial caps, full
free-wave cap, source cap, and full resolvent error; the complete proper
scalar clock/measure correction; all first quadratic volume terms; the
unweighted and linearly weighted scalar histories; both Green allowances;
and the final scalar, volume, and expansion intervals. Twelve reported
quantities match exactly, including all three final interval endpoints.
Normal and optimized independent reconstructions agree byte for byte in
`references/RESPONSE_NORMAL.json` and `references/RESPONSE_OPTIMIZED.json`.

The new free-wave bound is below `1.748e31`, and

```
X(U-U_free) <= (q X(U_free)+T FH+GH)/(1-q) < 6.563e27.
```

This includes every feedback iteration and all original bulk and boundary
forcing once. It is not an old-family cap or an error estimate for a
finite coupled iterate. At the shell the scalar free trace equals the new
initial scalar displacement evaluated at `r=x`, since its velocity is
zero and its reflection is even. The varying wedge estimate then gives
point error `(x-1)E_X`. The weighted and unweighted history allowances
are exactly `delta^3 E_X/6+delta^2 Q_M/2` and
`delta^2 E_X/2+delta Q_M`; no point-to-history extrapolation is used.

The exact Green kernel and its target derivative are independently checked
against the volume ODE, including the diagonal value zero and derivative
one. Their `10^-11` approximation allowances are verified with elementary
rational bounds. With `rmax=ell T=1/800000`,
`umax=ell delta=1/4000000`, and `Cmax=1/(1-rmax^2/2)`,

```
1/Cmax^2 <= k/(T-x) <= 1,
(1/Cmax^2)[1-2rmax umax/(1-2umax^2)]
  <= k_T <= 1/(1-2umax^2).
```

Each endpoint differs from one by less than `10^-11`. The exact upper
bound `k<=T-x` justifies the first-source volume allowance without any
missing kernel multiplier. The expansion allowance includes the derivative
kernel multiplier. The conversion to proper time is `(rho_b ell)^2`
for volume and `rho_b ell` for expansion.

Deliberately outward rounded consequences at the frozen target are

| Second parameter coefficient | Enclosure |
|---|---|
| Proper-volume scalar mean `M2` | `[1.1296e30, 1.9780e30]` |
| Relative volume `v2=V2/V0` | `[-5.4520e16, -2.6521e16]` |
| Mean expansion `Theta2=v2_s` | `[-9.021e21, -4.445e21]` |

No material error was found in the final formulas or error budget. These
are conditional second changes relative to the background on the same
anchored proper-time slice. The latter two signs do not mean negative
volume, negative total expansion, a contracting solution, or negative
volume acceleration. No time-interval sign beyond the single frozen
target is inferred. The result belongs to the newly specified explicit
preparation, and does not recover the old unspecified weighted-only family.

The conditioning argument only needs the second jet of the pullback of the original exact family. The smooth second-only flow is tangent to the shell and fixes the initial anchor. It preserves the physical observables; its second jet satisfies the linear relative-harmonic equations because the generator solves the homogeneous vector wave. No finite-parameter nonlinear harmonic property is asserted. Retained lower C4 regularity leaves the spatial differentiation buffer needed for second pullbacks into lower C2/C1, sufficient for the stated volume and clock functionals.


## Additional independent initial quadrature audit

The portable `reaggregate_initial.py` adds an independent standard-library
Fraction calculation. It does not accept the producer global I/J values as
unverified inputs. It re-encloses each bump value with a degree-40 alternating
exponential series, exact rational argument reduction, and outward squaring
on a finer 2^160 grid. The flat tails use the monotone bound with
`s>=4096/63`; both tail placements and the full partition are checked.

Using the retained raw rho, QH and derivative boxes as shared validated
premises, it independently reconstructs all 994 I/J quadrature contributions,
both global integrals, the normalization and metric moment, every suffix
primitive and xi interval, and both frozen point primitives. Every independent
unrounded result must lie inside its retained outward box. The separate sums
of the 994 recorded contributions must equal the recorded global integrals
exactly. Counts are 1,988 cell integral inclusions, 994 cell primitive
inclusions, two point primitive inclusions, and two exact global sum equalities.
Background, QH and f_y interval generation is not independently rerun here;
no actual oracle is called by any independent checker.

Five controls alter formulas while keeping the authenticated final witness:
omitted proper measure, halved scalar, omitted history error, wrong moment
factor and reversed primitive sign. Each must exit with the exact structured
scientific rejection reason and leave no output. A separate changed-report
control checks custody rejection. Positive input authentication and scientific
rejection counts are reported separately. Normal and optimized runs produce
identical canonical bytes for every check and control report.
