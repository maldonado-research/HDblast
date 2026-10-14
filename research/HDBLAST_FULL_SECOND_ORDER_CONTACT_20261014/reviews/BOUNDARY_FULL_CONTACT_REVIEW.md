# Independent boundary review of the complete finite-angular coefficient

14 October 2026. Internal AI-assisted analytic review, not external peer review.

Reviewed `FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`, SHA256
`69833a130d3e234c3ed304e1aa04671d869e5679254df4ee4aae1b306691dd8f`.

This binds to a fresh snapshot read after the root added its explicit reduced
quadratic source recipe, finite-angular initial gauge choices, and inherited
normal-force consequence. Those additions are included in this acceptance.

**Disposition: accepted under the explicitly stated model, input, finite-mode,
and controlled-collar hypotheses.** I find no remaining boundary/gauge or
constraint-source gap in this version's second-coefficient construction.
This acceptance does not assert an exact nonlinear family through contact,
an estimate uniform in all angular modes, or a numerically evaluated spacetime.

## Exact shell and gauge sources

I independently derived the full Eulerian scalar shell rows, their second
Regge–Wheeler reduction, and the relative wave-gauge exact laws in the
companion boundary notes. The root note's equations (8)–(10) match those
independent formulas, including the scalar source `−2(chi^w)^2`, the mixed
tension/lapse source, the full inverse-metric trace term, and the distinction
`k_yy=2 nu_2+2 nu_1^2`. The actual normal law follows from the relative
connection trace, rather than borrowing the incompatible absolute harmonic
condition of the proper-normal background.

The normal metric relation eliminates the physical tangential normal
derivatives without introducing unknown tangential derivatives. Thus its
claimed principal boundary structure is four mixed-component Dirichlet
waves and twelve normal-derivative waves. My subagent's independent full
and scalar stable-symbol calculation verifies the rank completion; simply
adding every harmonic constraint to Israel would instead leave four full
or two scalar tangential gauge modes.

All sources must be computed from the complete transformed first wave
fields. The note does this and supplies the correct inverse transformation
in (6). In particular, it does not retain the original purely conformal
tangential first field or homogeneous tensor source after a nontrivial first
shell tangential reparametrization.

## First gauge and initial compatibility

For covariant generator components, `xi_y=0` and
`xi_a,y−2 a_b xi_a=0` exactly preserve the boundary mixed first metric
components. The coefficient two is correct. The bulk pure-gauge harmonic
image is `(Box+Ric)xi`, so the forced generator in (5) sets the full first
wave-map vector to zero.

Zero initial generator and its time derivative make the pullback identity
on the initial surface with identity first jet. They preserve the physical
first canonical metric, extrinsic curvature, scalar and scalar momentum;
the known gauge wave equation adjusts lapse/shift velocities. The quadratic
first-generator corrections also have zero initial canonical contributions.
This is consistent with projecting the already established smooth localized
SECOND physical constraint correction, rather than differentiating the
separate exact fixed-finite-Sobolev initial family and assuming it is C-infinity.

The isotypic projection onto (11) commutes with the radial-background linear
constraint operator and preserves radial support. It acts on the complete
canonical tensor/density fields. The representation list is exact:
`Sym^2(1,1)=(0,0)+(1,1)+(2,2)+(2,0)+(0,2)`; the full SU(2) character identity
is independently checked, not only the dimension sum. This closes the
finite-angular initial correction issue. Removing the optional free initial
TT directions relies on the explicitly cited prior TT constraint-kernel
result; it does not remove any required scalar source.

For the second initial gauge conditions, the time-derivative coefficients
in (12) are correct. These five equations solve for the lapse and shift
velocities without changing the physical initial data. The constraint tensor
is the trace reversal of the trace-reversed bulk residual. In an initial
orthonormal Cauchy frame, (13) follows with the displayed positive one-half
factors. This gives the second initial subsidiary time derivative zero,
rather than assuming it independently. Since the physical corrections and
first fields vanish on an open shell neighborhood, the complete reduced
recursion gives the required corner hierarchy.

## Subsidiary constraint boundary system

The variable-tension cancellation in (18) is exact in intrinsic dimension
four: `−3 d(sigma/6)=(−sigma'/2)dphi=p dphi`. Both physical metric and
scalar junctions are required. With `C_n=0`, the reduced mixed Einstein row
gives `(nabla_n C)_a−K_a^b C_b=0`, which for background covector coordinate
components is `C_a,y−2 a_b C_a=0`. Literal coordinate Neumann would miss
this lower-order term.

The five-dimensional identity (17) has the correct potential coefficient:
its divergence is `−(Box phi−U') dphi`. Hence the subsidiary equation is
homogeneous when the scalar equation holds. At second order there is no
free constraint forcing once `C_0=C_1=0` and all reduced Taylor sources have
been retained. The mixed normal-Dirichlet/tangential-Robin system falls under
the same finite-mode uniqueness argument, so the claimed physical constraint
propagation follows. No independent normal-force condition is needed for
identical copies; the full physical equations imply its normal projection.

As an additional source/sign check, the independently derived longitudinal
mixed bulk residual at the shell is
`J_y=6u u_y−24a u^2−v v_y=3 S_K`, using the first shell rows. The general
RW linear mixed row evaluates to `−3 S_K` using the independently derived
second shell rows. Thus `DE_y=−J_y` holds at the shell with exactly the
normalization used in the source notes.

## Finite-mode characteristic construction and limit

I checked the signs in (15)–(16), the outgoing/incoming characteristic
directions, and the boundary trace ODE. Solving this finite-dimensional ODE
before the fixed-point map is essential: it converts dependence on the
outgoing trace into a small-time Volterra term. The unknown boundary trace
cannot instead be treated as arbitrary incoming data. The resulting map
has a small-step contraction for bounded lower-order matrices, including
nonsymmetric Robin coupling. The `U_z` integrability defect has zero time
derivative and vanishes initially. Smooth compatible data and the corner
recursion permit the differentiated construction and finite-interval
continuation. The claimed finite propagation and use of background extension
only outside causal support are consistent with this proof.

The all-angular unweighted Neumann stable determinant tends to zero at
glancing; the note correctly makes no uniform claim there. That limitation
does not invalidate a fixed finite angular matrix construction. The domain
and contact time are explicitly conditional on an existing first contact
solution and a controlled positive-cone background. The note does not use
this result to assert a nonlinear contact family or global extension across
the cone.

## Verification evidence

My standalone boundary checker passes 39 exact identities and 10 wrong-formula
controls per normal/optimized mode, with byte-identical scientific receipts.
It includes exact differentiation of all three relative-GH boundary sources.
Its actual extrinsic-curvature corruption exits one in both modes. The
independent principal checker passes 31 identities and four controls per mode,
with two actual trace-coefficient corruption rejections. These computations
verify the specified algebra; the analytic characteristic/constraint proof
was reviewed above and is not certified merely by check counts.

## Final frozen version readback

The final root manuscript has SHA256
`59c58b560b257a7311d4d00a72a1d4facb95c23377b88cf0323a4d57edc81949`.
I verified this hash and read back section 5's explicit requirement that
the selected second lapse/shift values vanish near the initial shell corner,
with zero values admissible. This clarifies the compatible choices used in
the accepted construction. The boundary acceptance extends to this final
version without changing the scientific checker sources or rerunning them.

The standalone boundary notes do not independently assert the complete
coefficient theorem. Their analytic proof obligations are met for this
selected finite-angular family by the combined root manuscript reviewed here.
