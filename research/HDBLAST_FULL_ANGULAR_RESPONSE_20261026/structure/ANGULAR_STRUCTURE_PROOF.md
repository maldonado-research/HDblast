# Exact angular blocks, selected TT parity, and the next construction

26 October 2026 UTC. This is a new bounded exact audit of retained data, not an
evolved response computation. Its checker authenticates all seven dependency
closures before any project import and again after the calculation. It imports
the frozen definitions and evaluates new exact algebraic queries; no frozen
producer or reproduction driver is run. The checker uses standard-library
rational arithmetic, an independent small Laurent algebra, and complete
polynomial identities on the sphere, without numerical angular sampling.

## Result

The full registered operator has invariant coordinate blocks of sizes
**5 + 8 + 8 + 1 + 1**. The selected actual bulk, shell and initial data have an
additional exact angular reflection symmetry: their two chiral TT coefficients
agree. Uniqueness therefore places the selected formal response in the
**22-dimensional** reflection-even space with blocks **5 + 8 + 8 + 1**.
Both TT slots should remain in an audit representation, with their difference
serving as a zero sentinel. This reduction does not set the surviving TT
coefficient to zero and does not apply to arbitrary independent second data.

## Exact invariant coordinate blocks

The basis order is the frozen `linear_operator.LABELS` order, zero indexed:

| Block | Indices | Fields | Dirichlet | Robin |
|---|---|---|---:|---:|
| scalar degree 0 | 0,3,6,9,12 | K,L,N,J,w | 1 | 4 |
| scalar degree 2 | 1,4,7,10,13,15,17,19 | K,L,N,J,w,v,u,H | 2 | 6 |
| scalar degree 4 | 2,5,8,11,14,16,18,20 | K,L,N,J,w,v,u,H | 2 | 6 |
| chiral TT plus | 21 | T+ | 0 | 1 |
| chiral TT minus | 22 | T− | 0 | 1 |

Every primitive matrix used by `assemble`, including both boundary projectors
and all four Robin generators, has zero entries between these blocks. This is
an exact support statement for independent formal background coefficients,
not an accidental numerical decoupling at one background point. The union of
the primitive directed support graphs is strongly connected within each listed
block. That establishes coordinate indecomposability; it does not rule out
further decompositions after a nontrivial simultaneous change of variables.
In particular no nonmean block feeds the mean block in this linear second
response problem. The first field enters through the prescribed quadratic
forcing; it does not make the response operator nonlinear in U.

The l=2 and l=4 blocks are scalar-type harmonic multiplets: K,L,N,J,w multiply
the corresponding scalar template; v,u multiply its gradient; H multiplies
its trace-free Hessian. The TT blocks are independent tensor-type harmonics,
not the trace-free scalar Hessian slots H. The two TT tensors obey, exactly,

    tr T± = 0,   div T± = 0,   rough-Δ T± = −6 T±,
    <T±,T±> = 1/6,   <T+,T−> = 0.

They are the lowest tensor-harmonic degree on the unit S3, represented here
by the self-dual and anti-self-dual Weyl parts, the SO(4) types (2,0) and (0,2).
These are tensor amplitudes in the fixed second-order gauge. The audit does
not establish a gauge-invariant radiative observable or a nonlinear wave.

## The two TT equations and two different parities

Write tau for either coefficient, t=ell*x, c=cosh(t), Ht=tanh(t),
h=Ht/rho, a=rho_y/rho, b=−(a_y+a²), d=rho^−2−a². Direct independent
Laurent evaluation of the assembled matrices gives

    tau_xx − tau_rr = A_T tau_x + B_T tau_r + D_T tau + F_T,
    A_T = −3 ell Ht,
    B_T = −3 ell rho a,
    D_T = ell² { −6/c² + rho²[−2a²+2h²−2b−8d+(4/3)V] },
    tau_r(x,0) = G_T(x).

All linear Robin coefficients in these rows vanish. The shell source is
generally nonzero. In the shell adapter's separate notation, with first
normal metric coefficient H=h_s Y, first scalar eta=v_s Y, and tangential
metric `m_ab=2(L_s Y kappa_ab+d_s Hess_ab Y)`, it is

    G_T+ = G_T− = −(8/3) ell rho (a h_s + sigma' v_s/3) d_s.

The shell variables h_s,d_s in this formula are not the bulk h,d above.
The factor comes from `YHESS_T±=2/3`, `QKab_T±=−(4/3)l0*d_s`,
`l0=a*h_s+sigma'*v_s/3`, and the physical forcing factor `2 ell rho`.
Thus the TT boundary is inhomogeneous Neumann; calling it a zero derivative
condition requires the additional fact G_T=0 for the specific source.

Normal reflection across the doubled shell and angular orientation reflection
are distinct operations. Under normal reflection, L=k_y0 and v=k_ya are odd
and satisfy the five Dirichlet rows; every other registered amplitude,
including T±, is even. An even extension can have an inhomogeneous one-sided
normal derivative (a cusp in the doubled coordinate). Normal parity does not
delete TT forcing.

For angular reflection, use the engine coordinates (n0,n1,n2,n3) and
R=diag(1,1,1,−1). Pullback includes a sign for each covector index equal to 3,
as well as the sign of the coordinate monomial. The new exact check finds

    R* T+ = T−,   R* T− = T+,

while every retained scalar, gradient and trace-free scalar Hessian template
is R-even. R fixes the chosen first scalar Y=n0²−n1². Consequently T++T−
is R-even and T+−T− is R-odd. This is not the antipodal parity of S3.

## Data containment and the strengthened selected parity

The inherited selected-space proof uses the multiplicity-free decomposition

    Sym²(1,1) = (0,0) + (1,1) + (2,2) + (2,0) + (0,2)

and multiplicities 5,8,8,1,1 in the relevant section space. Natural quadratic
expressions of the same aligned first state, including the specified initial
completion and physical shell laws, lie in its 23-dimensional selected image.
This is a hypothesis on the registered physical data; stabilizer symmetry
alone or an arbitrary choice of homogeneous second data would not suffice.

The later frozen G capsule supplies the previously missing complete geometric
coefficient maps and all-point polynomial reconstruction certificates. This
audit consumes those authenticated maps rather than rerunning their producer.
Every tensor map has the six field shapes 3,3,3,2,2,7. All eight stored tensor
maps have exactly identical symbolic term lists in their two TT slots:
QRic, QE, QE_matter, QEred, QEred_general, gauge_gradient, PC1 and metric_forcing.
The equality also holds in the general maps, although the actual simplified
QEred and Qphi maps still require C1=0 and its differentiated identity.

The generic registered shell adapter gives the equal TT source just displayed.
The registered initial adapter, tested with independent formal Laurent input
symbols, has zero TT initial values and zero ordinary-t velocities; multiplying
velocities by ell for x=t/ell preserves these zeros. Initial H and v slots
are also zero, but the l=2 and l=4 dynamics can generate them, so their omission
would not be an invariant reduction. Scalar/vector/tensor source projections
and the 23-slot packing preserve the registered basis conventions.

Subtract the two TT equations. Their difference has identical scalar bulk
coefficients, zero forcing, zero initial data and homogeneous Neumann boundary.
The inherited unique linear response theorem gives tau_plus−tau_minus=0.
Equivalently the R-even active space is 5+8+8+1 dimensions, with 5 Dirichlet and
17 Robin rows. With `tau=tau_plus=tau_minus` multiplying `Tplus+Tminus`, its
Gram entry is 1/3. One must not halve the source: each original TT coefficient
equals tau. If instead using coefficient tau_plus+tau_minus, both forcing and
normalization must be changed consistently.

The report exports the sparse 23-by-22 embedding E and 22-by-23 restriction R.
E duplicates the active TT coefficient into both original slots; R averages
the two original TT coefficients. Exact checks give `R E=I` and
`M E=E R M E` for all 25 primitive matrices M, and the report includes the
induced 22-coordinate Gram diagonal. The nonzero formal TT source maps do not
by themselves prove that the specific actual evolved TT forcing is nonzero.

## What is supported physically

The retained P construction already provides a conditional physical subsidiary
completion for the exact total formal response. In its notation,
`C2=DC[k]+QC` satisfies a homogeneous vector wave equation. The chosen physical
initial constraints and harmonic initial completion give zero initial C2 and
zero normal/time derivative. The primal normal gauge row gives the normal
covector component `(C2)_y=0` on the shell. The physical tensor and scalar shell
laws induce the tangential subsidiary Robin condition through Codazzi and the
matter equation. Zero-data uniqueness then gives C2=0 and the natural second
Einstein/scalar equations plus the copy-even physical shell laws.

Thus there is no missing independently prescribed tangential constraint row
to append to the square 23-equation primal system. The conclusion is conditional
on the retained actual smooth background/first solution, complete sources,
initial constraints, corner compatibility and stated regularity argument.
The mean calculation alone does not numerically evaluate this full coefficient.
`subsidiary_evidence.md` records the exact local proof citations and limitations.
Neither the inherited theorem nor the present finite-dimensional audit gives
a finite-amplitude nonlinear family, nonlinear error bound or instability.

## Precise next solvable construction

1. Use the existing first representative and its certified correlated jet
   errors to assemble all nonmean initial data, complete geometric/scalar
   forcing and physical/normal-gauge shell forcing. Retain the same second
   derivative convention, frame units, x-velocity factor ell, and forcing
   signs. Evolve the full 23 slots with equality sentinels, or use the exact
   22-active embedding and Gram entry 1/3 for the surviving TT coefficient.
2. Construct each independent block's total causal right-hand side
   `b=U_free+K F+L G` with the inherited reflected free propagator K and boundary
   lift L. Iterate `U=b+K(A U_x+B U_r+D U)+L(M U|shell)` or evaluate a candidate
   V with a certified integral residual. For the TT block U_free=0 and M=0.
   Both l=2 and l=4 retain their two Dirichlet and six Robin rows.
3. Bound the full value/x-derivative/r-derivative residual in a valid block or
   inherited full Gram norm, including input, interpolation and quadrature
   errors, and apply the retained contraction inverse bound. Existing global
   operator bounds can be inherited on these invariant subspaces; any sharper
   block bound requires its own calculation. No TT slot is dropped before
   this explicit parity argument and consistency check.
4. Report the full angular field and shell data with the certified C1 error.
   A numerical bound on the first-order gauge defect requires physical frame
   and ell conversion factors and QC input errors. Pointwise second-derivative
   Einstein residuals require a separate quantitative C2 certificate or a
   declared weak residual norm; a C1 certificate alone is insufficient.

## Evidence and replay

`ANGULAR_STRUCTURE_NORMAL.json` and `ANGULAR_STRUCTURE_OPTIMIZED.json` are
byte-identical. `STRUCTURE_CONTROLS.json` records rejection of introduced
cross-block coupling, an added TT Robin term, unequal TT bulk forcing, and
the wrong angular reflection. Checks use explicit exceptions and remain active
under Python optimization. All outputs are outside the frozen repository.

Run the new checker only, using a fresh output path:

    python -I -S -B exact_angular_structure.py --repo-root /path/to/repo --output /tmp/new-structure.json

Primary retained sources, relative to `research/`:

- N: `HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/angular/SELECTED_SPACE_CLOSURE.md`,
  `angular/EXACT_ANGULAR_OPERATOR.md`, `angular/angular_exact.py`, and
  `angular/source-audit/source_adapters.py`.
- G: `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json`,
  `projection/EXACT_GEOMETRIC_PROJECTION.md`, `projection/inputs/linear_operator.py`.
- P: `HDBLAST_COMMON_GAUGE_FORCING_20261016/evolution/FULL_SECOND_RESPONSE_CONSTRUCTION.md`,
  especially its sections on regularity and physical subsidiary completion.
- T: `HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md`.
