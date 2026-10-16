# Common-gauge second boundary completion: independent audit

16 October 2026. Separate working continuation, restricted to
`common-gauge-forcing-20261016/boundary-review/`. All frozen inputs are
read only. This is conditional coefficient mathematics. It supplies exact
boundary formulas and compatibility identities, not an evaluated second
evolution, a full nonlinear IBVP theorem, or an exact solution family.

## 1. Registered signs and the actual first representative

The authenticated retained contact and material-junction notes agree:
increasing proper radial `y` is **outward on each regular interior**.
Set `a=rho_y/rho`, `p=phi0_y`, and write `gamma` for the background
induced Lorentzian shell metric. Then

\[
\sigma_0=6a_b,\qquad \sigma'_0=-2p_b,\qquad
\sigma''_0=4\phi_b,\qquad\sigma'''_0=4.
\tag{1}
\]

The outward two-copy laws are
`sum K_out - sigma gamma/3 = 0` and
`sum n_out(phi) + sigma' = 0`. The inward convention
`n_in=-partial_y` therefore has precisely the equivalent signs
`sum K_in + sigma gamma/3 = 0` and
`sum n_in(phi) - sigma' = 0`. Reversing the normal does not change (1),
which defines `a,p` using the outward coordinate derivative. There is
one shell tension, not a tension separately added for each copy.

Use the actual pair `q=h_wave`, `eta=chi_wave`, not the old longitudinal
quadratic source. The certified generator has
`f|b=0`, `b_y|b=d_y|b=0`. Tangential differentiation of the first row
gives `f_t|b=0`; hence `q_ya|b=0` for all four shell-tangent indices.
It fixes the shell as a set and permits clock/angular transport along it.
In particular, the tangential tensor

\[
H=q_{yy}|_b,\qquad m_{ab}=q_{ab}|_b,\qquad
\eta=\chi_{\rm wave}|_b=vY,\qquad
A^a{}_b=\gamma^{ac}m_{cb},\quad\tau=\operatorname{tr}A
\tag{2}
\]

is generally anisotropic. `tr(A^2)` includes the time-angular and
anisotropic angular components with the Lorentzian index signs. It is
not a positive tensor norm and must not be replaced by an isotropic trace.

All sources below are second parameter derivatives. The coefficient
convention is `g=g0+epsilon q+epsilon^2 k/2`,
`phi=phi0+epsilon eta+epsilon^2 w/2`.

## 2. Physical affine Hessians, with normal normalization

For the copy-even fixed shell, define outward one-copy residuals

\[
R_{ab}=K_{ab}-\sigma(\phi)g_{ab}/6,\qquad
P=n_{\rm out}\phi+\sigma'(\phi)/2.
\]

At the shell `g_ya=0`, so `N=sqrt(g_yy)`,
`K_ab=g_ab,y/(2N)`, `n_out phi=phi_y/N`. For the affine first-pair
path, `N'=H/2`, `N''=-H^2/4`, and `(1/N)''=3H^2/4`.
This yields the direct Hessians

\[
Q_R=\frac{3aH^2}{4}\gamma-\frac H2m_y
 -\frac{\sigma''\eta^2\gamma+2\sigma'\eta m}{6},
\]
\[
Q_P=\frac{3pH^2}{4}-H\eta_y+\frac{\sigma'''}2\eta^2.
\tag{3}
\]

The actual first physical rows imply

\[
m_y=2am+(aH+\sigma'\eta/3)\gamma,\qquad
\eta_y=pH/2-\sigma''\eta/2.
\tag{4}
\]

Thus their values can be computed without differentiating a first
boundary identity:

\[
\boxed{Q_R=(aH^2/4-\sigma'H\eta/6-\sigma''\eta^2/6)\gamma
                  -(aH+\sigma'\eta/3)m,}
\]
\[
\boxed{Q_P=pH^2/4+\sigma''H\eta/2+\sigma'''\eta^2/2.}
\tag{5}
\]

The two-copy outward residual Hessians are `2Q_R,2Q_P`; the two-copy
inward Hessians are their negatives. The equations for a second
coefficient always have forcing equal to minus the chosen residual
Hessian. In particular, the one-copy physical equations are
`DR[k,w]=-Q_R`, `DP[k,w]=-Q_P`. No source has been divided or multiplied
by two solely because it is called a Taylor coefficient.

## 3. Actual normal harmonic row and all sixteen primal rows

The registered reduction uses the relative wave-map vector

\[
C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\Gamma[g_0]^A{}_{BC}),\qquad
E^{\rm red}_{AB}=E_{AB}-\nabla_{(A}C_{B)},
\]

with symmetrization including `1/2`. The correct second subsidiary
field is `C2=DC[k]+Q_C`, with `Q_C=D^2C[q,q]`. Since the actual
first field has `DC[q]=0`, the first lowering cross term vanishes and
the lower and upper normal Hessian components agree at the background.

For arbitrary boundary-orthogonal metrics the direct normal expression is

\[
C^y=\frac{(g_{yy})_y}{2g_{yy}^2}
 -\frac{\operatorname{tr}_{g_{ab}}(g_{ab,y})}{2g_{yy}}
 +a\operatorname{tr}_{g_{ab}}\gamma.
\tag{6}
\]

The last term contracts the fixed reference connection with the
**perturbed** tangential inverse metric. It cannot be discarded.
Its first row, using (4), gives

\[
H_y=-4aH+2a\tau+\frac{4\sigma'}3\eta.
\tag{7}
\]

Twice differentiating (6) directly on the affine path and using (4),(7)
gives the useful actual shell gauge value

\[
(Q_C)_y=4aH^2-3aH\tau+2a\operatorname{tr}(A^2)
 -\frac{4\sigma'}3H\eta+\frac{\sigma'}3\eta\tau.
\tag{8}
\]

Write `K=k_yy`, `T=tr_gamma k_ab`. The full primal completion is:

\[
\boxed{k_{ya}|_b=0\quad(4\hbox{ Dirichlet rows}),}
\]
\[
\boxed{k_{ab,y}-2ak_{ab}-(aK+\sigma'w/3)\gamma_{ab}
            =F_{ab}:=-2Q_{R,ab}\quad(10\hbox{ rows}),}
\]
\[
\boxed{w_y+\sigma''w/2+\sigma'K/4=F_\phi:=-Q_P,}
\]
\[
\boxed{K_y+4aK-2aT-4\sigma'w/3=G_N,}
\tag{9}
\]

where

\[
\boxed{G_N=\frac{4\sigma''}3\eta^2+4\sigma'H\eta-10aH^2
                 +8aH\tau-4a\operatorname{tr}(A^2).}
\tag{10}
\]

An independent useful join is

\[
G_N=\operatorname{tr}_\gamma F-2(Q_C)_y
             =-2\operatorname{tr}_\gamma Q_R-2(Q_C)_y.
\tag{11}
\]

Indeed after inserting the physical metric row,
`(C2)_y=(L_N[k,w]-tr_gamma F)/2+(Q_C)_y`, where `L_N`
is the left side of the last row in (9). Thus (10) is exactly `C2_y=0`,
not an additional physical source or an arbitrary incoming constraint.

The same result follows from the retained lapse rows with
`N=1+epsilon n1+epsilon^2 n2/2`. The conversion is
`n1=H/2`, **`K=2n2+2n1^2`**. Omitting the quadratic term changes
both physical and normal-gauge forcing. Equation (10) uses (7) to
eliminate `n1_y`; that elimination uses values of actual first rows
and needs no new source derivative estimate.

For the actual first harmonic the contraction in (10) can also be formed
without a componentwise isotropic bound. Define
`U0=u+b_t`, `L=u+tanh(t)b`, `J=cosh(t)d_t-b/cosh(t)`,
and `h=-4u+2f_y`. Then

\[
H=hY,\quad \eta=vY,\quad
\tau=(8u+2b_t+6\tanh(t)b-16d)Y,
\]
\[
\operatorname{tr}(A^2)=(4U0^2+12L^2-64Ld)Y^2
                  +4d^2|\operatorname{Hess}Y|^2-2J^2|DY|^2.
\tag{10a}
\]

With `S=n1^2+n2^2`, the ambient quadratic identities give
`|DY|^2=4(S-Y^2)<=4` and
`|Hess Y|^2=8-8S+24Y^2<=24`. In particular the time-angular term in
(10a) has a minus sign. This shows directly that the scalar trace
sources close in the span of `1,S,Y^2`; it gives a simple safe bound
`|tr(A^2)| <= 4|U0|cap^2+12|L|cap^2+64|L|cap|d|cap
+96|d|cap^2+8|J|cap^2` without assigning it a positivity interpretation.

In the scaled inward coordinate, `partial_y=-(ell rho)^(-1)partial_r`.
All twelve non-Dirichlet rows therefore become a finite matrix Robin
law with known inhomogeneous trace data. If tangent tensor components
are expressed in the radially parallel background orthonormal frame,
`partial_y(k_ab)_frame` corresponds to
`(k_ab,y-2ak_ab)_frame`; one must include this warp conversion.

## 4. Natural tangential and subsidiary compatibility

There is a useful exact off-shell one-copy identity, independent of
the wave reduction:

\[
E_{na}=D_bR^b{}_a-D_a\operatorname{tr}R-PD_a\phi.
\tag{12}
\]

It follows by substituting `K=sigma g/6+R` and
`n phi=-sigma'/2+P` into Codazzi. It includes the nonconstant-tension
matter cancellation. At the background `D_a phi0=0`; for the actual
first solution `R1=P1=E1=0`. Therefore the quadratic *natural* bulk
source obeys the boundary compatibility relation

\[
\boxed{Q^E_{ya}=D^bQ_{R,ba}-D_a\operatorname{tr}_\gamma Q_R.}
\tag{13}
\]

The background first residuals vanish, so variations of the normal,
raising, boundary connection and pullback multiply zero first or
background residuals in this identity. Those cross terms have not
been omitted for a nonzero residual. The scalar physical row is
essential to the exact cancellation, even though no separate `Q_P`
term survives in (13) because `D_a phi0=0`.

If (9) and the second reduced bulk equation hold, (12) gives `E2_na=0`.
Using the reduced equation, the normal primal row `C2_y=0`, and its
zero tangential derivatives gives the four induced boundary conditions

\[
\boxed{\partial_y(C2)_a-2a(C2)_a=0.}
\tag{14}
\]

These are Robin conditions on *covector coordinate components*, not
literal coordinate Neumann conditions. They are not four additional
primal boundary rows. Matching is automatic for identical copies; the
normal-force equation and tangential shell conservation follow from the
physical equations and matching, as in the retained exact identities.

The bulk trace-reversed Bianchi identity is

\[
B_0(Q^E):=\nabla_0^AQ^E_{AB}
 -\tfrac12\nabla_B^0\operatorname{tr}_0Q^E=-Q^\phi\,d\phi_0.
\]

For `Qred=QE-nabla_(A(Q_C)_B)`, it becomes

\[
B_0(Q^{\rm red})=-Q^\phi d\phi_0
 -\tfrac12(\Box_0+\operatorname{Ric}_0)Q_C.
\tag{15}
\]

Thus `C2=DC[k]+Q_C` has the homogeneous subsidiary wave equation,
normal Dirichlet row, and tangential Robin rows. This is a compatibility
identity, not a new proof of the needed subsidiary uniqueness or a bound
for a second evolution. The certified first gauge C3 gives `q` through
C2, `Q_C` through C1, and direct bulk sources through C0. The divergence
and `Box Q_C` in (15) can be interpreted distributionally; no new C4
gauge certificate or source C1 bound is asserted by differentiating (15).
The tangential boundary identity (13) only needs available lower jets.

## 5. Initial data and corner audit

The published quantitative initial-data checkpoint already supplies the
actual physical second Cauchy jet and a complete relative-harmonic
second spacetime initial one-jet. Its first correction has precisely
the initial gauge-wave acceleration

`f_tt=-6rho^2(u_y+4au)`, `b_tt=d_tt=0`, with
`xi=xi_t=0` initially. The current actual first gauge wave therefore
has exactly the same first spacetime initial one-jet used in that
completion, including the changed mixed metric velocity.

Consequently the initial data need not be replaced by an uncomputed
zero-scalar prescription or reconstructed from an old longitudinal
quadratic source. They can be adopted directly as the chosen harmonic
second data. Their scalar and normal scalar-velocity corrections are
nonzero. Their `g2_tt=-2rho^2 U^2` and mixed metric-velocity conversion
are retained. If instead a different second gauge convention is chosen,
the full second transformation, including `2 L_xi h1+L_xi^2 g0` and a
possible second generator, must be matched to these data.

The inherited completion establishes `C2|t=0=0` exactly. Its physical
Hamiltonian and four momentum constraints, together with the reduced
acceleration equations, imply `partial_t C2|t=0=0`. The accelerations
are determined by the reduced equations, not extra prescribed Cauchy data.

The master packet and both physical second data and initial gauge
completion are supported in `ell<=z<=2ell` with flat endpoint zero
extensions. The first gauge wave has zero initial vector and velocity.
The shell corner therefore has an open zero collar. Finite propagation
keeps the first fields and these quadratic shell sources zero before
the first contact front. This supplies the abstract corner identities
for (9), including all orders for which the inherited smooth solution
and differentiated recursion are justified. A new quantitative high-jet
bound does not follow merely from flatness.

The concrete remaining initial work is implementation: authenticate and
assemble these already computed data into the chosen finite-component
U2 solver, apply its actual coordinate/frame/modal conversions, and
check the initial and corner defects of its reconstructed solution.
Initial harmonic one-jets are available; a validated second evolution
and its higher target regularity have not thereby been computed.

## 6. Minimal closed problem and remaining computation

At the algebraic level the boundary specification is closed: the fifteen
metric components and scalar have four homogeneous mixed Dirichlet rows
and twelve Robin rows (9). No physical, normal-gauge, or independent
incoming tangential constraint datum is missing. The natural bulk source,
registered reduced source, initial one-jet, physical trace sources,
normal-gauge source, and induced subsidiary rows must all use this same
first representative.

For the selected single quadrupole the full equivariant symmetric square is

`(0,0) + (1,1) + (2,2) + (2,0) + (0,2)`.

It includes scalar degrees 0,2,4 and both tensor degree-two orientations.
No independent forced coexact-vector representation is generated by this
single state. This is a statement about full tensor/covector group
representations, not a scalar componentwise Laplacian cutoff; the solver
must still construct every necessary section multiplicity and coupling.

For any proposed second Volterra construction, the remaining concrete
steps are to assemble the actual finite bulk/boundary matrices and initial
profiles, prove a bound for that coupled inverse in a stated norm, and
validate a continuous reconstruction with all bulk, boundary, gauge,
initial, corner and rounding defects. A source C0 cap alone does not
certify a classical C2 solution, high Sobolev regularity, or a shell
time/normal derivative observation bound. Such a bound must match the
desired second response functional. These are requirements on the
construction, not extra free shell rows.

## 7. Independent exact checker and provenance

`check_boundary_completion.py` uses only the Python standard library.
It independently expands the exact normalized physical laws and direct
relative normal gauge law on three anisotropic rational Lorentzian
boundary jets, including arbitrary first intrinsic metric/scalar jets.
It verifies all ten metric components, scalar and normal gauge Hessians,
the metric/lapse conversion, the trace join (11), and a direct Codazzi
contraction giving (13). These are finite exact witnesses supplementing
the displayed all-point derivations, not numerical evaluations of the
actual master/gauge fields.

The retained normal and optimized reports each pass 192 comparisons and
six named controls: wrong scalar sign, missing factorial, isotropic source,
wrong scalar third-derivative term, omitted full matrix square, and omitted
metric/lapse quadratic conversion. The reports are byte identical.
Twelve separate fresh intended-error subprocesses also fail for their
specified reasons in normal and optimized modes, leaving no report.

Primary local references read for this audit:

- `gauge-wave-solution-20261016/COMPLETED_HANDOFF.md`.
- `gauge-wave-solution-20261016/future_work/first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md`.
- `scalar-contact-20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`.
- `scalar-contact-20261014/boundary/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md`.
- `quadratic-shell-forcing-20261015/assembled-source/inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md`.
- `second-initial-norms-20261015/gauge/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`.
- `second-initial-norms-20261015/COMPLETED_HANDOFF.md`.

`BOUNDARY_AUDIT_RECEIPT.json` pins these references and records the fresh
rejection processes. Hash authentication is not a fresh replay of every
inherited PDE certificate.

## 8. Independent quantitative audit of the actual shell majorants

`check_actual_shell_caps.py` is a separate independent implementation. It
reads the fixed retained master, gauge, coefficient and coordinate-shell
reports; authenticates their hashes; reconstructs the shell field caps
and all ten physical source component bounds with exact fractions and
**no intermediate rounding**; and compares them to the shell producer's
outward rounded bounds. It imports no producer code. Its 168 comparisons
include source values, proper clock derivatives and induced-shell angular
and tangential covariant derivatives, including the de Sitter boost terms.
For the mixed `0a` angular derivative, metric compatibility cancels the
large isotropic `A gamma` term before bounding.

The actual quantitative shell trace inputs give the convenient strict caps

| Quantity | Absolute component cap |
| --- | ---: |
| `H` | 0.034137 |
| `eta` | 46658 |
| `tau` | 0.0000030541 |
| `tr(A^2)` | 0.0000000000089232 |
| one-copy physical `Q_R` maximum | 1.4512e9 |
| one-copy physical `Q_P` | 4.3539e9 |
| coordinate metric forcing `-2Q_R` maximum | 2.9024e9 |
| normal gauge forcing `G_N` | 1.2e10 |
| lower normal affine gauge component `(Q_C)_y` | 1 |

The last two new caps follow directly from (8),(10), using the producer's
authenticated full field contractions, not a longitudinal replacement:

\[
|G_N|\le\tfrac43s_2 e^2+4s_1he+10a_*h^2+8a_*h\tau_*+4a_*A_{2,*},
\]
\[
|(Q_C)_y|\le4a_*h^2+3a_*h\tau_*+2a_*A_{2,*}
                +\tfrac43s_1he+\tfrac13s_1e\tau_*.
\tag{16}
\]

Here `h,e,tau_*,A2,*` bound `H,eta,tau,tr(A^2)` in absolute value,
and `a_*,s1,s2` bound `a,sigma',sigma''`. The exact rational majorants
are approximately `1.1609283815e10` and `0.582936557`, respectively.
They hold for `0<=x<=3`, `r=0`, all of the unit three-sphere, and refer
to this fixed background normal/component convention. The audit also
records proper clock derivatives of these two component/trace sources;
it claims no source normal derivative, full ambient covariant derivative
of `Q_C`, or bulk source C1 estimate.

`ACTUAL_SHELL_MAJORANTS_NORMAL.json` and its optimized companion are
byte identical. They pin the finalized physical shell report
`559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a`.
The independent fresh normal/optimized producer replays also reproduce
that complete report byte for byte, including 153 finite algebra comparisons
and six reason-specific rejection subprocesses in each mode. The review
receipt records the final source, input and report pins. This gives an
independent code/arithmetic audit of those physical source caps and the
new normal gauge cap; it does not replay all inherited PDE certificates.
