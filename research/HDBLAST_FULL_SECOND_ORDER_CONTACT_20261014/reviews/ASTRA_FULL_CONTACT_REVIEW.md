# Independent Astra review of the full finite-angular contact coefficient

Date: 14 October 2026. Internal mathematical review, not external peer review.

Reviewed artifact:
`FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md`

Exact reviewed SHA-256:
`59c58b560b257a7311d4d00a72a1d4facb95c23377b88cf0323a4d57edc81949`

## Verdict

**Accepted in its stated conditional scope.** I found no remaining fatal
mathematical gap in the construction of a smooth complete second Taylor
coefficient on the specified finite positive-cone collar, assuming the
frozen first-order contact solution and localized second-order physical
constraint correction used as inputs. The proof constructs all reduced
fields and then proves the full physical residual vanishes; it does not
merely count boundary rows or solve a tensor projection.

This acceptance does not assert an exact nonlinear contact family, a
uniform estimate over all angular frequencies, continuation through the
cone, or independence from the localized initial-constraint hypotheses.
Those stronger statements are correctly excluded by the reviewed note.

## 1. Gauge conversion and the normally hyperbolic system

The relative wave-map vector is appropriate with the reference connection
held fixed. Its linearized metric operator has principal part
`-Box_0/2`; the scalar principal part is `Box_0`. Coupling to the scalar
background introduces only first derivatives of the unknown metric or
scalar coefficient. Thus the finite-mode reduction has a common scalar
wave principal operator on every component.

I independently checked the first longitudinal gauge residual:

\[
DC_y[h_1]=-6(\Psi_y+4a\Psi),\qquad DC_\mu[h_1]=0.
\]

For lower-index gauge vectors,
\(DC[\mathcal L_\xi g_0]=(\Box_0+\operatorname{Ric}_0)\xi\).
The proposed first vector wave problem therefore sets the complete first
subsidiary field to zero. At the shell,

\[
(\mathcal L_\xi g_0)_{ya}
=\partial_y\xi_a+\partial_a\xi_y-2a_b\xi_a=0
\]

under its stated Dirichlet/Robin data. The vanishing normal component
makes the flow boundary-preserving, including when its tangential time
component is nonzero. Identical vectors on the two copies preserve their
matching.

The inverse quadratic conversion is correct. Inverting the flow gives
\(k^{\rm long}=k^w-2\mathcal L_\xi h_1^w+\mathcal L_\xi^2g_0\);
substituting \(h_1^w=h_1+\mathcal L_\xi g_0\) yields exactly the
printed \(k^w-2\mathcal L_\xi h_1-\mathcal L_\xi^2g_0\).
The scalar conversion has the same factors and signs. The note explicitly
recomputes the reduced quadratic bulk source from the transformed first
fields and retains \(Q_C=D^2C[(h_1^w)^2]\). It does not reuse a bare
longitudinal quadratic source in the new first gauge.

Zero initial \(\xi,\xi_t\) leave the physical initial metric,
extrinsic curvature and scalar Cauchy data unchanged; the generated
second time derivatives of \(\xi\) change coordinate gauge jets. The
quadratic pullback terms likewise do not change the spatial canonical
second data on that initial slice. This makes the use of the inherited
physical constraint correction consistent with the first gauge map.

## 2. Exact boundary completion and normal signs

The boundary is orthogonal only at the shell; the proof does not impose
Gaussian normal gauge throughout the collar. Consequently the metric
wave system retains its lapse and shift components.

With outward normal \(n=N^{-1}\partial_y\), \(g_{ya}=0\), and
\(g_{yy}=N^2\), the physical metric law gives
\(h_{ab,y}=\sigma Nh_{ab}/3\) and trace \(K=2\sigma/3\).
Direct contraction with the fixed reference connection gives

\[
C^y=\frac{N_y}{N^3}-\frac{2\sigma}{3N}
          +a_b\operatorname{tr}_h h_0.
\]

Thus the printed lapse boundary equation and its sign are correct. The
four mixed Dirichlet rows, ten physical metric rows, one scalar row and
one normal gauge row constitute four Dirichlet and twelve independently
solvable normal-derivative rows. Once the physical rows have been used,
the normal gauge row has no remaining tangential derivative or unknown
normal derivative on its right-hand side.

I expanded all three nonlinear boundary expressions to second order.
The displayed \(R_{ab}\), \(R_\phi\), and \(R_N\) agree, including
\(-16a_b\nu_1^2\), \(-2a_b\operatorname{tr}(A_1^2)\), and
\(6a_b\nu_1\operatorname{tr}A_1\) in \(R_N\). The relation
\(k_{yy}=2\nu_2+2\nu_1^2\) is essential and is retained. These equations
therefore give a genuine lower-order matrix Robin system for the twelve
remaining metric/scalar wave unknowns.

No symmetry of this Robin matrix is established or needed by the
characteristic argument. The proof does not infer an energy cancellation
from a nonsymmetric matrix.

## 3. Physical initial constraints and gauge completion

The finite SO(4) argument is a projection of a **linear, already sourced**
second-order constraint equation. It does not claim that a finite angular
subspace is closed under the nonlinear initial-data map. The background
constraint operator commutes with the compact-group isotypic projection,
so \(P\Pi q=\Pi Pq=S\). Rotations preserve radial support, and the
projection on the actual tensor/density bundles has finite multiplicity
for this finite collection of SO(4) representations.

The symmetric-square decomposition is correct for the chosen fixed
quadrupole family. All linear gauge maps are equivariant, and the full
quadratic source is a symmetric bilinear image of that same first
representation. The selected finite system therefore includes the
exceptional scalar monopole without dividing by its eigenvalue. The
argument excludes a forced coexact-vector representation for this state;
it need not exclude arbitrary independently added homogeneous vector data.

The initial gauge equations solve for five coordinate time derivatives,
with coefficients \(-1/(2\rho^2)\) and \(-1/\rho^2\) as printed.
They do not alter the chosen physical spatial metric, extrinsic curvature
or scalar Cauchy data. The use of the actual trace-adjusted Einstein
constraint is correct: when \(C^{(2)}=0\) on the initial slice,

\[
\mathscr H^{(2)}_{00}=\tfrac12\partial_0 C^{(2)}_0,
\qquad
\mathscr H^{(2)}_{0I}=\tfrac12\partial_0 C^{(2)}_I.
\]

The sourced initial physical constraints therefore give the required zero
initial time derivative of the subsidiary vector. Replacing this step by
\(\mathscr E_{00}^{(2)}=0\) alone would have been wrong; the reviewed
version does not do so.

One precision raised during this audit was resolved in the pinned final
version: the freely selected initial lapse/shift values must themselves
vanish near the initial shell corner. The note now explicitly requires
this and states that zero values are admissible. With that choice, every
first field, physical correction, gauge completion and boundary source
has zero corner jets. The reduced wave recursion then supplies the full
smooth corner hierarchy.

## 4. Finite-mode characteristic construction

In inward \(z=-x\), the signs are consistent:
\(p=U_t+U_z\) is outgoing at the shell and
\(m=U_t-U_z\) is incoming. The Robin law gives

\[
q'=p_{o,N}-Mq+g,\qquad m_N=p_{o,N}-2Mq+2g.
\]

The order of the construction matters and is correct in the note:
compute the outgoing characteristic functional, solve the finite trace
ODE, substitute it into the incoming characteristic functional, and then
integrate the time derivative to reconstruct \(U\).

For fixed initial data, the outgoing functional has difference
\(O(\delta)\|\Delta(p,m,U)\|_\infty\). Solving the trace ODE adds
one time integral, so the trace difference is \(O(\delta^2)\) in this
comparison. Both the incoming boundary functional and interior
characteristic integrals are therefore \(O(\delta)\). This is a genuine
contraction for sufficiently small \(\delta\), also with first-derivative
matrix couplings in the bulk. It does not require a favorable Robin sign.

At its fixed point, the defect in
\(U_z=(p-m)/2\) has zero time derivative and zero initial value; the
boundary trace defect likewise vanishes by its ODE. Smooth compatible
data bootstrap to a smooth solution. Later short intervals use the
compatible jets of the preceding solution, so the initial corner does
not need to be artificially made zero again at each restart. A finite
number of bounded steps covers the specified finite slab.

## 5. Subsidiary boundary problem

The exact Bianchi/scalar identity produces
\((\Box_g+\operatorname{Ric}_g)C=0\) for the reduced system. Since the
background and first subsidiary fields vanish, its second coefficient is
homogeneous with the background operator; no quadratic forcing survives.

The matter-dependent Codazzi cancellation is central and has the correct
outward sign:

\[
R_{na}=D_bK^b{}_a-D_aK
       =-\tfrac12\sigma' D_a\phi
       =(n\phi)D_a\phi.
\]

Thus the complete mixed Einstein–scalar residual vanishes at the shell,
not merely a vacuum Ricci term. On a reduced solution, \(C_n=0\) and
this identity give

\[
(\nabla_n C)_a-K_a{}^bC_b=0.
\]

In boundary-orthogonal coordinates,
\((\nabla_nC)_a=N^{-1}\partial_yC_a-K_a{}^bC_b\), so the coordinate
condition is \(\partial_yC_a-2NK_a{}^bC_b=0\). Its second coefficient
is exactly \(\partial_yC_a^{(2)}-2a_bC_a^{(2)}=0\), together with
\(C_y^{(2)}=0\). These are one Dirichlet and four Robin subsidiary
conditions, induced by the physical primal rows rather than appended as
four additional primal conditions.

The same finite-mode wave lemma gives subsidiary uniqueness. Its zero
Cauchy data then force \(C^{(2)}=0\), converting every reduced metric
row into the corresponding physical Einstein equation. This step closes
the incoming gauge/constraint issue that a primal characteristic count
alone would leave unresolved.

## 6. Domain and final scope

Every gauge wave source and every second-order source is supported within
the causal radial support of the first fields and chosen initial
correction. The characteristic construction preserves
\(z\le L+t\). Consequently, on any stated finite \(t_*\), all
calculations remain on the actual smooth positive-cone collar, where the
principal coefficients are nonzero and the matrices are bounded. The
coefficient extension beyond that support cannot affect the solution.
The inherited conformal-time bound controls the first source on this
whole collar; the chosen inequality \(t_c<t_*\) gives a nonempty interval
past its causal contact time. This is not an ambient slab through the
cone or a claim of nonzero radiation at the support edge.

I accept the complete finite-angular **second-order coefficient** claim
under the listed inputs. The independently frozen Astra bulk note keeps
its own narrower scope; the additional initial-data, gauge, boundary and
subsidiary construction in the reviewed root note is what closes the
combined coefficient theorem. No new heavy checker execution, repository
mutation, or publication was performed for this review.
