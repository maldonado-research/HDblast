# Exact finite-regularity evolution across a prescribed shell-interaction interval

15 October 2026 UTC. Internal mathematical checkpoint for the registered
HDBLAST Einstein–scalar model. This is a conditional analytic construction,
with its hypotheses stated below. Finite symbolic checks verify specified
identities; they do not computer-certify the analytic proof. No external
novelty or external peer review is claimed.

## 1. The result and its hypotheses

Retain the five-dimensional equations and identical-copy physical shell laws

\[
 E_{AB}=\operatorname{Ric}_{AB}-\phi_A\phi_B-\tfrac23U(\phi)g_{AB}=0,
 \qquad S=\Box_g\phi-U'(\phi)=0,
\]
\[
 K_{ab}=\tfrac16\sigma(\phi)h_{ab},\qquad
 n\phi=-\tfrac12\sigma'(\phi).
 \tag{1}
\]

Each copy uses its own outward unit normal toward the common shell. The
registered functions are
\(W=1-\phi+\phi^3/3\), \(U=(W')^2/2-2W^2/3\), and
\(\sigma=2W+\delta(1+c\phi)\), with the fixed registered constant
\(c=96714321425692/161839258930731\). Fix one admissible background member;
there is no uniform detuning-limit assertion. The argument uses its smooth
coefficients and the margins below, not a transfer of tail bounds from a
different scalar branch or tension model.

The background on the controlled positive-cone region is

\[
 g_0=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,\kappa_{S^3}],
 \qquad \phi_0=\phi_0(y),\qquad z=\log(b/r),\quad dy=-\rho\,dz.
 \tag{2}
\]

The hypotheses are these.

1. Fix a finite time slab \([0,T]\) and a spatial collar
   \([0,Z]\times S^3\) on which this background exists smoothly and
   \(\rho>0\). Fix a sufficiently small state neighborhood on which the
   time and spatial quadratic forms used in the positive auxiliary wave
   energy have uniform positive margins, the shell is timelike, and the
   exact boundary state flow remains defined. These are local hypotheses
   on one background and one finite slab.
2. Retain the support, reflection, adjoint-kernel and weighted-inverse
   hypotheses of the previously established exact localized initial-data
   theorem. Use that theorem at the single fixed value \(k=12\): it
   supplies a parameter-smooth exact canonical constraint curve with metric
   in \(H^{14}\), metric momentum in \(H^{13}\), and the prescribed scalar
   sectors. It equals the background on fixed open initial neighborhoods
   of the shell and vertex. Corrections have support in a fixed closed
   annulus. This input is reproduced without alteration in the capsule.
3. Select its aligned, reflection-invariant quadrupole tangent
   \(A=B=\operatorname{diag}(1,-1,0,0)\), with the previously stated smooth
   annular displacement and velocity profiles. Thus the ten necessary
   second-order charges vanish. The exact initial family, rather than
   vanishing charges alone, is the input to this theorem.
4. Fix a strict radial causal margin. If the initial support has
   \(z\le L\), choose a common bound \(c_*\) for the actual perturbed
   radial characteristic speeds and \(Z>L+c_*T\), including the fixed
   transformation collar. Near the background one can take
   \(c_*=1+\eta\) for fixed \(\eta>0\) after shrinking the state
   neighborhood. Smooth coefficient extensions or an artificial outer
   wall are construction devices outside this domain of dependence.

Write
\[
 \mathcal X_s=\bigcap_{j=0}^{s}C^j([0,T];H^{s-j}([0,Z]\times S^3)).
 \tag{3}
\]
The component norms use a fixed positive product bundle metric. Subtract the
background when taking these norms.

**Theorem.** Under these hypotheses there is a common \(\varepsilon_0>0\)
such that, for \(|\varepsilon|<\varepsilon_0\), the selected exact canonical
initial family has an exact Einstein–scalar evolution on both copies of
the controlled positive-cone region for the whole prescribed interval \([0,T]\), satisfying induced-field matching and both
physical laws (1) at the shell. In the fixed relative-harmonic gauge its
perturbation belongs to \(\mathcal X_{14}\). For the fixed completed reduced Cauchy data, including the chosen initial
lapse and shift values, the solution is unique in the stated small
state/regularity class in that gauge and has
continuous dependence in the top finite Sobolev norm. The parameter curve
is \(C^2\) into \(\mathcal X_{10}\), and hence into a lower piecewise
spacetime \(C^4\) norm. The initial gauge can be chosen so that its full
first spacetime tangent is the previously selected harmonic quadrupole.

The asserted spacetime domain is the controlled region specified above.
Where an actual regular interior background is separately available, the
solution patches to it on an open inner overlap on which the perturbation
vanishes. This argument does not infer a global regular foliation through
the cone vertex from positivity of the collar alone.

The parameter radius may shrink with the fixed slab, collar, support gap,
background and analytic constants. This is a finite-time, finite-regularity
statement. If the background and collar hypotheses permit a prescribed
\(T\) beyond the reference packet contact time \(t_c=\log(b/R_0)\), the
constructed interval includes that reference interaction window. The actual
first arrival need not equal \(t_c\). The theorem does not, by itself,
assert a nonzero gauge-invariant shell response.

## 2. Exact boundary reduction and initial completion

Use the relative wave-map vector
\[
 C^A=g^{BC}(\Gamma[g]^A{}_{BC}-\Gamma[g_0]^A{}_{BC})
\]
and the reduced equations \(E_{AB}=\nabla_{(A}C_{B)}\), \(S=0\), with
one half in symmetrization. There are sixteen common scalar-wave principal
unknowns: fifteen metric components and the scalar. At the shell impose
four coordinate rows \(g_{ya}=0\), the ten Israel rows, the scalar law,
and the single normal gauge row \(C_y=0\). No additional tangential gauge
rows are imposed on the primal system.

The exact boundary calculation, reproduced in
[the retained state-flow proof](inputs/NONLINEAR_BOUNDARY_STRAIGHTENING.md),
uses \(q_D=(g_{ya})\), \(q_N=(N,h_{ab},\phi)\), \(N=\sqrt{g_{yy}}>0\),
and
\[
 f_N=\tfrac23\sigma N^2-a_bN^3\operatorname{tr}_h h_{0,b},\qquad
 f_h=\tfrac13\sigma Nh,\qquad f_\phi=-\tfrac12\sigma'N.
 \tag{4}
\]
For \(U=q-q_0\), the Neumann-block boundary vector field is
\(\mathcal B(u)=-\rho_b[f(q_{0,b}+u)-f(q_{0,b})]\). Its local state flow,
with flow time equal to the inward coordinate near the shell and cut off
inside a fixed collar, gives an invertible point map \(U=\Psi(x,V)\).
It converts all boundary rows exactly to
\[
 V_D=0,\qquad \partial_zV_N=0,\qquad \dim D=4,\quad\dim N=12.
 \tag{5}
\]
The transformed equations retain the actual scalar principal metric:
\[
 g^{\mu\nu}(x,V)D_\mu D_\nu V=F(x,V,DV),\quad F(x,0,0)=0,
 \qquad g^{za}|_{z=0}=0\quad(a=t,\Omega).
 \tag{6}
\]
The last identity holds for every state with \(V_D|_0=0\), including every
iterate; it is not an assumption of angular symmetry of the evolved fields.

[The gauge companion](gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md)
chooses initial lapse and shift by
`alpha_epsilon=alpha_0+epsilon alpha_1` and
`beta_epsilon=beta_0+epsilon beta_1`, with the desired first tangent,
and solves the five initial harmonic equations for their time velocities.
The canonical fields remain the prescribed exact constraint data. The local
formulas map \(H^{14}\times H^{13}\) smoothly to reduced Cauchy data at
those orders and preserve the open background shell collar. In particular
initial \(C=0\). The physical Hamiltonian and momentum constraints are
projections of \(H_{AB}=E_{AB}-\tfrac12g_{AB}\operatorname{tr}E\), and
with the reduced equation they imply \(\partial_t C=0\) initially.

All finite initial time jets of \(V\) vanish on the open shell collar:
solve the equation recursively for its next time derivative, using
\(g^{tt}\ne0\), locality and \(F(x,0,0)=0\). Every required corner trace
therefore vanishes. A protected short reduced Cauchy evolution, cut off
smoothly in time within its precontact interval, is a regular trial seed
with these exact jets. This supplies the iteration seed without presuming
the desired through-shell evolution.

## 3. Variable-coefficient existence on the entire prescribed interval

The complete analytic construction is given in
[the variable mixed-wave proof](astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md)
and in [the independent derivation](sol/INDEPENDENT_QUASILINEAR_MIXED_PROOF.md),
with a separate [frozen linear construction](sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md).
Here are the estimates that connect those proofs to this model.

Set \(a=-g^{tt}>0\), \(b^i=g^{ti}\), \(c^{ij}=g^{ij}\). The frozen
principal operator, after an overall sign change, is
\(aV_{tt}-2b^iD_iV_t-c^{ij}D_iD_jV\). Its positive auxiliary energy is
\[
 \mathcal E[V]=\tfrac12\int
 [a|V_t|^2+c^{ij}\langle D_iV,D_jV\rangle+|V|^2].
 \tag{7}
\]
The exact principal boundary flux at the shell is
\(-c^{zz}\langle V_z,V_t\rangle\), because \(b^z=c^{zA}=0\).
It vanishes by (5). All coefficient commutators are retained. This energy
is analytic control; its positivity is not a statement about physical
gravitational energy or stability.

Time and angular commutation preserve (5). Normal derivatives are recovered
from the equation using \(c^{zz}>0\), rather than assigned fictitious
homogeneous boundary conditions. At order fourteen, the finite Leibniz
counts give tame estimates with low jet order ten. If
\(X_s(V;t)=\sum_{j=0}^s\|\partial_t^jV(t)\|_{H^{s-j}}\), the exact
lower-jet identity
\[
 X_{13}(V;t)\le X_{13}(V;0)+\int_0^t X_{14}(V;s)\,ds
 \tag{8}
\]
accounts for the pointwise source costs in normal recovery. Thus the
compatible frozen Picard successors obey
\[
 X_{14}(V^{n+1};t)\le C\varepsilon+
 C\int_0^t[X_{14}(V^{n+1};s)+X_{14}(V^n;s)]\,ds.
 \tag{9}
\]
Here \(\varepsilon\) denotes the size of the initial perturbation, bounded
by a constant times the absolute family parameter. The same exact initial
jets are preserved at every iterate. The underlying smooth frozen solver
is constructed variationally before taking limits; it does not assume
that a variable coefficient commutes with an angular projection.

An exponential majorant \(C_0\varepsilon e^{C_1t}\) closes (9) on the
whole fixed interval when \(C_0\varepsilon e^{C_1T}<R\), where \(R\) is
inside the chosen state and energy margins. Differences have zero initial
and boundary data, and the energy norm satisfies a Volterra estimate
\(e_{n+1}(t)\le C\int_0^t e_n(s)\,ds\). Iteration gives factorial
convergence. This proves existence on the prescribed interval directly.
It does not rely on a general finite-Sobolev restart theorem at arbitrary
later nonzero corners.

Smooth approximations preserve a common smaller open initial shell collar.
The tame difference estimate loses one derivative; the additional
Bona–Smith comparison places the sole extra derivative on the coarser
smoothed solution. Its \(O(N)\) high norm multiplies an \(o(N^{-1})\)
lower initial-data tail. Consequently the approximations converge strongly
in \(C_tH^{14}\cap C_t^1H^{13}\), and the equation supplies all the time
jets in (3). Top continuity is obtained by this comparison, not inferred
from weak compactness. The smoothing approximants are a device for the
reduced equations; the limiting data are the original exact constraint data.

Localized energy bounds propagation by the actual permitted characteristic
speeds. The strict margin in hypothesis 4 removes any influence of the
artificial outer construction on the asserted physical solution.

## 4. Exact physical constraints, both copies and parameter derivatives

The exact Einstein–scalar identity is
\[
 \nabla^AE_{AB}-\tfrac12\nabla_B\operatorname{tr}E=-S\phi_B.
 \tag{10}
\]
For a reduced solution and its exact scalar equation it gives
\((\Box_g+\operatorname{Ric}_g)C=0\). Codazzi and (1) imply
\(E(n,e_a)=0\), including the derivative of the field-dependent tension.
Together with \(C(n)=0\), the subsidiary boundary conditions are
\[
 C(n)=0,\qquad(\nabla_nC)_a-K_a{}^bC_b=0.
 \tag{11}
\]
In the actual shell coordinates these are
\(C_y=0\), \(\partial_yC_a-N\sigma C_a/3=0\). The coordinate Robin
coefficient contains two copies of the extrinsic-curvature factor, one
from the connection. A bounded collar conjugation changes these five
homogeneous rows to a fixed Dirichlet/Neumann split, preserving the scalar
metric principal symbol. The same variable energy proves subsidiary
uniqueness. The exact initial data \(C=\partial_tC=0\) therefore give
\(C=0\) throughout the interval. The reduced solution solves (1) exactly.

Use the same initial data and construction on both interior copies.
Uniqueness makes their identified induced fields equal; both physical
outward-normal shell laws hold. The normal-force equation follows from
the exact identity and reflection symmetry, as derived in the companion.
It is not an extra incoming boundary condition.

For parameter regularity, linearization of (6) and then a second
linearization have exactly homogeneous rows (5). The worst first-variation
coefficient is \(A_V[W]D^2V\); the worst second-variation source is
\(W D^2W\). The model's state-only principal coefficients and the proved
mixed estimates allow the conservative difference-quotient argument in
[the companion's section 5](gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md)
and its [finite-regularity variation completion](gauge-parameter/FINITE_REGULARITY_VARIATION_COMPLETION.md).
The latter constructs the variation limits from actual smoothed nonlinear
solutions and requires strong continuity only at the lower orders used
in the parameter conclusion. Four Sobolev orders are reserved: the solution curve is \(C^2\) in
\(\mathcal X_{10}\). Smooth composition with the inverse state flow
preserves this conclusion. In four spatial dimensions that class embeds
in the stated piecewise spacetime \(C^4\) norm. This completes the theorem
at one fixed finite regularity.

## 5. Relation to the earlier perturbative results

The first derivative of this family has the selected physical initial
quadrupole and the chosen harmonic gauge initial data. Linear uniqueness
identifies it with that first perturbation on the controlled interval.
Its second parameter derivative exists in the stated lower norm and obeys
the full second-order bulk, shell and gauge equations with the initial
second derivative supplied by the exact canonical family.

That second initial derivative need not equal the particular correction
selected in the earlier finite-angular Taylor construction. This theorem
does not silently identify those two second coefficients. Their respective
initial data remain part of their definitions.

The earlier noncommuting rotating-quadrupole obstruction is unchanged.
Within the restricted two-matrix family, passing the ten tests still does
not establish existence for every state with \(AB=BA\). The construction
here applies to the selected reflection-invariant aligned seed together
with the additional exact-data, background and support hypotheses above.

The next scientific target is to calculate a gauge-invariant shell response
for this family, prove or measure its nonvanishing with controlled errors,
and identify which coefficients depend on the choice of nonlinear initial
completion. Existence across a reference interaction interval supplies a
foundation for that calculation; it does not replace it. No global
stability, common spatial \(C^\infty\) parameter interval, cosmological
origin, observational agreement or APS readiness is established here.

## 6. Literature and verification scope

The construction uses established wave-equation methods and model-specific
boundary identities. The [primary-source review](literature/QUASILINEAR_CONTACT_PRIMARY_LITERATURE.md)
checks Shibata–Kikuchi's variable-coefficient Neumann theory, related
quasilinear results and reference-lifespan arguments against their actual
hypotheses. Their all-Neumann or scalar theorems are not quoted as a direct
mixed Einstein–scalar theorem. The mixed proof and model interfaces are
written here and internally cross-reviewed.

The reproducible capsule checks energy-divergence signs, the complete
Dirichlet/Neumann flux cancellation, lapse/shift harmonic completion,
Codazzi cancellation and the subsidiary Robin factor. Deliberately wrong
formulas are rejected in normal and optimized Python. These tests establish
the indicated finite algebraic results. The analytic existence and
parameter arguments remain written mathematical proofs with internal
review, not formal proof-assistant certificates or external peer review.
