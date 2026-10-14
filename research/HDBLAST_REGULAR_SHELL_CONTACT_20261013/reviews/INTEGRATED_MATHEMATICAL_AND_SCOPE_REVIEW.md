# Independent integrated mathematical and scope review

14 October 2026 UTC. Internal AI mathematical review, not external peer
review. The companion JSON pins each reviewed author proof and manuscript.
Authors' proofs and the draft were not modified. No scientific checker,
network operation or earlier checkpoint replay was run for this review.

## Verdict

No blocking mathematical defect was found in the reviewed proof chain or
public draft. The result supports actual smooth linearized coupled solutions
after first shell contact on the stated full ambient moving-ball slab,
restricted to the specified finite angular family. The physical pairing and
non-gauge claims have separate proofs and are not inferred from auxiliary
master energy or nonzero boundary trace.

One optional scope clarification remains: distinguish existence for arbitrary
admissible collar displacement/velocity data from the proved activation of
the selected pure-displacement and pure-velocity states. In particular, the
contact infimum a0 in the activation argument is the actual support infimum
of the selected nonzero profile, rather than merely a lower bound on its
support. Later manuscript paragraphs already use that stronger hypothesis,
so this is a wording refinement, not a defect in the construction.

The integrated-check-count placeholder is intentional and still pending.
This review makes no machine-replay completion or publication claim.

## 1. Domain and full ambient-slab coverage

For D=cos(tau)+cos(chi), |tau|<pi/2 and 0<=chi<=pi/2 give D>0.
The exact invariant is b²(cos(tau)-cos(chi))/D. Its chi derivative is
2b² cos(tau)sin(chi)/D², so its range is exactly
[-b²tan²(tau/2),b²]. The strict choice of tau_* keeps the closed rectangle
within the actual nonvanishing-f background interval.

Inverting the conformal map gives
tau=arctan((T+R)/b)+arctan((T-R)/b). For T>0,

    partial_R tau=-4bTR/[(b²+(T+R)²)(b²+(T-R)²)]<=0.

Its maximum on a complete moving-ball slice is therefore at R=0, with
tau=2arctan(T/b). Symmetry handles T<0. Thus the entire ambient slab
|T|<=b tan(tau_*/2), R<=sqrt(b²+T²), is covered, not only its shell or
an exterior wedge. The shell reaches b tan(tau_*), which is a different
bound and is not incorrectly used as whole-slice coverage.

The collar support x<=d0+|tau| is strictly away from the axis because
d0+tau_*<pi/2. Extension by zero near that axis is smooth in Cartesian
coordinates, even for l>=2 harmonics. The auxiliary coefficient extension
outside the support's causal collar does not alter the solution. The optional
d0+2tau_*<pi/2 also implies chi>|tau| and s>0 on the full support.
The construction does not claim that the older near-cone pair has reached
the shell within its prior coefficient domain.

## 2. Boundary conversion, energy and smooth existence

The conformal weight3/2 and cylinder curvature12 give the bulk potential
9/4+Q²V_c. Normal differentiation contributes -3/(2b) to the shifted shell
Robin coefficient. Tangential conjugation gives the surface term
(3/2)sec²(tau)+(3/4)tan²(tau), the positive inertia mu=c/Q_b and precisely
the first derivative mu' u_tau. The manuscript and both independent energy
notes agree on these coefficients and orientations.

The even energy identity retains boundary kinetic and gradient energy;
mu_min>0 controls the cross term (k-beta)z z_tau. The coefficient bounds
give a finite-slab Gronwall estimate, with no contraction or infinite-time
claim. For odd parity the actual alpha z² potential cancels the normal flux.
The epsilon trace inequality and added bulk mass give coercivity. Adding an
unrelated boundary potential and assuming an L² trace of bulk velocity
would be invalid, and the manuscript explicitly avoids that step.

For fixed l the radial Liouville transformation produces
v_tt-v_xx+[Q²V_c+(lambda+3/4)sec²x]v=0. The incoming/outgoing equations
p_t-p_x=-Wv and m_t+m_x=-Wv, together with the boundary ODE, supply exactly
one characteristic datum. Even parity gives
mu z''+(mu'+1)z'+(mu lambda+beta)z=p(t,0), m=2z'-p;
odd parity gives z'=p-alpha z, m=p-2alpha z.

After inserting the p integral, the unsmoothed 2z' term in m is integrated
once when reconstructing v, yielding the required small interval factor.
The z and z' updates likewise have a factor Delta t; coefficient bounds
and mu_min make the supremum-norm map contractive on sufficiently short
intervals. Linear continuation uses the same coefficient bounds independent
of amplitude. Highest derivatives satisfy successive linear Volterra
equations. Initial data vanish on a full boundary neighborhood, so all formal
corner jets vanish and no matching jump appears on the corner characteristic.
The defects v_x-(p-m)/2 and v(t,0)-z(t) have zero time derivative and zero
initial value. This supplies a smooth solution, rather than inferring
regularity from a formal Galerkin energy estimate.

Only a fixed finite harmonic set is used; no uniform-l estimates or general
angular Sobolev completion are asserted. Scalar uniqueness is sufficient
for the reconstructed solution family, not uniqueness in every gravity gauge.

## 3. Both activation arguments

All undifferentiated coefficients are time even and mu' is odd. Uniqueness
therefore gives time-even solutions for pure displacement data, and time-odd
solutions for pure velocity data. If the shell trace were identically zero
on (0,tau1), time symmetry extends that zero trace to (-tau1,tau1).
Either boundary law then forces v_x=0 there.

After separation the sideways equation is an ordinary1+1 wave equation
with x as evolution time and bounded lower-order potential. Its zero Cauchy
data at x=0 imply zero solution in x+|tau|<tau1. Evaluating v(0,x), or
v_tau(0,x) for the time-odd case, contradicts a nonzero initial profile with
support infimum a0<tau1. No higher-dimensional timelike Cauchy theorem or
analytic unique continuation is used.

Consequently some nonzero trace occurs in each interval immediately after
a0; it need not occur at a0 itself or at every later instant. At the shell
T=b tan(tau), so T_c=b tan(a0)=(b²-R0²)/(2R0)<T0.
The positive Cauchy weights Q0^(3/2)/w and Q0^(5/2)/w preserve support
infimum in the two master states (q,0) and (0,q). The extra Q0 in velocity
is correctly retained. Both selected solutions therefore activate within
the full ambient slab.

## 4. Full reconstructed junctions and constraints

The normal graph formula gives delta K mixed=(Psi'+2aPsi+a'zeta)Id
-rho^-2 Hess(zeta), and delta p_n=chi'+2pPsi+p'zeta.
Induced metric/scalar continuity for the odd part gives
zeta_o=-Psi_o/a and P_o'=-p²P_o/(3a). The latter is the odd master law.
The even embedding choice zeta_e=0 is sufficient and makes the full Israel
tensor residual vanish, including its trace. The even scalar residual
vanishes by B P_e'=rho^-2(Box_gamma+4)P_e. No tangential eigenvalue or
trace-free harmonic tensor was divided out.

Independently, the difference of Hamiltonian constraints gives
3a delta K_o=p delta p_n,o. Direct substitution gives
delta K_o=-H(Box_gamma+4)zeta_o and
delta p_n,o=-3aH(Box_gamma+4)zeta_o/p. Hence
sigma delta K_o+sigma' delta p_n,o=0 using sigma=6a, sigma'=-2p.
The normal force is not an additional odd wave equation. Common induced
fields and bulk Codazzi give tangential conservation.

The five-dimensional ambient-coordinate kernel of the even embedding
Hessian equation is correctly kept as a necessity/completeness issue.
Selecting zeta_e=0 does not claim that every possible embedding is straight.
Inherited bulk residual identities supply all bulk linearized Einstein and
scalar equations for this sufficient family. Smooth master existence and
these graph/junction identities together establish the coupled result.

## 5. Non-gauge proof

On a polar annulus, projection of a smooth gauge covector onto H_l gives
xi_a=f_a H_l, xi_i=J D_i H_l. The reconstructed metric's angular trace-free
and mixed components vanish. The scalar angular trace-free Lie derivative
is2J HessTF(H_l), whose squared norm is
(2/3)lambda(lambda-3)||H_l||². For l>=2 it is strictly positive, forcing
J=0; the mixed equation then forces f_a=0. Other angular sectors are
orthogonal to this scalar projection. A nonzero reconstructed metric in
the annulus cannot be a smooth Lie derivative. A pure gauge bulk-plus-shell
pair would in particular be pure gauge in the bulk, so adding embedding
variation cannot defeat this obstruction. No l=0,1, nonsmooth gauge group
or finite-cut charge classification is inferred.

## 6. Thin-annulus pairing and actual current conservation

The physical T=0 pullback has the scalar weight9r³/(4Omega f²), metric
weight6r^7/Omega and factors(D_r+4),(D_r+5). The extra unit is the
commutator[partial_T,E]=partial_T. With r partial_r log(m_r)=3-2sL,
radial integration by parts gives the exact diagonal coefficient2+9sL.

In I=(r1,r2), the gradient term is at least m_* integral|q_r|², and the
negative zeroth-order term is bounded by M_* integral q². Dirichlet
Poincare yields the stated lower bound
(pi²m_*/w_I²-M_*) integral q². It applies for each angular point before
integration. A regular positive-radius shell collar has m_* bounded away
from zero and M_* bounded while the interval shrinks. A sufficiently thin
annulus can therefore also lie inside the admissible activation collar.

The same-parity displacement/velocity pair has doubled physical pairing
2K(q,q)>0. Testing a kernel element against(0,q), then(v,0), proves weak
nondegeneracy on precisely the smooth compact thin-annulus image with a
fixed finite harmonic set. Symmetry of K is not required. It is neither a
physical energy theorem nor a statement about every completed state.

For arbitrary packet pairings, even matching makes the physical radial
current zero pointwise. Odd matching makes the completed shell divergence
equal the sum of radial bulk currents. No common separated spectral value
is used. The shell term includes both gravitational copies and one tension.

Coordinate-density transformation gives
omega^R-R_b'omega^T=omega^y/R_b, while omega_Sigma^T=omega_Sigma^t and
angular shell densities acquire dt/dT=1/R_b. The moving-domain derivative
therefore cancels the shell divergence with the sign in manuscript(28).
The future-oriented minus sign agrees with manuscript(29)'s positive
Cauchy pairing. Compact spheres have no angular boundary; the constructed
solution is smooth and zero near the origin, so there is no extra endpoint
flux. The cone is internal, not an unaccounted boundary.

Conservation transports the nonzero total pairing and nondegeneracy through
contact. It does not prove a nonzero local shell Wronskian or energy transfer;
the manuscript explicitly says so. Even boundary activation coexists with
zero physical radial symplectic current and adds no physical shell kinetic
field. The auxiliary dynamic-boundary energy is kept separate.

## 7. Public scope and remaining work

The draft accurately states a local linearized finite-harmonic solution
family, full junction sufficiency, two selected activating states and a
conserved nonzero total physical pairing. It keeps unresolved unrestricted
IBVP regularity, uniform angular estimates, completed phase spaces, full
gauge quotient, nonlinear integrability/stability and observational origin.
Prior scalar boundary methods are acknowledged; no general boundary theorem
is claimed new. Citation scopes agree with the frozen literature audit.
No external originality, peer review, APS acceptance or cosmological result
is established by this internal review.
