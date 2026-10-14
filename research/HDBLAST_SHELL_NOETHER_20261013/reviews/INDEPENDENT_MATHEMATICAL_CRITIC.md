# Independent critic: nonlinear internal-shell identity and rotational obstruction

14 October 2026 UTC. **PASS, with no remaining blocking mathematical finding**
for the reviewed classical statements. This is an internal AI mathematical
review, not external peer review, a formal proof assistant result, or an
originality determination. Exact reviewed files and SHA-256 pins are recorded
in the companion JSON receipt. No author proof, frozen source or repository was
edited, and no scientific checker or prior PDE calculation was executed by
this reviewer.

## 1. What passes and the domain of the result

The Noether note derives the completed internal cut identity for unrestricted
matched bulk field variations and common shell-tangent reference gauge
generators. The separate ADM note derives a spatial-rotation constraint
functional with its exact tilted-cut term. They independently give the same
necessary quadratic zero condition, with the same action normalization and
the appropriate Cauchy-orientation conversion. For compact first-order
initial data, the constructed rotating degree-two state has a nonzero
constraint Hessian and therefore cannot be the tangent to an exact family in
the stated differentiability class.

The hypotheses matter: the action and its functions U and sigma are fixed;
both outward GHY completions and one scalar-dependent tension are retained;
induced metric and scalar match exactly; two compact Cauchy balls have one
common closed cut and no external physical boundary; the timelike shell and
uniformly spacelike section remain nondegenerate. The claimed exact family
is parameter-C² in a classical topology controlling the needed spatial and
time derivatives, traces, pullbacks and products, with a common smooth
reference identification. Piecewise spacetime C⁴ and the explicit collar
hypotheses are sufficient for the displayed calculations. A weaker finite
energy family is not automatically included.

The obstruction allows arbitrary second-order corrections in every
gravitational and scalar sector and at the shell. It assumes compact support
only for the **first-order** initial perturbation and embedding jet. Neither
parity preservation of the exact family nor nonlinear closure of the master
image is needed. Changing the model parameters, adding an external cut charge,
or adding shell fields carrying independent angular stress changes the problem.

## 2. Unrestricted GHY variation and the relative cut

I checked the local boundary variation using a Gaussian chart for the
background while keeping q_nn and q_na unrestricted. With the R/2
normalization, covariant metric variation q_AB, outward spacelike normal and
K_ab=h_a^A h_b^B nabla_A n_B, its sum with the GHY variation is

```text
theta|W + delta(K epsilon_W)
 = -[ (K^ab-K h^ab) delta h_ab/2 + p delta phi ] epsilon_W + d beta,
beta = - i_v epsilon_W/2,   v^a = h^ab n^A q_Ab.
```

The minus sign and half factor in beta pass. The scalar tension variation
supplies sigma' delta phi once. Summing the two copies gives the full Israel
residual I and scalar residual P=p_plus+p_minus+sigma', not two tension terms.

For a common shell-tangent fixed vector xi, let A=h(nabla_n xi). Its normal
extension enters the Komar pullback as -i_(A+K.xi)epsilon_W/2 and the corner
evaluation as -i_(A-K.xi)epsilon_W/2. Their difference cancels A exactly:

```text
B_xi = sum Q_xi + i_xi ell - beta(G_xi) = -i_(I.xi) epsilon_W.
```

This identity needs no Killing assumption. It would fail if the corner were
dropped or its sign reversed. The exact fixed-parameter Cartan calculation is

```text
k_xi - (delta beta)(delta,G_xi)
 = delta B_xi - i_xi alpha(delta) + d(i_xi beta(delta)).
```

The explicit delta I, volume and index-raising terms in the note are correct.
Terms multiplying the nonlinear I vanish only after imposing I=0. Full
linearized matching/junctions remove the remaining delta I term. The last
exact differential integrates to zero on the **closed** cut; it remains an
edge term for a cut with boundary.

Outward-normal Israel sums and oriented forms cannot be mixed casually. The
notes consistently convert every boundary form, including theta, Q and beta,
to the same outward co-oriented density convention. Raw pullbacks using one
global doubled orientation would instead require the indicated side signs.
The future end of the timelike boundary has the opposite orientation to the
outward boundary of a future Cauchy section. This fixes
Theta_rel=sum int_C theta-int_S beta and its curl, rather than a second
additional corner. Curling the boundary relation and using Stokes verifies
the stated conservation sign.

## 3. Off-shell identity, moving embeddings and reducibility

For the chosen Lee-Wald potential and Komar normalization, I independently
checked the bulk identity by varying J_xi=theta(G_xi)-i_xi L:

```text
C_xi = E^A_B xi^B epsilon_A,  E=G-T,
omega(delta,G_xi) = delta C_xi + d(delta Q_xi-i_xi theta)
                    + i_xi(Euler-Lagrange form . delta fields).
```

Combining this with the previous cut identity gives the derived relative
functional H_xi=sum int_C C_xi+int_S B_xi and the displayed bulk/shell residual
terms. On a nonlinear solution and a full linearized solution, those terms
and delta H vanish. This proves relative gauge degeneracy; conservation of a
restricted master current alone would not prove it.

The moving-domain step is an actual change of variables in the full matched
action. The formula delta(Y*Phi)=Y*(delta Phi+L_chi Phi) retains embedding and
evaluation terms. Applying the bulk Noether decomposition to theta(L_chi Phi)
gives i_chi L+C_chi+dQ_chi. The displayed moving potential therefore retains
the domain displacement, constraint and gravitational cut terms. The barred
GHY corner is also retained, and the common shell is not added again.
Field-space differentiation must vary chi and the pullbacks consistently;
the barred definition of the curl avoids an ambiguous component shortcut.

An inverse ambient embedding transformation alone leaves the barred fields
unchanged. That auxiliary fact is **not** the proof of the nontrivial common
reference diffeomorphism's degeneracy: the latter is established by the
relative residual identity. The common, fixed worldvolume compensation is
essential. A background Lorentz Killing field alone moves a parametrized
embedding tangentially; its fixed common compensation makes the entire
background generator zero. On graph perturbations its differential is the
usual Lie action on the fields and graph displacement. No field-dependent
compensator formula is silently used.

At a reducible background G_0=0, the off-shell identity gives DH_0[w]=0 for
every smooth variation in the exact matching chart. Differentiating once in
a full linearized direction makes the derivatives of all displayed residual
terms zero. Hence D²H_0[v,v]=Omega_rel(v,DG_0[v]). A parameter-C² exact family
with H=0 must satisfy this necessary zero. The halved Taylor coefficient is
half the Hessian; neither its sign nor that factor changes the zero condition.

The scalar derivative f=Phi_s nonzero implies Ks=0 for any common
metric-scalar isometry, and then L_K eta=0. A flat translation cannot preserve
s on an open set. Thus the ten common Lorentz generators identified in the
note are correct. This does not classify an adjoint kernel in an unspecified
weak constraint domain.

The ordinary moving-shell normal and tangential residual relations also
check independently. From the spacelike-shell Gauss constraint,

```text
sigma Delta K + sigma' Delta p
 = 2 Delta E_nn + I_ab Delta K^ab + P Delta p,
sum E_na = D_b I^b_a - P D_a phi.
```

The four-dimensional trace reversal and all scalar signs agree. These exact
relations retain their quadratic products and explain the on-shell force
consistency; a first-order radial normal-force equation would not suffice.

## 4. Independent ADM route and arbitrary second-order cut data

I rederived the tilted normal relation from
n=cosh(eta)s+sinh(eta)u and tau=cosh(eta)u+sinh(eta)s. Differentiating along a
common cut-tangent J gives

```text
k(s,J) = -K(tau,J) - J eta.
```

The minus sign in front of J eta is fixed by the Lorentzian inner product.
With pi=sqrt(q)(k^ab-k q^ab)/2, integration of the ADM momentum constraint
gives

```text
C_J = sum int_C H_a J^a - int_S sqrt(c) I(tau,J)
    = sum int_C [pi^ab L_J q_ab + varpi J phi]
      + int_S sqrt(c) J(eta_plus+eta_minus).
```

The interface subtraction is derived from integration by parts, not guessed.
Pure scalar-dependent tension has zero exact mixed projection on (tau,J),
and introduces no independent intrinsic scalar kinetic momentum. This
functional vanishes on every exact bulk-constraint and Israel datum.

At the time-symmetric background, pi=varpi=eta=0 and J annihilates q_0,
phi_0 and the cut density. The bulk first derivative is zero. The cut first
derivative for **arbitrary** w is int sqrt(c_0)J(delta eta_sum[w])=0 by
closed-cut integration. Thus arbitrary second-order shell boundary values
cannot cancel the Hessian. Compact first-order support makes the remaining
corner Hessian 2 int delta sqrt(c) J delta eta_sum zero. Evaluation and
matching corrections are in the pulled-back w and have been removed by the
zero first derivative, rather than assumed absent.

The ambient reconstruction directly gives at T=0

```text
delta pi^ij = -3 p x_i x_j/(2 Omega),
chi = -3(D+4)q/(2 Omega² f),
delta varpi = Omega³ chi_T,
chi_T = -3(D+5)p/(2 Omega² f).
```

I checked the metric, lapse and vanishing initial shift in this calculation.
The additional unit in D+5 is [partial_T,E]=partial_T. Rotations preserve the
tensor basis s eta_AB-3X_AX_B, so its Lie derivative acts on F as claimed.
The physical one-copy pairing is K(q_1,p_2)-K(q_2,p_1), while the chosen
canonical/relative convention has the opposite sign. Summing both copies
therefore gives D²C_J=-Omega_phys(v,Jv) for these compact initial states.
The revised note distinguishes the one-copy pairing from the doubled total.

## 5. The nonzero witness and what it proves

Take q=QH and p=omega QJH, with real degree-l angular harmonic, JH nonzero,
omega nonzero, and a smooth profile Q supported in a sufficiently thin annulus
strictly inside the shell. Spherical integration by parts gives
int H J²H=-int(JH)². Both radial terms have the same diagonal K_rad(Q,Q);
no unproved symmetry of the full radial bilinear is needed. Thus

```text
Omega_one(v,Jv) = -2 omega K_rad(Q,Q) ||JH||²,
D²C_J[v,v] = 4 omega K_rad(Q,Q) ||JH||²
```

for the doubled state of either reflection parity. Odd-copy sign reversal of
both entries leaves its quadratic contribution unchanged.

The diagonal radial integration by parts produces 2+9r²L, not an assumed
constant sign. On a positive-radius compact annulus, m_r r² has a positive
minimum and the negative potential part is bounded. Dirichlet Poincare gives
the stated bound pi²m_min/width²-M_max. Shrinking the annulus makes this
strictly positive. This establishes positivity of this diagonal bilinear for
the witness, not positivity of gravitational energy.

For H=n_1²-n_2² and J=x_1 partial_2-x_2 partial_1, I checked
JH=-4n_1n_2, J²H=-4H, int H²=pi²/3 and int(JH)²=4pi²/3. The doubled Hessian
is therefore 16pi² omega K_rad/3, which is nonzero. These are real, smooth,
finite-angular first-order solutions admitted by the authenticated master
and junction construction; the degree-two metric projection excludes smooth
pure gauge. The support can be chosen in the thin controlled collar required
for the prior finite-slab theorem. The obstruction itself already follows on
the initial compact slice and does not require a universal contact claim.

Consequently the result is a genuine **nonintegrability of that tangent state
in the specified exact-family class**, rather than merely a nonzero auxiliary
scalar radiation current. It is not a claim that no nonlinear solutions
exist, that every first-order state is obstructed, or that weaker uncontrolled
families are excluded. Single-real-harmonic data with aligned position and
velocity can pass the ten symmetry tests; those zeros prove no converse.

The separate quadrupole note also passes. For one fixed Q, omega and reflection
parity, q=Q n^T A n and p=omega Q n^T B n give an eighteen-dimensional space
with real traceless symmetric A and B. The exact fourth moments give
int H_A H_B=Vol(S³)tr(AB)/12 and J_M H_A=H_[A,M]. The six Hessians are
therefore -2Vol(S³)omega K_rad[A,B]_ij/3 in the stated J_ij convention.
Their simultaneous zero is exactly [A,B]=0. All four boost pairings vanish
by spatial antipodal parity, even though boosts leave this fixed quadratic
subspace. Thus the ten-test commuting criterion is correct **within that
restricted common-profile family**. It is not an all-sector admissibility
classification or a sufficient nonlinear constraint condition. The revised
note properly calls A=aC, B=bC a degree-two specialization of the older
same-real-harmonic family. Orthogonality of the displacement and velocity
harmonics prevents their separately established contact traces from cancelling
as a complete angular trace.

## 6. Corrections incorporated and verification limits

Two precision issues in the first Noether draft were corrected before this
review was sealed. Vanishing of a straightening variation merely near the cut
does not remove the bulk i_chi L term. The revised comparison chooses the
free interior straightening variation and needed jet to vanish near the
**entire** initial Cauchy section, possible because the initial shell jet is
zero. Alternatively an interior-supported extension difference could be
removed using the proved gauge identity. The revised notation also defines
Omega_future=-Omega_rel and restricts its equality with the frozen physical
pairing to the compact-support states actually compared. No new arbitrary
after-contact component rederivation of the old odd cut current is claimed.

The independent constraint draft's one-copy/doubled notation was clarified.
Its desired nonzero obstruction witness was moved out of the deliberately
wrong-formula roster. Existing normal/optimized receipts now report 46
positive checks and 9 wrong controls; the Noether checker reports 6 identities
and 8 wrong controls. The finite symbolic tests support the explicit algebra,
including unrestricted normal/shift variations and the tilted-cut control.
They do not machine-certify the differential-form proof, exact-matching
configuration space or nonlinear PDE existence. I read the scripts and
checked the existing receipt equivalence/pins; I did not rerun them.

The selective literature note properly treats smooth external-boundary and
edge-mode methods as structural prior art, includes the physical-shell
Hamiltonian precedent and its erratum, and distinguishes primary readings
from historical metadata. This critic did not reread all third-party source
bodies. No external novelty or exhaustive coverage claim is supported.

No unresolved blocking inequality, orientation, factor, embedding or
differentiability defect was found after the documented corrections. The
result remains a necessary classical obstruction, with global evolution,
constraint-range sufficiency, all-sector completion, higher orders and
cosmological interpretation unresolved.
