# A rotational second-order obstruction for the two-interior shell

14 October 2026 UTC. Independent local derivation, confined to the authoring
directory. The public capsules and private archive are unchanged. This note
derives a spatial-rotation necessary constraint condition. It does not assume
the unrestricted nonlinear covariant shell Noether identity, prove nonlinear
existence, or infer integrability from zero quadratic candidates.

## 1. Statement and differentiability assumptions

Take the actual registered Einstein–scalar background

\[
 g=\Omega(s)^2\eta,\quad \phi=\Phi(s),\quad
 s=R^2-T^2,\quad \Omega>0,\quad f=\Phi_s\ne0,
\]

with two interiors and their single common timelike shell. The action is

\[
 S=\sum_{i=\pm}\left[\int_{M_i}\epsilon_g
 (R/2-(\nabla\phi)^2/2-U)+\int_\Sigma\epsilon_h K_i\right]
 -\int_\Sigma\epsilon_h\sigma(\phi).
\]

There are two GHY terms and one shell action. Each shell normal is the
outward normal of its own interior. No identification of these normals as
one equal directed normal on a smooth manifold is used.

Assume a hypothetical parameter-C² family of exact solutions through this
background admits the following initial-data description:

1. Its complete Cauchy section consists of two compact four-balls, with
   one common closed three-sphere cut of the timelike shell, no external
   physical boundary, and the smooth cone/vertex retained as internal geometry.
2. Smooth parameter-dependent identification maps pull the two moving balls
   and common cut back to fixed reference balls and a fixed reference cut.
   They preserve induced metric/scalar matching and the common cut chart.
   Bulk and shell embeddings are part of the family, rather than fixed
   external boundaries. Independent tangential maps on the two cut copies
   are excluded by matching.
   For the compact first-order states below, the identification maps can
   have identity value and zero first parameter derivative at the initial
   slice: the first-order shell position, velocity and geometry already
   vanish near the cut. Thus the stated tangent v is preserved rather than
   altered by an unspecified first-order coordinate change.
3. On each reference ball the spatial metric, scalar, lapse, shift and their
   first time derivatives, and their shell traces, are C² in the parameter
   in a topology controlling the products, traces and integrations below.
   A sufficient classical assumption is piecewise spacetime C⁴ regularity,
   C² parameter dependence in that norm on a compact slab, C⁴ embeddings,
   timelike shell and uniformly spacelike Cauchy section. Less regular
   Sobolev assumptions can replace this only after establishing the needed
   composition, trace and product estimates.
4. All exact bulk initial constraints and all Israel/scalar shell junctions
   hold; the common induced metric/scalar match exactly. Arbitrary
   second-order corrections are allowed in every gravitational/scalar
   sector and at the shell. No first-order support assumption is imposed
   on these corrections.

Under these hypotheses, every tangent state v whose first-order Cauchy data
and embeddings vanish near the shell must obey

\[
 D^2\mathcal C_J|_0[v,v]=0
 \quad\Longleftrightarrow\quad
 \Omega^{\rm phys}_0(v,\mathcal L_Jv)=0                 \tag{1}
\]

for any background spatial rotation J, using the physical symplectic
orientation of the published Cauchy formula. The equivalence and its sign
are derived below. There are explicit smooth finite-angular master states
violating (1), so they cannot be tangents to such exact families. This is a
linearization-instability example in the specified differentiability class.
It is not a statement about exponential growth, a scalar nonlinear closure,
or families in an uncontrolled weaker topology.

## 2. Exact ADM momentum and shell geometry

On each Cauchy ball let q_ab be the spatial metric, u the future unit normal,
and k_ab=q_a{}^A q_b{}^B nabla_A u_B. The R/2 normalization gives

\[
 \pi^{ab}=\frac{\sqrt q}{2}(k^{ab}-kq^{ab}),\qquad
 \varpi=\sqrt q\,u\phi,\qquad
 \mathcal H_a=-2D_b\pi^b{}_a+\varpi D_a\phi.             \tag{2}
\]

Here pi and varpi are spatial densities. Gauss–Codazzi gives
H_a=-sqrt(q)(G-T)_{ua}; thus the full normal Einstein constraint, rather
than its trace-reversed replacement, is used.

Let s_i be the outward spatial unit normal to the cut within ball i, n_i
the outward spacetime shell normal, and tau the common future unit normal
to the cut within the shell. Their boost angle eta_i is defined by

\[
 n_i=\cosh\eta_i\,s_i+\sinh\eta_i\,u_i,\qquad
 \tau=\cosh\eta_i\,u_i+\sinh\eta_i\,s_i.              \tag{3}
\]

The cut metric c and tau are common because of exact induced matching.
For a common vector J tangent to the cut, write K_i for shell curvature
with outward n_i. Symmetry of both extrinsic curvatures, together with
derivatives along J of the orthonormality relations, gives the exact identity

\[
 \boxed{k_i(s_i,J)=-K_i(\tau,J)-J\eta_i.}                \tag{4}
\]

For example invert (3) to u=cosh(eta)tau-sinh(eta)n and
s=cosh(eta)n-sinh(eta)tau. Then
k(s,J)=g(nabla_J u,s). The differentiated boost coefficients contribute
-J eta; the remaining terms give -K(tau,J). This derivation permits arbitrary
shell tilt and does not set the shell velocity or second-order correction
to zero. At an orthogonal cut eta=0 it reduces to k(s,J)=-K(tau,J).

The complete Israel residual is

\[
 \mathcal I_{ab}=\sum_i(K_{i,ab}-K_i h_{ab})+\sigma h_{ab}.
\]

Consequently, for every configuration with matching,

\[
 \mathcal I(\tau,J)=\sum_iK_i(\tau,J),                 \tag{5}
\]

because h(tau,J)=0. The trace terms and the single pure-tension shell stress
have zero mixed projection exactly. The scalar-dependent tension has no
independent scalar kinetic momentum; its force and scalar junction are
separate residuals, not extra tangential angular stress.

## 3. Full spatial angular constraint, including the interface

Choose J fixed on both reference balls, with the same tangential restriction
on their common cut. It is a background rotation but need not be a symmetry
of any perturbed datum. Integration by parts on each ball gives

\[
 \sum_i\int_{C_i}\mathcal H_{i,a}J^a
 =M_J-\sum_i\int_{\mathcal S}\sqrt c\,k_i(s_i,J),
\]
\[
 M_J=\sum_i\int_{C_i}
       (\pi_i^{ab}\mathcal L_Jq_{i,ab}+\varpi_iJ\phi_i).
\]

Inserting (4)–(5) yields

\[
 \sum_i\int_{C_i}\mathcal H_{i,a}J^a
 =M_J+\int_{\mathcal S}\sqrt c
      [\mathcal I(\tau,J)+J(\eta_++\eta_-)].            \tag{6}
\]

Thus the explicitly derived combined spatial constraint functional is

\[
 \boxed{\mathcal C_J:=
 \sum_i\int_{C_i}\mathcal H_{i,a}J^a
 -\int_{\mathcal S}\sqrt c\,\mathcal I(\tau,J)
 =M_J+\int_{\mathcal S}\sqrt c\,J\eta_{\rm sum}.}       \tag{7}
\]

The interface subtraction is fixed by the integration and Israel projection;
it is not an ad hoc extra charge. Both outward GHY completions supply their
own normalization in (2), and their corner angles occur once each in
eta_sum. The one common shell contributes its one Israel stress, whose
mixed contraction is zero. The last term in (7) is the geometrical corner
angular momentum. It must be retained when the cut is tilted. No unrestricted
identification of the full covariant corner potential is needed to prove
(7) or the compact-supported Hessian calculation below.

Every exact constraint-and-junction datum obeys C_J=0, whether or not its
surface data vanish. For an orthogonal foliation the corner term is zero;
(7) also covers a general C² foliation. A cut with an external boundary, an
artificial cone boundary, unmatched fields or a source carrying extra shell
angular stress would change (6) and this conclusion.

The bulk diffeomorphism and common worldvolume reparametrization are combined
in the pullback chart. If X_i are parametrized shell embeddings, the generator
is G_{J,j0}X_i=-J composed with X_i+(X_i)_*j0, where j0 is the common
background restriction of J. It fixes the background embeddings. On normal
graphs its first-order action is J zeta_i. Keeping j0 fixed avoids any
field-dependent compensator terms. This is the symmetry action meant in (1).

## 4. Why an arbitrary second-order shell correction cannot cancel

At the background T=0 slice,

\[
 q_{0,ij}=\Omega(r^2)^2\delta_{ij},\quad
 \pi_0=\varpi_0=0,\quad \eta_{i0}=0,
\]

and Jq0=Jphi0=J sqrt(c0)=0. Therefore DM_J|0[w]=0 for every
second-order correction w. The remaining corner derivative is

\[
 D\mathcal C_J|_0[w]_{
 \rm cut}=\int_{\mathcal S}\sqrt{c_0}\,J\eta_{\rm sum}[w]=0,
                                                               \tag{8}
\]

by integration on the closed cut and rotational invariance of its background
density. Equation (8) permits an arbitrary, nonzero second-order angle and
induced geometry. It is the explicit adjoint cancellation at the shell.

For data A(epsilon)=A0+epsilon v+epsilon²w/2+o(epsilon²), twice
differentiating the exact equation C_J(A(epsilon))=0 gives

\[
 0=D\mathcal C_J|_0[w]+D^2\mathcal C_J|_0[v,v]
   =D^2\mathcal C_J|_0[v,v].                              \tag{9}
\]

The corner Hessian is 2 int delta sqrt(c) J delta eta_sum. It is zero for
the admitted v because v and its normal/time derivatives vanish on an open
initial shell neighborhood. Moving-domain, graph-evaluation and second-order
matching effects are represented by the pulled-back w and by (7); their
linear contribution has been cancelled by (8), not omitted by support.
The remaining Hessian is

\[
 D^2\mathcal C_J|_0[v,v]=2\sum_i\int_{C_i}
 [\delta\pi_i^{ab}\mathcal L_J\delta q_{i,ab}
                   +\delta\varpi_i J\delta\phi_i].       \tag{10}
\]

With Omega_can=delta pi wedge delta q+delta varpi wedge delta phi,
integration of J on compact-supported data gives
D²C_J[v,v]=Omega_can(v,L_Jv). This fixes the factor two in (10).
The Taub convention using the coefficient of epsilon² instead of its second
derivative would use half of (10).

## 5. Direct ADM reconstruction and the physical orientation

Use the proved ambient reconstruction

\[
 h_{AB}=2F(s\eta_{AB}-3X_AX_B),\qquad
 \chi=-\frac{3(E+4)F}{2\Omega^2f},\qquad E=X^A\partial_A.
\]

Write master initial data q=F|0 and p=F_T|0, and D=r partial_r.
At T=0, directly from the metric and lapse/shift,

\[
 \delta q_{ij}=2q(r^2\delta_{ij}-3x_ix_j),\quad
 \delta N=r^2q/\Omega,\quad\delta\beta_i=0,
\]
\[
 \delta k_{ij}=\frac p\Omega(r^2\delta_{ij}-3x_ix_j),\quad
 \delta k=\frac{r^2p}{\Omega^3},\quad
 \boxed{\delta\pi^{ij}=-\frac{3p}{2\Omega}x_ix_j.}       \tag{11}
\]

Also

\[
 \chi=-\frac{3(D+4)q}{2\Omega^2f},\quad
 \chi_T=-\frac{3(D+5)p}{2\Omega^2f},\quad
 \delta\varpi=\Omega^3\chi_T.                            \tag{12}
\]

The extra unit in (D+5) is [partial_T,E]=partial_T; the background
coefficients have zero T derivative on this slice.
Equations (11)–(12) independently give

\[
 \int_{C_i}[\delta\pi^{ij}\mathcal L_J\delta q_{ij}
                 +\delta\varpi J\chi]=\mathcal K(Jq,p),
\]
\[
 \mathcal K(q,p)=\int dr\,d\Omega_3
 [m_r(D+4)q(D+5)p+n_rqp],
 \quad m_r=\frac{9r^3}{4\Omega f^2},\quad n_r=\frac{6r^7}\Omega.
                                                               \tag{13}
\]

The published one-copy future-Cauchy pairing is
Omega_1(v1,v2)=K(q1,p2)-K(q2,p1). The full Omega_phys is the sum of
the two copy pairings on this compact initial slice, hence 2 Omega_1 for
equal reflection parity. Direct calculation shows
Omega_can=-Omega_phys for that stated orientation. Thus, for the compact
initial states considered here,

\[
 \boxed{D^2\mathcal C_J|_0[v,v]
       =-\Omega^{\rm phys}_0(v,\mathcal L_Jv).}           \tag{14}
\]

This is an explicit constraint-Hessian identification with its action
normalization, not just a proposed symplectic candidate. Reversing the global
symplectic orientation reverses (14) but leaves its necessary zero condition.
The prior current formula passes this ADM check.

## 6. A real finite-angular nonintegrable state

Let H be a real S³ scalar harmonic with JH nonzero, Q a nonzero smooth real
profile compactly supported in a sufficiently thin positive-radius annulus
strictly inside the shell, and omega a nonzero real constant with the units
of inverse ambient time. Set

\[
 q=QH,\qquad p=\omega QJH.                               \tag{15}
\]

Since J is divergence-free on the sphere,
int H J²H=-int(JH)². No symmetry of the radial bilinear K_rad is assumed:
both terms use the same diagonal K_rad(Q,Q). Therefore

\[
 \Omega_1(v,Jv)=-2\omega K_{\rm rad}(Q,Q)\int(JH)^2,
\]
\[
 \boxed{D^2\mathcal C_J[v,v]
 =4\omega K_{\rm rad}(Q,Q)\int(JH)^2\ne0}              \tag{16}
\]

for either fixed reflection parity on the two interiors. Each contributes
the same value: the simultaneous sign reversal of q and p in the odd copy
does not reverse a quadratic form. Mixed reflection pairings cancel, which
does not affect this same-parity witness.

Compact-support integration by parts gives

\[
 K_{\rm rad}(Q,Q)=\int_I
 [m_r\{(DQ)^2+(2+9r^2L)Q^2\}+n_rQ^2]dr,
 \quad L=\Omega_s/\Omega+2f_s/f.
\]

Define m*=min_I(m_rr²)>0, M*=max_I m_r(-2-9r²L)_+, and w_I the
annulus width. Dirichlet Poincare gives

\[
 K_{\rm rad}(Q,Q)\ge(\pi^2m_*/w_I^2-M_*)\int_IQ^2dr>0
                                                               \tag{17}
\]

when the annulus is sufficiently thin. Regularity and f nonzero guarantee
such annuli in every controlled collar. This is positivity of a diagonal
auxiliary bilinear, not a positive gravitational energy assertion.

An explicit degree-two example is

\[
 H=n_1^2-n_2^2,\quad J=x_1\partial_2-x_2\partial_1,
 \quad JH=-4n_1n_2,\quad J^2H=-4H,
\]
\[
 \int H^2=\pi^2/3,\qquad\int(JH)^2=4\pi^2/3.
\]

Then the doubled Hessian is 16 pi² omega K_rad(Q,Q)/3. Normalizing H
to unit L² norm gives doubled Hessian 16 omega K_rad(Q,Q).

These smooth data belong to the proved finite-angular master construction.
They satisfy all first-order initial constraints and shell compatibility;
the reconstruction and admitted boundary laws produce the corresponding
linearized shell solutions on the controlled slab. At least a pre-contact
exact linearized family follows directly from finite propagation. The
registered finite-slab boundary theorem supplies after-contact continuation
when its support/time hypotheses hold. No assertion that every arbitrary
profile makes contact is needed for the nonintegrability result. The l=2
metric projection argument in the public junction proof independently
excludes smooth pure gauge for the nonzero state.

By contrast q=QH,p=PH has Omega(v,Jv)=0 for every rotation, and its boost
candidates vanish by antipodal parity. Those zeros remain necessary tests
only; they provide no nonlinear existence result. The larger finite-angular
linear data image contains the explicit obstructed state (15), and therefore
cannot be treated as an automatically integrable nonlinear tangent space.

## 7. Additional nonlinear junction consistency control

Let E_i=G_i-T_i, Delta K_ab=K_+ab-K_-ab, Delta p=p_+-p_-, and
S=sum p_i+sigma'. Exact matched h,phi give

\[
 \sum_i E_{i,na}=D_b I^b{}_a-S D_a\phi,
\]
\[
 \boxed{\sigma\,\Delta K+\sigma'\Delta p
 =2(E_{+,nn}-E_{-,nn})+I_{ab}\Delta K^{ab}+S\Delta p.}    \tag{18}
\]

The second identity follows by subtracting the two spacelike-shell Gauss
constraints and substituting sum K_ab=I_ab+(sigma-tr I)h_ab/3 and
sum p=S-sigma'. It retains the full quadratic corrections. These identities
explain why the embedding force/tangential conservation residuals add no
independent on-shell correction capable of invalidating (7). They are not
used to discard any off-shell term without checking its residual.

## 8. Verification and limits

The independently written `check_rotational_constraints.py` imports no
scientific implementation from the frozen inputs. It checks the Lorentz
normal geometry, all 16 ADM spatial matrix entries, gravitational and scalar
Cauchy normalization, commutator shift, moment-map Hessian factor, exact
degree-two harmonic/sphere moments, radial integration-by-parts coefficient,
nonlinear normal-force identity and a nonzero tilted-cut control. Deliberate
wrong controls detect an omitted corner angle, wrong boost sign, wrong
gravitational sign, lost Hessian factor two, lost D+5 shift and wrong force
residual. Normal and optimized runs write portable JSON; explicit injected
failure must exit nonzero in both modes. These checks establish algebraic
identities. The functional differentiability and PDE arguments are analytic
proof assumptions, not machine-certified existence theorems.

The rotational constraint proof establishes a real second-order obstruction
under the stated exact-family hypotheses, even if unrestricted boost/normal
shell covariant identities remain open. It neither classifies all KIDs nor
proves that the vanishing of the six rotation conditions, or of all ten
Lorentz candidates, suffices for integration. No originality or cosmological
claim accompanies this local result.

## Inputs

- Public `junctions/FULL_PACKET_JUNCTION_SUFFICIENCY_AND_CURRENT.md` in
  HDBLAST_REGULAR_SHELL_CONTACT_20261013, especially Sections 5–7.
- Public regular-cone ambient reconstruction and independent physical
  Cauchy-pullback review.
- Private archived `SECOND_ORDER_SHELL_CONSTRAINT_SCOPE.md` and
  `SECOND_ORDER_CHARGE_ASSESSMENT.md`, both read as unchanged inputs.
- Separate concurrent noether derivation supplied (18), cross-checked here
  from the shell Gauss constraint; its unrestricted covariant cut argument
  is outside the independent rotational theorem.
