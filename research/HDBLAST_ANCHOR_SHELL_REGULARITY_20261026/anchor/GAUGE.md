# Conditional second-order gauge and moving-target corollary

This is a formal, local statement about a smooth one-parameter metric family
and auxiliary data. It is not an existence theorem for that family. Assume a
fixed reference manifold and compatible pullback identifications. As in the
merged `HDBLAST_OBSERVER_FRAME_STRUCTURE_20261026` proof, let C0 = 0 and let
the exact registered first metric have the correlated master-plus-Lie form
`g1 = -4uY dy² + 2uY rho² h + L_zeta g0` on an open region. This gives the
special first-Weyl structure and vanishing three-screen magnetic projection
H1 for the background pair U0 tangent to a y-leaf and n0 normal to it. Let
all perturbations be sufficiently differentiable for the following expansions.

## Pullback covariance of the frame component field

Let A_lambda denote the entire augmented physical data, including g, an
oriented orthonormal frame (U,n,e_a), any shell, section, clock, scalar,
transport prescription, and labels used to choose it. With P the projector
onto the common screen and epsilon^(3) its consistently oriented three-form,
define the covariant spacetime tensor (using normalized symmetrization)

`H_mn = (1/2) P_m^a P_n^b epsilon^(3)_(a^{cd} C_{b)Ucd}`.

The three-form projects c,d into the screen. For a smooth oriented screen
frame e_a, define the local scalar-valued component functions
`q_ab(A_lambda)=e_a^m e_b^n H_mn`. The H_mn are spacetime tensor components;
the q_ab are scalar-valued functions only once the frame fields and their
labels have been specified. A different frame gives different q_ab in general. If all
geometric fields *and frame vectors* are pulled back by the same local
family of diffeomorphisms Phi_lambda, tensor naturality gives

`H(Phi_lambda* A_lambda) = Phi_lambda* H(A_lambda)`,
`q_ab(Phi_lambda* A_lambda) = Phi_lambda* q_ab(A_lambda)`.

Suppose Phi_0 = id and expand pullback on any scalar field as

`Phi_lambda* = 1 + lambda L_xi + (lambda²/2)(L_eta + L_xi²) + o(lambda²)`.

Here xi is the first-order generator. In a local coordinate chart define
`Phi_lambda^mu(x) = x^mu + lambda xi^mu(x)
 + (lambda²/2)[eta^mu(x)+xi^nu partial_nu xi^mu(x)] + o(lambda²)`;
this specifies eta and fixes the convention. If q = q0 + lambda q1 + (lambda²/2)q2 + o(lambda²), the
second derivative in the new identification is

`q2' = q2 + 2 L_xi q1 + (L_eta + L_xi²) q0`.

In the registered background and its exact first sector, q0=q1=0 as fields
on the comparison neighborhood, so `q2'=q2`. The factor 2 is required by
the second-derivative convention. A change in first-frame choice or a moving
frame must not be silently neglected; it is either pulled back here or handled
by the next paragraph. The sign of xi depends on which direction Phi is
defined, but the zero conclusion is unaffected.

An independently chosen smooth oriented orthonormal frame for the transformed
metric need not equal the pulled-back original frame. If their background frames
agree at each event, their first-order difference is an arbitrary infinitesimal
Lorentz transformation. The merged same-event structure lemma shows that its
contribution to q2 vanishes. Thus, *for this particular first metric*, the
second-order result survives both a perturbative diffeomorphism and arbitrary
first-order changes of the frame prescription, provided the background frame
is the same. This does not follow from q1=0 alone: a generic algebraic Weyl
counterexample with q1=0 but nonzero first frame derivative is in the merged
control. It also does not apply to an unrelated background observer.

## Moving points and target spheres

Let p_lambda be a smooth physical event curve with p_0=p, expressed in one
common reference identification. For the same smooth frame-component field,

`d²/dlambda² [q_ab(lambda,p_lambda)]|0
 = q2_ab(p) + 2 v^mu partial_mu q1_ab(p)
   + a^mu partial_mu q0_ab(p) + v^mu v^nu partial_mu partial_nu q0_ab(p)`.

If q0=q1=0 on a neighborhood, all displacement terms vanish. For a fixed
smooth frame with the registered background observer and normal, the merged
lemma gives H0=H1=0 as projected spacetime tensors, hence also q0=q1=0. In particular,
changing the clock, scalar level, shell embedding, angular label or radial
transport may move the target by O(lambda), but cannot alter q2 at the common
background point. This statement requires smoothness and a common limit p;
it does not assert that the finite-lambda targets coincide.

More generally, for smooth embeddings i_lambda:S³ -> M with a fixed identified
background embedding i0, the geometric pullback of the covariant spacetime
tensor, `i_lambda* H_lambda`, has second parameter derivative `i0* H2`:
every variation
of i_lambda multiplies H0 or H1. If a scalar projection uses a smooth family
of basis functions/tensors T_lambda, measures mu_lambda, and a smooth nonzero
normalization with the **same** zeroth-order choices, each variation of T,
mu, or normalization likewise multiplies H0 or H1. Its leading lambda²
coefficient is therefore unchanged. A sphere-wide comparison requires this
pointwise background alignment, orientation, basis, and measure on the whole
sphere, not merely at one point. The abstract S³ needs an atlas or globally
defined geometric tensor fields; no single angular coordinate chart covers the
sphere. Frame-component arguments can be made chart by chart with consistent
transitions, or with a separately established global frame. The projection basis
and measure must be genuine globally defined objects (or compatible local ones)
for the sphere-wide integral. Changes of zeroth-order angular convention or
orientation change components/signs and must be accounted for.

## Conditional application to the proposed clock

The earlier anchoring note defines a candidate clock from a marked shell
section and intrinsic shell proper time, transported in a collar along a
nonvanishing spacelike scalar-gradient normal, provided the flow intersects
the shell transversely and uniquely and the resulting clock gradient is
timelike. When the augmented marked data and both transport constructions are
specified and transform covariantly, this clock's unit observer and scalar
normal have the registered background pair. If a smooth screen frame and
smooth target sphere have the registered background identification, the two
arguments above imply that the leading magnetic q2 is the registered one,
without computing the first variations of that observer or of the target
location. The earlier numerical interval can transfer only under all its
original response/error premises as well as these geometric premises, and
with the same orientation, normalization, and projection. In particular the
small amplitude is not a physical force or observed radiation.

This is a conditional comparison, not a proof that the collar exists in the
full nonlinear problem, that a scalar-gradient lower bound holds, that a
physical observer has been measured, that the selected initial section is
unique, or that any global gauge-independent observable has been constructed.
A frame or event rule singular at lambda=0 need not satisfy the hypotheses.
The complete compatible second-order Einstein–matter bulk/shell solution and
finite-amplitude exact spacetime remain open.

## Negative sanity checks

1. If q1 is not identically zero, moving events can matter: let
   `q_lambda(x)=lambda*x` and `p_lambda=p+lambda*a`. Then the second derivative
   changes by `2a`, matching `2 L_xi q1`.
2. If only q1=0 but first Weyl lacks the special structure, an independently
   chosen first-order frame can change q2 (see merged exact counterexample).
3. If two prescriptions disagree already at zeroth order, the infinitesimal
   same-frame lemma does not compare them. Any statements about distinct
   background observers require a separate calculation.
4. A geometric physical claim requires specifying the extra marked initial
   section, clock, target and transport; a perturbative pullback identity alone
   does not choose or observe any of them.
