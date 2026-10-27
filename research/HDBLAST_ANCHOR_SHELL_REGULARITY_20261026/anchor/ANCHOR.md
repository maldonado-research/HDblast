# Conditional shell-clock event anchor and the leading screen tensor

This is an explicitly conditional continuation of the published registered
screen-curvature and frame-structure calculations. It is a kinematic construction,
not a proof of a second-order Einstein–matter solution, a preferred cosmic
observer, global gauge invariance, or the higher-dimensional Big Bang hypothesis.

## Data and local construction

Let a sufficiently smooth one-parameter family of oriented, time-oriented
five-dimensional metrics, scalar fields and timelike shells be supplied. In
addition, supply a smooth marked spacelike section of each shell, a consistent
orientation, and three independent scalar angular labels on a local patch of
that section; these markings are physical preparation data in the sense that
they transform along with the fields under a diffeomorphism. The background
is `g0=dy²+rho(y)²[-dt²+cosh²(t) kappa]`, `phi0=phi0(y)`, the shell is at
`y=ys`, and the marked section is its `t=0` slice. The background labels
are a chosen local coordinate chart on the `S³` section.

On the shell use the future unit timelike intrinsic normal to the marked
section to generate intrinsic Gaussian normal coordinates until their first
caustic. Let `sigma` be the shell proper time and `alpha_s^A` be the three
labels transported constantly along that congruence. Set `tau_s=sigma/rho_s`
with the fixed background normalization `rho_s=rho(ys)>0`. Thus
`tau_s,0=t` and `alpha_s,0^A=Omega^A` on the shell.

In a shell-to-target tube, suppose the physical scalar gradient is spacelike
and nonzero and let `n_phi` be its unit vector oriented to agree with `+∂y`
at the background. Suppose each integral curve through the tube meets the
shell transversely exactly once, smoothly in the data, without flow crossings.
Extend `tau_s` and the `alpha_s^A` by

```
n_phi(tau)=0,    tau|shell=tau_s,
n_phi(alpha^A)=0, alpha^A|shell=alpha_s^A.
```

Assume `d tau` is timelike. Then

```
U_geo = -grad(tau) / sqrt(-g^-1(d tau,d tau)),
n_geo = n_phi,
```

are smooth orthogonal unit vectors. Because `n_phi(tau)=0`, the common screen
is tangent to the surfaces of constant `phi` and `tau`; the labels give local
coordinates on those surfaces. A smoothly oriented orthonormal screen frame
may be constructed from the label-coordinate tangent vectors and the screen
metric, fixing its background orientation to match the registered frame. A
local orthonormalization convention can be used in each patch. An intrinsic
shell-geodesic congruence and `U_geo` need not agree at the perturbed shell
unless further compatibility is imposed.

The local existence assertions are assumptions about the family, not verified
bounds on the physical solution. In particular, they require the shell chart,
scalar flow and their inverse to exist through the target, transversality,
`g^-1(dphi,dphi)>0`, `d tau` timelike, and enough differentiability. The
registered frame theorem also requires the *exact correlated first metric*
`g1=-4uY dy²+2uY rho² h_dS4+L_xi g0`; it is not an arbitrary first metric.

## Inherited background check at the registered target

The published retained input `ACTUAL_POTENTIAL_JET_BOUNDS.json` gives, conditional
on its authenticated regular-branch theorem, `-7/50000 <= p=phi0,y <= -3/25000`
and `129 <= rho <= 131` on `0<=t<=3 ell`, `z>=0`, `t+z<=8 ell`, with
`ell=10^-6` and `dy=-rho dz`. The target `(x,r)=(1/4,5/4)` and its fixed-time
background shell-to-target segment lie within that domain: `t=ell/4`,
`0<=z<=5 ell/4`, so `t+z<=3 ell/2 < 8 ell`. Thus the inherited input gives
`|p|>=3/25000` on that background segment and the background radial map is
nonsingular. This is an application of a previously certified *conditional*
enclosure, not a new verification of the regular branch. It does not bound
the perturbed gradient, prove a noncaustic transport tube, or establish that
a smooth metric family exists. `checks/check_inherited_background.py` pins
the cited input and verifies the domain arithmetic; it does not re-prove its
numerical bounds.

## A geometrically identified event

Define the five scalar labels `Z=(phi,tau,alpha^1,alpha^2,alpha^3)` and fix
numerical label values `(phi0(y*), t*, a*^1,a*^2,a*^3)`. If
`p=partial_y phi0(y*) != 0`, then at the background

```
Z0=(phi0(y),t,Omega^1,Omega^2,Omega^3),
D Z0 = diag(p,1,1,1,1)
```

in the chosen angular chart. The implicit function theorem therefore gives a
unique smooth local event `X_lambda` with `Z_lambda(X_lambda)=Z0(X0)` for all
sufficiently small lambda, subject to the tube assumptions above. This local
statement is not a bound on the size of the admissible interval and does not
provide a globally unique target sphere.

With `Z_lambda=Z0+lambda Z1+O(lambda²)` and
`X_lambda=X0+lambda X1+O(lambda²)` in the registered coordinate chart,
implicit differentiation gives

```
X1^y = -phi1/p,       X1^t = -tau1,       X1^A = -alpha1^A,
```

all evaluated at `X0`. Here `phi1`, `tau1` and `alpha1` are Eulerian
first variations of the already constructed scalar fields. The shell
embedding displacement is included in their boundary conditions; it must not
be silently discarded. For example, if the boundary-clock perturbation is
pulled back by the shell embedding, the effective fixed-coordinate datum is
`beta=beta_pull-zeta^t`, since `tau0=t`; see the prior anchoring note.

## Leading screen-curvature response

Distinguish the spacetime tensor from its frame-labelled functions. With
`P^mu_nu=delta^mu_nu+U^mu U_nu-n^mu n_nu`, let `epsilon^(3)` be the
oriented three-form induced on the screen from the spacetime volume form by
`U,n`, with the prior screen orientation. Define the covariant spacetime
screen tensor by projecting its free indices,

```
H_mn = (1/2) P_m^a P_n^b epsilon^(3)_(a^{cd} C_{b) U c d}.
```

The indices `c,d` are screen indices; the three-form already projects them.
For a chosen oriented orthonormal screen frame `e_a`, its local scalar-valued
functions are `q_ab=e_a^m e_b^n H_mn`. A different frame gives different
component functions in general. Sphere embedding pullbacks below apply to
`H_mn`, not to an individual frame component treated as a spacetime tensor.

Write `C_lambda=lambda C1+lambda² C2/2+...`. The registered background Weyl
tensor and complete first screen tensor vanish throughout the smooth
background region, so `H0=H1=0` and also `q0=q1=0` as local functions.
At each `x`,
`U_geo,0=U_reg,0=rho^-1 ∂t` and `n_geo,0=n_reg,0=∂y`; their first relative
frame change is a general infinitesimal Lorentz transformation. The published
frame-structure theorem, under the exact first-metric hypothesis above, then
gives equality of the complete local second derivatives in matched frames at
the same fixed event:

```
H_geo,2(x) = H_reg,2(x).
```

The physical target is moving in the registered chart. For any smooth
component function `q_lambda(x)=q0(x)+lambda q1(x)+lambda² q2(x)/2+...`,
Taylor expansion gives

```
d²/dlambda² q_lambda(X_lambda)|0
  = q2(X0) + 2 X1^mu partial_mu q1(X0)
    + X2^mu partial_mu q0(X0)
    + X1^mu X1^nu partial_mu partial_nu q0(X0).
```

All additional terms vanish because `q0` and `q1` vanish *throughout a
neighborhood*, not merely at the point. The same calculation holds in the
smooth matched orthonormal frame: frame variations multiplying lower-order
zero screen tensors vanish, and the special first-Weyl structure controls
relative frame changes. Consequently the anchored event and frame yield the
same leading local `H2` components as the registered calculation, with the
same background frame and orientation.

For a smooth family of sphere embeddings `i_lambda`, the geometric pullback
`i_lambda^* H_lambda` has second derivative `i0^* H2`: variations of the
embedding multiply `H0` or `H1`. To transfer a sphere-wide projected
coefficient or its previously bounded sign, these local hypotheses and
comparisons must hold at every point of the entire target sphere. Its atlas must cover the sphere; label overlaps,
background angular identification, oriented frame transformations, projection
basis, and integration measure must be consistently matched. First-order
changes in a smooth basis or measure multiply `H1=0`, so they do not contribute
to the leading second-order coefficient. These facts do not prove the needed
global perturbed collar, establish new perturbed estimates for `p`, or choose
the physical markings. The inherited coefficient's interval remains conditional on its own
source and response hypotheses.

## Covariance under changes of identification

Let `psi_lambda` be a smooth local diffeomorphism with `psi_0=id`, and define a
second description by pullback of the *entire augmented configuration*: metric,
scalar, shell, marked section, labels, orientation and transport data. By
uniqueness of geodesic and flow constructions, the newly built scalar labels
satisfy `Z'_lambda=psi_lambda^* Z_lambda`. The new event with the same numerical
labels is therefore

```
X'_lambda = psi_lambda^(-1)(X_lambda).
```

Writing `xi=(d psi_lambda/dlambda)|0`, the first variations obey
`Z1'=Z1+xi^mu partial_mu Z0` and `X1'=X1-xi(X0)`. The latter also follows by
substituting the transformed labels into the implicit-derivative formula
above. An equivariantly constructed oriented frame transforms with the pullback,
so its evaluated curvature components at the labelled event agree in both
descriptions. This is ordinary covariance for the augmented data and valid
tube; it is not invariance after changing, omitting or fixing those physical
markings by hand.

Equivalently, for the full spacetime screen tensor constructed covariantly
from the augmented `g,phi,tau` and their normalizations, the perturbative
pullback formula has `H2' = H2 + 2 L_xi H1 + (terms acting on H0)`.
Here `H0=H1=0` as tensor fields in the background neighborhood, so `H2'=H2`
under those changes of identification. This conclusion concerns this specific
augmented observable and the stated background data; it neither makes the
full Weyl tensor's `C2` gauge invariant nor chooses a unique observer from the
unmarked symmetric spacetime. The construction and conclusion do not establish
an exact nonlinear solution, matter evolution or observational relevance.
