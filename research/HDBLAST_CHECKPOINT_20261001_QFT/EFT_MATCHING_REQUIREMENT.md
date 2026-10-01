# Exploratory EFT matching requirement

This analytical extension was developed after registration. It is not a new quantum-mode result or an external novelty claim.

## The original shell action is not closed under a chi loop

The original tuned tension has a fixed cubic W and a quadratic detuning:
sigma=2(1−phi+phi³/3)+delta(1+c phi+d phi²/2).

For one real shell scalar with M²=m0²+gbar²(phi−phi_star)², the four-dimensional heat-kernel coefficient, modulo total derivatives, contains

a2 = 1/2 [M²+(xi−1/6)R]² +(Riemann²−Ricci²)/180.

At minimal chi coupling xi=0, this includes M⁴/2−M²R/6+R²/72 and the curvature combination above. Pole and counterterm signs require a specified regularization convention; the operator content is the relevant conclusion here.

M⁴ generates an independent (phi−phi_star)^4 shell potential. Its phi^4 term is absent from the original sigma basis. M²R generates intrinsic R, phi R and phi² R terms, and curvature-squared terms also occur. These change shell metric/scalar variations and junction conditions; changing sigma alone is insufficient. Absorbing lower powers into W would additionally change its imposed relation to the bulk potential.

An effective theory can include and match these operators. This is a requirement for specifying the quantum extension, not an impossibility theorem. Some curvature combinations simplify on conformally flat backgrounds or are topological, without restoring the original operator basis.

A separate divergent scalar kinetic term is not automatically forced by this one-loop scalar determinant. Finite derivative terms and other loops require separate analysis; the Yukawa condensate wavefunction counterterm cannot simply be transplanted.

## A quantitative illustration, not a vacuum prediction

For a constant background in flat-space MS-bar,

V1=M⁴/(64pi²)[ln(M²/mu²)−3/2],
J1=partial_phi V1=M² M²_phi/(32pi²)[ln(M²/mu²)−1].

Take the derivative at fixed renormalization scale mu. Evaluating that fixed scale numerically at an endpoint mass is different from differentiating a field-dependent mu=M(phi).

The local potential contributes rho=V1, p=−V1 and J=V1_phi. It is vacuum/potential energy, not produced radiation. Once absorbed into the renormalized tension it must not also be counted as a second matter source.

In the checkpoint units, Delta sigma_hat=b V1_hat and Delta sigma_model=(b/rb)V1_hat, with b=kappa5² H0³ and sigma_hat=rb sigma_model. The same relation applies to the phi derivative.

The executed code/matching_diagnostic.js uses the OLD fixed A1 endpoint phi=1.0361021913225144 and G=100, phi_star=0.5. It uses the previously saved full-history characteristic mass screen Lambda_hat=53.61022411468734 and b=Lambda_hat^(-3), rather than recomputing b from the slightly different endpoint mass. The new mode experiment ends earlier at s=6.9.

At constant mu_hat=M_hat_endpoint, the illustrative shifts are Delta sigma_hat≈−0.1273090 and Delta sigma_hat_phi≈−0.6332572. At constant mu_hat=1, they are approximately +0.5485726 and +4.4096736. Both are unretuned MS-bar pieces; changing mu also requires running/matching the renormalized coefficients, so these are not two physical predictions. Finite conditions can fix the vacuum value and force at a chosen point. They cannot remove the need for the additional operator basis.

This screen places the largest mass near the assumed gravity cutoff and supplies no hierarchy or UV certification. G=gbar/H0 is a mass ratio for dimensionless phi, not by itself a conventional dimensionless perturbativity bound.

## Tension and force enter the geometry together

For the old geometric decomposition v=rb²[sigma²/36−sigma_phi²/48+U/6], a frozen-phi, fixed-U, potential-only replacement gives

Delta v =
sigma_hat Delta sigma_hat/18
−sigma_hat_phi Delta sigma_hat_phi/24
+(Delta sigma_hat)²/36−(Delta sigma_hat_phi)²/48.

The executable illustration reports this algebraic shift. It excludes induced curvature terms, the revised junction conditions and scalar/bulk/Weyl readjustment; it is not a new cosmological solution or a predicted vacuum shift.

Finally, a constant-background Coleman–Weinberg expression cannot replace causal mode evolution during production. A large-mass derivative expansion fails at this experiment's mass-zero crossing. Absolute source calculations must specify the common effective action, finite matching, state and causal evolution; the finite in-minus-out difference leaves the reference vacuum contribution undetermined.

See LITERATURE_AND_RENORMALIZATION.md for the retrieved curved-space counterterm and state references.
