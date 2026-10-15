# Independent referee review of the conditional initial checkpoint

15 October 2026. **Accepted within the stated initial-data scope; no blocking
mathematical defect found.** This is an independent internal AI review, not
external specialist peer review. Exact candidate-file hashes are in
`REFEREE_ACCEPTANCE_PINS.json`; the review does not authenticate a subsequent
publication, ZIP, replay manifest or later file change.

The reviewed claim is well-defined: bounds for the specified physical second
parameter derivatives `(gamma2,K2,phi2,Pi2)`, their full first spatial covariant
derivatives, the local initial relative-harmonic gauge velocities through
order two, and `(k,nabla[g0]k,phi2,dphi2)` on `t=0` after that initialization.
The norms are maxima of absolute individual components in the fixed adapted
background frame, including the derivative slot. They are not a Lorentzian
contraction, an arbitrary boosted-frame supremum, or a Sobolev norm.

## Independent mathematical assessment

I rederived the inherited Hamiltonian Hessian from the radial-lapse identity
`R[gamma]=R[h]-12j^2-6L^-1 j_y-2L^-1 Delta_h L`, with
`L^2=1-4 epsilon uY` and `h=rho^2(1+2 epsilon uY) kappa`.
Its angular contribution is
`rho^-2[-336u^2Y^2+36u^2|DY|^2]`; the radial coefficients agree with the
copied constraint theorem. Expanding the scalar gradient, potential, K and
Pi terms gives its remaining source terms. Independently expanding
`D_J(K^J_I-delta^J_I tr K)-Pi d_I phi` gives the displayed F and G,
including the radial `-6u_y u_t/rho` term. Thus the source being bounded uses
physical covariant K and unit-normal Pi, not a momentum-density or
spacetime-affine source. The local inverses have the correct sector
denominators `-6,10,42` and scalar compensator `-3f_y/(2p)`.

The Hamiltonian and scalar-momentum reductions are exact identities for the
actual smooth data. In particular `v v_y` cancels from `G-F_y`. Differentiating
these identities before interval evaluation is legitimate. Normalized radial
Taylor coefficients and second parameter derivatives are treated separately;
the latter multiply `epsilon^2/2`. Both interval implementations round
outward using exact rationals, retain the fourth bump derivative needed for
`phi2_y`, and apply `D_y=-rho^-1 D_z` as a differentiated operator. Its square
is `rho^-2 D_z^2+(a/rho)D_z`, not merely a constant rescaling.

The angular weights and gradient estimates bound the whole S3. The spatial
metric derivative includes the easily missed components
`(nabla_a gamma2)_{yb}=-af delta_ab`. In contrast, metric compatibility makes
`nabla K2=d(q) tensor gamma0` for the pure-trace physical K2. The smooth compact
support follows from finite differentiation of a flat bump and division by
p only on a neighborhood where its nonzero gap is certified.

I independently checked the zero-shift contracted-Christoffel formulas. The
fixed reference connection must, and does, use the perturbed inverse metric
in `R^A=g^{BC}Gamma[g0]^A_BC`. The all-point ADM expansion then gives the stated
lapse/shift coefficients. The first harmonic representative changes its
mixed metric velocity while preserving the physical Cauchy tangent; it is
not the unchanged longitudinal spacetime one-jet.

The second conversion retains `g2_tt=-2rho^2U^2`,
`(dot g_ty)_2=dot beta2^y-8U dot beta1^y`, and the lapse cross term in
`(dot g_tt)_2`. In the fixed background frame,
`nabla_0 e0=a ey` and `nabla_0 ey=a e0` contribute `+2aU^2` to the mixed
derivative. Consequently
`(nabla_0 k)_{0y}=-26UU_y-(3/2)f_y-3af`; the residual coordinate
`-2aU^2` cancels. The scalar conversion is
`nabla_0 phi2=Pi2+2U Pi1`. Together with the spatial connection terms, the
listed component families exhaust the initial first spacetime derivative.

These calculations support the four physical value/full spatial derivative
caps `(1.06e24,4.25e29)`, `(7.96e7,1.96e13)`, `(5.13e33,3.06e39)` and
`(2.00e11,5.02e16)`, in the order gamma2, K2, phi2, Pi2. A separate coarse
exact calculation using only strict physical caps and elementary first-data
bounds gives `|nabla[g0]k|<6.4e29` and `|dphi2|<3.06e39` initially.
The dominant new metric term is the mixed time-radial derivative.

## Explicit check of the retained collar domain

The copied older background prose proves a collar of length `8e-10`, while
the later retained potential report specifies `Z=8ell=8e-6`. I independently
verified the larger first-exit enclosure from the authenticated shell box.
For the outer rectangle I, the exact interval images are those of

`F=(-rho^2 a, -rho p, -rho(U'-4ap), rho^-1+rho p^2/3)`.

Writing `J_i=I_shell,i+[0,1/125000]F_i(I)`, the independently evaluated
intervals satisfy `J_i subset interior(I_i)`. The lower and upper distances
from J to the respective faces of I are strictly larger than these positive
rational numbers:

| Component | Lower-face margin exceeds | Upper-face margin exceeds |
|---|---:|---:|
| rho | 9/10 | 107/100 |
| eta | 27/10^7 | 69/10^8 |
| p | 28/10^7 | 41/10^7 |
| a | 377/10^6 | 622/10^6 |

The exact margins and independent Horner evaluation of U' are implemented in
`check_referee_acceptance.py`. Smoothness of this ODE near the rectangle and
the retained actual shell data give the usual first-exit/continuation
argument. This resolves the apparent difference in collar lengths without
requiring a numerical trajectory or extending the original short-note claim.
The original existence/branch theorem and shell-box authentication remain
retained hypotheses of this checkpoint, not newly proved here.

## Conditional gauge consequence, literature and boundaries

The normal-gauge implication is correctly conditional. With the reduced
equation `E_AB=nabla_(A C_B)` and vanishing initial C, the Einstein residual
`H_AB=E_AB-g_AB tr(E)/2` satisfies `H_nn=X_n/2`, `H_nI=X_I/2`, where
`X=nabla_n C`. The physical constraints therefore force its first and second
parameter coefficients to vanish only when the reduced acceleration equations
are imposed. Lower-order normal/connection variations multiply already
vanishing lower-order residuals. This does not construct an exact nonlinear
constraint family or an evolved reduced solution.

The comparison with Lindblom et al. and Pretorius is appropriately limited to
the established GH initialization/subsidiary architecture. The scalar Cauchy
and conservation comparisons with Choquet-Bruhat, Isenberg and Pollack, and
Delay's localization comparison, supply no packet constants or shell theorem.
I reviewed the supplied source comparison and ledger; I did not independently
repeat their external literature search. The dimension and scalar-potential
extension needed for the initial normal identity follows from the local
algebra above. Matter conservation for subsequent propagation needs the scalar
equation; incoming subsidiary boundary conditions remain an additional issue.

The enormous upper bounds are conservative finite envelopes, not estimated
responses, lower bounds, admissible finite epsilon values or physical forces.
Flatness on an initial shell collar supplies initial compatibility, not later
shell constraint preservation. No spatial C1 bound of all gauge time
velocities, positive-time spacetime C1 bound, first gauge-wave solution,
transformed evolving source, second evolution, shell response, nonlinear
remainder or observational conclusion is accepted here.

## Evidence and portability

I inspected the derivations and relevant arithmetic/source implementations,
and replayed the background, bump, direct gauge, formal ADM, formal source,
independent norm and frame checks from the assembled candidate. Every complete
JSON result equals its registered candidate report. Point-jet witnesses and
intentional-error controls supplement the all-point proofs; their number is
not a proof of a field identity or evolution theorem.

The separate referee checker authenticates the 40 reviewed candidate files
before recomputing the exact collar and coarse frame checks. It excludes its
own directory, the final manifest, replay roster and common replay runner
from its pin set. Normal and optimized runs produce identical output:

```sh
python -I -B fresh-referee/check_referee_acceptance.py --root . --output /tmp/fresh-referee-normal.json
python -I -O -B fresh-referee/check_referee_acceptance.py --root . --output /tmp/fresh-referee-optimized.json
```

Those output paths must be new. Acceptance applies to the pinned scientific
files and the bounded claim above; the separate final capsule/publication
replay remains responsible for publication integrity.
