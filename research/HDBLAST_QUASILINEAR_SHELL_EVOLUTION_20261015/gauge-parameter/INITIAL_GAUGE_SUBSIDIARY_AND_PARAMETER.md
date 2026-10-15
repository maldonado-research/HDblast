# Initial harmonic completion, exact subsidiary constraints and a finite parameter budget

15 October 2026 UTC. This note proves the local gauge and subsidiary identities.
The parameter conclusion is conditional on the explicitly declared variable
coefficient mixed boundary existence and energy theorem below. It does not
independently supply that evolution theorem or a lifespan through contact.

## 1. Inputs, conventions and scope

Retain the exact canonical initial family of
`HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`.
Fix one integer `k>=12` and put `s=k+2>=14`. The spatial metric has class
`H^s`, the contravariant metric momentum density has class `H^{s-1}`, and
the prescribed scalar sectors have the same sufficient orders. The parameter
curve is smooth in these fixed spaces, has the selected physical first
tangent, and equals the background on one fixed open initial shell collar.
Nothing here constructs a common spatially smooth family by intersecting
infinitely many parameter intervals.

Write `u` for the future unit normal to the initial spacelike section and
`n` for the outward spacelike normal to the timelike shell. These are
different normals. Set

\[
 g=-\alpha^2dt^2+\gamma_{ij}(dx^i+\beta^i dt)(dx^j+\beta^jdt),\quad
 u=\alpha^{-1}(\partial_t-\beta^i\partial_i),\quad
 \kappa=\tfrac12\mathcal L_u\gamma .
\]

For the action normalization `R-|d phi|^2-2U`, the canonical definitions are
`pi^{ij}=sqrt(det gamma)(kappa^{ij}-tr(kappa) gamma^{ij})` and
`varpi=2 sqrt(det gamma) u(phi)`. Consequently, in four spatial dimensions,

\[
 \kappa_{ij}=\frac{\pi_{ij}-\frac13\pi\,\gamma_{ij}}{\sqrt{\det\gamma}},
 \qquad v_\phi:=u\phi=\frac{\varpi}{2\sqrt{\det\gamma}}.                 \tag{1}
\]

Indices on the density `pi` in (1) are moved with `gamma`, and
`pi=gamma_{ij}pi^{ij}`. Equivalently one may use the checkpoint's fixed
canonical-to-geometric normalization throughout; changing a canonical
normalization must change both (1) and the constraints consistently. Only
the smooth algebraic invertibility and the exact physical constraints are
needed below.

The physical Einstein residual, scalar residual and relative harmonic vector
are

\[
 E_{AB}=\operatorname{Ric}_{AB}-\phi_A\phi_B-\tfrac23Ug_{AB},\quad
 S=\Box_g\phi-U',\quad
 C^A=g^{BC}(\Gamma[g]^A_{BC}-\Gamma[g_0]^A_{BC}).                    \tag{2}
\]

The reduced equations are `E_AB=nabla_(A C_B)` and `S=0`; parentheses include
the factor one half. All statements use this convention and the Ricci
commutator convention for which the subsidiary equation below is
`Box C+Ric(C)=0`.

## 2. Exact local harmonic Cauchy completion with the chosen first gauge tangent

The lapse and shift on the initial section are freely chosen positive/small
gauge fields. To preserve the already chosen **full spacetime** first tangent,
let `alpha_1,beta_1` be the initial lapse and shift variations induced by the
previous harmonic first perturbation `q_1^wave`. Choose, for example,

\[
 \alpha_\varepsilon=\alpha_0+\varepsilon\alpha_1,\qquad
 \beta_\varepsilon=\beta_0+\varepsilon\beta_1 .                    \tag{3}
\]

These variations have support inside the initial annulus and are zero in
the protected shell collar. Shrink the parameter interval to keep the lapse
positive. At zero background shift, explicitly
`beta_1^i=gamma_0^{ij}(delta g_{tj})` and
`alpha_1=-delta g_{tt}/(2 alpha_0)`.
For a nonzero background shift, differentiate the displayed ADM metric to
obtain the unique analogous algebraic expressions. For the longitudinal
initial representation used here, `alpha_1=rho Psi` and `beta_1=0`; the
zero initial gauge vector and its first time derivative leave that initial
metric unchanged when passing to `q_1^wave`.

This differs from simply prescribing background lapse for every parameter.
The latter retains the physical canonical tangent but can change its full
spacetime gauge representative. Equation (3) removes that ambiguity.

Construct the initial metric algebraically from `(gamma,alpha,beta)` and put

\[
 \partial_t\gamma=2\alpha\kappa+\mathcal L_\beta\gamma,\qquad
 \partial_t\phi=\alpha v_\phi+\beta^i\partial_i\phi.                \tag{4}
\]

Set `R^A=g^{BC}Gamma[g_0]^A_{BC}`, evaluated on this initial metric, and
`Gamma_gamma^i=gamma^{jk}Gamma[gamma]^i_{jk}`. The remaining time derivatives
are specified explicitly by

\[
\begin{split}
 (\partial_t-\beta^j\partial_j)\alpha
    &=\alpha^2\operatorname{tr}_\gamma\kappa-\alpha^3R^0,\\
 (\partial_t-\beta^j\partial_j)\beta^i
    &=\alpha^2(\Gamma_\gamma^i-R^i-\beta^iR^0)
       -\alpha\gamma^{ij}\partial_j\alpha .                     \tag{5}
\end{split}
\]

These are pointwise formulas in the canonical fields, first spatial
derivatives and fixed reference connection. To verify (5), use
`Gamma[g]^A=-|g|^{-1/2}partial_B(|g|^{1/2}g^{AB})`, which gives

\[
 \Gamma[g]^0=\frac{\operatorname{tr}\kappa}{\alpha}
             -\frac{(\partial_t-\beta\cdot\partial)\alpha}{\alpha^3},
 \quad
 \Gamma[g]^i=-\beta^i\Gamma[g]^0
       -\alpha^{-2}(\partial_t-\beta\cdot\partial)\beta^i
       +\Gamma_\gamma^i-\alpha^{-1}\gamma^{ij}\partial_j\alpha.
\]

Thus `C=0` exactly initially. Formulas (1), (3)--(5) give a smooth map from
the fixed canonical spaces to metric/scalar fields in `H^s` and their time
velocities in `H^{s-1}`. Products, inverse metrics and square roots are smooth
in these spaces on the uniformly positive neighborhood (`s-1>4/2`). No
elliptic gauge solve and no change in the physical canonical data occurs.
All formulas are local and agree with the background in the open shell
collar. Differentiating them at zero gives exactly the desired harmonic
first tangent: its canonical and lapse/shift values agree, and `DC=0`
uniquely fixes the five remaining lapse/shift time velocities.

The physical Cauchy constraints are projections of

\[
 H_{AB}:=E_{AB}-\tfrac12g_{AB}\operatorname{tr}_gE,\qquad
 H(u,u)=H(u,e_i)=0.                                              \tag{6}
\]

They are not the condition `E(u,u)=0` alone. Given a reduced solution through
the completed data, let `X_A=(nabla_u C)_A` on the initial section. Because
`C=0` there, all tangential derivatives of `C` vanish. In an orthonormal
`u,e_i` frame the reduced equations give

\[
 E(u,u)=X_u,\quad E(u,e_i)=\tfrac12X_i,\quad E(e_i,e_j)=0,
 \quad \operatorname{tr}E=-X_u.
\]

Equation (6) therefore says `X_u=X_i=0`. Coordinate `partial_t C` also
vanishes because `partial_t=alpha u+beta` and tangential derivatives vanish.
This is a conditional consequence of the reduced equations at time zero;
the Cauchy data themselves prescribe only the usual field values and first
time velocities, not independent second metric velocities.

## 3. Initial corner identities after exact boundary straightening

Use the exact state flow of `NONLINEAR_BOUNDARY_STRAIGHTENING.md` to write
the perturbation as `U_D=V_D`, `U_N=Psi(z,t,Omega,V_N)`. Its zero state is the
background and the transformed boundary laws are exactly

\[
 V_D|_{z=0}=0,\qquad \partial_z V_N|_{z=0}=0 .                    \tag{7}
\]

Both completed initial perturbation and its first velocity vanish on one
open shell collar. The transformed equation has zero-state solution, smooth
coefficients and nonzero `g^{tt}`. Recursively solve that equation for each
initial `partial_t^j V`: if all preceding jets vanish on this open collar,
their spatial derivatives and the right side vanish there as well. Hence
every time jet defined at the chosen finite regularity vanishes there.
All Dirichlet and Neumann corner identities whose traces are required
therefore hold exactly. Parameter variations and difference initial data
inherit this same open collar. This argument supplies the initial jets;
the quasilinear iteration must still preserve its corresponding compatible
jets, as addressed in the separate evolution argument.

## 4. Exact subsidiary bulk and shell system

The exact identity

\[
 \nabla^A E_{AB}-\tfrac12\nabla_B\operatorname{tr}E=-S\phi_B       \tag{8}
\]

follows by contracted Bianchi and differentiating the scalar terms. With
the exact scalar equation and reduced Einstein equations it gives

\[
 \Box_g C_B+\operatorname{Ric}_B{}^A C_A=0.                       \tag{9}
\]

There is no independent forcing in (9). Arbitrary analytic forcing of a
reduced wave equation does not in general have this property.

At the shell the physical laws are
`K_ab=sigma(phi) h_ab/6` and `n(phi)=-sigma'(phi)/2`; its dimension is four.
Codazzi, with `K=1/2 L_n h`, gives

\[
 \operatorname{Ric}(n,e_a)=D_bK^b{}_a-D_a\operatorname{tr}K
      =-\tfrac12\sigma'(\phi)D_a\phi=(n\phi)D_a\phi.
\]

Thus `E(n,e_a)=0`, including field-dependent tension. The imposed normal
harmonic primal boundary row is `C(n)=0`. Differentiating this scalar along
the shell gives `(nabla_a C)(n)=-K_a{}^b C_b`, hence the reduced mixed
Einstein row gives the homogeneous boundary system

\[
 C(n)=0,\qquad (\nabla_n C)_a-K_a{}^b C_b=0.                      \tag{10}
\]

In the shell coordinates `g_ya=0`, `g_yy=N^2`, this is precisely

\[
 C_y=0,\qquad \partial_y C_a-\frac{N\sigma(\phi)}3 C_a=0.         \tag{11}
\]

Indeed `Gamma^b_{ya}=N K_a{}^b`; its contribution to `nabla_n C_a` supplies
the second copy of `K`. Omitting it would incorrectly halve the coordinate
Robin coefficient. Tangential derivatives of `C_y|_b=0` vanish; terms
containing `C_y` vanish there. No four tangential harmonic rows are added
to the primal sixteen-row boundary problem.

For completeness, uniqueness of (9)--(11) needs only the elementary
variable coefficient **linear** energy estimate at this already obtained
solution, not a new quasilinear construction. Use the fixed positive
component control metric and a product connection preserving the split
`R dy` and tangential one-forms. In inward `z`, write (11) as
`C_y=0`, `-partial_z C_a+M C_a=0`, where
`M=-rho_b N sigma/3` on the tangential block. Extend this bounded smooth
coefficient into a thin fixed collar and put
`C_y=W_y`, `C_a=exp(chi(z)M) W_a`, with `chi(0)=0`, `chi'(0)=1`.
Then `W_y=0`, `partial_z W_a=0` at the boundary. The scalar principal
metric remains `g^{mu nu}`, whose cross terms `g^{za}` vanish there.
Lower-order bundle/connection terms are bounded because the solution has
class `s>=14`. On the fixed near-background hyperbolic neighborhood, the
positive time energy is equivalent to
`||W_t||_2^2+||D_x W||_2^2+||W||_2^2`. Its boundary flux is a scalar
multiple of `<W_t,W_z>` and vanishes by the homogeneous orthogonal split.
Integration by parts and Cauchy--Schwarz give
`E_C'(t)<=A(t)E_C(t)` with integrable bounded `A` depending on this fixed
solution and interval. Zero initial `C` and `partial_t C` imply `C=0` by
Gronwall. Thus an exact reduced solution through the completed canonical
data satisfies the full physical equations.

The subsidiary argument needs neither positivity of the physical
gravitational energy nor symmetry of an arbitrary primal Robin matrix.
The positive energy here is an auxiliary analytic norm.

## 5. Explicit conditional two-parameter-derivative budget

Let `E_r={V: partial_t^j V in C([0,T];H^{r-j}), 0<=j<=r}` with the sum of
these derivative norms, on the controlled four-dimensional spatial region,
and let `V_epsilon` be the
solutions of the exact transformed homogeneous boundary system (7). Assume
the separate quasilinear theorem gives all of the following on a common
finite interval: uniform `E_s` bounds; exact (7); the stated hyperbolic
margins; uniqueness; and the variable coefficient linearized energy
estimates at all orders `s-4<=r<=s-1`, with source order `r-1`, including
the time-commuted source norms and initial jets used for normal recovery.
The estimates must allow the bounded first/zero-order terms produced by
linearization. This is a precise dependency, not an imported general
mixed-boundary theorem.

Write the transformed equation as

\[
 A^{\mu\nu}(x,V)\partial_\mu\partial_\nu V=F(x,V,\partial V).
\]

`A` depends on the state, not its derivatives. Define the linearization

\[
 L_V W=A^{\mu\nu}(V)\partial_{\mu\nu}W
       +A_V^{\mu\nu}[W]\partial_{\mu\nu}V
       -F_V W-F_P\partial W.                                   \tag{12}
\]

The first and second variations have the homogeneous boundary conditions
(7), initial values obtained by differentiating the completed initial data,
and equations

\[
 L_V W=0,\qquad
\begin{split}
 L_V Z={}&-2A_V^{\mu\nu}[W]\partial_{\mu\nu}W
          -A_{VV}^{\mu\nu}[W,W]\partial_{\mu\nu}V\\
        &+F_{VV}[W,W]+2F_{VP}[W,\partial W]
                       +F_{PP}[\partial W,\partial W].          \tag{13}
\end{split}
\]

The term `A_V[W] partial^2 V` in (12) is a zero-order matrix coefficient
in `W`; it is not a second derivative of `W`. Sobolev multiplication in
dimension four and the uniform `E_s` bound control it at first-variation
order `s-1`. Equations (12)--(13) give uniform `W` in `E_{s-1}` and `Z`
in `E_{s-2}`: the worst displayed second-variation source is
`W partial^2 W` in order `s-3`. The same counts apply to the finite
time-commuted spaces, with one Sobolev order charged for each time
derivative and the highest two normal derivatives recovered by the PDE.
All relevant orders are safely above the multiplication threshold.

Formal variation values do not establish actual parameter differentiability.
For clarity the following are the sufficient comparison estimates, obtained
from the declared linearized energy estimates and the integral Taylor
formula for smooth `A,F`:

| Quantity, with `delta=epsilon'-epsilon` | Sufficient comparison norm |
| --- | --- |
| `V_epsilon'-V_epsilon` | `O(abs(delta))` in `E_{s-1}` |
| `(V_epsilon'-V_epsilon)/delta-W_epsilon` | `O(abs(delta))` in `E_{s-2}` |
| `W_epsilon'-W_epsilon` | `O(abs(delta))` in `E_{s-2}` |
| `(W_epsilon'-W_epsilon)/delta-Z_epsilon` | tends to zero in `E_{s-4}` |
| `Z_epsilon'-Z_epsilon` | tends to zero in `E_{s-4}` |

Here all constants are uniform on a smaller compact parameter interval.
To see the derivative costs explicitly, a coefficient difference times
`partial^2(V_epsilon'-V_epsilon)` in the first quotient has source order
`s-3` after division by `delta`, giving the second row. In the quotient
of first variations, the new worst terms contain
`partial^2(W_epsilon'-W_epsilon)` or
`partial^2((V_epsilon'-V_epsilon)/delta-W_epsilon)`, hence have at least
source order `s-4`. Reserving order `s-5` for that comparison permits
`E_{s-4}` convergence with a spare derivative. Taylor coefficient remainders
are small in the same product spaces because the first three rows are
bounded and convergent there. Subtracting (13) at the two parameter values
proves continuity of `Z` in the same lower norm. Smooth initial data curves
give the corresponding convergent initial quotient remainders. Boundary
remainders are exactly zero because (7) is parameter-independent. The
initial open background collar supplies every variation corner identity.

Consequently the common solution curve is `C^2` into `E_{s-4}`. This
deliberately conservative four-loss statement does not assert sharp
variation orders or top-norm differentiability. It remains conditional on
the separate evolution/energy theorem declared above; equations (12)--(13)
alone are not enough to claim it.

The smooth pointwise inverse boundary flow and ADM changes preserve this
lower-space `C^2` conclusion. For `s>=14`, every spacetime derivative of
total order at most four is controlled in `H^{s-8}`, whose exponent is at
least six, greater than `4/2`. Spatial Sobolev embedding and time continuity
therefore give `C^2` parameter dependence into the piecewise spacetime
`C^4` norm on the compact controlled region. No claim concerns a common
infinite-order spatial family.

## 6. Identical copies and the normal-force consequence

Use exactly the same completed initial data and reduced boundary construction
on the two identical interior copies, with each copy's own outward normal
pointing toward its shell. Uniqueness makes the copy fields identical under
the copy identification. The induced metric and scalar therefore match,
and the prescribed identical-copy physical extrinsic curvature and scalar
normal derivative laws hold on both sides. This is the symmetric branch
only, not a theorem for independent nonsymmetric bulk data.

The normal shell displacement equation is not an additional incoming
boundary row. Once the full bulk equations, matching, Israel law and scalar
jump law hold, the normal projection of the distributional
Einstein--scalar identity gives the shell force identity. Equivalently,
in an oriented jump convention the normal-force balance involves the
average extrinsic curvature and average scalar normal derivative multiplied
by their jumps, together with the jump of normal bulk stress. Reflection
interchanges the two sides, so all averages odd under that reflection and
the jump of identical normal bulk stresses vanish. The force residual is
therefore zero. A convention change from oriented jumps to two outward
normals changes individual signs, not this zero consequence. This argument
requires the full physical equations established by the subsidiary step;
it would not justify a force claim for arbitrary forced reduced solutions.

## 7. What this note closes and what it does not

Closed here: a local exact initial harmonic completion that preserves both
canonical data and the selected first gauge tangent; initial zero subsidiary
velocities from the correct physical constraints; initial open-collar corner
jets; exact bulk/physical shell subsidiary identities; auxiliary-energy
subsidiary uniqueness conditional on an obtained regular solution; and a
precise four-loss parameter argument conditional on declared mixed boundary
evolution estimates. These do not replace construction of that evolution,
its compatible iteration or its prescribed lifespan. In particular this
note alone is not an exact shell-contact existence theorem or evidence for
a higher-dimensional origin of the Big Bang.
