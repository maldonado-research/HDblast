# Identifying the observed first derivative with the registered master wave

At the fixed observer clock $q=-\tau/\rho_b$, the derivative of the shell scalar for the selected canonical family is

$$
 \left.\partial_\epsilon\Phi_\epsilon(q,n)\right|_0
     = S(-q)Y(n),\qquad Y(n)=n_1^2-n_2^2,\qquad 0\le q\le25\ell/16.
$$

Here $S$ is the embedding-included physical scalar coefficient in the registered past linear waveform. This identification combines the previous linear reconstruction and its waveform enclosures with the actual differentiated family and the common observer domain. It retains the model, background, regularity, reconstruction, and linear construction premises of those results. All source paths below refer to commit `609fd60518150d23c0995631aeffa6a4958e0088`.

## Matching the model and all Cauchy data

Use the two identical interiors on the regular branch near $\phi=+1$, with $\delta=1/1000$, $c=96714321425692/161839258930731$, and $\ell=10^{-6}$. The two radial coordinates satisfy $dy=-\rho\,dz$. In each copy the normal points outward and the physical rows are $K_{ab}=\sigma h_{ab}/6$ and $n_{\rm out}(\phi)=-\sigma'/2$. The reference connection is that of $g_0=dy^2+\rho^2(-dt^2+\cosh^2t\,\kappa)$, with $\phi_0=\phi_0(y)$. The scalar field, oriented axes and affine perturbation parameter are fixed. [Registered model and full seed, §§1–2](../HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/inputs/REGISTERED_EXACT_FAMILY_BINDING.md) · [Actual differentiated physical equations, §§3–4](../HDBLAST_EVOLVED_PARAMETER_20261116/EVOLVED_PARAMETER_AND_PHYSICAL_CHAIN.md)

Let $B(z)=\exp[4-\ell^2/((z-\ell)(2\ell-z))]$ inside $(\ell,2\ell)$ and zero elsewhere. The master has $w(0)=B$ and $w_t(0)=e^{-z}B$. With $p=\phi_{0,y}$, $\alpha=(\rho/\rho_b)^{3/2}p/p_b$, and $k=\alpha_z/\alpha$, reconstruction at the initial section gives

$$
 u=\alpha B/\rho^2,\quad u_t=e^{-z}u,\quad
 v=3\alpha(B_z+kB)/(\rho^3p),\quad
 v_t=3\alpha e^{-z}(B_z+(k-1)B)/(\rho^3p).
$$

In a background orthonormal spatial frame with radial slot first, put $s_1=u_tY/\rho$. The first geometric data are

$$
 \gamma_1=uY\,\mathrm{diag}(-4,2,2,2),\quad
 K_1=s_1\,\mathrm{diag}(-2,1,1,1),\quad
 \phi_1=vY,\quad \Pi_1=v_tY/\rho.
$$

The current canonical metric and momentum correction has zero value and derivative at parameter zero. The half-action momentum density gives $\pi_1/\rho_0=\mathrm{diag}(-3s_1/2,0,0,0)$ and $\varpi_1/\rho_0=v_tY/\rho$. Applying the actual inverse $K_{ij}=2\lambda(\gamma_{ik}\gamma_{jl}\tau^{kl}-\gamma_{ij}\operatorname{tr}_\gamma\tau/3)$, $\Pi=\lambda\omega$, recovers the displayed geometric pair, since $\lambda_0=1$ and background momenta vanish. [Source reconstruction and half-action, §2](../HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/inputs/REGISTERED_EXACT_FAMILY_BINDING.md) · [Canonical correction, equation (28)](../HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/NONLINEAR_CONSTRAINT_BALL.md) · [Selected seed](../HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113/SEED_ENVELOPE.md) · [Current Cauchy map, §§1–3](../HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md)

Both constructions use initial lapse $\rho(1+\epsilon uY)$ and zero shift. The full physical relative-harmonic one-jet consists of the displayed spatial and scalar values, $h_{tt}=-2\rho^2uY$, $h_{ti}=0$, and these coordinate-time velocities.

$$
 \dot h_{tt}=-2\rho^2u_tY,\quad
 \dot h_{ij}=2\rho K_{1,ij},\quad
 \dot h_{ty}=-6\rho^2(u_y+4au)Y,\quad
 \dot h_{tA}=0,\quad \dot\phi_1=v_tY,
 \qquad a=\rho_y/\rho .
$$

The dot denotes the initial coordinate-time derivative. The mixed radial velocity follows from the reference-connection term in the exact shift completion. [Current Cauchy map, equations (5)–(7)](../HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_CAUCHY_DATA.md) · [Full first jet in product variables, §3](../HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/SYMMETRY.md)

## A harmonic representative on the past slab

The past master on $-3\ell\le t\le0$ satisfies the dynamic shell law, and its support lies in $z\le 2\ell+|t|$. Retain the longitudinal reconstruction for this master, including its physical scalar junction. A harmonic representative on the *past* slab can be constructed directly. Write $\xi^y=fY,\ \xi^t=bY,\ \xi^A=dD^A Y$. The gauge defect is cancelled by $(\Box_0+\mathrm{Ric}_0)\xi=-DC[h_{\rm long}]$. The radial system for $U=(f,b,d)$ is

$$
 U_{tt}-U_{zz}=A U_t+B_1 U_z+D U+F,\qquad
 f|_{z=0}=0,\quad b_z|_{z=0}=d_z|_{z=0}=0,\quad U(0)=U_t(0)=0.
$$

The symbols $B_1,D$ denote coefficient matrices. Equations (9)–(10) of the gauge derivation give their smooth coefficients and $F=(6\rho u_z-24a\rho^2u,0,0)$ from the evolved longitudinal metric amplitude. With $s=-t$ solve for $U(-s,z)$ in its original component basis. The matrix multiplying the elapsed first derivative becomes $-A(-s,z)$. The other coefficients and the source are evaluated at $(-s,z)$. This construction uses the master with elapsed initial velocity $-e^{-z}B$. [Past master construction, §§3–4 and 6](../HDBLAST_PAST_SHELL_CONTACT_20261104/primary/PROOF.md) · [Gauge operator, equations (1), (5)–(12)](../HDBLAST_FIRST_GAUGE_WAVE_20261016/operator/GAUGE_OPERATOR_AND_CONSTRUCTION.md)

For existence take $0\le z\le6\ell$ with homogeneous Dirichlet data for all three gauge components at the outer endpoint. The free mixed wave with these separated componentwise rows has the usual energy estimate on $H^1\times L^2$ with the Dirichlet traces imposed. In the energy norm augmented by $\|U\|_2$, Duhamel iteration of the smooth bounded first-order and zeroth-order terms on $0\le s\le T=13\ell/8$ converges by the iterated Volterra bound $(C T)^j/j!$. The source is smooth by the retained master theorem. All corner compatibility conditions hold since the source and its initial jets vanish near both endpoints at $s=0$. Differentiating the system in $s$ gives the same energy estimate at each order with known lower-order terms. Using the equation to recover $U_{zz}$ then yields all required mixed regularity up to the walls. The exterior cone energy yields $U=0$ for $z>2\ell+s$ because the free principal radial speed is one, the source vanishes there, and the outer boundary flux is zero. In particular $2\ell+T=29\ell/8<6\ell$. This construction uses the direct energy estimate on the past interval.

The transformed pair $(h_{\rm long}+\mathcal L_\xi g_0,\ \phi_{1,\rm long}+\mathcal L_\xi\phi_0)$ satisfies the linearized bulk and relative-harmonic equations. The vector equation gives $DC=0$ throughout, including the normal-gauge row on the shell. At the shell $f=0$ preserves the embedding and $b_z=d_z=0$ preserves the mixed metric rows. Covariance preserves the two physical junction laws with their outward normals. The transformed first derivatives agree with the full harmonic one-jet above, including $f_{tt}(0)=-6\rho^2(u_y+4au)$. Apply the same fixed point-chart differential and elapsed tensor reflection to this pair and to the actual derivative. A covariant component acquires one minus sign for each temporal index under the reflection, and initial elapsed differentiation supplies another minus sign. Consequently even the mixed metric velocity matches. Both fields satisfy the same homogeneous four Dirichlet and twelve Neumann rows at the shell, and the candidate vanishes near the outer wall. Both solve the same linearized normalized bulk equation at the zero background state, since the point chart is a smooth invertible representation of these reduced equations and boundary rows. The transformed smooth candidate lies in the actual linear uniqueness class. Uniqueness with the common initial one-jet gives equality on $0\le s\le T$. [Gauge and boundary match, §4](../HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/inputs/REGISTERED_EXACT_FAMILY_BINDING.md) · [Point chart and its twelve physical/normal-gauge rows](../HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md) · [Reflection and walls, §1](../HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/REDUCED_EVOLUTION_AND_CAUSAL_SUPPORT.md) · [Linearized uniqueness, §8](../HDBLAST_EVOLVED_PARAMETER_20261116/MIXED_DIFFERENCE_AND_LINEARIZED_ESTIMATES.md)

## From the fixed shell to the moving cut

On the shell the gauge scalar term is $\xi^y p_b=0$. The background shell scalar is constant, so the first proper-clock and angular-label displacements also contribute zero. At parameter zero $t=-q$ and the initial sphere labels are $n$. Thus the scalar derivative on each moved observer equals the longitudinal physical trace

$$
 \partial_\epsilon\Phi_0(q,n)= S(-q)Y(n),\quad
 S(t)=\frac{3\cosh^{-3/2}t}{\rho_b^3p_b}\,[w_z(t,0)-\tfrac12 C_b w(t,0)] .
$$

Here $\partial_\epsilon\Phi_0$ denotes differentiation of the family at zero. The actual family and all cuts used here have the common $q$-domain $0\le q\le25\ell/16<T$, with positive density and the strong parameter regularity required for integration. The earlier signed wave estimates now refer to this specific actual derivative. Their hypotheses on the linear physical reconstruction and waveform enclosures remain in force. [Embedding-included trace, §§1–2](../HDBLAST_PAST_FUTURE_SCALAR_TRANSFER_20261104/primary/PROOF.md) · [Common effective observer domain, §§1–2](../HDBLAST_EFFECTIVE_OBSERVER_20261116/EFFECTIVE_OBSERVER_DOMAIN_AND_RESPONSE_BOUNDS.md)
