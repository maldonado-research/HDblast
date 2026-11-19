# How often can the emitted particles collide after the source switches off?

On the quiet shell, a particle reaches a probability of at least one collision of 50% on number average at $y^2\kappa H_0^4=0.02823$ for the extrapolated contact law. Applying a single elastic partial-wave unitarity limit across the initial spectrum requires $\kappa H_0^4\le 4.2131\times 10^{-7}$. At illustrative production coupling $y=0.1$, at most **1.074 in ten million** of the initial particles collide by the end of the measured window in the contact scenario specified below. This result concerns opposite-charge scattering among particles supplied by the modeled source. [Numerical results](CONTACT.json) · [Calculation](transport.py)

## Collision law and units

Consider one massless Dirac species produced in equal particle and antiparticle populations. Assume a homogeneous isotropic classical dilute gas obeying the Boltzmann equation, initially with the one-body spectrum reconstructed from scalar pair emission. Take only elastic $\psi\bar\psi$ collisions, isotropic in their center-of-momentum frame. With $\sigma_{\rm tot}=\sigma_{\rm tr}=\kappa\mathsf{s}$, the coefficient $\kappa$ has units of energy to the power $-4$ in natural units. The transport cross section uses the weight $1-\cos\theta_{\rm cm}$. The invariant $\mathsf{s}=2pp'(1-\mu)$ uses physical momenta and the incoming directional cosine $\mu$. The massless Møller flux is $1-\mu$. The proper collision rate of one tag is

$$
\Gamma(p)=\int dn_{\bar\psi}(p')\int_{-1}^{1}\frac{d\mu}{2}
(1-\mu)\,2\kappa pp'(1-\mu)
=\frac{8}{3}\kappa p\rho_{\bar\psi}.
$$

Only the opposite species contributes to this target density. The pair event count equals the count of particles of either one sign. Write $k=ap/H_0$, $x=H_0\tau$, $h=H/H_0$ and $\lambda=y^2\kappa H_0^4$. Define the initial spectrum $dN$ per comoving volume in units $y^2H_0^3$ for one sign and $M_1=\int k\,dN$. After the production gate closes, massless expansion and elastic collisions conserve total conformal energy. Charge symmetry keeps the two species' energies equal. Thus $\rho_{\bar\psi}=y^2 H_0^4 M_1/a^4$ even if the target spectrum changes shape. Before its first collision a tag has constant $k$. These facts yield its instantaneous rate, integrated depth and first-collision probability

$$
\frac{\Gamma(k,x)}{H(x)}=\frac{8\lambda k M_1}{3a(x)^5 h(x)},\qquad
J(x)=\int_{x_g}^{x}\frac{dx'}{a(x')^5},\qquad
\tau(k,x)=\frac{8}{3}\lambda k M_1 J(x),\qquad
P_{\ge1}=1-e^{-\tau}.
$$

The uncollided-tag population solves the homogeneous Boltzmann loss equation with this rate. For the exact contact law the first-collision probability therefore remains valid during nonlinear elastic redistribution under the stated assumptions. It does not count subsequent collisions. Initial pair correlations, Pauli factors, additional targets, further emission, inelastic processes, and a response of the geometry to collisions require additional modeling.

## Number averages and pulse controls

For each computed pair with total comoving energy $\omega$ and total momentum magnitude $q$, scalar rest-frame isotropy makes either daughter's $k$ uniform between $(\omega-q)/2$ and $(\omega+q)/2$. Set $u=(8/3)\lambda M_1J$. The survival probability at fixed $(q,\omega)$ is $\exp[-u(\omega-q)/2](1-e^{-uq})/(uq)$. Weighting it by emitted pair counts gives the fraction of the initial tags that survive without collision. A direct angular quadrature of $1-e^{-uk}$ provides an independent check. First collisions of any sign-preserving tag do not select a thermal final spectrum.

| Wall | Gate closing time | $J(2.16)$ | Mean initial-tag depth per $\lambda$ at 2.16 | $\lambda$ for 50% colliding at 2.10 | $\lambda$ for 50% colliding at 2.16 |
|---|---:|---:|---:|---:|---:|
| Quiet | 2.050000 | 0.03804155 | 25.4908 | 0.0605041 | 0.0282303 |
| Plus | 2.050089 | 0.03801496 | 25.3854 | 0.0608162 | 0.0283481 |
| Minus | 2.049257 | 0.03832044 | 26.0734 | 0.0586709 | 0.0275956 |
| Far | 2.049887 | 0.03809255 | 25.5346 | 0.0603248 | 0.0281813 |

On the quiet shell $\ln(a_{\rm end}/a_{\rm gate})=0.0193866$. The full-population instantaneous mean $\Gamma/H$ at the gate is $1304.84\lambda$ using the initial spectrum. Its crossing of unity at $\lambda=0.00076638$ has little bearing on a first collision during this short interval. The endpoint instantaneous values in the JSON average the rates along original freely redshifting tags with their initial weights. At large depth, conditioning on tag survival changes those weights.

For quiet upper-band tags in the entire source population, the mean depth at 2.16 is $45.9442\lambda$. Removing the lower-band targets reduces it to $1.90360\lambda$. Here the two bands label the energy of the *original pair* that made the tag. The two populations' one-particle momentum spectra overlap. The results include the complementary lower-band tests. The numerical integration uses the native scalar-clock gate and each wall's measured scale factor. Eight-point Gaussian integration in each wall interval and a separate Simpson integration on the measured clock agree in $J$ within $7.5\times10^{-11}$ relative. Repeating the reconstruction with 72 instead of 96 frequency nodes and three instead of five integration nodes per native time interval changes the mean depths and the 50% and 90% thresholds by at most $2.6\times10^{-11}$ relative.

## Where the contact law can apply

To illustrate a restricted physical interpretation, impose one elastic partial-wave bound $\sigma_{\rm tot}\le 16\pi/\mathsf{s}$, with the normalization for distinguishable particles in a single angular-momentum-zero channel. This follows from $\sigma_0=4\pi\sin^2\delta_0/p_{\rm cm}^2$ and $p_{\rm cm}=\sqrt{\mathsf{s}}/2$. The assumed cross section incorporates any spin weighting into the one effective channel. On the initial continuous frequency band with the retained parent-momentum grid, $k_{\rm max}=(120+\max q)/2$. Two head-on opposite-species particles can approach $\mathsf{s}_{\max}/H_0^2=4k_{\rm max}^2/a_g^2$. Requiring contact scaling across this original support implies $\kappa H_0^4\le 16\pi/(\mathsf{s}_{\max}/H_0^2)^2$. For the quiet wall at this coefficient and $y=0.1$, the formal first-collision fraction is $1.07394\times10^{-7}$. The corresponding values across all four walls range from $1.06944\times10^{-7}$ to $1.09767\times10^{-7}$.

An elastic collision can move a daughter above the initial single-particle energy range. One admissible rate-level extension takes $\sigma=\min(\kappa\mathsf{s},16\pi/\mathsf{s})$. At subsequent energies the extension obeys $\sigma\le\kappa\mathsf{s}$. Conservation of target conformal energy then keeps the contact-law rate an upper bound for every surviving tag. The fraction just reported also bounds that extension from above. Other energy dependences, extra angular channels, stronger infrared scattering, additional species and a longer history are separate physical cases. The large formal $\lambda$ thresholds quantify what the extrapolated rate law would require. Evaluating them as physical predictions would require a microscopic interaction valid at those energies and a production calculation valid at the chosen $y$.

## Replay

From the study root, with the repository or archive root supplied to the command, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB transport/transport.py --repo ../.. --output /tmp/CONTACT.json
~~~

The per-spectrum rates use the sibling daughter calculation and the prior archived scalar gate and wall data. [Daughter spectrum](../daughters/DISTRIBUTIONS.md) · [Quantum-memory edition](https://doi.org/10.5281/zenodo.24374531)
