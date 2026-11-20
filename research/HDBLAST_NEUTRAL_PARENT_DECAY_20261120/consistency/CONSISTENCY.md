# Which particles could supply the model's decay and drag rates?

The fluid calculation inputs a Yukawa decay width and a separate drag rate of six per wall time. At the simulated peak, one possible electromagnetic normalization gives an upper limit of **0.0378 per wall time** for elastic scattering through a single isotropic partial wave off the opposite charge species. The limit already assumes the largest scattering amplitude allowed by unitarity in that channel. A screened, weak-vector scattering model gives an even lower rate. Other particle masses, wall tensions and interactions can change the conclusion. A second example below explicitly realizes the requested drag in a local kinetic calculation at one instant.

## Normalize the abundances and couplings

Use $c=\hbar=k_B=1$, curvature length $L$, wall tension $\sigma$, dimensionless mass $\mu=mL$, wall action $S_0=\sigma L^4$ and wall clock $\tau=t/L$. The classical force parameter for a daughter that carries one of the three independent gauge charges is $Q=q/(\sqrt3\mu)$. For a local, slowly varying, surface-dominated probe the canonically normalized coupling is

$$
\alpha(K)=\frac{q^2}{4\pi S_0K}=\frac{3Q^2\mu^2}{4\pi S_0K}.
$$

The simulation fixes $Q=10^8$. It leaves $S_0$ and $\mu$ unspecified. The value $\alpha=1/137.036$ below is an illustrative **late** normalization. It does not identify the fields with the visible photon. The surface coefficient grows by about 13.2 million from its initial value. Extending a fixed-charge weak-coupling approximation back to that initial epoch is a separate issue. [Earlier normalization and state counts](https://doi.org/10.5281/zenodo.24405438)

Let the neutral spin-zero parent have $M=2mR(\psi)$ and equal Yukawa couplings to three independent Dirac flavors. With six daughter species counting both signs, the *sum* over the open channels gives the proper vacuum width

$$
L\Gamma=\kappa R(1-R^{-2})^{3/2},\qquad
\kappa=\frac{\mu}{4\pi}\sum_{i=1}^3y_i^2,
$$

above $R=1$. Holding $\kappa=100$ therefore fixes $y_i^2=4\pi\kappa/(3\mu)$ in the equal-branching interpretation. At $R=1.05$, $L\Gamma=2.9765$. Finite temperature, final-state blocking, inverse decay, stimulated effects, finite mode spacing and field-induced changes to the spectral function can modify the vacuum law.

In dimensionful proper time let $\mathcal D=n_P\Gamma$ be the number of decays per wall four-volume. For parents at rest relative to the wall, the number sources are $-\mathcal D$ for the parents and $+2\mathcal D$ for the daughters, added to the respective $\dot n+3Hn$ equations. The parent energy equation has sources $n_P\dot M-M\mathcal D$. The daughter's energy and pressure satisfy $\dot\rho_d+3H(\rho_d+p_d)=J E+M\mathcal D$. The scalar equation supplies the force $-n_PM_{,\psi}$, canceling $n_P\dot M$ in the summed energy. Isotropic birth injects zero mean momentum, so a drift equation written per daughter has dilution $-(2\mathcal D/n_d)hu$ for stream momentum per mass $hu$, in addition to electric, expansion and drag terms. Evolving total daughter momentum and energy density keeps the empty initial state regular. The two-body momentum deposited at birth is $p_b=m\sqrt{R^2-1}$ per particle and is generally nonzero. The decay itself supplies the kinetic energy $M-2m$ per pair. The much larger subsequent electrical work is separately accounted in the fluid calculation. These sources conserve number with the stated stoichiometry and total stress energy. A thermodynamic entropy theorem for the reacting mixture would also need a consistent inverse reaction and chemical potentials.

## Pauli exclusion at birth

At maximum intrinsic gas energy in the primary run, $a=9.76159$, $D=n_dm/\sigma=9.02384\times10^{-7}$, $T/m=1760.77$, the bulk drift speed is 0.00147 and 99.675% of parents have converted. Each sign/flavor has two spin states. Putting the entire daughter population at that density into cold nearly stationary Fermi spheres would require

$$
x_F^3=(p_F/m)^3=\frac{\pi^2D S_0}{2\mu^4}
$$

at zero drift. This is a **batch-packing comparison**. During production, accelerating or scattering earlier pairs clears some birth states. The maximum injection momentum is $0.32016m$. With $\mu=100$ and late $\alpha=1/137.036$, $p_F=22.26m$ in the batch comparison. With $\mu=10^5$ at the same $\alpha$, it falls to $0.223m$. Neutral bosonic parents have no analogous exclusion limit on their occupation of a single state. That does not give their fermionic daughters additional states.

The field-free round sphere supplies an independent discrete check. For one four-component Dirac flavor, the positive-energy momenta are $(l+3/2)/(aL)$ with multiplicity $2(l+1)(l+2)$ for one charge sign. Filling both spatial-Dirac branches through $l_{\rm max}=\lfloor a\mu p_b/m-3/2\rfloor$ gives $2(l_{\rm max}+1)(l_{\rm max}+2)(l_{\rm max}+3)/3$ states per sign and flavor, or zero when $l_{\rm max}<0$. At this radius a field-free two-body decay with $R\le1.05$ cannot even reach the lowest free angular mode if $\mu<0.480$. The small phase charge $q=1$ at fixed $Q=10^8$ gives $\mu=5.77\times10^{-9}$ and zero available states in that free comparison. A background gauge connection can shift angular levels, and a time-dependent field can perform additional work. The comparison is conditional on the indicated free, instantaneous geometry. [Round-sphere Dirac spectrum](https://arxiv.org/abs/gr-qc/9505009) · [Cosmological Dirac mode count](https://arxiv.org/abs/1603.00037)

A complementary diagnostic uses the actual injection history. Write $x_b=p_b/m$, dimensionless rate $\gamma_d=L\Gamma$ and dimensionless clock derivatives with a prime. At the birth shell the changing shell momentum, redshift and electric force sweep radial momentum states at signed speed $x_b'+Hx_b-QX\cos\vartheta$ in units of $m/L$. For an infinitesimally narrow shell, collisionless propagation smoothly to either side and Pauli occupation at most one, the **largest total sweep rate** in dimensionless number-density units is $12 x_b^2\langle|x_b'+Hx_b-QX\cos\vartheta|\rangle/(2\pi^2)$. The requested production rate in those units is $\gamma_d(\rho_P/(\sigma R))S_0/\mu^4$. The factor twelve counts six species with two spins. This necessary condition permits empty incoming states and maximal filling on the departing side in every direction. The actual angular occupation can be more restrictive. For $\mu=100$ with the reference late $\alpha$, requested rate divided by that generous capacity rises from 0.790 at $\tau=10.36$ to 1.36 at 10.4, 1.98 at 10.5 and 8.46 at 11.0. These are sampled local checks on the computed field. Finite birth bandwidth, collisional broadening, alternative quantum states and backreaction require their own kinetic treatment. The actual probability of blocking is not inferred from a density ratio alone.

## Derive the linear drag of a specified collision channel

The fluid drag relaxes a small relative **momentum**. Counting all collisions can overestimate this relaxation because forward scattering transfers little momentum. Near the peak $T\gg m$ and the classically modeled streams drift slowly. Consider an isotropic dilute ultrarelativistic Boltzmann distribution with a small opposite drift of the two charge signs. The target density is $n_t=n_d/6$ for each sign in its own gauge sector. Two-body elastic scattering with an azimuthally symmetric centre-of-momentum differential cross section has transfer cross section $\sigma_T=\int(1-\cos\theta_*)\,d\sigma$. To first order in opposite bulk speeds $\pm v$, the incoming joint weight is perturbed by $v(p_{1z}-p_{2z})/T$, while the directed momentum per particle is $4Tv$. In the collision average the change in the first particle momentum is $-(p_{1z}-p_{2z})\sigma_T/(2\sigma)$. The relaxation coefficient is consequently

$$
\nu_{\rm drag}=\frac{n_t}{8T^2}\left\langle v_{\rm M}\sigma_T(p_{1z}-p_{2z})^2\right\rangle,
$$

where $v_{\rm M}=1-\cos\vartheta$ is the massless Møller flux velocity. The average samples independent equilibrium momenta of each sign. With isotropic elastic $s$-wave unitarity $\sigma_T\le4\pi/k^2$, $k^2=p_1p_2(1-\cos\vartheta)/2$. Angular averaging the squared momentum difference and using $\langle p\rangle=3T$, $\langle1/p\rangle=1/(2T)$ gives

$$
L\nu_{\rm drag}\le \frac{\pi n_t L}{T^2}
=\frac{\pi D S_0}{6\mu^3(T/m)^2}.
$$

The maximum thermally averaged collision frequency in this model is twice that drag ceiling. With late $\alpha=1/137.036$, the ceiling is $0.03777/L$ for $\mu=100$ and $0.003777/L$ for $\mu=1000$. The requested $6/L$ exceeds even a saturated $s$ wave by a factor 159 at $\mu=100$. The ceiling decreases as $1/\mu$ at fixed $\alpha,Q,K,D,T/m$. An interaction of longer range can excite multiple partial waves. Other targets, radiation and strongly anisotropic distributions require a different calculation. [Partial-wave normalization](https://arxiv.org/abs/2405.02222) · [Relativistic gauge transport](https://arxiv.org/abs/hep-ph/0302165)

For a separate weak-vector benchmark in the screened plasma, put $m_D^2 L^2=DQ^2/[K(T/m)]$ for Boltzmann carriers with both signs of one gauged flavor. This gives $m_D L\simeq19.7$. The spin-averaged massless distinguishable vector $t$ channel has its forward envelope $d\sigma/d\Omega\le\alpha^2/[4k^2(\sin^2(\theta_*/2)+m_D^2/(4k^2))^2]$ in the simple static screening model. Integrating the transfer weight, then averaging over the gas, gives the conservative envelope
$$
L\nu_{\rm vector}\le\frac{\pi n_t L}{T^2}\frac{\alpha^2}{2}\ln[1+20T^2/m_D^2].
$$
One way to check the logarithm is to set $x=p_1/T,y=p_2/T,c=\cos\vartheta$ and $W=x/y+y/x-2c$. For the two thermal momenta $\langle W\rangle=3$ and $\langle Wxy(1-c)\rangle=30$. Dropping the negative term in the transport logarithm and applying concavity of $\ln$ with weight $W$ gives the displayed bound. At $\mu=100$ and the reference $\alpha$ it is $2.13\times10^{-5}/L$. This is an elastic, electrostatically screened Born-gas calculation. It does not bound all weak-QED transport processes, collective instabilities, multiple targets or non-Abelian sectors. Screening is evaluated within the same local surface-dominated Boltzmann approximation. The mean free path required for $\nu L=6$ is about $0.017$ of the wall radius at the peak.

The spin-zero parent provides another explicit channel through $i\bar i\to j\bar j$, $j\ne i$. At a representative head-on momentum $k=3T$ the program uses the spin-averaged scalar $s$-channel cross section
$$
\sigma_{i\to j}(s)=\frac{y_i^2y_j^2}{16\pi s}
\frac{(s-4m^2)^2}{(s-M^2)^2+M^2\Gamma^2}.
$$
There are two distinguishable final flavors. For $\mu=100$ and the reference $\alpha$ their combined representative rate is $4.66\times10^{-4}/L$. A full quantum collision operator would involve additional channels, statistical factors, inverse reactions and averaging the resonance over actual distributions. The program suppresses the local rate estimates for the $q=1$ ultra-light example, where their scale and continuum assumptions fail.

## A constructive rate benchmark

Imposing both cold batch packing under $0.32016m$ and the hot, isotropic $s$-wave rate at the specified peak requires $\mu\gtrsim5.34\times10^9$ in this comparison. This follows by combining the cold-capacity upper bound $S_0\le 2x_b^3\mu^4/(\pi^2D)$ with the hot-rate lower bound $S_0\ge6(\nu L)\mu^3(T/m)^2/(\pi D)$. Cold batch packing is stronger than gradual injection and acceleration require.

For example, choose $\mu=10^{11}$ and $S_0=10^{47}$. The free cold capacity is 7.37 times the required count and the hot single-channel drag ceiling is $15.24/L$. The local gauge coupling is $1.81\times10^{-14}$ and the summed Yukawa coupling divided by $4\pi$ is $10^{-9}$. Introduce an additional elastic interaction between opposite signs with a unitary isotropic phase shift $\sin^2\delta_0(k)=k^2/(k^2+k_0^2)$. Its transfer cross section is $4\pi/(k^2+k_0^2)$. Choosing $k_0/T=2.25396$ and doing the full thermal momentum and angular integral gives **$L\nu_{\rm drag}=6.00000$**. Varying the quadrature order from 24 to 40 changes the fitted $k_0/T$ by $1.7\times10^{-5}$. This explicitly calibrated local scattering law satisfies elastic unitarity and provides an example of the extra microscopic input needed. The simulation has yet to evolve its temperature- and density-dependent rate, establish that the distribution reaches local equilibrium, or derive this phase shift and scale normalization from a fundamental action.

Keep that same physical $k_0$ and evaluate the law along six hot, low-drift states of the continuously decaying solution. Each target density is $n_t=D\sigma/(6m\gamma)$, using the one opposite-charge population for that field, and the ratio $k_0/T(\tau)$ changes as the gas cools. At wall times 11.1, 11.4, 11.5, 11.7, 12.0 and 12.4 the inferred $L\nu_{\rm drag}$ is respectively **4.306, 5.295, 5.494, 5.736, 5.912 and 5.997**. The temperature-to-mass ratio at these points falls from 2082 to 1761 while the converted fraction rises from 83.3% to 99.65%. The kernel therefore approaches its endpoint calibration during this portion of expansion. These values re-evaluate the scattering law on the stored fluid history. Evolving that law back into the Einstein–Maxwell–matter system could change the history itself. At earlier times the two streams can have large opposite velocities or be cold, outside the hot small-drift reduction.

Run `python -I -B consistency/consistency.py --output /tmp/CONSISTENCY.json` from the study directory. The only dynamical input is `coupled/DECAY.json`. [Calculation](consistency.py) · [Results and parameter choices](CONSISTENCY.json)
