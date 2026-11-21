# How much ordinary radiation could the wall supply?

In the higher-mass example, the hot wall grows in radius by a factor of 2.25 under prescribed drag and 2.16 with the explicit elastic collision rate. A separate question is whether its whole closed volume has enough energy to account for the photons observed in our part of the universe. At the selected wall normalization, even converting every available unit of surface energy at turnaround into equilibrium radiation leaves its entropy **at least $6.34\times10^{43}$ times too small** to supply the photons in one present Hubble sphere under adiabatic evolution. “Hot” in the fluid simulation measures an energy fraction. The astronomical comparison also requires a physical energy density and a volume. [Budget and parameters](BUDGET.json) · [Particle feedback](../feedback/FEEDBACK.md)

## The photon target

For a two-polarization blackbody, the entropy density in units of Boltzmann's constant is $s_\gamma/k_B=(4\pi^2/45)(k_BT_\gamma/\hbar c)^3$. The measured temperature $T_\gamma=2.72548$ K gives $1.4792\times10^9$ per cubic metre. A present Hubble sphere has radius $c/H_0$ and, in the spatially flat reference cosmology, volume $4\pi(c/H_0)^3/3$. Taking the Planck base-ΛCDM inference $H_0=67.4$ km s⁻¹ Mpc⁻¹ gives $S_{\gamma,H}/k_B=1.6020\times10^{88}$. This counts photon entropy inside one specified sphere today. The full particle horizon encloses a larger volume. Varying the two inputs to opposite ends of the quoted individual uncertainties gives $1.566$–$1.639\times10^{88}$ as a sensitivity check. [Observational references](../LITERATURE.md)

## Energy constraints, stored stress and expansion history

Let $R_{\rm w}=aL$ be the physical radius, $e=\rho/\sigma$ the added surface energy density, $S_0=\sigma L^4$ the normalization and $m_b$ the wall-side bulk-mass parameter. The closed-wall Friedmann constraint in natural units reads
$$
a^2[1+(aH)^2]=a^4(2e+e^2)+m_b .
$$
For nonnegative $e$ and $m_b$, $a^4e\le a^2[1+(aH)^2]/2$. At the first expansion maximum $H=0$ this reduces to $a^4e\le a^2/2$. The energy density has the tighter ceiling $e\le[ a^2(1+\sqrt{1+a^{-2}})]^{-1}$. Live runs check the constraint using interpolated $e$ and $m_b$ at the turnaround, and track the mass throughout the preceding evolution. Wave flux and pressure work vary during the pulse, so the turnaround ceiling by itself need not bound earlier energy. [Wall evolution](../feedback/feedback.py)

Suppose a homogeneous ideal bath of relativistic particles with zero chemical potentials can use the entire added energy $e\sigma$. For $g_s=g_\rho=g$, the entropy of the complete three-sphere of volume $2\pi^2R_{\rm w}^3$ is
$$
\frac{S_{\rm eq}}{k_B}
=\frac{8\pi^2}{3}\left(\frac{\pi^2 g}{30}\right)^{1/4}
 (S_0 a^4 e)^{3/4}.
$$
The coefficient is 64.0701 at $g=106.75$, the high-temperature Standard Model benchmark. This is a deliberately generous conversion of total added energy, including electric, kinetic and coherent energy. A degenerate or dilute heavy gas, or a model with other relativistic degrees of freedom, has a different relationship between energy and entropy. The calculation concerns ordinary radiation entropy and omits bulk, horizon and gravitational entropies.

| Example action, half-step joint-vacancy evolution | $4\times10^{10}$ model | $4\times10^{12}$ model |
|---|---:|---:|
| $S_0$ | $3.609\times10^{48}$ | $1.042\times10^{54}$ |
| Turnaround radius in units of $L$ | 34.5982 | 34.5983 |
| Turnaround upper $S_{\rm eq}/k_B$ | $2.03\times10^{40}$ | $2.53\times10^{44}$ |
| Present Hubble-sphere photons / turnaround upper | $7.89\times10^{47}$ | $6.34\times10^{43}$ |

For the higher-mass half-step joint-vacancy run, converting the *measured* turnaround density gives $2.2031\times10^{44}$ in entropy units. Before the turnaround the calculation records the largest value of $a^4e$ on accepted steps; it gives a larger hypothetical equilibrium entropy of $2.9731\times10^{44}$ at $\tau=6.23855$. The largest recorded $a^2[1+(aH)^2]$ gives $3.0504\times10^{44}$ after the same conversion. These are sampled pre-turn stress and curvature inventories, including the stored energy before the daughters carry it. The actual plasma entropy has yet to be calculated. Even setting $S_0$ to the maximum allowed by the earlier angle-integrated sharp-shell screen on its fixed higher-mass history increases its turnaround upper bound only to $2.9260\times10^{44}$. Reaching today's photon target at the same dimensionless turnaround radius under this ideal conversion would require at least $S_0=2.63\times10^{112}$. Changing the normalization also changes particle couplings and state-count constraints. [Recorded inventories](BUDGET.json)

## Restoring measured gravity and testing a 5 MeV bath

Matching the low-density Friedmann term $2\rho/(\sigma L^2)$ to $\rho/(3\bar M_4^2)$ fixes $S_0=6(\bar M_4 L)^2$. Taking $\bar M_4=2.435\times10^{18}$ GeV, the higher-mass example gives $L=3.3764\times10^{-8}$ m, a physical daughter mass of 23.38 TeV and a physical turnaround radius of $1.1682\times10^{-6}$ m. The finite evolution from the initial slice to turnaround lasts $4.36\times10^{-15}$ seconds in wall proper time. The most favorable physical turnaround radius consistent with the target radiation entropy under the same gravity and gas assumptions would be about **6.02 Mpc**. At the actual hot-radius doubling, a hypothetical full-energy equilibrium conversion gives $T_{\rm eq}\simeq1.18\times10^4$ GeV; the fluid kinetic parameter $\theta m\simeq5.81\times10^7$ GeV has a different meaning. These conversions do not establish thermalization or a measured cosmic temperature.

One can remove the adjustable $S_0$ from the simultaneous energy and volume question. On any expanding slice with measured $G$, nonnegative $m_b$, and a bath drawing at most the wall energy,
$$
\rho_{\rm rad} R_{\rm w}^2\le 3\bar M_4^2[1+(aH)^2],\qquad
\frac{S_{\rm rad}}{k_B}\le
\frac{4\pi\,90^{3/2}}{45}\frac{g_s}{g_\rho^{3/2}}
\left(\frac{\bar M_4}{T}\right)^3[1+(aH)^2]^{3/2}.
$$
Take $T=5$ MeV and $g_\rho=g_s=10.75$ as an illustrative ordinary early-universe bath. Requiring $S_{\rm rad}\ge S_{\gamma,H}$ on the same slice would need $aH\ge1.2401\times10^8$ at this temperature. The archived #192 higher-mass reference has $aH=2.11122$ at the start of its 90%-hot interval. Its slice therefore permits at most $1.071\times10^{65}\,k_B$ in a conventional 5 MeV bath, a shortfall of at least $1.50\times10^{23}$. This hot-onset value comes from the preceding archived trajectory. At hot-radius doubling, the new live half-step joint-vacancy trajectory gives $aH=0.52375$ and a corresponding ceiling of $1.2083\times10^{64}\,k_B$. These are tests of what could fit on the specified slices at 5 MeV. Neither epoch has been identified with nucleosynthesis.

## Entropy carried by the explicit elastic model

At the same higher mass and normalization, a fixed elastic scattering model gives a 90%-hot radius growth of 2.15625 on the native step and 2.15583 on the half step under directional capping. The peak intrinsic random-energy share is 97.62%. The corresponding joint-vacancy runs give the same values to the displayed precision. The table below uses energies measured along each trajectory. Each value is the equilibrium entropy *if* the entire measured added energy were converted into the specified ordinary radiation bath. The fluid solver follows directed and random energy and does not evolve that radiation entropy. [Elastic evolution](../boost/RESULTS.md) · [Trajectory measurements](../boost/EVOLUTION.json)

| Source and step | At hot-radius doubling, $S_{\rm eq}/k_B$ | At turnaround, $S_{\rm eq}/k_B$ | Most favorable turnaround upper, $S/k_B$ | $aH$ at first hot / doubling |
|---|---:|---:|---:|---:|
| Directional cap, native | $2.07919\times10^{44}$ | $2.07333\times10^{44}$ | $2.50115\times10^{44}$ | 1.99253 / 0.40719 |
| Directional cap, half | $2.07916\times10^{44}$ | $2.07332\times10^{44}$ | $2.50112\times10^{44}$ | 1.99202 / 0.40662 |
| Joint vacancy, native | $2.07919\times10^{44}$ | $2.07333\times10^{44}$ | $2.50115\times10^{44}$ | 1.99253 / 0.40719 |
| Joint vacancy, half | $2.07916\times10^{44}$ | $2.07331\times10^{44}$ | $2.50112\times10^{44}$ | 1.99202 / 0.40661 |

The higher-resolution directional cap gives $R_{\rm w}=1.1600\times10^{-6}$ m at the turnaround and $T_{\rm eq}=1.1152\times10^4$ GeV at hot-radius doubling. Today's Hubble-sphere photons carry $7.73\times10^{43}$ times the equilibrium entropy from the modeled turnaround energy. Even the generous turnaround upper requires a subsequent total entropy multiplication of at least $6.41\times10^{43}$. Adiabatic expansion alone conserves total entropy. Before turnaround, the largest sampled equilibrium conversion from the stored stress is $2.9731\times10^{44}\,k_B$ at $\tau=6.23855$, a factor $5.39\times10^{43}$ below the target. [Complete four-run budgets](BUDGET.json)

Independently of the adjustable normalization, keeping a hypothetical 5 MeV bath at the *measured* first-hot value $aH=1.99202$ limits that bath to $9.30\times10^{64}\,k_B$. The minimum entropy multiplication from that limit to today's target is $1.72\times10^{23}$. At doubling the same test gives $1.057\times10^{64}\,k_B$ and $1.52\times10^{24}$ respectively. The geometric requirement to fit the entire target on a single 5 MeV slice remains $aH\ge1.24\times10^8$. These limits compare different physical normalizations at a specified temperature and expansion rate. The actual equilibrium temperature inferred for the chosen normalization is given above.

An interpretation in which a single finite pulse already supplied today's photon entropy followed by adiabatic evolution must clear these limits. This particular wall turns around and contracts, so it has no modeled path to our present Hubble sphere. A later entropy-producing episode, further bulk work, a much larger hot domain, or a different cosmological solution would require its own energy and expansion calculation. [Budget code](budget.py)

From the study root run `python -I -B cosmology/budget.py --base prior --feedback feedback/FEEDBACK.json --evolution boost/EVOLUTION.json --output /tmp/BUDGET.json`.
