# Radiation from collisions through the neutral vector resonance

A charged particle and its antiparticle can radiate a photon while they scatter. We calculated the photon energy from zero up to the maximum allowed by four-momentum conservation, including emission both before and after the pair encounters the model's neutral vector resonance. On the 16,384-pair trajectory through wall time 7.05, the resulting energy is **$1.1011\times10^{-14}$ of the endpoint gas energy**. The four particle and time resolutions gave $1.1011$–$1.1120\times10^{-14}$. This is a calculation for the inherited charge, isolated resonance and recorded binary collision histories. [Numerical results and controls](HARD.json)

| Time | Collisions in the largest sample | Total emitted energy / endpoint gas energy | Photons above 1% of the kinematic endpoint per charged particle |
|---:|---:|---:|---:|
| 6.67 | 51,015 | $2.4222\times10^{-15}$ | $2.2217\times10^{-14}$ |
| 7.05 | 154,876 | $1.1011\times10^{-14}$ | $7.0750\times10^{-14}$ |

At 7.05, photons above a quarter of their event's kinematic endpoint carry about **83.45% of the single-photon energy**. Photons emitted before the collision can lower the pair's energy enough to bring a higher-energy pair back onto the resonance. On the finest late run, 34.30% of the total emitted energy comes from initial-state photons above a quarter endpoint that take a pair from $s>1.30M^2$ to $0.70M^2\le s'\le1.30M^2$. These limits define the quoted resonance window. At an incoming center-of-mass energy of 8,300 daughter masses, the calculated photon count in the 10%–25% band is 77.60% of the value from extending the low-energy formula into that band. At higher pair energies, radiative return can enhance the hard spectrum. Summed over all photon energies and these actual encounters, the calculation gives 1.436 times the energy from extending the angle-averaged low-energy formula to the endpoint. [Spectrum and sampled energies](HARD.json)

## What interaction is being tested?

Work in units of daughter mass $m=1$. Let $s=(p_1+p_2)^2$, with $p_1,p_2$ the incoming fermion and antifermion, and $p_3,p_4$ the outgoing fermion and antifermion. We treat the particles as point Dirac fermions with vector couplings, on shell between interactions. The added photon couples minimally to the same conserved charge on the four external lines. The mediator has zero photon charge. At the real-emission vertices we use the tree diagrams. We retain the phenomenological running width of the elastic calculation,

$$
D(u)=\left[u-M^2+i\,0.30(u+2)\sqrt{1-4/u}\right]^{-1},
\qquad M=8300.
$$

The elastic vector coupling obeys $g_V^2=12\pi(0.30)$. With $b^2=1-4/s$, the spin-averaged squared Born amplitude is $g_V^4 s^2 |D(s)|^2[3-2b^2+b^4\cos^2\theta]$. Both the angular shape and the integrated cross section then reproduce the elastic kernel that generated the event records. A complete same-species vector or gauge theory may also carry crossed-channel exchange. The isolated $s$ channel here keeps the interaction used in those records. Quantities are reported at first order in $\alpha=3Q^2(mL)^2/(4\pi S_0K)$. It equals $1.44382\times10^{-16}$ on these histories. The charged-particle trajectories and their effective charge follow the predecessor model. The physical spin and photon channel are model assumptions. [Collision record and charge convention](../events/README.md)

## Amplitudes, angle integration and energy

Write $k$ for the photon momentum and $\epsilon$ for either physical polarization. Define $x=2k^0/\sqrt{s}$, $s'=s(1-x)$ and $\xi=x/(1-4/s)$. The archived bands use $\xi$, which runs exactly from zero to one for two massive final particles. Incoming radiation is proportional to $D(s')$, while outgoing radiation is proportional to $D(s)$. For example, the two initial diagrams have the numerator

$$
\bar v(p_2)\left[
 \gamma^\mu\frac{\not p_1-\not k+1}{-2p_1\cdot k}\not\epsilon
 +\not\epsilon\frac{\not k-\not p_2+1}{-2p_2\cdot k}\gamma^\mu
\right]u(p_1).
$$

The two final diagrams instead use

$$
\bar u(p_3)\left[
 \not\epsilon\frac{\not p_3+\not k+1}{2p_3\cdot k}\gamma^\mu
 +\gamma^\mu\frac{-\not p_4-\not k+1}{2p_4\cdot k}\not\epsilon
\right]v(p_4).
$$

The two amplitudes interfere at specific outgoing momenta. Exchanging the final fermion and antifermion leaves the squared initial and final terms even and reverses the interference term by charge conjugation. The interference therefore integrates to zero for the charge-blind photon count and energy used here. We evaluate it explicitly before integration in the numerical controls. This cancellation uses complete final-angle acceptance.

The unobserved pair's angle-integrated spin tensor is proportional to the projector transverse to its total momentum. Its spatial coefficient in the pair rest frame is $4(s'+2)/3$ for the final pair, or $(s+2)/3$ for the initial pair including the spin average. The remaining spin sums are evaluated with explicit four-component spinors. The massive collinear regions are resolved with the substitution $\tanh y=\beta\cos\theta$. The upper band's recoil integral is split across the resonance and uses $\sqrt{\log(s'/4)}$ to handle the final pair threshold. With the spin and angle averages so defined, the phase-space factor is

$$
\frac{1}{\sigma_B}\frac{d\sigma_\gamma}{dx}
=\frac{\alpha s x}{4\pi}\frac{\beta(s')}{\beta(s)}
\frac{\langle|\widetilde{\mathcal M}_\gamma|^2\rangle}
     {\langle|\widetilde{\mathcal M}_B|^2\rangle}.
$$

The tildes remove $g_V^2$ from both amplitudes and the photon charge from the radiative amplitude. The energy fraction at each $x$ in the pair center-of-mass frame is $x/2$. After integrating over both charge labels and all outgoing angles, parity and initial charge interchange make the mean photon three-momentum zero in that frame. The emitted four-momentum is consequently proportional to the incoming pair four-momentum. The conformal shell-energy contribution of each recorded encounter is

$$
w a (E_1+E_2)\int_0^{1-4/s}\frac{x}{2}\,\frac{d\sigma_\gamma}{\sigma_B}.
$$

Here $w$ is the weight of one charged simulation particle, $a$ is the shell radius when the collision occurs, and $E_i$ are the incoming energies in the shell frame. Dividing their sum by the measured gas energy at the endpoint gives the table. A simulation record is used once per encounter. The cross section averages over the allowed outgoing angles at the actual incoming $s$ of that encounter. The result estimates radiative energy statistically on the elastic histories. The photon count requires a positive lower energy cut, while the energy integral is finite down to zero.

## Checks and limits

- Contracting the primitive amplitudes with $k^\mu$ gives the photon Ward identities separately for each emitting pair. Applying the external Dirac equations first reproduces the same spin amplitudes.
- An independent integration over all incoming and outgoing angles, with full spin and interference sums at $s=12$ and $x=0.17$, agrees with the reduced initial and final results to $8\times10^{-13}$ relative. The charge-odd interference cancels to rounding error. The soft coefficient of one emitting pair approaches $(2/\pi)[(1+b^2)\operatorname{atanh}(b)/b-1]$ per unit $\alpha$ and per unit $\log x$. The ultrarelativistic initial-state radiator also recovers the familiar $(1+(1-x)^2)(\log s-1)/(\pi x)$ times the shifted Born cross section. [Matrix-element checks](emission.py)
- We compare 20- and 40-node quadrature on seven invariant energies, including the fine late tape's approximate median, 90th percentile, 99th percentile and maximum. For the 25%–100% band, those calculations agree to $3\times10^{-8}$ relative in both energy and photon count. The inclusive rates are tabulated on 257 points in $\log s$ across the complete eight-tape domain and interpolated locally with four points. Removing every other tabulation point changes the folded energy and count by less than two parts per million on every tape. Separate outgoing and incoming energy grids for $\xi>0.25$ are archived for reuse. [Numerical controls and grids](HARD.json)
- The first-order single-collision calculation treats radiation as incoherent from one encounter to the next and the photon as freely propagating afterward. Finite time in the evolving fields, a plasma dispersion relation, absorption, other scattering channels, secondary cascades, multiphoton probabilities and gravitational response require additional physics. Virtual radiation enters an infrared-safe inclusive *event rate* at the same perturbative order. This first-order energy observable counts the energy carried by a real photon and is infrared finite. [Independent analytic radiator comparison](https://arxiv.org/abs/hep-ph/0505236)

Run the following command after the study has regenerated its binary collision records.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B hard/emission.py \
    --artifacts /path/to/events --output /tmp/HARD.json --jobs 2
~~~
