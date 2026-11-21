# Can a specified elastic interaction relax the current?

The earlier fluid evolves its charged-particle momentum with a prescribed relaxation rate $L\nu=18$. In its $mL=4\times10^{12}$ history an isolated spin-zero elastic resonance with coupling $y^2/(4\pi)\le1$ falls short even at the largest action allowed by the earlier angle-integrated, sharp-birth screen. A separate broad vector resonance, with one fixed mass and width, supplies enough momentum relaxation throughout the sampled hot, low-drift expansion in a dilute thermal calculation. The restricted calculation below tests the particles after they have gained random motion. A companion calculation evolves a massive vector collision rate from the first births through the expanding hot interval, assuming two internally thermal particle streams. [Full-history evolution](../boost/EVOLUTION.json) [Earlier hot-particle calculation](https://github.com/maldonado-research/HDblast/pull/192)

## The thermal average

Consider one charged Dirac flavor and its opposite charge. Let the two isotropic Maxwell–Boltzmann populations have a common temperature $T$, infinitesimal opposite drift speeds, and proper number density $n_t$ per sign. For an azimuthally symmetric elastic channel, $\sigma_T=\int(1-\cos\theta_*)\,d\sigma$ with the scattering angle in the center-of-momentum frame. The linear relaxation rate of each sign's mean momentum in the wall frame is

$$
\nu=\frac{n_t}{8T^2}\left\langle (1-\cos\vartheta)\sigma_T(s)(p_{1z}-p_{2z})^2\right\rangle
=\frac{\pi n_t}{T^2}\int_0^\infty w(z)F(zT)\,dz .
$$

Here $s=2p_1p_2(1-\cos\vartheta)$, $z=\sqrt{s}/T$, and $F(\sqrt{s})=s\sigma_T(s)/(16\pi)$. Averaging the squared longitudinal momentum difference over overall orientation supplies a factor $1/3$. The two incoming radial weights are $x^2 e^{-x}/2$, $y^2 e^{-y}/2$, with $x=p_1/T$, $y=p_2/T$ and a uniform relative cosine. The ratio to the spin-independent $16\pi/s$ ceiling has integrand $(x/y+y/x-2\cos\vartheta)F/3$ against that joint probability. Changing integration variables to $v=x+y$, $d=x-y$ and $z$ gives $v\ge z$ and $|d|\le\sqrt{v^2-z^2}$. Integrating $d$ and then $v$ yields

$$
 w(z)=\frac{z^3 K_2(z)+z^4 K_1(z)}{24},\qquad
 \int_z^\infty w(t)\,dt=\frac{z^4 K_2(z)+3z^3 K_3(z)}{24}.
$$

The kernel integrates to one and its second moment is 20. Direct quadrature in the two momenta and their relative angle reproduces the reduced integral within $1.3\times10^{-5}$ in three tests. This weights the loss of directed momentum per particle. Counting encounters uses a different thermal kernel. The derivation assumes locally field-free scattering, random spin, dilute particles, $T\gg m$ and a small relative drift. In the reference hot interval $2205<T/m<5553$ with $|v|<1.90\times10^{-4}$. The inferred massless Maxwell–Boltzmann fugacity per sign, $\pi^2 n_t/(2T^3)$, stays below $1.62\times10^{-14}$ on this thermal ansatz. Occupation of individual newly created, cold states remains a separate issue. [Drift moment and normalization](https://github.com/maldonado-research/HDblast/blob/dc99ff5c5894abcff3e11e02e4b97e19da31de00/research/HDBLAST_NEUTRAL_PARENT_DECAY_20261120/consistency/CONSISTENCY.md) · [Thermal checks](RESONANCE.json)

## Two elastic poles

Take one pole per flavor, with the same parameters for all three flavors and one open elastic decay channel per pole. The running imaginary part of its propagator is $\gamma s$, with $\gamma=\Gamma_R/M_R$ on shell. The dimensionless resonance factor is

$$
 B(s)=\frac{\gamma^2s^2}{(s-M_R^2)^2+\gamma^2s^2},\qquad 0\le B\le1 .
$$

For a neutral scalar $y\phi\bar f f$, a spin average over two incoming spin states per sign and a sum over the final spins gives $\sigma_T=y^4s/[16\pi((s-M_R^2)^2+\gamma^2s^2)]$, with $\gamma=y^2/(8\pi)$ in the massless, single-channel width. Thus $F=B/4$. For a massive neutral vector $g_VV_\mu\bar f\gamma^\mu f$ the spin-averaged $s$-channel pole gives $\sigma_T=g_V^4s/[12\pi((s-M_R^2)^2+\gamma^2s^2)]$, $\gamma=g_V^2/(12\pi)$, and $F=3B/4$. The vector differential distribution is proportional to $1+\cos^2\theta_*$. Its first angular moment vanishes, making its transfer and total cross sections equal. The spin factors reflect the accessible resonance polarizations. Each example adds a new interaction; the weak classical gauge coupling that supplies the electrical force is fixed separately. A spin-one pole and an isotropic spin-independent $s$ wave have different angular distributions. These expressions specify an elastic, single-channel pole model; additional decay channels lower the elastic branching fraction. A literal field theory also has crossed diagrams, medium effects, radiative corrections and other possible reactions. The calculation retains just the stated pole. Couplings near unity and a 30%-wide pole warrant that qualification. These new mediators are independent of the neutral parent that creates the charged daughters. [Particle Data Group review](https://doi.org/10.1103/PhysRevD.110.030001)

For the scalar coupling range $y^2/(4\pi)\le1$, increasing $\gamma$ raises $B$ everywhere, so it suffices to test $\gamma=1/2$. A search over $M_R/T$ gives a maximum thermal $\langle B\rangle$ near 0.42979 at $M_R/T=3.166$. A separate upper enclosure covers **every** $0\le M_R/T<\infty$. It partitions the thermal measure using the Bessel tail above and bounds $B$ on each rectangle in $(z,M_R/T)$. The tails in $z$ are allowed to contribute at full strength. The bounds use double-precision Bessel values. Doubling both resolutions lowers the numerical upper from 0.43332 to **0.43156**, before the factor $1/4$. At the previous trial action $S_0=1.04156\times10^{54}$ the scalar pole can supply at most **53.95%** of the required rate at the most demanding state. Even raising $S_0$ to $1.26590\times10^{54}$, the previous collisionless, angle-integrated sharp-shell upper action on that history, raises this fraction only to **65.57%**. This excludes the isolated scalar-pole realization of the assumed drag under these conditions. It gives no exclusion for a complete theory with other collision channels or for other histories.

For the vector example, a scan guided the selection $M_V/m=8300$, $\Gamma_V/M_V=0.30$. The corresponding $g_V^2/(4\pi)=0.90$. These values remain fixed as the model cools. A read-only rerun bins every accepted hot, low-drift state in 1% temperature intervals, recording the extrema of the required action in each. The selected action is five times the reference drag floor. Bounding the thermal integral over each entire temperature bin encloses $\nu_{\rm vector}/(18/L)$ between **1.0505 and 1.5505** across the two resolutions. Ratios sampled at the greatest required action in each bin lie between 1.0717 and 1.5133. The archive reports both to distinguish a conservative interval from sampled extrema. The excess over the target varies with temperature. Matching the original constant drag exactly is a different requirement. [Thermal-history bins](THERMAL.json) · [Elastic-pole analysis](resonance.py)

## Putting the rate back into the wall equations

The live test replaces the relaxation coefficient consistently in the directed-momentum equation, the gas-energy reconstruction and the drag-heating ledger. At every Runge–Kutta stage the vector rate uses the evolved temperature and density. The spin-weighted integral is tabulated on 1,025 logarithmic temperature nodes and linearly interpolated in log temperature. The finer comparison halves the time step and doubles the number of table intervals. All runs start with zero neutral and charged excitation and keep the original finite bulk and gravitational equations.

| Charged birth law | Drag protocol | Peak intrinsic energy fraction | Radius growth while at least 90% hot |
|---|---|---:|---:|
| Original unrestricted | Switch to the vector at first hot, slow-drift state | 97.743% | 2.25286 |
| Directional hard cap | Same hot-state switch | 97.743% | 2.25286 |
| Independent pair vacancies | Same hot-state switch | 97.743% | 2.25286 |
| Independent pair vacancies | Turn drag off at that state | 97.741% | 2.25286 |
| Independent pair vacancies | Use the vector in its thermal domain and the old drag outside | 97.742% | 2.25206 |
| Independent pair vacancies | Use the vector in its thermal domain and zero drag outside | 55.651% | No 90%-hot interval |
| Same zero fallback, finer integration | Same | 55.642% | No 90%-hot interval |

The first hot-state switch occurs near $\tau=8.147$. With the vacancy source, the live vector rate ranges from $19.28/L$ to $27.24/L$ during the remaining expansion. Turning the drag off at the same time gives nearly identical radius growth. By this switch the wall had already accumulated 99.999% of its eventual signed electrical work. Sustaining a previously heated gas gives little information about the scattering that first heated it.

The last three rows allow the vector to act before the gas becomes hot. Write $s(x)=3x^2-2x^3$ for $0<x<1$, extended by 0 and 1 beyond the endpoints. The vector weight is $w=s((\Theta-100)/100)s((0.1-|v|)/0.05)$. The local drag is $w\nu_{\rm vector}+(1-w)\nu_{\rm outside}$. This smooth choice keeps every evaluation of the linear ultrarelativistic law inside its tested domain. With $\nu_{\rm outside}=18/L$ the wall heats. With $\nu_{\rm outside}=0$ the gas misses the hot threshold. In the latter run the ratio $\int w|a^3 JX|d\tau/\int|a^3 JX|d\tau$ is only **0.0181%**. The absolute-work denominator counts repeated forward and reverse electric energy exchange. This quantity identifies where the rate approximation acted; it does not assign electrical work to scattering. The zero fallback tests a restricted closure. The actual physical cross section of a massive vector at cold or fast drift requires the massive, relative-stream collision integral. Substituting zero there cannot rule out that interaction. [All continuation runs](LIVE.json) · [Joint wall and gas evolution](live.py)

The companion full-temperature calculation uses a separate reduction of the moving-stream integral. At $T/m=10$ and stream rapidity 4, twelve independently scrambled sequences of $2^{21}$ momentum-pair samples each give $5.19678\times10^{-7}$ for $m^2\nu/n_t$. The reduced integral gives $5.19726\times10^{-7}$. Their difference is 0.0093%, with a relative standard error across scrambles of 0.0099%. [Nested sampling and uncertainty](CONVERGENCE.json)

From the capsule root, reproduce the products with

```sh
python -I -B scattering/history.py --base prior --output /tmp/THERMAL.json
python -I -B scattering/resonance.py --history /tmp/THERMAL.json --sweep prior/kinetics/SWEEP.json --output /tmp/RESONANCE.json
python -I -B scattering/live.py --base prior --output /tmp/LIVE.json --jobs 2
python -I -B scattering/convergence.py --output /tmp/CONVERGENCE.json --jobs 2
```

The history runner checks the pinned prior solver signature, the daughter rest, the hot interval, and the original drag floor. The live runner verifies pair-source accounting and energy conservation. The thermal fit, the choice of mediator, and the switch protocol are exploratory. The collision calculation uses masses and lengths in model units. A separate estimate matches the low-density gravitational coupling to its measured value. [Conditional physical scale](../cosmology/BUDGET.md)
