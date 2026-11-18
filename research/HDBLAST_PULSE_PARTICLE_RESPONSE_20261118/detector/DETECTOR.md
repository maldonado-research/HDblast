# Can a comoving probe distinguish the pulse and quiet histories?

A switched detector offers a way to ask the question without defining particles at a final cosmological time. For the chosen detector, the positive pulse increases its excitation probability by 0.3030% relative to the quiet control. The negative pulse decreases it by 0.3025%. Both changes include the effect of the time-dependent spectator mass on the field encountered by the detector.

## Experiment

Use dimensionless wall proper time $s=H_0\tau$, with $ds=e^{B_{\rm wall}}dT/\rho_b$, and identify comoving momenta by fixing $a(0)=1$ in every history. The shell is spatially flat. The spectator is a real, conformally coupled scalar with mass $m/H_0=100|\phi-0.5|$. Its canonically normalized mode $v_\kappa=a\chi_\kappa$ satisfies

$$
v_\kappa''+\omega_\kappa^2 v_\kappa=0,\qquad
\omega_\kappa^2=\kappa^2+a^2(m/H_0)^2.
$$

Primes denote derivatives in dimensionless conformal time $d\eta=ds/a$. All primary histories start in the same Gaussian state, with $v_0=(2\omega_0)^{-1/2}$ and $v'_0=-i\omega_0v_0$. For a comoving two-level detector with dimensionless gap $\delta=20$ and switching $S(s)=\sin^4[\pi(s-s_l)/(s_r-s_l)]$ inside its window and zero outside, the leading excitation probability from band $B$ is

$$
\frac{P_B}{\lambda^2|\mu_{eg}|^2}
=\int_B \frac{\kappa^2 d\kappa}{2\pi^2}
\left|\int_{s_l}^{s_r} ds\, S(s)e^{-i\delta s}\,\frac{v_\kappa(s)}{a(s)}\right|^2 .
$$

Here $\lambda$ is the weak detector coupling and $\mu_{eg}$ its transition matrix element. This expression follows by inserting the vacuum mode expansion into the two-time Wightman response. It integrates the complex amplitude before taking its modulus, retaining phase interference. The primary experiment selects $5\le\kappa\le18$ using the same comoving filter on each history. For the two principal histories, the corresponding physical mode frequencies through $s=0.15$ have an upper bound of $53.084H_0$, below the specimen scale $M_5=53.611H_0$. This screen applies to the field modes. The finite switch itself has a broad frequency response. Such a filter describes an engineered spatially extended probe.

## Response and source controls

For $0.10\le s\le0.15$, the quiet normalized excitation probability is $1.7550463\times10^{-4}$. Values below use the finest available pulse history and subtract that quiet value. The early window $0.04\le s\le0.09$ gives zero contrast at the reported floating-point precision on this grid.

| Pulse applied to the bulk | Change in normalized probability | Change relative to quiet |
|---|---:|---:|
| Positive | $+5.317768\times10^{-7}$ | +0.302999% |
| Negative | $-5.308462\times10^{-7}$ | −0.302469% |
| Half positive | $+2.657722\times10^{-7}$ | +0.151433% |
| Half negative | $-2.655396\times10^{-7}$ | −0.151301% |

With the positive-pulse mass history on the quiet geometry, the change is $5.314367\times10^{-7}$. With the pulse geometry and quiet mass, it is $3.388400\times10^{-10}$. Their mixed remainder is $1.2522\times10^{-12}$. These switches of inputs are diagnostic counterfactuals, evaluated at equal proper time. The response at this setting is dominated by the mass change.

Transporting the initial conditions with leading WKB modes on each frequency history produces $3.969865\times10^{-7}$, about 74.7% of the exact contrast. This illustrates how the evolving frequency alone changes what the detector sees. Using the WKB slope in the initial state for the full evolution gives +0.30464%. Neither the detector probability nor its change specifies an asymptotic occupation, a temperature, or deposited heat. The state, switch, filter and coupling are choices, and the spectator's backreaction is omitted.

## Verification and replay

Across the four underlying history grids the positive contrast in units of $10^{-7}$ is 5.30200, 5.31401, 5.31701 and 5.31777. Doubling the mode integration resolution gives 5.3177714. An independent standard-library Runge–Kutta integrator with Simpson momentum quadrature gives 5.3177736. A constant-scale-factor oscillator matches the analytic switched amplitude within $2.0\times10^{-16}$. For an expanding $a=e^s$ control with constant mass at $\kappa=0$, doubling the step count reduces the absolute amplitude error from $1.53\times10^{-10}$ to $3.84\times10^{-11}$.

Use the prior incident-pulse package's authenticated evolution/wall_histories.zip and B1_summary.json as the two positional arguments below. [Prior data release](https://doi.org/10.5281/zenodo.24354840)

    python -I -B detector.py PATH/TO/evolution/wall_histories.zip PATH/TO/B1_summary.json > detector-output.json
    python -I -B replay_detector.py PATH/TO/evolution/wall_histories.zip PATH/TO/B1_summary.json > replay-output.json

The first program requires NumPy. Its checked output is in DETECTOR.json. Arrays of detector probabilities index initial preparation, switching window and energy gap in the order given at the top of the JSON. It also records exploratory wider-band and spatially Gaussian responses extending above the conservative frequency screen. Those figures depend on continuing the free spectator theory beyond that scale. The second program reproduces the screened primary result and relative squeezing at $\kappa=10$. Its checked output is in REPLAY.json.
