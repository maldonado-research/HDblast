# How much energy can the matched pulse put into this spectator?

The quiet B1 wall and its two signed pulse histories keep the spectator mass near 50 times the initial Hubble scale through proper time $s=H_0\tau=0.15$. For a single freely propagating, conformally coupled real scalar in the specified initial states, all initially admitted modes together have an instantaneous excitation-energy ceiling of $27.60H_0^4$ at the final cut on the fine grid. This bound covers the expansion and pulse together. It holds for the recorded, interpolated backgrounds and a fixed set of comoving modes. The actual mode calculation retains their phases and gives a much sharper estimate.

## Physical frequency and the B1 cutoff

Write $\mu=m/H_0=G|\phi-\phi_*|$, $G=100$, $\phi_*=0.5$, $a(0)=1$, and express comoving $k$ in units of $H_0$. B1 registers a cutoff on the **maximum mass and characteristic production momentum** $\sqrt{q}$, with $q=G|d\phi/ds|$ at a crossing. Its gas uses quadrature nodes across a momentum distribution. For this audit we impose an additional per-mode energy screen $E_k=\sqrt{\mu^2+(k/a)^2}\leq\Lambda$, where $\Lambda=M_5/H_0=b^{-1/3}=53.6107358$ at $b=6.49\times10^{-6}$. The distance to this scale supplies no estimate of unknown higher-order operators. [B1 registration, §§1.2–1.4](https://github.com/maldonado-research/HDblast/blob/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/REGISTRATION.md)

The initial mass is $49.9382910$. The largest eligible $k$ is $K_*=\sqrt{\Lambda^2-\mu(0)^2}=19.5007200$ at the initial slice. The replay bounds $a\sqrt{\Lambda^2-\mu^2}$ over every segment of the same cubic Hermite interpolant used for the modes. For all three histories on both grids its minimum is the initial value. Thus the full $5\leq k\leq18$ band satisfies the screen. At $k=18$ the initial energy uses 99.016% of the scale, while $k=20$ already starts at $53.79436H_0$. A wider-band continuum calculation needs separate interpretation. [Matched pulse histories](https://github.com/maldonado-research/HDblast/tree/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_INCIDENT_PULSE_20261118/evolution)

## Frequency variation and preparation

In dimensionless conformal time $d\zeta=ds/a$, the canonical variable $v_k=a\chi_k$ obeys $v_k''+\Omega_k^2v_k=0$, with $\Omega_k^2=k^2+(a\mu)^2$. Conformal coupling cancels the $a''/a$ term in this spatially flat geometry. With $h=d\log a/ds$ and $\dot\mu=d\mu/ds$, the first frequency-variation indicator is

$$
\epsilon_k=\frac{|\Omega'_k|}{\Omega_k^2}=\frac{|\mu(h\mu+\dot\mu)|}{[\mu^2+(k/a)^2]^{3/2}}.
$$

The table reports extrema over recorded fine-grid nodes through $s=0.15$. Expansion contributes to the indicator even on the quiet history.

| Fine-grid history | Sampled mass range, $H_0$ | Max $\lvert\dot\mu\rvert/\mu^2$ | Max $\epsilon_5$ | Max $\epsilon_{18}$ |
|---|---:|---:|---:|---:|
| Quiet | 49.93829–49.93834 | $2.03\times10^{-7}$ | 0.01980 | 0.01745 |
| Positive pulse | 49.83339–49.93832 | 0.00465 | 0.02354 | 0.02066 |
| Negative pulse | 49.93829–50.04321 | 0.00464 | 0.02436 | 0.02122 |

An instantaneous oscillator vacuum has $v'_0=-i\Omega_0v_0$. The local zeroth-order WKB mode uses $v'_0=(-i\Omega_0-\Omega'_0/(2\Omega_0))v_0$. Both take $v_0=1/\sqrt{2\Omega_0}$. Projected into the instantaneous convention, the WKB tangent already has $n_{k,0}=\epsilon_{k,0}^2/16$, from $2.43\times10^{-5}$ at $k=5$ to $1.74\times10^{-5}$ at $k=18$. This difference arises at preparation on the shared initial slice. These finite-band conventions alone give no renormalized stress tensor.

## Phase-independent upper bound

For the instantaneous basis, the oscillator evolution gives $|d|\beta|/d\zeta|\leq |\Omega'/(2\Omega)|\sqrt{1+|\beta|^2}$. Integrating bounds $\operatorname{arsinh}\sqrt{n_k}$ by its initial value plus half the total variation of $\log\Omega_k$. The cubic-interpolation check bounds $d\log(a\mu)/dT>0.755$ on every interval of all six histories. Thus at $s_f=0.15$,

$$
n_k(s_f)\leq\sinh^2\!\left[\operatorname{arsinh}\sqrt{n_{k,0}}+\tfrac12\log\frac{\Omega_k(s_f)}{\Omega_k(0)}\right].
$$

The replay bounds a cubic and its derivative by their Bernstein control values on every saved time step. It also bounds the derivative of the allowable momentum ceiling on each interval. Where the derivative bound has either sign, it uses interval lower bounds on $a$ and $\sqrt{\Lambda^2-\mu^2}$ instead. These are floating-point bounds on the declared interpolation. Evolution between samples in the underlying continuum bulk theory remains subject to the spatial and temporal discretization errors.

The right side decreases with $k$ for both initial preparations here. For a band $[k_0,k_1]$, let $n_{\max}$ be its value at $k_0$. A conservative bound for the finite-band energy above the instantaneous oscillator ground state is $\hat\rho_{\rm exc}\leq \Lambda(k_1^3-k_0^3)n_{\max}/(6\pi^2a_f^3)$ for one real degree of freedom, with $\hat\rho=\rho/H_0^4$. Using the larger bound from the WKB tangent covers both preparations.

| Fine-grid history | Ceiling $\hat\rho_{\rm exc}$, $5\leq k\leq18$ | Ceiling $\hat\rho_{\rm exc}$, $0\leq k\leq K_*$ | B1 radiation term / recorded $h_f^2$ if all the latter energy converted at the cut |
|---|---:|---:|---:|
| Quiet | 20.755 | 27.460 | $1.629\times10^{-4}$ |
| Positive pulse | 20.654 | 27.328 | $1.621\times10^{-4}$ |
| Negative pulse | 20.854 | 27.592 | $1.637\times10^{-4}$ |

Each tabulated ceiling rounds upward. The last column applies the existing B1 ledger $[r_b\sigma(\phi_f)Y/18+Y^2/36]/h_f^2$, with $r_b=7.838285538073243$ in the code's units and $Y=b\hat\rho$. It measures the size of a hypothetical complete conversion at that instant while holding the geometry fixed. The calculation covers initially admitted modes and the specified preparations, under free evolution. Additional species, modes entering a moving cutoff later, a different initial state, continual dissipative transfer, or gravitational backreaction require separate dynamics. [B1 Friedmann and radiation ledger](https://github.com/maldonado-research/HDblast/blob/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py)

## Production and thermal interpretation

The original mass-zero trigger at $\phi_*=0.5$ has no crossing on these histories. The corresponding B1 gas remains empty. In the separate stable free-field model, a time-dependent mass can still change the state. B1's $y=1$ fermion-decay prescription would give an optical depth $\int\mu^2/(8\pi E_k)\,ds$ of 0.282–0.297 over this interval for an already-present excitation in the primary band (nodal trapezoidal estimate). Applying this decay probability to a coherent vacuum amplitude needs an interacting model. The B1 massless daughter ledger itself contains no thermalization dynamics. [B1 gas and decay equations](https://github.com/maldonado-research/HDblast/blob/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/shell_matter.py)

As an explicitly post-hoc sensitivity check, moving the mass zero to $\phi_*=0.0012$ would produce a pair of near-wall passes in the fixed positive-pulse history. Fine-grid nodal interpolation gives $s=0.10777,0.13304$, $q=8.98,8.29$, and nominal production durations $1/\sqrt q=0.334,0.347$. The separation is only 0.0757 of the first duration. The outward field excursion reaches 0.0156 of the first nonadiabatic width $\sqrt q/G$. An estimate that treats these as independent asymptotic single passages lacks the required separation. This sensitivity changes the B1 model parameter. The rolling B1 control with its original threshold is a different experiment.

## Replay

Use Python 3 with its standard library and the two files from the published pulse packet. The script checks both full input SHA-256 digests and runs both grids and all three histories. [B1 summary and archived histories](https://github.com/maldonado-research/HDblast/tree/fe35172f7c94dbd616c677045d59c6616c194d9b/research/HDBLAST_INCIDENT_PULSE_20261118)

```sh
python -I -B validity.py /path/to/evolution/wall_histories.zip /path/to/B1_summary.json > /tmp/validity-replay.json
cmp VALIDITY.json /tmp/validity-replay.json
```
