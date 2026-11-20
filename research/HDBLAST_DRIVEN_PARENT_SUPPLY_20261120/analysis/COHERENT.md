# What the homogeneous field can supply

A driven, real scalar can store a macroscopic amount of energy in a mode with zero spatial momentum. Whether that mode acts like cold matter depends on its pressure and on how quickly its mass and the geometry change. A count made at a single time and an averaged cold-matter history answer different questions.

## Whole-sphere normalization

Use $c=\hbar=1$, dimensionless time $\tau=t/L$, dimensionless radius $a$, wall tension $\sigma$, and $S_0=\sigma L^4$. The dimensional parent field is $\Phi=\sqrt{\sigma}L\chi$ and its dimensionless mass is $\omega=ML$. On a unit round three-sphere of volume $2\pi^2$, the homogeneous, minimally coupled action reads

$$
S_\chi=2\pi^2 S_0\int d\tau\,a^3\left[\frac{\dot\chi^2}{2}-\frac{\omega^2}{2}(\chi-v(\psi))^2\right].
$$

Here dots are derivatives with respect to $\tau$. Let $q=\chi-v$, $e_\chi=(\dot\chi^2+\omega^2q^2)/2$, and $p_\chi=(\dot\chi^2-\omega^2q^2)/2$. The canonical momentum is $2\pi^2S_0a^3\dot\chi$. Diagonalizing the untransformed Hamiltonian on a specified time slice gives a coherent displacement occupation

$$
N_{\rm coh}=2\pi^2 S_0 I,\qquad I=\frac{a^3 e_\chi}{\omega}.
$$

This expression is exact for the mean-field contribution in that instantaneous oscillator convention. It uses $\dot\chi$, including during a moving-minimum drive. For a quantum state, split its oscillator coordinate and momentum into their means and centered fluctuations. The quadratic expectation of the number operator adds the coherent term to the fluctuation occupation, with the vacuum half-quantum subtracted from the latter. The driven calculation follows the means and classical stress. It never assumes a measurement returned an integer particle count.

The comparison fluid starts at $a_i=2$ with rest density $D_0=10^{-4}$ in wall-tension units, $mL=1$ and $M_iL=1.9$. Its target is $I_{\rm target}=a_i^3D_0/(M_iL)=4.2105263\times10^{-4}$ or $N_{\rm target}=0.0083112458\,S_0$ on the whole sphere. At the opening the driven calculation gives $N_{\rm coh}=0.0139090\,S_0$, a ratio of **1.67352**.

The source-only calculation sets the charge response to zero. Absolute number estimates based on a hypothetical charged descendant carry an additional normalization assumption. In the local surface-dominated convention $\alpha=3Q^2(mL)^2/(4\pi S_0 K)$, $\alpha<1$ on the late plateau requires the following values. Earlier portions of the trajectory can have a different $K$.

| Daughter model used for the normalization | $Q$ | $K$ | Minimum $S_0$ (strict lower bound) | Minimum $N_{\rm target}$ |
| --- | ---: | ---: | ---: | ---: |
| Earlier single-stage reference | $10^8$ | 13200.04 | $1.80857\times10^{11}$ | $1.50315\times10^9$ |
| Two-stage hot-growth reference | $2\times10^8$ | 26.04 | $3.66716\times10^{14}$ | $3.04787\times10^{12}$ |

The ratio of coherent action to the benchmark cancels $S_0$. The earlier vacuum-production calculation propagated a free **conformally coupled** real scalar about zero mean on an already seeded trajectory. Its $n=0$ canonical mode $u\propto a\chi$ has conformal frequency $\sqrt{1+a^2\omega^2}$. Minimal coupling, a moving minimum and gravitational response change the system. The free-vacuum count on that earlier trace cannot be added to the coherent count here. A separate forced conformal calculation tests the relation on the earlier fixed geometry. [Earlier vacuum-production study](https://github.com/maldonado-research/HDblast/pull/189)

## Energy source and pressure

The instantaneous equation $\ddot\chi+3H\dot\chi+\omega^2q=0$ gives

$$
\dot e_\chi+3H(e_\chi+p_\chi)=\omega\dot\omega q^2-\omega^2q\dot v
   = U_{,\psi}\dot\psi .
$$

The driver equation carries the force $-U_{,\psi}$, so its work cancels this exchange in the combined matter energy. The term $\omega^2v^2/2$ is part of the potential; it contributes driver force while the minimum moves. At the last primary sample, a short solver step past the opening, moving-minimum work is **0.00123858**, varying-mass work is **0.0000165980**, and the integrated pressure work $3\int H a^3 p_\chi\,d\tau$ is **−0.000154125**. These give **0.00140930** in comoving parent energy $a^3 e_\chi$. Multiply by $2\pi^2S_0/L$ for dimensional whole-sphere energy. The comoving residual is below $7\times10^{-17}$ in this run.

A direct differentiation gives two finite-time identities,

$$
\dot I=-\frac{a^3}{\omega}\left[(3H+\dot\omega/\omega)p_\chi+\omega^2 q\dot v\right],\qquad
2a^3p_\chi=\frac{d}{d\tau}(a^3q\dot\chi)+a^3\dot v\dot\chi .
$$

Once $v$ stops, the second identity makes the comoving-volume-weighted pressure integral exactly zero between two extrema with $\dot\chi=0$. Expansion work also weights pressure by $H$. Changing $H$ can therefore change the action across an entire period even when the volume-weighted pressure integral vanishes. To express the dependence explicitly, define $G=a^3q\dot\chi$ and $w=(3H+\dot\omega/\omega)/\omega$. Then $|G|\le I$, $\dot I=-w\dot G/2$, and the corrected action $J=I+wG/2$ satisfies $\dot J=\dot w G/2$. For $|w|<2$ one has $|d\ln J/d\tau|\le|\dot w|/(2-|w|)$. Right at the gate, where $w$ is large, this estimate ceases to control the evolution.

For the primary preparation the available unforced interval before the mass starts moving spans **1.339** cycles of the phase $\int\omega\,d\tau/(2\pi)$. It contains one full period bounded by field extrema, from $\tau=1.57233$ to $4.94170$. The comoving-volume-weighted pressure divided by similarly weighted energy is $-7.1\times10^{-9}$ there. The ordinary proper-time average ratio is **−6.36%**, the expansion-weighted one is **+3.74%**, and $I$ drops **12.03%** across the period. Across the entire source-off, constant-mass interval the sampled $I/I_{\rm target}$ ranges from **1.082** to **2.167**.

At the gate $p_\chi/e_\chi=0.7652$, $|H|/\omega=0.1715$ and $|\dot\omega|/\omega^2=5.1852$. In the canonically normalized quadratures $(\sqrt{a^3\omega}\,\chi,\sqrt{a^3/\omega}\,\dot\chi)$ the off-diagonal creation/annihilation mixing relative to oscillator frequency is $|3H/\omega+\dot\omega/\omega^2|/2=2.8499$ at this time. The instantaneous occupation remains defined. The slowly varying, pressureless particle approximation is uncontrolled on that slice. A gas of decaying, on-shell parents and its daughter distribution require additional microphysics.

## Following the source with decays disabled

When the parent is held stable after its channel opens, the same $mL=1$ trajectory can be followed to its first maximum radius. $I/I_{\rm target}$ is **1.60040** at the end of the mass change. From the first saved constant-mass state at $\tau=5.78$ to the last saved expanding state at $\tau=23.70$, its proper-time mean is **1.35060**, and its sampled extremes are **1.14566** and **1.59356**. The field remains coherent and its instantaneous pressure alternates sign. The 11 completed half-periods include five full periods between like extrema. Across those five periods the proper-time pressure-to-energy integral ratio is **−1.215%**, the expansion-weighted ratio is **−1.145%**, and the action rises **2.514%**. The comoving-volume-weighted ratio is zero by the virial identity. Reconstructing from saved rows gives $1.4\times10^{-6}$ for that last ratio, with magnitude below $1.6\times10^{-5}$ after doubling the sample spacing. Direct half-period integrals accumulated by the solver give absolute ratios below $2\times10^{-11}$.

The action at turnaround can be bounded without resolving its phase at that instant. In a constant-mass source-free segment, $d\ln I/d\ln a=-3p_\chi/e_\chi$, and $|p_\chi|\le e_\chi$. For either stored state $j$ on a monotonic segment adjoining the maximum,
$$
 I_j(a_j/a_{\max})^3\le I_{\rm turn}\le I_j(a_{\max}/a_j)^3.
$$
Intersecting the bounds from the neighboring expansion and contraction samples gives the intervals below, conditional on the simulated radii and energies. The heavier case has a narrower drive and a smaller amplitude. It is a separately prepared example. Its target action uses the same reference rest density and is smaller by a factor of 16.

| Stable parent continuation | Opening $I/I_{\rm target}$ | $\lvert\dot\omega\rvert/\omega^2$ at opening | Complete full periods before turnaround | Largest half-period $\lvert\int p\,d\tau/\int e\,d\tau\rvert$ | $I_{\rm turn}/I_{\rm target}$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| $mL=1$, amplitude 0.01, width 0.15 | 1.674 | 5.185 | 5 | 3.03% | 1.351515–1.351525 |
| $mL=16$, amplitude 0.0013, width 0.075 | 1.855 | 0.483 | 82 | 0.0601% | 1.822724–1.822727 |

The heavier example completes 164 half-periods before the wall turns. At their endpoints its action spans **1.82223–1.82272** times the target. Those extrema and signed integrals use the solver's own fine time steps. Their precision does not rely on resolving each half-period in the saved samples. This supports a cold-matter-like *averaged stress* for the stable field across many oscillations. Allowing real decays can change the field, pressure, geometry, and all of these numbers. [Coupled histories](../coupled/PARENT.json)

Run `python -I -B analysis/coherent.py --output /tmp/coherent.json` from the study root. Radius and field between saved samples are reconstructed in cubic Hermite form using their supplied velocities. The phase averages use four-point Gauss–Legendre quadrature on each saved time interval. Velocity zeros define complete periods. [Code](coherent.py) · [Results](COHERENT.json)
