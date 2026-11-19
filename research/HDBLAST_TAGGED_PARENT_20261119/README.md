# Count fermion pairs from one extra scalar on the rolling shell

An explicit scalar–fermion interaction offers a test that reads the daughter particles directly. We put one extra scalar at rest on each of four previously simulated shell histories and computed the additional massless fermion–antifermion pairs at weak coupling. Against the quiet shell, the nearby positive front **reduces** the pair count by 0.0674%, the nearby negative front **increases** it by 0.4703%, and the distant-compensated positive front **increases** it by 0.0616% for the wider interaction window. The same signs occur for both switching windows, both bulk meshes and a common daughter energy cutoff. The result is an exploratory follow-up to the energy-accounting study. [Archived study](https://doi.org/10.5281/zenodo.24356257) · [Project discussion](https://github.com/maldonado-research/HDblast/issues/157)

| Wider switch, fine mesh | Excess pairs divided by $y^2$ | Final daughter energy divided by $y^2 H_0$ | Pair count change | Energy change |
|---|---:|---:|---:|---:|
| Quiet | 0.494140 | 10.52853 | — | — |
| Nearby positive | 0.493807 | 10.52145 | −0.06742% | −0.06723% |
| Nearby negative | 0.496464 | 10.59277 | +0.47033% | +0.61017% |
| Distant positive | 0.494445 | 10.53966 | +0.06164% | +0.10570% |

Each pair is counted once. The numbers are the coefficients of $y^2$ per added parent; quantitative use requires a small enough $y$ for higher orders to be negligible. The quiet count contains 0.33919 from annihilation and 0.15495 from excess stimulated creation. The latter is about 31% of the inclusive count.

The interaction switch is $\sin^4[\pi(s-s_l)/(s_u-s_l)]$ on its support and zero outside. Shortening the proper-time window from $1.15$–$2.05$ to $1.25$–$1.95$ changes the quiet coefficient from 0.49414 to 0.40835. For the shorter switch the changes in pair count are −0.07120%, +0.47150% and +0.05721%, respectively. At the common readout $s=2.10$ the ratios compare energies measured by the comoving shell observer. For the wider quiet run a screen at total comoving pair energy $\omega=60$ excludes 0.49% of the full count and 1.20% of the final energy. Each retained daughter has energy below the chosen gravitational scale throughout the switch. These are measured sensitivities to the switch and frequency range. [Spectral calculation and resolution tests](primary/TAGGED.md) · [Independent integration](independent/CHECK.md)

## Interaction and observable

Consider an explicit Yukawa interaction $-y f\,\chi\bar\psi\psi$ on the **prescribed** B1 shell. Start with one extra parent scalar in the zero-comoving-momentum mode of a normalization box and subtract a run with the same geometry, initial mode functions, coupling switch and daughter vacuum, but with no extra parent. The observable is the excess number and energy of Dirac particle–antiparticle pairs after the coupling returns to zero. This differences out vacuum production common to both runs. A single tagged rest parent is a test particle; the result is per parent and supplies no cosmological energy density.

With $v=a\chi_{\bf 0}$ and $\Psi=a^{3/2}\psi$, the conformal interaction is $-y f v\bar\Psi\Psi$ for the rest mode with its usual box normalization. Massless Dirac daughters are conformal plane waves on the spatially flat homogeneous shell. The unit-Wronskian parent obeys $v''+(a^2m^2-a''/a)v=0$ and its initial data are those of the physical oscillator ground. The extra particle occupies this same mode. Define $F(\omega)=\int d\eta\,f(s(\eta))v(\eta)e^{i\omega\eta}$. At leading order in $y^2$ the spectrum for total conformal energy $\omega=2|\mathbf p|>0$ is

$$
\frac{1}{y^2}\frac{d\Delta N_{\rm pair}}{d\omega}
=\frac{\omega^2}{8\pi^2}\left(|F(\omega)|^2+|F(-\omega)|^2\right).
$$

The terms describe parent annihilation and excess stimulated creation. In the latter channel, creating another scalar in the occupied mode has matrix-element squared 2, compared with 1 in the reference vacuum. The net weight after subtraction is 1. Creating scalars in other modes cancels in the subtraction. The spin sum is $4p\cdot q=2\omega^2$ in a mostly-minus convention with opposite massless momenta. The two daughters are distinguishable. Integrating and applying Parseval gives two independent ways to compute the inclusive number. The observed pair energy redshifts by $a_f^{-1}$ at the common final cut.

$$
\frac{\Delta N_{\rm pair}}{y^2}=\frac{1}{4\pi}\int |(fv)'|^2\,d\eta,\qquad
\frac{\Delta E_{\rm pair}}{y^2}=\frac{1}{8\pi^2 a_f}\int_{\mathbb R}|\omega|^3|F(\omega)|^2\,d\omega .
$$

The specified switch contributes to the energy budget. Even in flat space with a constant mass $m$ and an initially resting parent, one obtains $\Delta N_{\rm pair}=\Gamma_0\int[f^2+\dot f^2/m^2]\,dt$, with $\Gamma_0=y^2m/(8\pi)$. Daughter energy counts both channels positively. Parent depletion counts annihilation minus excess creation. With both endpoint interactions off, external work on the excess state is $\Delta E_{\rm pair}-m(P_{\rm ann}-P_{\rm stim})$. In flat space for a real smooth switch this equals $2\Gamma_0\int\dot f^2dt/m+2\int(\omega+m)\,dP_{\rm stim}$. An external control can supply energy even in a model with explicit daughter degrees of freedom. The flat result provides a calibration. The cosmological histories have additional mass and gravitational work.


## Scope and reproduction

The positive daughter energies arise from an explicit tree-level interaction between fields. The parent is a single prepared rest mode on stored wall histories. The full distribution of produced parents, nonzero parent momenta, interactions among daughters, and feedback of this new source on the geometry are open. The finite coupling windows also carry work. A cosmological radiation density or thermal history would require those additional calculations. [Independent flat-space calibration](flat/NOTES.md) · [Microscopic precedent for scalar decay on an expanding background](https://arxiv.org/abs/1904.12343v2)

From this directory in the HDblast repository, run the replay with NumPy and SciPy installed. The programs authenticate the original shell archive and the quiet crossing normalization. The workflow runs on Python 3.12 with NumPy 2.3.5 and SciPy 1.17.0.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B reproduce.py
~~~
