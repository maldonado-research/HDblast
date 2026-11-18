# Can expansion look radiation-dominated when little ordinary radiation is present?

The same cutoff run used for the scalar-work calculation has a later reliable interval. Use the dimensionless proper time $s=H_0\tau$ and the saved expansion rate $h=H/H_0$. The registered run has $\phi_*=0.5$, $G=100$, $y=1$ and detuning 0.1. We chose the later window $4\le s\le6$ for this retrospective comparison. Both original reliability intervals end after 7.43. The window spans only 0.165 e-folds. Its duration matters for interpreting a good fit. [Earlier energy-transfer calculation](https://github.com/maldonado-research/HDblast/tree/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_REHEATING_TRANSFER_20261118)

## Shape versus the contribution to expansion

A flat, expanding four-dimensional radiation universe in general relativity with constant gravitational coupling has $H\propto\sqrt{\rho_r}$. If radiation is separately conserved it also has $H\propto a^{-2}$. In the simulation $R$ is proportional to the material radiation density, with a fixed conversion factor within each run. To allow a free normalization in the first comparison, set $q_i=h_i/\sqrt{R_i}$. The largest fractional discrepancy at positive calibration $C$ is $\max_i |C/q_i-1|$. The two extreme values determine its minimum. Equating the positive and negative discrepancies yields
$$
 C_{\rm opt}=\frac{2q_{\min}q_{\max}}{q_{\min}+q_{\max}},\qquad
 e_{\rm min}=\frac{q_{\max}-q_{\min}}{q_{\max}+q_{\min}} .
$$
For the second comparison use $q_i=h_i a_i^2$. These are separate one-parameter fits evaluated at 48 nodes on each radial grid, including linearly interpolated endpoint columns.

| Diagnostic on $4\le s\le6$ | Coarse | Fine |
| --- | ---: | ---: |
| Smallest maximum fractional error fitting $h=C\sqrt R$ | 1.1560% | 1.1565% |
| Smallest maximum fractional error fitting $h=C a^{-2}$ | 0.8331% | 0.8335% |
| $[h(6)^{-1}-h(4)^{-1}]/2$ | 2.1031 | 2.103 (rounded) |

The last row requires no free normalization. For flat, separately conserved pure radiation, $d(1/H)/d\tau=2$. During the modeled window $a^4R$ still increases by about 1.30%. These residuals measure departures of the numerical histories from the comparison laws. They are not observational uncertainties or rigorous time-continuum bounds.

On the fine grid ordinary radiation contributes **1.075%–1.126%** to the code's decomposition of $H^2$. The Weyl term contributes about 102.4%–105.9%, while the vacuum term contributes about −7.0% to −3.5%. The components are signed; a term can exceed the total. The two fitted shapes can therefore look close to radiation expansion over this window even while the modeled material radiation supplies a small part of the expansion. A normalization inferred by fitting the expansion itself does not independently calibrate the gravitational response to ordinary radiation. [Source, grid and exact sample definitions](INPUTS.json) · [Numerical replay](RESULTS.json)

## A conditional energy scale

In natural units the B1 convention is $b=\kappa_5^2H_0^3$ and the stored background radius is $r_b=1/H_{0,\rm code}$. The physical equivalent of the model radiation density obeys $\rho_r/H_{0,\rm phys}^4=r_b R/b$. Define a reference temperature by placing that energy density in a single equilibrium bath of massless species. For illustration take $g_*=2+(7/8)(4+6)=10.75$ for photons, electrons, positrons and three neutrino species with the same temperature and negligible chemical potentials. The associated scale is
$$
 \frac{T_{\rm eq}}{H_{0,\rm phys}}
 =\left(\frac{30 r_b R}{\pi^2 g_* b}\right)^{1/4}.
$$
Among the saved nodes on $1\le s\le4$, the largest value on either grid rounds to **2.2043**. At $s=4$ it is **1.9884**. The physical value of $H_0$ in an energy unit, the identity of the decay products, and their interaction and equilibration rates must be supplied to obtain a plasma temperature in MeV. A radiation density at a selected time and a temperature defined via a decay rate crossing the expansion rate are different quantities. This comparison assigns no BBN limit to HDblast. [B1 registration and units](https://github.com/maldonado-research/HDblast/blob/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/REGISTRATION.md) · [Particle Data Group, Big-Bang Cosmology, relativistic energy densities](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-bbang-cosmology.pdf)

Run `python -I -B replay.py` in this directory. The CSVs contain the needed clock, material and geometric columns. The optional NumPy extractor accepts the original archived B1 directory and a grid name. It authenticates the source files against the recorded digests.
