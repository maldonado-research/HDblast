# How much light would direct particle–antiparticle annihilation produce?

A charged particle can annihilate with its antiparticle into two real photons. The two photons share the full incoming energy, including both rest masses. We applied this channel to the archived dedicated-vector shell simulation, assuming unpolarized spin-½ Dirac daughters coupled to their own local four-dimensional electromagnetic field. The calculation includes the history of both charge populations from birth. It uses the electromagnetic coupling fixed by the primary model.

At wall time 7.05, the finest run yielded **$9.94\times10^{-29}$ of the charged gas's final conformal energy** in these photons. The count was **$1.90\times10^{-26}$ photons per daughter** born by that time. Uniformly extrapolating the simulation to the complete round shell with $6.69\times10^{39}$ daughters gives $1.27\times10^{14}$ emitted photons. This large absolute count represents a very small fraction of the available particles and energy. The photon energy is about $3.00\times10^{-24}$ of the final charged rest-energy budget. [All eight realizations and quadrature checks](ANNIHILATION.json)

| Wall time | Elastic collisions in inherited trajectory | Positive particles | Time steps | Photons per final daughter | Photon energy / final gas energy |
|---|---|---:|---:|---:|---:|
| 6.45 | Active | 4,096 | 1,024 | $4.76\times10^{-31}$ | $6.58\times10^{-32}$ |
| 6.67 | Active | 16,384 | 4,096 | $8.56\times10^{-30}$ | $6.82\times10^{-31}$ |
| 7.05 | Active | 4,096 | 4,096 | $1.46\times10^{-26}$ | $8.36\times10^{-29}$ |
| 7.05 | Active | 8,192 | 4,096 | $2.42\times10^{-26}$ | $1.16\times10^{-28}$ |
| 7.05 | Active | 8,192 | 8,192 | $2.01\times10^{-26}$ | $1.02\times10^{-28}$ |
| 7.05 | Active | 16,384 | 8,192 | $1.90\times10^{-26}$ | $9.94\times10^{-29}$ |
| 7.05 | Off | 8,192 | 4,096 | $8.37\times10^{-29}$ | $6.71\times10^{-30}$ |

Each sample also has the same number of negative particles. The second-to-last column counts both photons from each event. The final column includes their production-time expansion factor, so emitted energy at different wall times can be compared on a common final slice if the photons propagate freely on the shell. The four interacting late samples span a factor of 1.39 in energy fraction. Their spread includes the sensitivity of the underlying charged-particle histories.

## Cross section and upper bound

Put $m=c=\hbar=1$ within the two-body formulas. Let $p,q$ be wall-frame three-momenta, $E_1=\sqrt{1+|p|^2}$, $E_2=\sqrt{1+|q|^2}$, $g=E_1E_2-p\cdot q$ and $s=2(1+g)$. The invariant Møller flux is $v_{\rm M}=\sqrt{g^2-1}/(E_1E_2)$. Multiplying the standard tree-level two-photon annihilation cross section by the full flux gives the finite rate kernel

$$
K=m^2\sigma_{2\gamma}v_{\rm M}
=\frac{\pi\alpha^2}{E_1 E_2(g+1)}
\left[(g^2+4g+1)\frac{\operatorname{acosh}g}{\sqrt{g^2-1}}-(g+3)\right].
$$

The ratio at $g=1$ has continuous limit one. The code calculates $g-1$ from positive relative-rapidity and squared direction-difference contributions. Its controls include the co-moving threshold for different boosts and the center-of-mass expression from detailed balance with the inverse Breit–Wheeler process. The tree-level diagrams exchange the charged Dirac particle. The separate neutral vector resonance sets the elastic transport history. An on-shell spin-one resonance cannot decay to two identical real photons. [Heitler cross section and reference, Geant4 Physics Reference Manual](https://geant4.web.cern.ch/documentation/dev/prm_html/PhysicsReferenceManual/electromagnetic/electron_incident/annil.html)

Writing $g=\cosh y$ gives $\operatorname{acosh}(g)/\sqrt{g^2-1}=y/\sinh y\le1$. Also $E_1 E_2\ge(g+1)/2$. Together,

$$
0\le\frac{K}{\pi\alpha^2}\le
2\frac{g^2+3g-2}{(g+1)^2}=2+\frac{2(g-3)}{(g+1)^2}\le\frac{17}{8}.
$$

The last maximum occurs at $g=7$. This ceiling applies to **every incoming pair**, including slow pairs unresolved by sampling. The lower sign follows from the positive total cross section.

## Two charges, three flavors and the energy ledger

The plateau coupling is $\alpha=3Q^2\mu^2/(4\pi K_{\rm wall} S_0)=1.44382\times10^{-16}$, with $Q=10^7$, $\mu=mL=4\times10^{12}$, $K_{\rm wall}=2.54$ and $S_0=1.0415592372880226\times10^{54}$. Each positively charged particle encounters the negative particles of its own species. Define $C=S_0/(6\mu^3)$, the scale factor $a(t)$ and the accepted total comoving rest budget $D(t)$, summed over both signs and all flavors. The per-daughter exposure per unit model time from the target density and dimensionless kernel is $CD(t)\langle K\rangle/a^3$. A charge pair has two participants and makes two photons. With $D_f$ the final rest budget, the expected count relative to all final daughters and photon conformal energy are

$$
\frac{N_\gamma}{N_{\chi,f}}=\frac{C}{D_f}\int \frac{D^2}{a^3}\langle K\rangle dt,\qquad
E_{\gamma,\mathrm{conf}}=\frac{C}{2}\int \frac{D^2}{a^3}\left\langle K\,a(E_1+E_2)\right\rangle dt.
$$

The full unit round three-sphere has physical volume $2\pi^2a^3 L^3$ and $N_{\chi,f}=2\pi^2S_0D_f/\mu$, including all three species and both charges. The final charged rest ledger in conformal units is $a_fD_f$. [Inherited shell density, coupling conventions and pulse](../../HDBLAST_AXIAL_CURRENT_20261121/live/README.md)

Using $D\le D_f$ and the smallest radius on the recorded interval bounds the integrated annihilation probability per daughter by $(17\pi\alpha^2/8) C D_f \Delta t/a_{\min}^3$. At 7.05 this is below $3.60\times10^{-22}$ for any daughter. A sharper count ceiling uses the number born by each time step. An energy ceiling additionally sums, for every tag, the larger conformal energy at the two ends of the particle push, both evaluated at the larger step radius. Convexity of $\sqrt{a^2+|u|^2}$ bounds the interpolated comoving momentum $u$ through the push. The resulting late ceilings are below $9.24\times10^{-23}$ photons per final daughter and $1.28\times10^{-22}$ of the final gas energy across the four interacting grids. The energy ceiling is conditional on these particle histories and time discretizations. It does not require the sampled pair weights to find the largest kernel.

## Integration and reproducibility

At every live time step the estimator samples 256 pairs in each of four independently seeded batches. The positive charge's cumulative proposal probability is stratified, and the negative charge is independently sampled. Each proposal uses a half-uniform, half-inverse-energy mixture. Dividing by its exact discrete probability restores the equal-weight pair average. Each pair contributes only during the overlap of its two lifetimes. An independent uniform substep samples this interval, with comoving momenta interpolated through the force push before the split elastic collision. The estimator uses its own random generator and leaves the elastic collision stream and fields untouched. The replay verifies the archived elastic encounter count and final gas energy. In the finest late case, the standard deviations across the four quadrature values are 0.27% of their mean for photon counts and 0.19% for energy. These are numerical diagnostic spreads, not confidence intervals for the physical model. [Estimator and bound](pairs.hpp) · [Replay](integrate.py)

With the primary studies in `research/` and the extended dedicated history generated through 7.0505, run the following from the repository root.

```sh
python -I -B research/HDBLAST_AXIAL_CURRENT_20261121/fields/extract.py \
  --study research/HDBLAST_POLARIZATION_20261121 --cases dedicated_half --end 7.0505 \
  --output-dir /tmp/hdblast-annihilation --output /tmp/hdblast-annihilation/FIELDS.json
OPENBLAS_NUM_THREADS=1 python -I -B research/HDBLAST_COLLISION_PHOTONS_20261121/annihilation/integrate.py \
  --study research/HDBLAST_AXIAL_CURRENT_20261121 --polarization research/HDBLAST_POLARIZATION_20261121 \
  --scattering research/HDBLAST_SCATTERING_20261121 --history /tmp/hdblast-annihilation/dedicated_half.npz \
  --output /tmp/hdblast-annihilation/ANNIHILATION.json --jobs 2
```

The calculation is a first-order emission diagnostic on an already populated shell. It assumes unpolarized free Dirac particles and vacuum two-photon kinematics. Bound states, Sommerfeld corrections for sufficiently slow relative motion, spin correlations, in-medium masses, inverse reactions, reabsorption, five-dimensional photon transport and further particle species each change the problem. This process has no infrared energy cutoff and conserves four-momentum at each annihilation. A photon temperature or hot-universe origin would require additional physics.
