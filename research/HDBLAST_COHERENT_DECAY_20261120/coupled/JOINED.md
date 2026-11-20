# From an initially unexcited wall to charged-particle heating

A single finite-source evolution now follows the production of a neutral field, loss of its energy to oppositely charged particles, and electrical heating of those particles. The initially empty particle sector reaches 90% of the added wall energy in intrinsic random motion at radius **12.26756**. It remains above that fraction until the expanding wall reaches radius **35.45899**, a growth factor of **2.89047**. At exactly double the starting hot radius, the internal share is **97.2877%**. This is a conditional result for a prescribed effective model. The microscopic decay law, charge coupling, collision rate and spectrum still require an independent physical derivation.

## What is being joined

The closed wall begins at radius 2, with zero expansion rate, bulk mass 4, vanishing scalar and electromagnetic fields and zero particle density. The two finite proper-time pulses propagate through dynamical, finite-depth bulk collars. The incident dose parameter is 0.123. A pulse excites the wall scalar \(\psi\). A compactly supported shift \(v(\psi)\) of a second field’s oscillator minimum excites the initially vanishing parent \(\chi\). Its stored potential is \(\tfrac12\omega(\psi)^2(\chi-v)^2\) and the wall scalar receives the opposite reaction. This uses the source of the previous production study and the parity-even gauge-energy transfer of the earlier heating calculation. The full wall/bulk initial-value calculation is run again for every case. [Parent production](https://github.com/maldonado-research/HDblast/pull/191) · [Earlier heating model](https://github.com/maldonado-research/HDblast/pull/190)

The primary uses \(\mu=mL=400\), source amplitude \(2\times10^{-5}\), and field support \([0.03,0.0306]\). The mass \(\omega=ML=2\mu R(\psi)\) moves from \(R=0.95\) to 1.05 across an interval of width 0.003 around \(\psi=0.53\). The end of the parent drive precedes this interval. All frequencies use wall proper time. A local damping closure converts energy into daughters when their kinematic threshold opens,
\[
 \ddot\chi+(3H+\Gamma)\dot\chi+\omega^2(\chi-v)=0,
 \qquad \Gamma=\kappa R(1-R^{-2})_+^{3/2},
 \qquad P=\Gamma\dot\chi^2.
\]
The dimensionless quantity \(P\) is the local power removed from the oscillator. The primary has \(\kappa=40\), \(Q=q/(\sqrt{3}m)=2\times10^8\) in project units, and drag coefficient \(\nu L=18\). The maximum prescribed width-to-parent-mass ratio is 0.0014174. Matching the rate scale to a single Yukawa channel would imply \(y^2/(4\pi)=0.1\). The large charge-to-mass coefficient and the parent-field dependence of a microscopic daughter mass still matter. These dimensionless values alone do not establish a consistent microscopic model.

## Rest mass, birth energy and gravity

The three isotropic pairs of opposite charge streams start with zero rest mass, momentum and internal energy. Let \(D\) denote their total comoving rest energy, \(S\) their comoving directed momentum summed with a common sign convention for the two charge classes, and \(E=a^3\rho_g\) their laboratory-frame comoving energy. We evolve \(U=E-\sqrt{D^2+S^2}\) to resolve initially cold states. An isotropic decay at rest supplies opposite charges with zero mean directed momentum. For the comoving pair count \(N\) and total birth energy \(W_b\), the source is
\[
 \dot W_b=a^3P,\qquad \dot D=\dot W_b/R,\qquad
 \dot N=\dot W_b/\omega,\qquad D=2\mu N.
\]
The remaining \((1-1/R)\dot W_b\) is kinetic energy at birth. The comoving rest variable grows as decays occur. There is no assigned initial inventory.

The Taub–Mathews equation of state recovers each stream’s proper pressure and internal energy from \(D,S,U\). With \(v_g\) its directed speed, \(J=QD v_g/a^3\) is the electromagnetic current. The equations include \(\dot S=QDX-(H+\nu)S\), an opposing current in the Maxwell equation, and the full gravitational pressure of kinetic motion, electromagnetic fields, scalar fields and gas. The symmetric drag transfers directed momentum into random motion while conserving total gas energy. For \(W_X=\int a^3JX\,d\tau\), the two energy ledgers are
\[
 a^3\rho_\chi +\int 3Ha^3p_\chi\,d\tau + W_b=W_v+W_M,\qquad
 E+\int 3Ha^3p_g\,d\tau=W_X+W_b .
\]
Here \(W_v,W_M\) are the moving-minimum and changing-mass work debited from \(\psi\). A separate surface field ledger uses \(a^4\) weighting. The coherent fields, the gas and both gravitational collars evolve simultaneously. The auxiliary independent radiation fluid stays zero.

## Why changing the final gauge coefficient matters

In the previous narrow-source production case, the final gauge coefficient was 26.04. With all sectors active this choice gives only about 28% intrinsic particle energy at the reference decay scale. At the start of the final scalar release, the measured incoming electric displacement suggests a matched coefficient \((KX)^2/(2\Delta V)=2.4863\) for the stored scalar potential \(\Delta V=0.061\). We chose a final value of **2.54**. The first gauge step and the reciprocal interpolation during the last release have the same form as in the earlier models. This parameter was chosen using the actual incident-wall trajectory. It is an engineered match between two parts of the model.

The primary bounces at radius 1.88140. The parent’s instantaneous action is 2.522 times the previous supply target when the source switches off. The final scalar drop completes at time 6.23797 with the surface field 99.99984% electric and \(aH=3.6460\). Decay opens at time 6.24180 and radius 10.24955. The continuously hot interval starts at time 6.87263 and ends at the first expansion maximum, time 39.75226. At the doubling time 14.35652 the gas pressure-to-density ratio is 0.333333331. The intrinsic random share peaks at **98.6728%**. The separate effective hot-fluid pressure closure does not establish a measured temperature or microscopic equilibration.

At the intrinsic peak, the electrical field has supplied 56.5103 in the same comoving proper-volume units as the direct parent transfer of 0.00223092, a factor of 25,331. There are \(2.65598\times10^{-6}\) pairs per dimensionless coordinate volume by the end of the run. The summed rest energy is 0.00212478, with birth kinetic energy 0.000106133. The initial burst and the stored wall scalar provide energy. Most of the daughter energy after birth comes through the electric field. Pressure work during expansion subsequently reduces the comoving gas energy.

At 48 initial radial elements and six late elements the absolute residuals in the Friedmann constraint, Maxwell ledger, daughter ledger and parent ledger are respectively \(1.28\times10^{-10}\), \(6.42\times10^{-9}\), \(1.91\times10^{-9}\), and \(4.60\times10^{-10}\). The particle-rest and pair-number check is within \(9.5\times10^{-16}\) in comoving rest-energy units. The minimum sampled bulk mass is 0.10013 and lapse is 0.99778.

## Resolution and comparison runs

| Case | Peak intrinsic share | Expanding e-folds continuously above 90% |
|---|---:|---:|
| Matched primary, 48 × 6, \(\kappa=40\) | 98.673% | 1.06142 |
| 24 × 6 with primary parameters | 98.612% | 1.03857 |
| 64 × 6 with primary parameters | 98.678% | 1.06349 |
| 96 × 9 with primary parameters | 98.68% | 1.06412 |
| 48 × 6 with slower decay, \(\kappa=4\) | 98.461% | 0.99967 |
| \(\kappa=0.4\) | 96.570% | 0.67913 |
| \(\kappa=0.04\) | 82.160% | 0 |
| No decay, or no incident pulse, or no parent drive | 0% | 0 |
| Uncharged daughters | 0.0021% | 0 |
| Zero explicit drag | 62.192% | 0 |
| Final gauge coefficient 2.34 or 2.84 | 99.362% or 97.176% | 1.21807 or 0.67296 |
| Original narrow-source gauge coefficient 26.04 | 28.033% | 0 |
| 26.04 with drive amplitude reduced to \(1.3\times10^{-5}\) | 98.522% | 0.91818 |
| Historical \(mL=1\) broad drive, \(\kappa=2000\), final coefficient 26.04 | 96.589% | 0.57805 |

Unless specified, parameter-control runs use 18 × 3 elements. Changing dose to 0.12295 lowers the peak to 22.02%. At 0.12305 the gate does not open before radius 120. The width-scale comparison on the old gauge trajectory spans \(\kappa=4,40,100,200,2000\) and peaks between 27.85% and 28.09%. Increasing the number of carriers alone cannot fix an inadequate supply of electric energy. In the zero-drag run the continuous single-fluid closure mixes birth cohorts of different velocities, so some modeled internal energy can remain without explicit drag. The historical light-mass case has \(\Gamma/M\simeq28.35\) in its mass plateau and is a nonperturbative comparison.

These controls map a selected range of model behavior. The positive trajectory has no microscopic derivation of the moving minimum, gauge profile, width or drag, and assumes a fixed daughter mass even as \(\chi\) oscillates. A literal Yukawa interaction, radiative corrections, width memory, Pauli blocking, and a collision model require additional tests. A finite hot interval in a closed wall that eventually recollapses does not identify the cause of the observed Big Bang.

## Reproduction

```sh
python joined.py --jobs 3 --output replay.json
python verify.py replay.json
```

Use Python 3.9 or later and a C++17 compiler. The first command starts each case from the initial time and records its solver SHA-256. The verifier reconstructs the proper internal energy, total pressure, particle counts, current speed limit, all energy transfers and the doubling interval from the stored states. `JOINED.json` supplies the complete parameter sets, three sampled trajectories and summaries for all 35 cases. High-frequency parent oscillations can be under-resolved in plotted samples. The production power and all accumulated source, current and pressure work use the native integration stages. Scalar events and hot thresholds use interpolation within native steps.
