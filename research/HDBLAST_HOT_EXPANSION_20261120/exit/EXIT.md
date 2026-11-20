# Sustained particle heating during expansion in a two-stage wall model

A finite incoming pulse drives an initially cold, closed wall through contraction, a bounce and a scalar-potential plateau. In the selected calculation the radius grows from 10.75233 to 35.90603 while at least 90% of its added energy stays in the particles’ internal motion. This is a factor of **3.33937** in radius, or **1.205782 e-folds**. At the first doubling, radius 21.50467 and wall time 12.20762, the internal share is **98.7258%** and the particle pressure is one third of its total energy density to nine decimal places. At that doubling the total added wall pressure divided by added energy is 0.34181. The bulk mass term contributes 0.5655% of the positive gravitational source \(2e+e^2+\mu/a^4\) in the Friedmann equation. The wall subsequently turns around. The earlier one-stage carrier-decay case gave 0.0077 e-folds at the same threshold. Both numbers describe dimensionless trajectories of a specified model. [Previous carrier-decay study](https://github.com/maldonado-research/HDblast/pull/188)

## Where the particle energy comes from

The incoming higher-dimensional electromagnetic pulse initially raises the wall scalar onto a potential of height 0.111. A magnetic field produced by the pulse remains on the sphere. At the first downward step the scalar gives up 0.050 of its potential and the parity-even gauge coefficient rises from 0.04 to 3400.04. Its subsequent evolution generates electric displacement before the remaining 0.061 is released. At the second step the gauge coefficient drops to 26.04. For approximately frozen displacement, a decreasing gauge coefficient raises electric energy. Newly produced charges draw this energy from the electric field. An imposed inter-stream drag converts directed motion to isotropic random motion. All energy transfers enter the wall and bulk equations, including the work from changing the parent mass. The separate radiation-fluid variable stays zero.

For positive gauge coefficient \(K\), surface electric field \(X\), and magnetic amplitude \(B\), write the surface field energy as
\[
\rho_{\rm EM}=\tfrac12 K B^2+\frac{(KX)^2}{2K}.
\]
The direct change due to the gauge coefficient is \(\dot K(B^2-X^2)/2\), with opposite work on the scalar. The measured displacement at the start of the second step is \(KX\simeq-1.80\). Matching its electric energy to the remaining potential gives the design estimate \(K_{\rm end}\sim(KX)^2/(2\times0.061)=26.59\). Expansion, radiation into the bulk, magnetic induction and changing scalar speed remain active in the integration. At the end of the step the electric field carries 99.99805% of the coherent surface field energy and \(aH=3.48861\). At the first hot state \(aH=3.31038\). There is sufficient expansion left to double the radius.

The shape of the gauge coefficient and its two release scales were chosen to achieve this ordering. They have no independently established microscopic derivation. The positive result is a numerical existence example within the effective equations. The effect requires an initial wall with a cold population of unstable parents, a particular scalar potential, and a prescribed collision law. The calculation does not predict the continuous radiation era of observed cosmology, an observed particle spectrum, primordial perturbation amplitudes, or the origin of the initial wall.

## Finite initial data and reproducible model

At time zero \(a=2\), \(H=0\), the scalar and all gauge fields vanish and the parent rest density is \(D_0=10^{-4}\). The initial bulk mass 3.99679984 solves the wall constraint. The sources are the same two proper-time limited inner-tube pulses as in the previous study, with dose parameter 0.123 and cutoff at dimensionless source time 0.0039783904. Both Gaussian-normal collars have depth 0.002. The remesh at wall time 0.08 occurs after the pulse. Its outgoing boundary remains open throughout the run. The bulk mass and wall geometry continue to evolve. No external energy is injected after the source cutoff.

With \(S(x)\) the same compact smooth step as in the earlier solver, define \(r=S((\psi-.004)/.008)\), \(s=S((\psi-.36)/.01)\), \(z=S((\psi-.5265)/.003)\). The potential is
\[
 V(\psi)=0.111r-0.050s-0.061z .
\]
Before the final step \(K_{\rm pre}=K_{\rm launch}(\psi)+3400s\), where the launch coupling of the original model has reached 0.04 by the preparation interval. During the decreasing step the interpolation is
\[
 K^{-1}=\frac{1-z}{K_{\rm pre}}+\frac{z}{K_{\rm launch}+26}.
\]
Both endpoints are positive, and the solver uses the analytic derivative in the scalar and Maxwell equations. Controls with an increasing final gauge step use the original interpolation linear in \(K\).

The neutral parents’ mass ratio \(R=M/(2m)\) grows smoothly from 0.95 to 1.05 across the interval of width 0.003 centered at \(\psi_*=0.530\). The local two-body decay width is \(\Gamma L=2000R(1-R^{-2})_+^{3/2}\). The two daughters have equal and opposite charges. Three equal stream pairs give isotropic total stress. The dimensionless charge to mass coefficient is \(Q=2\times10^8\) and the phenomenological inter-stream drag is \(\nu L=18\). The Taub–Mathews equation of state supplies a local relationship between pressure and internal energy. For dimensionless temperature variable \(\Theta=P/\rho_0\), the internal fraction uses \(D_{\rm daughters}[\tfrac32\Theta+\sqrt{1+9\Theta^2/4}-1]/e\), where \(D_{\rm daughters}\) is the wall-frame daughter rest-mass density and \(e\) includes all added wall material, scalar and gauge energy. This definition excludes the boost of internal heat, directed stream energy, rest mass and the fixed wall tension from its numerator. The dimensional temperature and microscopic equilibration rate remain undetermined. [Taub–Mathews equation of state](https://arxiv.org/abs/0704.1679)

## Control and resolution checks

The 48-element launch and six-element late primary bounces at radius 1.881040. Preparation begins near wall time 4.0601 and the final drop completes at 6.65085. Decay opens at 6.65713, half of the parents have decayed by 6.67506 and 90% by 6.70210. The continuous hot interval starts at 6.78411 and ends at the radius maximum near time 40.86225. The largest intrinsic share reaches 99.5355%. Electrical work on the gas at its peak is 61.27135 in the same local-comoving-volume units as direct decay injection 0.00088097, or about 69,550 times as much. The sum of daughter energy and pressure work matches these inputs. The maximum absolute residuals of the Friedmann constraint, surface energy ledger, particle energy ledger and gauge energy ledger are respectively 1.28×10⁻¹⁰, 4.55×10⁻¹⁰, 6.55×10⁻¹⁰ and 3.12×10⁻⁸. The minimum sampled local bulk mass is 0.09654 and the minimum lapse is 0.99778. The gauge field’s separate fourth-power-of-radius ledger gives net scalar work 640.839, net bulk work −2.97648 and electric work on the particles 637.862.

| Change at fixed model inputs unless stated | Peak intrinsic share during post-exit expansion | Continuous hot e-folds |
|---|---:|---:|
| Primary, 48 early elements and six late | 99.535% | 1.205782 |
| 24 early elements and six late | 99.568% | 1.205926 |
| 64 early elements and six late | 99.533% | 1.205763 |
| No parent decay | 0% | 0 |
| Uncharged daughters | 0.000080% | 0 |
| No explicit inter-stream drag | 10.63% | 0 |
| Omit the preparation gauge rise | 0.00047% | 0 |
| Keep the high gauge coefficient at final release | 2.64% | 0 |
| Open decay early at \(\psi_*=0.48\) | 0.85% | 0 |
| Open decay late at \(\psi_*=0.58\) | 96.41% | 0.954 |
| Reduce \(Q\) by ten | 99.57% | 1.179 |
| Increase final gauge increment from 26 to 31 | 97.82% | 0.804 |
| Move final release from 0.5265 to 0.534 with matched gate shift | 91.75% | 0.106 |

The parameter controls use 18 early and three late elements. The matched reference at that resolution gives 1.205954 hot e-folds. A zero-incident-dose run contracts from radius 2 to the stopping radius 1.6 without bouncing or producing daughters. At dose 0.122995 and 0.123005 the hot intervals last 1.19037 and 0.95636 e-folds. At 0.122985 the release does not complete before radius 120. At 0.123015 the completed release produces only 0.06361 hot e-folds. With dose 0.12295 or 0.12305 the scalar fails to complete the second release before the expanding wall reaches the reported radius ceiling of 120. Lowering the final coefficient to 21, or moving the final release to 0.516 with the same relative decay onset, also fails to complete the release over the simulated interval. Missing post-exit keys in the stored output signify that the required event was never reached. The sharp dependence on the chosen interaction and incident dose limits inference about typical pulses.

### Can a longer plateau work by itself?

A comparison using the same finite source gave the one-stage magnetic model a longer plateau and retuned its decay threshold for the electric-field maximum. Each of the three runs starts from the cold initial state and uses 18 early and three late elements. For each plateau, a neutral calculation located the first post-exit trace point with an electric fraction of at least 98% on the expanding branch. The decay threshold used that point’s scalar value, rounded to four decimal places, for the charged rerun. All three use \(Q=2\times10^9\), \(\nu L=18\) and \(\kappa=2000\). 

| Final release field, width | Gauge increment | Decay field | \(aH\) at exit completion | Peak hot share | Hot e-folds |
|---|---:|---:|---:|---:|---:|
| 0.366, 0.010 | 42,000 | 0.7443 | 2.845 | 94.34% | 0.0910 |
| 0.3685, 0.005 | 130,000 | 0.7410 | 3.867 | 96.76% | 0.1443 |
| 0.36965, 0.002 | 510,000 | 0.7110 | 5.477 | 94.25% | 0.0651 |

The longer plateaus allow faster expansion at scalar exit. Waiting for the coherent magnetic energy to become electric costs much of that expansion. At the first 90% internal-energy crossing the measured \(aH\) values are only 0.458, 0.586 and 0.374. These selected one-stage controls miss the doubling target. They do not exhaust the space of one-stage potentials or gauge functions. The two-stage branch stores electric displacement before the final positive potential is exhausted.

## Reproduction

Run the following commands in this directory with Python 3.9 or later and a C++17 compiler. The first command replays the primary and all controls from the initial time. The verifier independently reconstructs the proper internal energy and gas energy from the trace, checks the hot interval and doubling, and closes the gas and gauge ledgers.

```sh
python growth.py --jobs 2 --output replay.json
python verify.py replay.json
```

The additional magnetic controls can be reproduced with `python growth.py --cases plateau_cases.json --output plateau_replay.json` and compared to `PLATEAU.json`. For the primary calculation alone, add `--only primary` to the first command. `EXIT.json` contains all launch parameters, trace column names, the primary time series, native-step events, controls and the SHA-256 of the solver. `cases.json` specifies shared defaults and per-case changes. All units are the dimensionless wall units used in the parent study.
