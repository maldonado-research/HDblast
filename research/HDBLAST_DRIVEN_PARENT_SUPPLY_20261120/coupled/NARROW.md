# Parent supply near a kinematically narrower decay threshold

A targeted five-run check uses the same initially unexcited wall, positive canonical oscillator and finite dose 0.123 as the principal study. We raise the nominal daughter mass parameter to \(mL=400\), so the initial parent frequency is 760. The compact source starts at \(\psi=0.03\). Its narrower field width of 0.0006 corresponds to a measured duration of 0.004914 in wall time, or 3.735 radians at the initial oscillator frequency. The same former fluid rest-energy budget now requires an action of \(0.0008/760=1.05263\times10^{-6}\). The target particle number changes with mass.

| Drive amplitude | Field width | Early × late elements | Ratio at source off | Ratio at opening |
|---|---|---:|---:|---:|
| 0.000015 | 0.0006 | 18 × 3 | 1.41325 | 1.40355 |
| 0.000020 | 0.0006 | 18 × 3 | 2.52150 | 2.57803 |
| 0.000020 | 0.0006 | 48 × 6 | 2.52177 | 2.57769 |
| 0.000030 | 0.0006 | 18 × 3 | 5.72830 | No opening before radius 120 |
| 0.000020 | 0.15 | 18 × 3 | Below resolved supply | Below resolved supply |

The broad-source run returned a residual on the order of \(10^{-23}\) times target. This small value has no established relative numerical accuracy. Slowly varying the minimum for many oscillation periods excites little residual motion. Stronger driving in the narrow field interval changes the self-consistent geometry and can prevent the later gate from opening.

For the refined \(g=0.000020\) run the kinematic opening occurs at wall time 6.22963 and radius 10.20504. It has \(I=2.71336\times10^{-6}\), \(|H|/\omega=0.00044754\), and \(|\dot\omega|/\omega^2=0.0290854\) at opening. Comoving drive work, mass work and expansion work to that event are 0.00201804, 0.000151679 and −0.00000099918. No daughter transfer or dissipative dynamics is included.

For interpretation only, insert the earlier two-body effective width coefficient \(\kappa=2000\) into
\[
 \Gamma/ M=\frac{\kappa}{2mL}(1-R^{-2})_+^{3/2}.
\]
It rises from zero at \(R=1\) to 0.0002231 at \(R=1.001\), 0.006915 at \(R=1.01\), 0.02644 at \(R=1.025\) and **0.07087** at \(R=1.05\). Locally the measured \(\dot R_*=23.2683\) gives \(\Gamma L\sim\kappa[2\dot R_*(\tau-\tau_*)]_+^{3/2}\) immediately above threshold if the trajectory continues. The absolute width and its first derivative vanish at the opening, while its fractional growth rate increases without bound as the opening is approached from above. The simulation stops at the opening.

This reduction in the kinematic width ratio follows from the mass and available phase space. Matching to a single-channel vacuum Yukawa expression requires \(y^2/(4\pi)=\kappa/(mL)=5\). Near-threshold phase-space suppression does not establish a perturbative coupling. A larger microscopic theory could have a sum of decay channels. Their masses, individual couplings and dynamics would have to be provided.

Reproduce all five cases with the production solver and compare the outputs with the archived `NARROW.json`.

```sh
python production.py --cases narrow_cases.json --output narrow_replay.json
```
