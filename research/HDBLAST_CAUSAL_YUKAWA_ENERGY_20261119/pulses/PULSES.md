# What the incoming pulses change

We applied the coupling to the same retained fork state on four archived wall histories. All arms use the same smooth proper-time switch from 1.35 to 2.05 and the same readout at 2.10. The quiet history carries the rolling background. The three other histories inherit the public simulation labels plus, minus and far. Each number below is the pulsed result minus the quiet result. The units are coefficients of $y^2H_0^4$ per comoving volume. Positive daughter energy in this table means that the pulse increased the energy emitted relative to quiet. [Protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-6998631142) · [Native histories](https://doi.org/10.5281/zenodo.24357388)

| Fine wall, comoving cutoff 120 | Daughter energy | Relative to quiet | Parent energy, bare | Parent energy, static | Switch work, bare | Switch work, static |
|---|---:|---:|---:|---:|---:|---:|
| Plus | −0.03315 | −0.1473% | +0.08559 | +0.07627 | −0.05675 | −0.00474 |
| Minus | +0.22716 | +1.0091% | −0.30018 | −0.42173 | +0.28031 | +0.02273 |
| Far | +0.02917 | +0.1296% | +0.01141 | −0.03517 | +0.01507 | +0.00106 |

The signs of the daughter-energy differences persist on both wall grids at cutoffs 60 and 120. A finer grid changes the three cutoff-120 differences from −0.03375, +0.22835 and +0.02875 to the tabulated values. The far pulse gives more outgoing energy than quiet in both calibrations, although the inferred change in the parent energy reverses sign between them. This is a specific counterexample to inferring the parent-energy response from the pair signal on these walls. The pulse arms and quiet history have the same prescribed switch. A causal inference that the five-dimensional pulse supplies a particular fraction of the emitted energy would require a self-consistent dynamical source for that switch and a matched stress tensor. [Numerical output](PULSES.json)

## Independent check of switching work

The primary calculation integrates the mixed parent–bath covariances. Here a second solver propagates only the uncoupled positive-frequency parent mode $u=\chi$ from the fork. It obtains the mixed covariance entering the switch work from two causal integrals in conformal time. Define $F=af$, $R(\tau)=2\int_q^\Lambda K_q(\omega)\sin(\omega\tau)d\omega$ and $D(\tau)=\int_q^\Lambda K_q(\omega)\cos(\omega\tau)d\omega$. Extend the drive by zero before the fork and write $R_u=R*(Fu)$ and $D_u=D*(Fu)$ with the integration restricted to earlier times. For occupation $n$ and anomalous covariance $c=\langle bb\rangle$, the bare mixed expectation divided by $y$ is

$$
2\,\mathrm{Im}(u D_u^*)-(2n+1)\,\mathrm{Re}(u R_u^*)-2\,\mathrm{Re}(c u R_u).
$$

Integrating this expression against $df/ds$ in proper time gives the bare switch work. In the static prescription add $\int ds\,a f(df/ds)\langle\chi^2\rangle_0\,2\int_q^\Lambda K_q(\omega)d\omega/\omega$. The frequency integration uses 96 Gauss nodes per parent momentum. The causal time integrals use trapezoidal integration with Richardson extrapolation from 512 and 1,024 intervals. Agreement with the covariance solver for all paired switch-work differences is within $2.1\times10^{-7}$ in the table’s energy units. This includes differences below 0.002. [Both calculations](pulses.py)
