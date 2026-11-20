# Can an initially unexcited wall supply the neutral parent?

In one fully coupled calculation the finite incident pulse leaves enough energy in a coherent neutral field to exceed the parent supply assumed in the earlier carrier model. Starting from zero, its comoving oscillator action reaches **1.514 times the target** when the forcing ends and **1.674 times the target** at the putative decay opening. A stable-parent continuation keeps enough action until the expansion stops. This is an explicit source of classical oscillations. Their conversion into independent particles, their decay, and subsequent heating require further physics.

## Initial conditions and source

The closed, gravitating wall begins at radius \(a_i=2\) with \(H=0\), bulk mass \(\mu_i=4\), and zero added wall energy. The gauge fields and both canonical wall fields have zero initial values and velocities. The background wall tension and the two external, finite-duration incident pulses are part of the initial-boundary model. We use the earlier finite Gaussian-normal collars, their physical time-dependent boundary conditions, the pulse dose 0.123 and the two-stage parity-even wall potential and gauge coefficient. All runs start before the pulses arrive. The geometry and bulk fields respond to the changing wall stress. [Earlier finite-source heating study](https://github.com/maldonado-research/HDblast/pull/190)

Let \(\psi\) be the wall scalar and \(\chi\) the new homogeneous, minimally coupled real field. Its mass frequency is \(\omega=ML=2(mL)R(\psi)\) and \(u=\chi-v\). The extra potential is
\[
 U_\chi=\tfrac12\omega^2 u^2.
\]
The moving minimum is a specified function of the dynamical field \(\psi\). For the primary,
\[
 v(\psi)=0.010\,[S((\psi-0.03)/0.075)-S((\psi-0.105)/0.075)].
\]
Here \(S(x)=0\) for \(x\le0\), \(S=1\) for \(x\ge1\), and \(S=[1+\exp(1/x-1/(1-x))]^{-1}\) in between. The forcing ends once \(\psi\) crosses 0.18. Both the linear driving term and the quadratic \(v^2\) term are included. The mass ratio follows the same smooth profile \(R=0.95+0.10 S((\psi-0.53)/0.003+0.5)\). We use \(mL=1\) in the reference case. The two-body opening is \(R=1\), reached at \(\psi=0.53\). Decay is disabled even after that point in the stable continuation.

The equations are
\[
 \ddot\chi+3H\dot\chi+\omega^2u=0,\qquad
 U_{\chi,\psi}=\omega\omega_{,\psi}u^2-\omega^2u v_{,\psi}.
\]
The \(\psi\) equation receives \(-U_{\chi,\psi}\). The wall energy and pressure receive \(\rho_\chi=(\dot\chi^2+\omega^2u^2)/2\) and \(p_\chi=(\dot\chi^2-\omega^2u^2)/2\). The electromagnetic field, the scalar, the oscillator and gravity therefore evolve together. The comoving source terms satisfy
\[
 a^3\rho_\chi + \int 3Ha^3p_\chi\,d\tau
 = \underbrace{-\int a^3\omega^2u v_{,\psi}\dot\psi\,d\tau}_{W_v}
  +\underbrace{\int a^3\omega\omega_{,\psi}u^2\dot\psi\,d\tau}_{W_M}.
\]
The driver-scalar ledger includes \(-W_v-W_M\). The model has no external clock forcing \(\chi\). Its chosen interaction \(v(\psi)\) is phenomenological.

## Meaning of the supply threshold

For an oscillator with fixed frequency, \(I=a^3\rho_\chi/\omega\) is its instantaneous action per unit coordinate volume about the current minimum. We compare to the former prescribed parent rest density \(D_0=10^{-4}\), for which
\[
 I_{\rm req}=\frac{a_i^3D_0}{2(mL)R_i}=0.00042105263\quad(mL=1).
\]
The full unit-three-sphere factor and common wall-action normalization cancel in the ratio. In a coherent-state interpretation the displaced mean field contributes this action in the appropriate quantum normalization. The calculation uses a classical homogeneous degree of freedom. It supplies neither vacuum fluctuation production nor a microscopic decay probability. Oscillatory pressure, changing frequency and the local Hubble rate determine whether a pressureless-particle description is justified. A rapid mass ramp can change the relationship between energy and an adiabatic particle number.

## Source, gate and late survival

The 48-element launch and six-element late solution bounces at \(a=1.88106\). The source turns on at wall time 0.19684 and ends at 1.32793. With \(v=0\) again, \(I=0.0006375323\), or 1.51414 times the target. At the beginning of the mass change the ratio is 1.75441. The opening occurs at time **5.76487**, radius **7.88555** and action **0.0007046384**. At that event the comoving work values are \(W_v=0.0012385805\), \(W_M=0.0000165590\), and \(\int 3Ha^3p_\chi\,d\tau=-0.0001541373\). For the numerical trajectory through opening, the maximum oscillator, driver-scalar and total-wall energy residuals are \(1.83\times10^{-16},1.06\times10^{-10},4.54\times10^{-10}\). The Friedmann residual is \(1.28\times10^{-10}\). The minimum sampled bulk mass and lapse are 0.10013 and 0.99778.

At opening \(|\dot\omega|/\omega^2=5.1852\) and \(|H|/\omega=0.17150\), so the 1.674 ratio is a particularly poor place to identify a conserved particle count. We therefore also let the stable field complete its mass change and evolve to the first expansion maximum. The final constant mass \(\omega=2.1\) begins at time 5.77456 with \(I/I_{\rm req}=1.60040\). The radius reaches **20.07573** at time **23.73048**. The 11 complete expanding half-cycles between adjacent zeros of \(\dot\chi\) give endpoint ratios from **1.33945 to 1.35045**. In the last half-cycle \(\int p_\chi d\tau/\int\rho_\chi d\tau=-0.0004560\). The volume-weighted ratio is zero to numerical precision. For free motion at fixed mass, \(d(a^3\chi\dot\chi)/d\tau=2a^3p_\chi\), so this last property is an exact virial check independent of how slowly the geometry changes.

A lower bound covers phases between the stored samples. At constant mass and \(v=0\), \(d\ln I/d\ln a=-3p_\chi/\rho_\chi\). Since \(|p_\chi|\le\rho_\chi\), throughout an expanding interval with radii \(a_i\) and \(a_j\),
\[
 I(a)\ge I(a_i)(a_i/a_j)^3.
\]
From mass completion to turnaround the smallest such bound is **1.13758 times the target**. This uses the numerical state at the start of each interval and the exact differential inequality within it. It tests whether oscillation phase alone can explain the excess. Physical decays during this continuation would remove parent energy.

## Controls

| Change | Ratio after source | Ratio at opening | Outcome |
|---|---:|---:|---|
| Reference, 48 early elements | 1.51414 | 1.67352 | Opens |
| 24 early elements | 1.51416 | 1.67403 | Opens |
| 64 early elements | 1.51414 | 1.67347 | Opens |
| No incident pulse | — | — | Contracts past radius 1.6 with both scalars at zero |
| No driving term | 0 | 0 | Opens with zero supply |
| Reverse the sign of \(v\) | 1.51414 | 1.67352 | Same stress and opening as reference |
| Driving amplitude 0.002 | 0.05972 | 0.06286 | Opens below target |
| Driving amplitude 0.005 | 0.37439 | — | No opening before radius 120 |
| Driving amplitude 0.008 | 0.96393 | — | No opening before radius 120 |
| Driving amplitude 0.015 | 3.47174 | 2.80119 | Opens |
| Driving amplitude 0.020 | 6.34335 | — | No opening before radius 120 |
| Pulse dose 0.12295 | 1.52441 | 1.69245 | Opens |
| Pulse dose 0.12305 | 1.50413 | — | No opening before radius 120 |

The parameter controls use 18 early and three late elements unless a grid is specified. The sign reversal uses the same mesh as the reference. Action after the excitation and reaching the later decay condition are different tests. Stronger driving can move the scalar trajectory far enough to miss the opening throughout the calculated interval.

A second selected case uses \(mL=16\), \(g=0.0013\) and support \([0.03,0.105]\). Its larger frequency produces 164 complete constant-mass half-cycles before turnaround. Their endpoint ratios lie between 1.82223 and 1.82272 using that mass’s own target \(0.00002631579\). The all-phase conservative bound is 1.77693. The mass-dependent targets refer to the same former input rest-energy density and different particle numbers.

The earlier rapid heating closure used an effective width scale \(\kappa=2000\). Interpreted as a microscopic two-body width, at \(R=1.05\) it gives \(\Gamma/M=[\kappa/(2mL)](1-R^{-2})^{3/2}\simeq28.35\) for \(mL=1\). A particle with that assigned width cannot be treated as a narrow resonance. Joining this production stage to daughter creation calls for a consistent decay model, including energy transfer, charge birth and feedback on the same geometry.

## Reproduction

```sh
python production.py --jobs 3 --output replay.json
python verify.py replay.json
```

The runner uses Python 3.9 or later and a C++17 compiler. To rerun just the opening and survival trajectories, add `--only primary stable heavy_stable` to the first command. The full verifier also checks null and sign-reversed controls. `PARENT.json` contains the measured states, native event values, each parameter set, the case-specific stopping reason, the solver hash and half-cycle averages. Blank opening results mean that event was not reached during the reported run.
