# What field survives from the early pulse at the wall?

A short tangential disturbance can seed a field that persists until particle birth in this expanding-wall example. The size depends strongly on its spatial wavelength and the prescribed history of the gauge coupling. This calculation takes the incoming waveform **at the wall** as its input. The companion bulk calculation instead propagates a perturbation from the outer source. [Results](transfer.json) · [Finite-collar calculation](../geometry/README.md)

## Field produced from an initially empty mode

Write the extra incident characteristic as $i_y(t,x)=\alpha i_{\rm axial}(t)\cos(kx)$, with $i_{\rm axial}$ measured directly in the reference five-dimensional simulation. The additional incoming mode has peak amplitude fraction $\alpha$ and spatially averaged squared amplitude fraction $\alpha^2/2$. Its absolute spatial phase is arbitrary because the background and particle births are translation invariant in the local patch. The time origin fixes the signed waveform's phase relative to the evolving wall. We also tested causal delays of 0.10 and 0.50 in wall units. These represent changed boundary inputs.

For the fine reference grid the radius is 10.302 when births start at wall time 6.251625. Both scattering prescriptions share that prebirth geometry. For unit $\alpha$, the following spatial RMS values include the electric and magnetic fields in quadrature. The final column gives the amplitude needed to meet the earlier necessary initial-field bound for 1 percentage point accuracy at the 0.1% electrical-work marker, with uniform density and no later incident drive on the dedicated-interaction history. Satisfying that bound leaves the distribution to be tested by particle evolution.

| Comoving wave number | Physical wave number at birth | RMS field at birth | Amplitude fraction for the early bound |
|---:|---:|---:|---:|
| 2 | 0.194 | 0.1464 | 0.0608 |
| 5 | 0.485 | 0.6141 | 0.0145 |
| 10 | 0.971 | 0.5788 | 0.0154 |
| 20 | 1.941 | 0.000712 | 12.51 |
| 50 | 4.854 | 0.002316 | 3.84 |

For $k=5$ the corresponding fraction in the shared-interaction history is 0.00993. Comparing a continuous local wave number on a finite logarithmic grid from 0.1 to 100 in birth-time physical units gives a largest sampled RMS response of 0.7620 at 0.6978. The discrete closed-wall modes and their incident amplitudes require separate treatment. The planar-patch approximation is especially restrictive for long wavelengths.

Truncating the waveform once it contains 99.008% of the prebirth time integral of $i^2$ leaves an electric field of $-0.8387$ in mode 5 at birth, compared with $-0.8684$ for the entire waveform. The truncation time was 0.00555, approximately 6.246 wall times before birth. The remainder alone gives $-0.02968$. These signed contributions sum. Under delays of 0.10 and 0.50 the opening RMS changes from 0.6141 to 0.6150 and 0.6241. At the dedicated 1% work marker the prebirth drive alone gives RMS 0.58996, and applying the source only after birth gives 0.00000424 in the passive field calculation.

## Where the energy comes from

For physical $E_y=E(t)\cos(kx)$ and $B_z=B(t)\sin(kx)$, define $P=Ka^2E$ and $b=a^2B$. The local Maxwell boundary model in proper time gives

$$
\dot P=-kKb/a-2P/K+4a^2 i,\qquad \dot b=kP/(aK).
$$

Here $K$ is the time-dependent coefficient of the gauge kinetic term and $a$ the recorded wall radius. Magnetic and electric curl terms exchange energy. The spatially averaged comoving field energy $U=(P^2/K+Kb^2)/4$ obeys

$$
\dot U =2a^2 iP/K-(P/K)^2+\tfrac14\dot K\,[b^2-(P/K)^2].
$$

The first two terms comprise the net boundary exchange. The third measures the work done on the mode by the evolving gauge coupling. With unit input in mode 5, integrated boundary driving supplies 1.9722, the radiative loss is 103.8959, and the coupling supplies 5,495.5122. The resulting comoving energy at birth is 5,393.5884. **The incident pulse seeds the mode. The evolving gauge coefficient supplies its net stored energy in this example.** Holding $K$ at its initial or final value on the same trajectory yields opening RMS fields of 0.001440 and 0.00003276, respectively. These controls isolate the coefficient change at fixed waveform and geometry. Recomputing the background with a different gauge evolution would change those histories.

The fields scale linearly and energies quadratically with $\alpha$ in this probe. At $\alpha=0.0145$ the spatially averaged opening energy in mode 5 is about 0.185% of the recorded homogeneous axial electromagnetic energy. This ratio only compares energies at birth. Testing the inhomogeneous scalar, metric and matter response during amplification still requires their coupled evolution.

## Numerical controls

The centered implicit integrator advances canonical fields using midpoint sources, retaining the work of each change of $K$ separately. An exact constant-coefficient matrix exponential checks amplitudes and signs. Refining the time step twice reduces the maximum field error by factors near four. The comoving $k=2$ equations with unit input reproduce the original axial prebirth equation on the recorded geometry. On the finer reference, sampled field errors relative to each peak are $1.51\times10^{-6}$ for electric and $2.31\times10^{-6}$ for magnetic. The energy-balance residual normalized by each mode's running peak energy is below $1.4\times10^{-11}$ for the four reference runs. In the finer history, keeping alternate source knots changes the mode-5 opening field vector by less than $7\times10^{-8}$ relative. Coarse and fine native background runs give 0.614086 and 0.614072 for the same RMS at birth. [Integrator](evolve.cpp) · [Capture](capture.py) · [Checks and field histories](run.py)

`run.py --prior PATH_TO_PREDECESSOR --output-dir OUT` rebuilds four short stretches of the predecessor's five-dimensional simulation, measures the incoming characteristic, compiles the mode solver and writes the results. The four generated `traces/*.npz` carry physical electric cosine and magnetic sine coefficients at every original particle-history knot, with and without continued driving after birth. The mode labels in `k` are physical at birth and `k_comoving` are comoving. After birth these externally driven traces omit newly induced daughter currents. They provide the reference input for the separate particle calculation. No extra transverse drive and zero initial transverse fields give the zero perturbation exactly.
