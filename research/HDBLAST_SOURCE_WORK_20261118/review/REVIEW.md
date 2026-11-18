# Independent checks of the finite-band matter sector

The run assigns a canonical Gaussian oscillator to each quadrature point. For physical frequency $E$ in units of $H_0$, let $b_k=a^3 E$, $U_k=\langle q_k^2\rangle$, $P_k=\langle p_k^2\rangle$ and $C_k=\langle\{q_k,p_k\}\rangle/2$. Its comoving excitation energy is $E n_k=(a^3 E^2U_k+P_k/a^3-E)/2$. At fixed canonical covariance, derivatives of that energy with respect to log scale factor and scalar background give $-3a^3p_k$ and $a^3j_k$. This also specifies the pressure and force contributed by the instantaneous ground-energy subtraction. The positive quadrature nodes make $E>0$ at the mass crossing. There the force and assigned damping rate vanish, while the covariance can retain excitation. The physical number density is the weighted sum of $n_k$ divided by $a^3$.

On 48 deterministic pure and thermal squeezed states, finite differences agree with the pressure and force derivatives to scaled errors below $5.3\times10^{-10}$ and $2.2\times10^{-10}$. The three time derivatives used by the junction solver agree with a separate directional difference to $1.5\times10^{-8}$. The native unit conversion $B/r_b$ agrees to $2\times10^{-16}$. For the passive control, all three additional junction values and their derivatives vanish while the local matter evolution is identical on a prescribed wall. Decay deposits $\sum w_k\Gamma_k E_kn_k/a^3$ into radiation. The physical total continuity equation closes to $1.3\times10^{-13}$ on the test states. These errors are normalized by the larger of one and the absolute reference.

At preparation the ground state adds zero density, pressure and force, and its total density derivative vanishes to roundoff. The native Hamiltonian and momentum constraints can therefore inherit their pre-existing discretization residuals. The initial pressure and force derivatives are finite. Switching on the band can introduce a corner in the time derivatives of the metric and scalar boundary conditions. The model makes a choice of initial surface and vacuum. The native sensitivity runs give its measured effect on the readout.

## Oscillator transfer on the numerical walls

For each frozen midpoint scale and mass the transfer calculation applies the exact Hamiltonian rotation
$q\mapsto q\cos(E\Delta s)+(p/b_k)\sin(E\Delta s)$,
$p\mapsto p\cos(E\Delta s)-b_k q\sin(E\Delta s)$.
The prescribed vacuum damping follows by mixing the rotated covariance with the instantaneous ground covariance with weights $\exp(-\Gamma_k\Delta s)$ and $1-\exp(-\Gamma_k\Delta s)$. The radiation accumulator receives the comoving excitation energy removed in the step times $a$. The method checks $U_kP_k-C_k^2\ge 1/4$ throughout the transfer, up to roundoff. It uses the same finite-band model as the native method and an independent propagator.

The compact wall archive supplies each history's proper clock, scale, scalar, and their first derivatives. Cubic Hermite interpolation supplies 2,200 and 4,400 proper-time midpoints from preparation to $s=2.10$. Second-order Richardson extrapolation gives the reported predictions. The native comparison uses eight neighboring archived samples for its endpoint interpolation. Oscillatory pressure and force are more sensitive to the sampled geometry and the endpoint than the integrated number and radiation.

| Maximum scaled error over 16 arms per mesh | Coarse | Fine |
|---|---:|---:|
| Comoving number | 5.36e−8 | 8.77e−9 |
| Radiation density | 4.26e−8 | 5.46e−9 |
| Scalar excitation density | 5.39e−8 | 8.80e−9 |
| Scalar excitation pressure | 8.67e−5 | 4.63e−6 |
| Scalar force | 1.71e−5 | 9.77e−7 |

The test reproduces the matter evolution on all 32 simulated walls. The bulk geometries themselves come from the original native solver. In the reciprocal and passive controls the matter sectors share initial data at preparation. Their evolution already differs at the pulse fork. On the fine grid for $y=1$, the passive and reciprocal forks have $\phi_b=0.292126$ and $0.292169$ respectively. The final quiet radiation shift of about −0.242% therefore measures the effect of coupling over the full history. Subtracting the quiet response within each history and then comparing those pulse contrasts gives a much smaller effect. With the eight-point endpoint interpolation, the reciprocal-minus-passive shifts of the normalized near-positive, near-negative, and distant-compensator radiation contrasts are −0.000309, −0.000224 and −0.000488 percentage points. These values belong to the constraint-matched pulse family, the chosen damping law, cutoff and initial state.

Run from the capsule root, with the previous pulse capsule in the same parent directory. Both commands write JSON to standard output.

```sh
python review/check_interface.py
python review/replay_native.py coupled/WALLS.zip ../HDBLAST_PULSE_PARTICLE_RESPONSE_20261118
```

The saved outputs are `INTERFACE_CHECK.json` and `NATIVE_TRANSFER.json`.
