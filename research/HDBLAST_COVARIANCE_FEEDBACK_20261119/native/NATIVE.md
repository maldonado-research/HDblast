# Four preparations in the evolving five-dimensional shell

A finite-band Gaussian scalar field can exert the same initial energy density, pressure and force on the B1 shell in four different preparations. Advancing that field and the gravitational and scalar equations together yields distinct later densities and expansion rates. In the quiet case, changing the sign of the field–momentum covariance at the fork increases the final fine-grid parent energy density by 1.37224% and increases the dimensionless expansion rate by $4.86705\times10^{-6}$. An energy-matched canonical impulse increases the density by 0.66119% and the expansion rate by $2.50875\times10^{-6}$. [Full numerical results](NATIVE.json) · [Wall time series](WALLS.zip)

## Model and fork

The preceding closed-shell run prepared an instantaneous oscillator vacuum at coordinate time $T=0.70$ for one minimally coupled real scalar with $m/H_0=100|\phi_b-0.5|$ and comoving momenta $0<q/H_0<8$. The matter and bulk evolved together to $T=1.11$ with the published three signed junction laws. The native calculation uses the original finite-band instantaneous-ground subtraction for density, pressure and scalar force, $H_0=1/r_b$, 24 Gauss–Legendre momentum nodes, constraint damping 10 and fourth-order spatial differences away from the boundary. The fine/coarse mesh spacings at the shell are 0.0005 and 0.001, with CFL 0.5. The physical scale factor is $a=\exp(A_b-A_{*,0})$ with the historical quiet-crossing normalization on each mesh. All wall runs here keep the Yukawa coupling at zero. The daughter readout uses these backgrounds perturbatively in a separate calculation. [Historical coupled solver](https://github.com/maldonado-research/HDblast/tree/838834ff58748b559a14726a3a16b8716051be35/research/HDBLAST_SOURCE_WORK_20261118/coupled)

The simulations branch from the same full native state at $T=1.11$ on each mesh. This is $s=H_0\tau=1.30371703$ on the fine grid. Let $X=\langle\chi^2\rangle$, $P=\langle\Pi^2\rangle$, and $C=\langle\{\chi,\Pi\}/2\rangle$ denote the three canonical covariance components for each momentum mode. The kinetic and field moments evolve through the native Klein–Gordon equations. With proper clock $s$ and $E=\sqrt{m^2+(q/a)^2}$ these read

$$
X'=2 C/a^3,\qquad P'=-2a^3 E^2 C,\qquad C'=P/a^3-a^3 E^2 X.
$$

In the retained, erased and reflected preparations, $X$ and $P$ are unchanged and $C$ is assigned the old value, zero and minus the old value, respectively. These preparations have identical oscillator occupations mode by mode. The erased covariance describes a mixed Gaussian state. The reflected covariance has the original determinant for each oscillator. The fourth state implements a uniform canonical momentum shear $\Pi\mapsto\Pi+\kappa\chi$. Its new $C$ and $P$ are $C+\kappa X$ and $P+2\kappa C+\kappa^2X$. We set $\kappa=-2\sum wC/\sum wX$ for $w=q^2dq/(2\pi^2)$, giving the same integrated $P$ and therefore the same initial density, pressure and force. This impulse changes occupation at individual momenta. On the fine mesh $\kappa=-0.31674094$. The finite-band canonical map leaves $XP-C^2$ unchanged. [Fork data](FORK.json)

All four use the same initial bulk geometry and bulk velocities. The normal derivatives in all three junction conditions agree at preparation. So do the native Hamiltonian and momentum constraint residuals and the derivative of the metric $A$ junction. The altered covariances generally change the first derivative of the pressure and force, hence the prescribed derivatives of the $B$ and $\phi$ junctions. Their new values apply for evolution forward from the fork. The original smooth history does not in general have these new derivatives through the preparation event. A separate pair of preparations controls for this difference by matching the first derivatives too. [Initial junction analysis](../initial) · [Matched derivative experiment](../native_smooth/SMOOTH.md)

At the same fork we introduce either no incoming perturbation, the registered positive or negative near-compensated perturbation, or the positive perturbation with distant compensation. Each pulse's compensator was solved on its own prepared Cauchy data. Its compact-support matching residual is below $9\times10^{-16}$ across all runs. The full simulations continue to $T=1.62$, which covers the common proper readout $s=2.10$. The statistics below use linear interpolation of the recorded time series to that common readout. [Evolver](native.py) · [Reduction](summarize.py)

## Native quiet readout on the fine mesh

| Preparation | $H/H_0$ | Parent energy $\rho/H_0^4$ | Parent pressure $p/H_0^4$ | Scalar source in native hatted units |
|---|---:|---:|---:|---:|
| Retained | 0.176688942 | 97.630281 | −31.722768 | 252.814498 |
| Erased | 0.176691376 | 98.300223 | −31.237945 | 253.184881 |
| Reflected | 0.176693809 | 98.970008 | −30.753373 | 253.555455 |
| Canonical impulse | 0.176691451 | 98.275806 | −30.150280 | 251.019714 |

The quiet coarse-grid increases in $H/H_0$ from erasure, reflection and the impulse are respectively $2.49281$, $4.98514$ and $2.58994$ in units of $10^{-6}$. For each preparation, the positive near-compensated pulse lowers the final parent energy density and the other two pulses raise it on both meshes. At fine resolution the respective pulse-minus-quiet density contrasts are $-0.18291,+0.98384,+0.07508$ for retained covariance, and $-0.18504,+0.99742,+0.07669$ for the impulse in units of $H_0^4$. The expansion-rate changes use a different order of signs: the nearby negative pulse increases $H/H_0$, and both positive pulses decrease it.

Summing the scalar work and pressure work gives the native comoving-energy balance. Across all 32 runs the largest drift of $a^3\rho-W_{\phi}+W_p$ relative to the fork is $1.91\times10^{-7}$ in the native hatted comoving energy units. The smallest recorded covariance determinant is 0.249999999999759. The largest normalized Hamiltonian and momentum residuals recorded near the shell node through the readout are $2.78\times10^{-5}$ and $8.24\times10^{-5}$. These are finite-resolution diagnostics. The retained trajectories also reproduce the original archived coupled wall series. The baseline geometry and the matter exchange are specific to this model.

## Reproduction

The two restart states are bundled with SHA-256 pins. With an HDblast repository containing the historical checkpoint and rolling-pulse scripts, the original $T=0$ to 1.11 evolution regenerates both archived restart states from its source initial data. All numerical arrays agreed bit for bit on the verification machine. Rebuilding the restart states is independent of recomputing the 32 post-fork branches. [Prefix replay](PREFORK.json)

~~~sh
python -I -B reproduce.py /path/to/HDblast
OPENBLAS_NUM_THREADS=1 python -I -B -O reproduce.py /path/to/HDblast --full --prefork --jobs 2
~~~
