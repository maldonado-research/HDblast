# Quantum source work, daughter radiation and their feedback on the B1 shell

The rolling B1 shell already produces excitations without an incident pulse. In the finite-band model studied here, energy supplied by the rolling shell scalar drives a Gaussian quantum field, and a specified decay transfers energy into massless daughters. The field and daughters supply density, pressure and reciprocal force to the original five-dimensional solver as it evolves. At the common readout $s=H_0\tau=2.10$, the near-compensated positive input decreases the daughter density by **0.1714%**, while the same positive front with distant compensation increases it by **0.0994%**. The near-compensated negative input gives **+1.0167%**. These percentages refer to the quiet radiation density on the same mesh and in the same model.

| Fine-grid coupled result at $s=2.10$, $y=1$ | Quiet | Near positive | Near negative | Distant positive |
|---|---:|---:|---:|---:|
| First $\phi_b=0.5$ crossing, $s$ | 1.486093 | 1.485455 | 1.486725 | 1.485455 |
| Daughter density $\rho_R/H_0^4$ | 37.15759 | 37.09390 | 37.53539 | 37.19453 |
| Remaining massive density $\rho_\chi/H_0^4$ | 44.52235 | 44.43760 | 44.94708 | 44.54733 |
| Comoving work supplied by the rolling scalar | 157.10694 | 156.80686 | 158.62604 | 157.20028 |
| Comoving work against expansion | 2.99566 | 2.98994 | 3.02267 | 2.99643 |
| Native daughter-radiation contribution to $H^2$ | 0.22512% | 0.22478% | 0.22729% | 0.22535% |

Comoving energies use $a_{*,0}=1$, the quiet crossing scale factor in the preceding B1 calculation on that mesh, and densities in units of $H_0^4$. Percentages in the final row refer to the measured expansion rate of each solution. The recorded fine-grid daughter terms remain below 0.228% of $H^2$ through this cut. Additional species, higher momenta, an altered preparation, or extended continued production can change that result. A massless-fluid equation of state and energy transfer do not demonstrate a thermal distribution. [Full numerical report](COUPLED.json) · [Original rolling calculation](https://doi.org/10.5281/zenodo.24355584)

## The finite-band shell model

Use one minimally coupled real scalar and $m^2/H_0^2=G^2(\phi_b-\phi_*)^2$ with $G=100$ and $\phi_*=0.5$. In the following expressions use $H_0$ units, $a=\exp(A_b-A_{*,0})$, $h=d\ln a/ds$, $p_k=k/a$, $E_k=\sqrt{p_k^2+m^2}$ and $d\nu=k^2 dk/(2\pi^2)$ on $0\le k\le8$. The production and decaying gas prescription used in the earlier B1 model is replaced. All three original junctions, the bulk Einstein–scalar equations, geometric reference and bulk RK4 remain in use.

For each real oscillator take the field $\chi$, canonical momentum $\Pi=a^3 d\chi/ds$ and covariance entries $X=\langle\chi^2\rangle$, $P=\langle\Pi^2\rangle$, $C=\langle\{\chi,\Pi\}\rangle/2$. At fixed $a,\phi_b$ the canonical ground has $X_0=(2a^3E_k)^{-1}$, $P_0=a^3E_k/2$, $C_0=0$. Use the explicitly **instantaneous, finite-band subtraction** of $E_k/2$ in the comoving Hamiltonian, including its derivatives with respect to scale factor and scalar. This choice fixes the reference stress and source throughout the calculation. It is a model convention with a preferred comoving cutoff, and supplies no covariantly renormalized vacuum or all-mode stress. Write $D_X=X-X_0$, $D_P=(P-P_0)/a^6$, and
$$
\rho_k=\tfrac12(D_P+E_k^2D_X),\qquad n_k=a^3\rho_k/E_k .
$$

The massive density, physical pressure and reciprocal scalar current are
$$
\rho_\chi=\int \rho_k\,d\nu,\quad
p_\chi=\tfrac12\int [D_P-(m^2+p_k^2/3)D_X]d\nu,\quad
j_\chi=G^2(\phi_b-\phi_*)\int D_X d\nu .
$$
In operator form $-3a^3p_\chi$ is the derivative of the subtracted comoving Hamiltonian with respect to $\ln a$ at fixed canonical covariance, and $a^3j_\chi$ its $\phi_b$ derivative. Pressure and current retain phase-sensitive correlations absent from a number-only replacement.

The closed massive-field control uses $y=0$. The dissipative comparison uses $y=1$ and the specified instantaneous-basis decay rate $\gamma_k=y^2m^2/(8\pi E_k)$ in proper time. With a dot denoting $d/ds$, the sector evolves by
$$
\dot X=2C/a^3-\gamma_k D_X,\quad
\dot P=-2a^3E_k^2 C-\gamma_k(P-P_0),\quad
\dot C=P/a^3-a^3E_k^2X-\gamma_k C .
$$
This Gaussian vacuum-attenuation channel is completely positive. The propagated covariance is a positive weighted combination of symplectically transported physical covariances, including the initial state and the instantaneous ground states inserted by attenuation. The weights sum to one. The microscopic decay kernel is a Markovian modeling assumption, especially across the short formation interval.

The massless daughters obey $\dot\rho_R+4h\rho_R=Q$, with $Q=\int\gamma_k\rho_k d\nu$. The loss from the massive sector cancels this gain in the total Ward identity,
$$
\dot\rho_\chi+3h(\rho_\chi+p_\chi)=j_\chi\dot\phi_b-Q,\qquad
a^3(\rho_\chi+\rho_R)=\int a^3 j_\chi\,d\phi_b-\int 3 a^3(p_\chi+\rho_R/3)\,d\ln a .
$$
Both work integrals advance inside the native RK4 stages from zero at preparation. The maximum recorded discrepancy divided by $\max(1, |\int a^3j_\chi d\phi_b|)$ is $4.2\times10^{-9}$ across 32 baseline runs through the cut. A small residual establishes internal numerical agreement under this closure. [Quantum sector](quantum_shell.py)

## Initial data, units and feedback

Both spatial meshes start with the original rolling B1 branch, $dc=0.01$. At conformal $T=0.70$ prepare the Gaussian ground of the canonical field/momentum Hamiltonian, common to all arms ($s\simeq0.769$). The number, subtracted density, pressure, current and daughter density all start at zero. The density derivative vanishes analytically because the mass and momentum drift of the reference energy exactly matches the drift of the ground-state Hamiltonian expectation. The three junction values and the density derivative entering the momentum constraint therefore match the empty initial matter sector. On the coarse mesh $|\partial_T\rho|<2\times10^{-20}$ in native junction units.

The ground preparation has a higher-derivative corner. The new native pressure and current derivatives at preparation are about $-0.00339$ and $-0.00705$ per conformal time on the coarse mesh. Their induced time derivatives in the $B$ and $\phi$ normal boundary laws differ from those of the preceding empty sector. The method imposes the new boundary data from this cut forward. Moving the preparation time provides a numerical check on its influence. The four incident/quiet arms in each group subsequently share the *same evolved bulk and matter states* up to $T=1.11$.

At $T=1.11$ solve the paired gravitational constraints on that new rolling background for the same pulse shapes as the preceding experiment. Nominal front supports are $[-0.14,-0.10]$ in the conformal radial coordinate. Near compensation starts at $Z=-0.145$ and distant compensation at $Z=-0.60$. All matched first crossings occur after the primary front has access, before either compensation signal can arrive. The distant scalar control can first arrive at $T=1.71$, after the complete native run ending at $T=1.62$. The incoming front includes the metric perturbation required by the constraints. The new negative branch uses near compensation. Each pulse shooting receipt and endpoint mismatch is in the compact dataset. Fine/coarse spatial resolutions are $0.0005/0.001$ at the wall, with CFL 0.5 and the historical constraint damping coefficient 10. The largest post-fork normalized shell Hamiltonian and momentum residuals over both grids and all arms through the cut are $2.8\times10^{-5}$ and $8.3\times10^{-5}$, respectively. These are finite-grid checks.

Set $b=6.49\times10^{-6}$ and $r_b=7.838285538073243$ in solver units. Multiply the hatted stress and source by $b/r_b$ to obtain the native $\kappa_5^2$ junction quantities. Feed total density into the $A$ junction, total density and pressure into $B$, and the reciprocal scalar source into $\phi$. For the passive comparison, propagate exactly the same field and daughter equations with these matter terms omitted from the bulk junctions. The passive walls are empty of gravitating shell matter in this comparison. At the cut, coupling lowers the quiet daughter density by $0.09015 H_0^4$ relative to passive. Its *additional* changes in the near positive, near negative and distant positive pulse-minus-quiet contrasts are $+0.0000395,-0.001002,-0.000261$ in the same units. The sign change between the two positive fronts persists under reciprocal feedback.

The native Friedmann split assigns $r_b^2[\sigma R/18+R^2/36]$ to daughters, where $R=(b/r_b)\rho_R$. It books the mixed term $r_b^2RX_{\rm mat}/18$ with massive matter, $X_{\rm mat}=(b/r_b)\rho_\chi$. Dividing by $h^2$ gives the tabulated fractions. Allocating the mixed term to radiation would change the fine quiet fraction from 0.225122% to 0.225135%. These are contributions to the modeled expansion, distinct from the daughter share of the selected matter band. The largest sampled physical frequency through $s=2.10$ is $51.39H_0$ across all baseline histories, below $M_5=53.61H_0$. This frequency screen gives no error estimate for excluded modes or unknown short-distance physics.

## Sensitivity and reproduction

On the coarse grid, moving the preparation to $T=0.65$ and $0.75$ moves its proper time to 0.7101 and 0.8290. Across all three preparations the quiet daughter density ranges from 36.7794 to 37.3674. The signed near-positive contrast ranges from −0.17443% to −0.17399%, near-negative from +1.01946% to +1.02485%, and distant-positive from +0.09675% to +0.09799%. Doubling the fine-grid Gauss–Legendre quadrature from 24 to 48 momentum nodes changes each contrast by less than 0.0000001 percentage points. These checks cover the tested preparation and quadrature choices. [Sensitivity data](SENSITIVITY.json)

Run the standard-library check from this directory. The full numerical path uses Python 3.12, NumPy 2.3.5 and SciPy 1.17.0. The bulk replays require the authenticated historical input paths listed in the capsule.

~~~sh
python -I -B verify_coupled.py
python -I -B verify_coupled.py --full
python -I -B matrix.py /path/to/repository /tmp/b1-quantum-work
python -I -B matrix.py /path/to/repository /tmp/b1-quantum-work --sensitivity
python -I -B summarize.py --build-from /tmp/b1-quantum-work --archive /tmp/replayed-walls.zip
python -I -B sensitivity.py /tmp/b1-quantum-work
~~~

The first matrix command regenerates eight common preparations and 32 post-fork histories. The sensitivity command adds twelve post-fork histories on three additional preparations or quadratures. The ZIP has one CSV prehistory per group and one wall CSV per arm. Its comoving occupation integrates $n_k$ over $d\nu$. A physical number density is obtained by dividing by $a^3$. Simulations that reconstruct the solver cache from scratch also solve the native B1 static background. [Compact histories](WALLS.zip) · [Native driver](run_coupled.py) · [Independent cut checks](verify_coupled.py)
