# Which angular parts of an early pulse can leave a field on the wall?

The original HDblast pulse lives in the lowest angular Maxwell pattern on a round three-sphere. Here an additional source sends a prescribed angular pattern from the outer edge of the five-dimensional region. The added field starts at zero throughout the region and on the wall. We follow the wave through the radial mesh, its reflection, and the electric and magnetic fields it leaves on the expanding wall. The geometry, original matter, and scalar-dependent wall gauge coefficient follow the archived benchmark. The additional field acts as a probe.

## Field present at the first particle births

The first column below is the magnitude of the curl eigenvalue on a unit sphere. Dividing by the wall radius gives the physical wavenumber in inverse model lengths. At birth the radius is 10.302. Each row represents a separate run with one normalized source and the same emitted time profile. The fields scale with its amplitude. The energies scale with its amplitude squared. Signed field amplitudes specify phase relative to the applied pulse.

| Angular eigenvalue | Wall electric field at first birth | Signed magnetic coefficient | Stored comoving energy | Net work from incoming and outgoing waves | Work from changing gauge coefficient |
|---:|---:|---:|---:|---:|---:|
| 2 | −0.20708 | 0.000172 | 613.39 | −9.83 | 623.22 |
| 5 | −0.83827 | −0.004703 | 10,051.22 | −189.63 | 10,240.85 |
| 7 | −0.99739 | −0.008838 | 14,229.86 | −269.23 | 14,499.09 |
| 10 | −0.69774 | −0.009371 | 6,964.66 | −130.82 | 7,095.48 |
| 20 | −0.00372 | −0.000105 | 0.20 | 1.85 | −1.65 |

The largest birth field among all integer eigenvalues 2 through 40 occurs at 7, on both the time-refined and radially refined runs. With that pattern at 1% of the reference launch amplitude, the electric coefficient at birth is about 0.00997. Its stored energy is 0.232% of the reference axial sheet energy at birth. The energy ledger shows that the varying wall coupling supplies the largest contribution. The initial wave supplies a seed to that process. For comparison, the pulse with eigenvalue 50 reaches the wall with peak inward characteristic 8.49, close to the 8.50 of the original wave. Its remaining electric field at first birth is below \(10^{-9}\) in each computed mesh. The tiny residuals for the high modes vary with the radial discretization. [Primary measurements](PROBE.json) · [Contiguous integer scan](SPECTRUM.json)

## Which change on the wall supplies the field?

We repeated the collar calculation on the same geometry with the same applied pulse. In one test only the probe's gauge coefficient stops evolving when the scalar reaches the start of the final release. In the original model that stage lowers the coefficient from 3400.04 to 2.54. For eigenvalue 7 the field at particle birth falls from magnitude 0.99739 to 0.0007670 when that decrease is suppressed. Its comoving energy falls from 14,229.86 to 24.835. The final decrease therefore accounts for a factor of about 1,300 in the electric amplitude and 573 in the stored energy in this comparison. On the doubled radial grid the corresponding amplitudes at its own birth are 0.99986 and 0.0007692. In a second control the probe coefficient stays at its original launch value of 0.001 throughout. Its birth electric amplitude is 0.00000440. The background scalar and geometry remain those of the benchmark, allowing the probe boundary response to be varied on its own. [Gauge evolution controls](GAUGE.json)

## Bulk mode and boundary conditions

Write the metric as \(h_{ab}(x)dx^a dx^b+r(x)^2d\Omega_3^2\). On the unit three-sphere take a divergence-free one-form \(Y\) with \(*dY=s\kappa Y\), helicity \(s=\pm1\), integer \(\kappa\ge2\) and \(\langle |Y|^2\rangle=1\). Substitution of the tangential Maxwell potential \(A=U(x)Y\) into the bulk action yields, per unit sphere volume,

\[
S_{\rm bulk}= -\frac12\int\sqrt{|h|}\left[r D_aU D^aU+\frac{\kappa^2}{r}U^2\right]d^2x,
\qquad D_a(r D^aU)-\frac{\kappa^2}{r}U=0 .
\]

The wall sheet contributes \(\frac12\int K(a\dot U^2-\kappa^2 U^2/a)dt\). On the two symmetric retained collars, in the archived Gaussian-normal sign convention, define electric coefficient \(X=\dot U/a\), signed magnetic coefficient \(b=\kappa U/a^2\), incoming characteristic \(i\) and outgoing characteristic \(o=X-i\). The mode evolves according to

\[
\frac{d}{dt}(a^2 K X)=4a^2i-2a^2 X-\kappa K a b,
\qquad \frac{d}{dt}(a^2 b)=\kappa a X .
\]

The factor four in the source follows from the two-sided junction and the characteristic normalization used for the background. In the interior, the angular potential contributes \(-\kappa^2 N U/r^3\) to the electric evolution. The lapse, areal radius and their gradients are evaluated at every node and Runge–Kutta stage of the original gravitational evolution. The perturbation uses the same upwind radial fluxes and interpolation at mesh transfer as its parent. [Angular harmonic identities, equations 13–20](https://arxiv.org/abs/1709.08020) · [Original coframe and boundary equations](https://github.com/maldonado-research/HDblast/tree/315e1edf085d601a93a616947db75d0929f9f8fe/research/HDBLAST_INCIDENT_WAVE_WALL_JUNCTION_20261120) · [Numerical implementation](mode_rhs.hpp)

The source follows the original waveform in the proper time of the outer tube. The spatial normalization makes an additional angular amplitude \(\alpha\) correspond to \(\alpha^2\) of the reference angular-mean squared incoming characteristic there. This refers to the launched component. The net work on the wall includes the outward characteristic as well. Orthogonality and the spatially uniform background coefficients make the mode equations diagonal. An initially absent higher harmonic with zero higher-harmonic source stays zero, including across a time-dependent \(K\). The calculation specifies the additional source at the outer tube. Generation of this angular structure in a parent bulk event, and the metric and matter changes it would produce, remain open. Tests with finite amplitude must account for the possible importance of the additional stress, including local cross terms with the original fields.

Under the conventional Lorentz force \(q(\mathbf E+\mathbf v\times\mathbf B)\), Faraday's law contains minus the spatial curl of \(\mathbf E\). For a single helicity eigenform use \(\mathbf E=XY\) and \(\mathbf B=-s bY\) in that same spatial orientation. For a local linearly polarized superposition of opposite helicities one may use \(\mathbf E=X\cos(kx)\hat y\) and \(\mathbf B=b\sin(kx)\hat z\). A plane cosine of unit peak amplitude has spatial RMS \(1/\sqrt2\). A local plane representation changes the spatial pattern sampled by particles relative to an individual global harmonic.

## Changing polarization while holding the launch power fixed

A concrete source can provide the additional angular component without increasing the power anywhere on the emitting tube. In complex coordinates on the unit sphere let \(\eta=-\mathrm{Im}(\bar z_1 dz_1+\bar z_2 dz_2)\) be the unit Hopf form and \(Y_\kappa=\sqrt{\kappa-1}\,\mathrm{Re}[z_1^{\kappa-2}(z_1 dz_2-z_2 dz_1)]\). The Hopf direction is orthogonal to this mode at every point. The squared norm is \( |Y_\kappa|^2=(\kappa-1)|z_1|^{2\kappa-4} \), with unit spatial average since \( |z_1|^2 \) is uniform on \([0,1]\) under the round volume measure. Multiply the original time profile at the tube by the unit one-form

\[
p_\alpha=\sqrt{1-\alpha^2 |Y_\kappa|^2}\,\eta+\alpha Y_\kappa,\qquad |\alpha|<1/\sqrt{\kappa-1}.
\]

The squared length is one at every point and the exact projection on \(Y_\kappa\) is \(\alpha\). The Hopf flow preserves \( |z_1| \), so the axial term also has zero divergence. At eigenvalue 7, a 1% transverse amplitude diverts 0.01% of the angular-mean incoming characteristic power into this mode. It requires a maximum local tilt of 1.40 degrees. The changed axial amplitude and its additional angular structure are second order in \(\alpha\). Their joint effects on matter and geometry would require a new evolution. An emitter that changes only the amplitude or timing of its axial direction, \(f(\Omega,t)\eta\), has zero direct projection on this particular transverse family. Propagation with broken angular symmetry, other harmonics and a physical source of the necessary rotation are separate questions.

## Energy and numerical checks

For an added source with no charged response, the exact sheet energy identity is

\[
\frac{d}{dt}\left[\frac12 a^4 K(X^2+b^2)\right]
 =2a^4(i^2-o^2)+\frac12a^4\dot K(b^2-X^2).
\]

The two work terms have separate integration variables throughout the simulation. The maximum sampled energy residual, divided by the maximum of the stored energy plus the absolute values of the two ledgers, is below \(1.4\times10^{-7}\) on the base runs and \(4.4\times10^{-9}\) on the refined runs for the main set of modes. These values characterize overall energy balance. At high eigenvalues the two work terms nearly cancel and that precision cannot resolve a tiny residual energy at the final instant.

The lowest mode reproduces the original incoming characteristic and electric and magnetic amplitudes to \(2.2\times10^{-14}\) in absolute units throughout the prebirth stage on all primary meshes. A second unforced mode at eigenvalue 50 stays zero. After births the added field continues as a source-free probe of the matter background with evolving boundary waves. The archived axial field includes the charged current. Downstream particle simulations specify their own treatment of the extra current. The measured axial electric and magnetic histories at the archived compact times agree within \(10^{-11}\) on the corresponding original meshes.

Halving the time step changes the eigenvalue-7 electric field at birth by about 0.0021%. Doubling all radial cell counts at fixed Courant factor shifts first birth from 6.251625 to 6.250052. At the refined mesh's own birth time its electric field is −0.99986, a 0.25% change. Comparing at the same time as the original run gives a 0.013% change. The primary studies follow both dedicated and shared late particle interactions. Up to first birth the backgrounds agree. [Reproduction script](probe.py)

The main run stores electric, magnetic and incident radial amplitudes at the archived particle-history times as arrays named `E_cos`, `B_sin`, and `I_cos`. Axis zero follows the explicitly stored integer `kappa`. The `k` array uses inverse physical length at the original birth radius. The `--spectrum` run does this for every integer from 2 to 40. Each run resolves the outer pulse internally from zero even though the exported orbit traces start at birth.
