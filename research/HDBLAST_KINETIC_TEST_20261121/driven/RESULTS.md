# Following newly born particles through the pulse

We let the charged particles keep their individual momenta from birth. At each encounter their outgoing momenta follow the differential cross section of the same vector mediator used in the fluid calculation. A wall electric field responds to the new current. This replay uses the earlier shell expansion, source of new particles and incoming bulk wave as supplied histories. The separate live calculation also allows the shell and the bulk wave to respond. [Fully coupled calculation](../live/RESULTS.md) · [Original fluid](https://doi.org/10.5281/zenodo.24408924)

## What happens over time?

Each entry below is a mean across four random seeds, with 4,096 representative particles at the end. The first two columns use the current-rest frame of a single sign and flavor. A thermal population would have a longitudinal pressure share of one third. The energy-spectrum discrepancy is the largest difference between the observed and fitted thermal fractions of particles below a common energy threshold. For example, a difference of 0.20 means that for some energy cut the two distributions differ by twenty percentage points in how many particles lie below it. The last column compares the rest-frame internal energy per particle, evaluated in the frame with vanishing energy flux, to total wall energy as defined below. [All runs and controls](DRIVEN.json)

| Channel and reference time | Longitudinal share of stream pressure | Energy-spectrum discrepancy | Fraction never scattered | Internal energy share |
|---|---:|---:|---:|---:|
| Dedicated, half early electrical work | 95.8% | 0.326 | 32.6% | 59.3% |
| Shared, half early electrical work | 91.3% | 0.276 | 52.3% | 78.0% |
| Dedicated, first fluid hot crossing | 45.2% | 0.198 | 20.8% | 90.8% |
| Shared, first fluid hot crossing | 45.2% | 0.256 | 31.7% | 90.1% |
| Dedicated, original hot radius doubled | 33.5% | 0.011 | 0.3% | 97.5% |
| Shared, original hot radius doubled | 33.7% | 0.048 | 4.7% | 97.7% |

“Half early electrical work” here uses half of the integrated absolute electrical work up to the first hot crossing of the **original fluid**. The other times are also determined before running the particles. The analytic counting test uses the original fluid's collision rate and a whole-run work clock. In this replay, the encounter rate is computed from the actual particle pairs, so it can differ even when the marker times are close. [Event-counting test](../history/RESULTS.md)

In both channels the collisionless controls have a time-averaged internal energy share above 92% near reference doubling, even as essentially all pressure in an individual stream remains axial. Particles born at different times can acquire opposite momenta from a reversing field. Summing those particles can yield a large rest-frame energy spread without filling transverse directions. Collisions in the tested model can later redistribute the momentum over directions. Agreement of the late pressure and of the one-dimensional energy spectrum separately tests two properties of equilibrium. Additional correlations, chemical equilibration and production of ordinary radiation require their own tests. [Angular and energy diagnostics](DRIVEN.json) · [Three-axis calculation](ORIENTATION.md)

## Birth, acceleration and collisions

Each representative has fixed comoving rest weight. We stratify the accepted accumulated births into equal weights and interpolate the time of each midpoint. The parent mass sets the daughter's birth momentum. The birth direction is sampled from the original angle-dependent weight. The comoving transverse momentum is constant between encounters and the comoving axial momentum changes by the electric impulse. Independent random directions in azimuth implement axial symmetry. In the shared-channel case, a collision can change flavor. Two of the three equally probable destinations have a perpendicular acceleration axis. The more detailed three-axis implementation retains the azimuth, follows its magnetic rotation, and explicitly selects the destination axis. [Axisymmetric implementation](ensemble.cpp) · [Three-axis implementation](ORIENTATION.md)

Use the supplied radius $a$, background electric field $X_0$, per-particle charge-to-mass ratio $Q$, and post-release gauge coefficient $K$. Define $Y=a^2 K(X-X_0)$ and $Z=a(C-C_0)$ for the wall gauge potential $C$. The inherited wall Maxwell equations with an unchanged incoming characteristic give

$$
\dot Y=-4KZ/a-2Y/K-Q(w\sum_i v_{z,i}-D_0 v_0)/a,\qquad
\dot Z=Y/(aK),\qquad \dot I=QaX_0+QY/(aK).
$$

Here $D_0 v_0$ is the original comoving current divided by $Q$, $w$ is the comoving rest represented by one sample, and $I$ is the signed impulse. A fixed-field control sets the response to zero. The density used for the binary hazard is the wall-frame number per volume of one sign and flavor, with the relative flux included for each individual pair. The collision sampler proposes candidate pairs at a bounding rate and accepts with probability proportional to $\sigma(s)v_{\rm M}$. It samples the massive-vector center-of-momentum angular law and boosts each outgoing pair back to the wall frame. The paired charge-reflected events maintain equal sign populations. [Pair selection and literature](../SOURCES.md)

The replacement total energy density is the reference total with its old gas and wall field energies removed and the re-evolved ones inserted. The internal energy share is $n_{\rm lab}(e_L/n_L-m)/\rho_{\rm repl}$ summed over streams, where the Landau frame diagonalizes the stream energy flow. The spectral comparison uses the particle-current rest frame, weights each wall sample by $E'/E$, and fits the temperature to its mean rest-frame energy. The tabulated instantaneous results can fluctuate when the axial current reverses. A separate integral averages over the preceding $0.25L$ of wall time at onset, doubling and turnaround. [Spectral estimator](spectra.hpp) · [Sampling windows and outputs](run.py)

At each split transport step the change of particle energy is decomposed into electrical work and expansion losses using differences of squared momenta and rationalized square roots. The ledger also transports each sampled birth from its own time to the collision midpoint. Comparing final particle energies with injection plus work minus redshift checks this accounting independently of the reported pressure fraction. The simulations vary the kinetic time step, the native wall step, sample size, and seeds. Small late unscattered fractions and instantaneous onset values can vary across these choices. At 8,192 and 16,384 particles with two seeds each, the shared unscattered fraction at doubling is 9.37% and 8.25%, versus 4.69% at 4,096. The dedicated fraction moves from 0.31% at 4,096 to 4.53% and 4.71% on the larger runs with half the kinetic step. At 16,384 particles, the mean late energy-spectrum discrepancies are 0.047 and 0.082 and the axial pressure shares are 0.336 and 0.337 for the dedicated and shared cases, respectively. These results give a consistent direction of pressure relaxation while leaving the exact late population uncertain. [Higher resolution](LARGE.json) · [Three-axis controls](ORIENTED_CHECK.json)

From the study root, run `python -I -B study.py --verify` to regenerate these and the accompanying calculations. The individual driver takes `--native` pointing to the output directory of `history/export.py` and `--output` pointing to a new JSON file. Add `--large` for the two larger particle counts.
