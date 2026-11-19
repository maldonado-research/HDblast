# A phase-sensitive test with unchanged spectral energies

Consider two descriptions of the outgoing high-frequency bath in the HDblast Gaussian model. Both place the same mean energy in every oscillator. In the second description we erase the relative phase between the bands 60–90 and 90–120. An energy measurement gives the same answer in both. A measurement that recombines the bands gives different results. This comparison concerns pair-channel modes in the surrogate model. Applying such a measurement to physical fermions or to the sky would require a detector and a propagation model. [Test design](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7015108648)

## Correlations at leading order

Fix a native parent momentum $q$. The free parent is $\chi(s)=h(s)\chi_0+k(s)\Pi_0$, with real canonical solutions $h,k$. The band starts in its vacuum. Write
$$
W=\left(\begin{matrix}X&C+i/2\\ C-i/2&P\end{matrix}\right),\quad
v_j=\sqrt{\frac{g_j^2}{2\omega_j}}\int ds\,F(s)e^{i\omega_j\eta(s)}
\begin{pmatrix}h(s)\\ k(s)\end{pmatrix},\qquad
N_{jk}=v_j^\dagger Wv_k .
$$
Here $X,C,P$ are the symmetrized parent moments at the fork, $[\chi_0,\Pi_0]=i$, and $g_j^2$ contains the integration weight. $N$ is the coefficient of $y^2$ in the normally ordered bath correlation $\langle b_j^\dagger b_k\rangle$ in the interaction picture. The conformal clock starts at the fork. The use of $ds=a\,d\eta$ includes the scale factor in the drive. Summing $\omega_j N_{jj}/a_f$ over frequencies and native momentum weights gives the physical energy per normalized comoving volume at the final readout. [Independent emission-energy derivation](../quantum/SPECTRUM.json)

The matrix factors through two parent coordinates, so its rank is at most two per momentum. For a centered, pure, single-mode Gaussian parent, $XP-C^2=1/4$ and $W$ has rank one. The archived retained state satisfies this purity condition to about $1.5\times10^{-12}$ on the fine grid. In the four tested histories the largest second eigenvalue of $N$ is below $3.1\times10^{-15}$ of its trace within each momentum. The bound concerns this leading one-excitation correlation. It neither assigns a rank to the many-body state of the universe nor involves higher-order scattering. A mixed Gaussian control, formed by doubling the symmetrized covariance while preserving the commutator, admits a second packet.

## Two frequency bands

Within each momentum we split the dominant quiet-shell eigenpacket at frequency 90 and normalize the two pieces to obtain filters $L,R$. We lock their phases on the quiet trajectory. Their normalized sum and difference give two hypothetical output ports with individual occupations $I_+$ and $I_-$. The visibility below is $(I_+-I_-)/(I_++I_-)$, after combining the momenta with the archived weights. Each pulsed trajectory uses the filters fixed on the quiet shell. For the control $N_{\rm block}=P_LNP_L+P_RNP_R$, the projections here select the two entire frequency bands. Removing cross-band phase this way leaves every diagonal entry and the total energy fixed. It makes the two port occupations identical.

| Incoming history | Proper high-band energy coefficient | Visibility | Visibility with timing jitter 0.10 |
|---|---:|---:|---:|
| Quiet | 0.9326764 | 7.218% | 0.361% |
| Near negative | 0.9596445 | 8.280% | 0.443% |
| Near positive | 0.9295983 | 6.797% | 0.323% |
| Distant compensated | 0.9368036 | 7.013% | 0.337% |
| Same histories with cross-band phase erased | Same in every row | 0 | 0 |

The energy column is in the same proper comoving normalization as the native quantum calculation. Timing jitter is the standard deviation of a common random time offset in dimensionless conformal units. Averaging the phase of each matrix element gives the exact damping factor $e^{-\sigma_\eta^2(\omega_j-\omega_k)^2/2}$. It preserves the energy. Without jitter, the two quiet-calibrated filters together capture over 99.94% of the modeled high-band occupation on each pulse. Refining the frequency node count per band from 24 to 48 and the Gauss order per native time interval from 6 to 10 changes all visibilities by less than $2\times10^{-12}$. [Results including individual momenta, jitter controls and refinements](PACKETS.json) · [Replay code](packets.py)

As an exploratory sensitivity check, normalize the excitation matrix to unit trace and treat the momentum channels as orthogonal. Its trace distance from the cross-band-dephased state is the sum of the singular values of the cross-band blocks, with the same momentum weights and occupation normalization. On the quiet trajectory this distance is 0.03609 without timing jitter. Even an ideal measurement on a single surrogate excitation, given the two states with equal prior probability, could distinguish them correctly only about 51.80% of the time. The two calibrated ports attain the corresponding trace-distance witness in that case. Visibility describes repeated measurements with a phase reference.

After the gate closes, the free band's total energy and homogeneous radiation pressure depend only on the diagonal occupations. Both comparison states therefore supply the same free-band density and pressure at this order. No new gravity solution for a physically dephasing interacting history is implied.

The exercise specifies an observable that distinguishes a coherent pair-channel signal from a phase-randomized signal with the same measured energies. The phase-randomized state is a comparison state. Claiming a thermal state would require a population law, interactions and relaxation times as well.
