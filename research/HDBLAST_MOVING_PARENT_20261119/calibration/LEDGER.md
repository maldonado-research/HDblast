# Moving parents, switching work and the endpoint Hamiltonian

## Exact momentum balance and a flat-space energy check

Prepare one extra traveling scalar quantum with momentum $\mathbf q$ in a normalization box, compared to a vacuum run with the same real spatially homogeneous switch $f$. At leading order in a Yukawa interaction $y f\chi\bar\psi\psi$, annihilating the tagged scalar gives the massless Dirac pair momentum $+\mathbf q$. Stimulated creation of a second scalar in the occupied traveling mode gives the pair momentum $-\mathbf q$. The squared creation matrix element has weight two in the tagged run and one in the reference. Pair production from other vacuum modes cancels in the difference. If the excess branch weights are $P_a$ and $P_c$, the inclusive pair count is $P_a+P_c$ and the excess daughter comoving momentum is $\mathbf q(P_a-P_c)$. The parent's change in comoving momentum is the negative of this result. Homogeneous switching supplies no spatial impulse. This momentum statement uses spatial translation symmetry and also holds on the prescribed homogeneous shell. For an isotropic parent distribution the vector sum is zero.

In Minkowski spacetime with constant $m>0$, let $E=\sqrt{q^2+m^2}$ and $F(\nu)=\int dt\,f(t)e^{i\nu t}$. In a massless two-body final state of total energy $\omega$ and momentum magnitude $q$, the Dirac spin sum is $2(\omega^2-q^2)$, and Lorentz-invariant two-body phase space has support $\omega\ge q$. The uniform-parent normalization contributes $1/(2E)$. With the standard rest width $\Gamma_0=y^2m/(8\pi)$ the differential excess probabilities are

$$
dP_{a,c}=\frac{y^2}{16\pi^2 E}(\omega^2-q^2)\,
             |F(\omega\mp E)|^2\,d\omega,\qquad \omega\ge q .
$$

The daughter energy is $D=\int\omega(dP_a+dP_c)$. The excess energy removed from the parent sector is $E(P_a-P_c)$ since the incoming mode diagonalizes the free Hamiltonian at both endpoints. With the interaction off at both endpoints the net external switching work is therefore $W=D-E(P_a-P_c)$. The same quantity follows directly from the explicit variation of the interaction. For either transition with signed energy gap $\Delta=\omega\mp E$, the first-order commutator gives its work kernel
$$
-2\int dt\,\dot f(t)\!\int_{-\infty}^t dt'\,f(t')\sin[\Delta(t-t')]=\Delta\,|F(\Delta)|^2 .
$$
Integration by parts and symmetrization in $t,t'$ give the equality. Multiplying by the transition spectral measure and integrating both channels recovers the stated energy balance.

**Sign of the switching work.** Write the annihilation detuning as $\nu=\omega-E$, and $d=E-q>0$. For $0<\nu<d$ both detunings $\pm\nu$ occur. Since $f$ is real, they have equal squared Fourier weights. Their summed work integrand, after dropping the common positive prefactor, is

$$
\nu\big[(E+\nu)^2-q^2-(E-\nu)^2+q^2\big]|F(\nu)|^2
       =4E\nu^2|F(\nu)|^2 .
$$

For $\nu>d$ only positive annihilation detuning contributes. The stimulated-work integrand $(\omega+E)(\omega^2-q^2)|F(\omega+E)|^2$ is nonnegative throughout its domain. Thus $W\ge0$ for any real switch with convergent indicated integrals. An interaction switched on forever at constant strength lies outside this finite-switch derivation.

**Light-cone correction to the time-local pair count.** Folding both branches onto the signed frequency axis and applying Parseval yields, for $C=y^2/(16\pi^2E)$,

$$
P_a+P_c=
\frac{y^2}{8\pi E}\left[m^2\int f^2\,dt+\int\dot f^2\,dt\right]
+C\int_{-q}^{q}(q^2-\omega^2)|F(\omega-E)|^2\,d\omega .
$$

The final term restores the interval where the polynomial continuation is negative but a physical pair is kinematically forbidden. It vanishes identically at zero parent momentum. A smoothly stretched switch recovers the time-dilated golden-rule rate $\Gamma_0m/E$ at fixed $q/m$. To use that approximation, its Fourier width should be small compared with the distance $E-q=m^2/(E+q)$ from the annihilation resonance to the pair threshold. The narrow-width approximation can fail for fast parents even when $mT$ is large.

For $f(t)=\exp[-t^2/(2T^2)]$, define $L=(\Gamma_0m/E)\sqrt{\pi}T$. The table uses $L$ for counts and $EL$ for energies. All entries are leading coefficients, which require weak coupling and small integrated transition probabilities for use as absolute yields.

| $q/m$ | $mT$ | $T(E-q)$ | Pair count | Signed depletion | Daughter energy | Switching work |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1.000 | 1.500000 | 1.471605 | 2.521859 | 1.050255 |
| 0 | 16 | 16.000 | 1.001953 | 1.001953 | 1.005859 | 0.003906 |
| 4 | 4 | 0.492 | 1.230082 | 1.230082 | 1.278227 | 0.048146 |
| 16 | 16 | 0.500 | 1.201554 | 1.201554 | 1.204527 | 0.002972 |
| 16 | 64 | 1.998 | 1.000369 | 1.000369 | 1.000612 | 0.000244 |

The calculation evaluates the shifted Gaussian moments analytically, then integrates the positive spectra, paired-detuning work and spacelike correction directly at two resolutions. A log-scaled integration resolves the stimulated branch when its probability underflows ordinary floating point. The equal displayed counts and depletions at larger boosts reflect rounding. The sixteen cases in the report span $q/m=0,1,4,16$ and $mT=1,4,16,64$. [All flat controls](FLAT.json) · [Earlier rest calculation](https://github.com/maldonado-research/HDblast/tree/0c79f123c00918f16647df9c868b4ca889336dc5/research/HDBLAST_TAGGED_PARENT_20261119/flat)

## Can the same energy subtraction be used on the shell?

The expanding background has no time-translation symmetry. In addition to switching work it changes the scale factor and scalar mass. Even the free physical Hamiltonian at the readout contains terms that the two branch probabilities alone cannot determine. For one propagated normalized mode $\chi_q$ let $p_q=a^3\dot\chi_q$ and $E_q=\sqrt{q^2/a^2+m^2}$. A dot denotes the dimensionless shell proper-time derivative. The readout coefficients are

$$
A_q=|p_q|^2/a^3+a^3E_q^2|\chi_q|^2,\qquad
B_q=p_q^2/a^3+a^3E_q^2\chi_q^2,\qquad
A_q^2-|B_q|^2=E_q^2 .
$$

For a nonzero pair of opposite momenta, the free physical energy expressed using the initial-mode operators is $A_q(N_{\bf q}+N_{-\bf q}+1)+B_q b_{\bf q}b_{-\bf q}+{\rm h.c.}$, where the conjugate applies to the last term alone. At zero momentum use $A_0(N_0+1/2)+(B_0 b_0^2+{\rm h.c.})/2$. The identity follows from the unit Wronskian. The initial physical vacuum makes $B_q=0$ at preparation. Propagation across the time-dependent mass changes these coefficients.

At $s=2.10$ on the fine quiet wall, $A_q/E_q$ is 2.97852 at $q=0$ and 1.38759 at $q=8$. The corresponding $|B_q|/A_q$ is 0.94196 and 0.69327. At $q=8$ the energy quantum is $51.37836H_0$, the diagonal coefficient is $71.29212H_0$, and the anomalous coefficient has magnitude $49.42500H_0$. The independent mode calculation halves its maximum time step on four arms of each bulk grid. The largest shift in any energy coefficient is below $2.35\times10^{-6}H_0$, with maximum scaled Hamiltonian-identity error below $5.31\times10^{-12}$. [Readout coefficients and checks](ENDPOINT.json)

These coefficients establish that the endpoint energy is sensitive to pair coherence in the transported basis. Measuring the interaction-induced change in that coherence requires another calculation. In perturbation theory, the off-diagonal energy operator can read virtual, zero-daughter amplitudes of order $y^2$ through interference with the initial state. A closed cosmological energy budget thus needs the parent two-point function through that order, its consistent subtraction, and the switching, mass and gravitational work. The directly calculated daughter energy and the change in parent momentum can already be tested without making that identification.

Run the flat check with the Python standard library. The endpoint check requires NumPy and SciPy and authenticates the wall archive and both clock normalizations before evolving.

~~~sh
python -I -B flat_moving.py
OPENBLAS_NUM_THREADS=1 python -I -B endpoint.py /path/to/HDblast
~~~
