# How much work does the rolling field do on a wave?

We compare two states of a free massless spectator on the same prescribed Einstein–scalar spacetime. At the initial time one has a single extra quantum in a chosen wave and the other is the specified vacuum. Their energy difference measures the cost of the extra excitation after both states evolve. The rolling parent field can add or remove energy from it. The changing metric can also do work. The two contributions can have opposite signs.

The samples below cover 0.005 model-time units on each side of the first minimum passage of the parent with $\phi_m=-0.5$. The polar and equatorial passages occur at 0.48839541 and 0.58913424 on their respective backgrounds. Each row uses one initially occupied mode, one $\ell=2$ angular state and the archived matter weight $Z=1/f(\Theta)$. Energies are in inverse model-time units for $\hbar=1$. Positive work increases the daughter's energy relative to the comparison state.

| First contact | Initial radial label | Work from changing metric | Work from changing $Z$ | Change in energy difference |
|---|---:|---:|---:|---:|
| Poles | 2 | −0.00879772 | +0.00121334 | −0.00758438 |
| Poles | 4 | −0.01682502 | −0.00056256 | −0.01738758 |
| Poles | 8 | −0.02993796 | +0.00033974 | −0.02959822 |
| Equator | 2 | −0.00661354 | −0.00171735 | −0.00833089 |
| Equator | 4 | −0.00435660 | +0.01961020 | +0.01525361 |
| Equator | 8 | +0.01854098 | +0.02785102 | +0.04639201 |

For example, the rolling field added about 0.40% of the initial energy of wave 4 in the equatorial run during this window. Work by the metric reduced its net gain. During that same interval, on exactly the same geometry, the rolling field removed energy from wave 2. The mode labels describe global spatial patterns. The integral spans the full spatial slice, including regions that have yet to cross. These numbers depend on the preparation of the waves and the time window. They leave spontaneous production, heating, and a local release of energy specifically caused by crossing unresolved.

## Work at fixed canonical state

The field variable, angular basis, initial state preparation and geodesic proper time $x$ follow the [parent study](https://doi.org/10.5281/zenodo.24379598). The projected real action for one fixed angular block is

$$
S=\tfrac12\int dx\,[\dot X^T M\dot X-X^TKX],\qquad P=M\dot X .
$$

Write $Q_i$ and $P_i$ for the complex solutions initialized as positive-frequency mode $i$ of the initial Hamiltonian. Their labels use the increasing eigenfrequencies on the perturbed initial slice. Replacing the launch vacuum by the state with one quantum in mode $i$ increases the field covariance by $2\operatorname{Re}(Q_iQ_i^\dagger)$, and similarly for the momentum. Other initial modes cancel from this difference. The difference in the integrated Hamiltonian is

$$
E_i=P_i^\dagger M^{-1}P_i+Q_i^\dagger KQ_i .
$$

At the bounce this is the mode frequency. At each endpoint we also evaluate it with the instantaneous generalized eigenvectors as $\sum_j\omega_j(|\alpha_{ji}|^2+|\beta_{ji}|^2)$. State evolution implies that the derivative of the energy equals the expectation of the explicit Hamiltonian derivative. If $v_i=M^{-1}P_i$, each contribution to the power has the form

$$
\mathcal P_{a,i}=Q_i^\dagger \dot K_a Q_i-v_i^\dagger\dot M_a v_i,
\qquad a\in\{g,Z\}.
$$

The first partial derivative varies the spatial metric at fixed $Z$ and the second varies $Z$ at fixed metric. The same power can be evaluated as an integral of local stresses. Let $T$, $G_r$, and $G_t$ denote the excess kinetic, radial-gradient and angular-gradient energy densities, including the proper measure, so $E_i=\int(T+G_r+G_t)$. Write $h_r$ and $h_t$ for the radial and tangential expansion rates. The weak forms give

$$
\begin{aligned}
\mathcal P_{g,i}&=\int[-(h_r+2h_t)T+(2h_t-h_r)G_r+h_rG_t],\\
\mathcal P_{Z,i}&=\int \frac{\dot Z}{Z}(G_r+G_t-T),\qquad
\frac{\dot Z}{Z}=-\frac{f_{,\Theta}}{f}\dot\Theta .
\end{aligned}
$$

The relative kinetic and gradient contributions explain why a positive instantaneous mode count need not give a positive coupling work. In a coupled system the parent target-space equation receives the opposite explicit exchange. Here the waves propagate as probes and the simulated spacetime supplies the work. The calculation includes no backreaction of the probe on the parent or gravity.

The time integral of each power is evaluated independently of the endpoint energies. Over whole evolutions, half windows and complete crossing windows, the largest absolute discrepancy is $5.11\times10^{-14}$ for the six one-quantum comparisons and their constant-$Z$ controls. At first passage the code independently varies the metric and matter times in the archived Hamiltonian, at fixed canonical phase-space state, and checks each partial derivative by centered differences and Richardson extrapolation. Angular quadrature, background spatial resolution, field time step and time quadrature are also refined. The largest work or energy-change difference is $7.05\times10^{-10}$ when changing the maximum radial label from 24 to 36, and below $1.2\times10^{-14}$ from 36 to 48. For $Z=1$ the coupling contribution vanishes in every case.

## What instantaneous subtraction gives for the launch vacuum

We also evolve the entire specified finite-band vacuum. At each later time we subtract the ground-state energy of the *instantaneous* finite Hamiltonian. This measures

$$
E_{\rm inst}=\langle H\rangle_{0,\,\mathrm{evolved}}-
             \tfrac12\sum_j\omega_j
             =\sum_{ji}\omega_j|\beta_{ji}|^2.
$$

The derivative of the instantaneous ground energy in each metric or coupling direction is the ground-state expectation of that partial Hamiltonian derivative. Subtracting that power from the evolved-vacuum power gives a work decomposition whose sum integrates to $\Delta E_{\rm inst}$. We check both the mode sum and the integrated balances. This subtraction depends on the choice of canonical field and the slicing. Calculating a renormalized stress tensor requires a prescription for its ultraviolet singularities that respects the stress-energy identities.

| Maximum even radial label | Equatorial bracket metric work | Equatorial bracket coupling work | Change in $E_{\rm inst}$ |
|---:|---:|---:|---:|
| 12 | +0.018814 | +0.049814 | +0.068628 |
| 24 | +0.039412 | +0.093687 | +0.133099 |
| 36 | +0.049515 | +0.116216 | +0.165731 |
| 48 | +0.056519 | +0.131924 | +0.188442 |

At a *fixed* maximum label of 36, simultaneous background, field-time, angular-quadrature and time-quadrature refinement changes the vacuum works and energy changes by less than $2.3\times10^{-12}$. Adding higher radial labels still changes the inferred energy substantially. The totals include one fixed angular state. Other angular sectors and their multiplicities are outside this sum.

### Explicit ultraviolet control in closed de Sitter space

The homogeneous expanding control admits a closed-form result even with constant $Z$. Let $h$ be the de Sitter expansion scale, $a=\cosh(hx)/h$, $H=\dot a/a$, $s=\sinh(hx)=aH$, $\eta=\arctan s$, and $L=n+1$. The exact mode solution with the same instantaneous bounce vacuum gives the occupation

$$
N_n=\tfrac14\left[
 \frac{(\sin(L\eta)/L-s\cos(L\eta))^2}{L^2-1}
 +\left(\frac{s\sin(L\eta)}{L}\right)^2\right].
$$

At any fixed $x>0$ we have $N_n=s^2/(4L^2)+O(L^{-3})$ and $\omega_n=\sqrt{L^2-1}/a$. In one fixed $\ell=2,m$ block with even $n$, the leading term in the energy sum is proportional to the harmonic series. Consequently

$$
\sum_{\substack{2\le n\le n_{\max}\\n\ \mathrm{even}}}\omega_n N_n
       = \frac{(aH)^2}{8a}\log n_{\max}+O(1).
$$

At the bounce the occupations vanish. At positive time the logarithm occurs even though the matter coefficient is constant. This analytic control demonstrates an explicit divergence in the instantaneous excitation-energy definition. It does not establish the ultraviolet coefficient for the inhomogeneous backgrounds, where different radial modes can mix. The code checks the formula against the exact field and momentum solution and computes partial sums to radial label 524288. Comparing two states that differ by one extra initial quantum removes this common vacuum contribution.

### Does the launch state itself have suitable short-wavelength behavior?

A distinct test compares two *initial states* at the bounce, using the same original-field Hamiltonian for both. Keep the exact homogeneous closed de Sitter background and restrict to the nonzero spatial harmonics $n\ge1$. Write $L=n+1$, $\lambda=L^2-1$, $a=(h\cos\eta)^{-1}$ and $d x=a\,d\eta$. For a real scalar with constant kinetic normalization, choose the following explicit reference solution for each spatial oscillator:

$$
u_n(\eta)=\frac{h}{\sqrt{2L\lambda}}
   (L\cos\eta+i\sin\eta)e^{-iL\eta}.
$$

The rescaled wave $a u_n$ solves $(a u_n)''+[L^2-2\sec^2\eta](a u_n)=0$. It has the positive-frequency short-wavelength branch. Its exact phase frequency and squared amplitude in conformal time are $W=L\lambda/(L^2+\tan^2\eta)$ and $1/(2W)$ respectively, with $W\sim L$ away from the conformal-time endpoints. The reference concerns these smooth, normalized nonzero modes. The homogeneous massless mode requires its own state. This comparison does not rely on a global de Sitter-invariant vacuum for that degree of freedom.

At the bounce $M=a^3=h^{-3}$ and $\omega=h\sqrt\lambda$. The reference data $Q=u_n(0)$, $P=a^2u'_n(0)$ obey $i(Q^*P-P^*Q)=1$ and $\dot Q/Q=-ih\lambda/L$. The instantaneous launch state uses $\dot Q/Q=-ih\sqrt\lambda$. Resolving the reference data in the instantaneous basis gives a per-oscillator squared mixing amplitude

$$
|\beta_n|^2=
 \frac{1}{4L\sqrt\lambda(L+\sqrt\lambda)^2}
 =\frac{1}{16L^4}+O(L^{-6}).
$$

The energy of the reference ground state in the launch Hamiltonian exceeds the instantaneous ground energy by $\omega|\beta_n|^2$ per real oscillator. With the full multiplicity $L^2$ on the three-sphere, each shell contributes $hL/[4(L+\sqrt\lambda)^2]=h/(16L)+O(L^{-3})$. Thus

$$
E_{\rm ref}(0)-E_{\rm inst.ground}(0)=\frac{h}{16}\log n_{\max}+O(1).
$$

Dividing by $2\pi^2 a(0)^3$ gives the energy-density coefficient $h^4/(32\pi^2)$. The code verifies the mode normalization and the direct Hamiltonian difference and recovers the logarithmic coefficient to 1.1 parts per million from the sum over labels 1,000,001 through 2,000,000. A subtraction common to both states cancels out of this difference.

The time-evolved energy above an instantaneous reference in one angular block and this initial full-sphere state comparison answer separate questions. Both diagnose limitations of extending the archived vacuum convention to arbitrarily short wavelengths. On this homogeneous control, choosing the instantaneous data for finitely many low harmonics and matching the explicit reference for higher harmonics gives a concrete finite-mode change of state. Extending the same idea to the distorted HDblast slices requires a short-distance state construction for their own spatial metric and kinetic factor.

From the study root, reproduce the complete calculation with

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B work/work.py --prior prior --output /tmp/WORK.json
~~~

[Calculation](work.py) · [Results, signed work and refinement checks](WORK.json)
