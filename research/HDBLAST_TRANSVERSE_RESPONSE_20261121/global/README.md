# A density requirement for fields made by the particles

How much spatial structure would the newborn particles need to create their own sideways fields? This test permits arbitrary spatial wavelengths and orientations. The particles may concentrate, and increased target density also increases the allowed binary collision rate. The upper limit on both effects is set by one quantity, $C$, the largest local number density divided by the spatial mean for the same charge sign and flavor over the interval. The phase-space distributions can otherwise vary freely.

With initially zero perturbing fields, approaching the reference thermal transverse momentum distribution to within one percentage point at **0.1% of the whole pulse's absolute electrical work** requires $C$ above **15** for the dedicated interaction and above **14** for the shared interaction. By the 1% work marker the corresponding lower bounds have fallen to **1.43** and **1.29**. A pointwise lower bound says that at least one position, time and population would have to reach that contrast. It leaves the spatial volume involved undetermined. The calculated thresholds are necessary conditions. The 5% work marker gives no exclusion at $C=1$ by this test. [Four numerical histories](CLUMPING.json)

| Work completed | Dedicated cap excluded, native / halved step | Shared cap excluded, native / halved step |
|---|---:|---:|
| 0.1% | 15.31 / 15.46 | 14.53 / 14.71 |
| 1% | 1.437 / 1.442 | 1.295 / 1.300 |

Here and in the table the target is a Kolmogorov distance of at most 0.01, meaning all cumulative transverse-momentum fractions agree to within one percentage point. Entries are rounded numerical boundaries using conservative endpoint bounds. The separate force requirement for a continuously emitting source resolves a much smaller tail of very recent births. [Cohort and force tests](../history/RESULTS.md)

## Bounding the field energy before knowing its spectrum

Take a periodic locally flat patch on the reference expanding shell with homogeneous gauge coefficient $K=2.54$. Use a spatial mean to define the field norm $\|\Psi\|_2^2=\langle|\delta E|^2+|\delta B|^2\rangle$, with $\Psi$ including spatial variations of the axial field. An imposed uniform parallel field and its parallel guide are accounted for in the background. Assume the inhomogeneous Maxwell evolution is energy conserving in conformal variables, or dissipative, between sources. The calculation below includes the daughter current as its sole new source. Transverse radiation from the bulk and amplification by additional wall degrees of freedom would require source terms of their own.

The effective current normalization for one of the three symmetric gauge sectors is the one used in the archived wall equation, $Q D/a^3$ if all charges stream maximally, with $Q=10^7$ and $D$ the tracked comoving daughter rest. Positive and negative charges are included. The same convention gives the longitudinal electric work in the parent solver. If the local density of each sign is bounded by $C$ times its mean, the normalized sum of both sign densities has mean one and second moment at most $C$. All particle speeds are at most one. Hence $\|J\|_2\le Q D\sqrt C/a^3$. Removing the prescribed uniform current is an orthogonal projection when calculating nonzero modes. Maxwell's energy estimate gives

$$
\|a(t)^2\Psi(t)\|_2 \le a_0^2\|\Psi(0)\|_2
   + \frac{Q\sqrt C}{K}\int_{t_0}^t\frac{D(s)}{a(s)}\,ds.
$$

This uses the integral in proper wall time and the comoving Maxwell source $a^2 J/K$. The bound holds for curved wavefronts represented as superpositions in the local patch. There is no summation penalty for admitting additional Fourier modes. The full five-dimensional dynamics of spatial variations on the wall have yet to be evolved.

## Bounding what the particles encounter

Take a final time $T$ and use $D(T)$ to normalize probabilities over all particles then present. Consider only surviving unhit particles born after a candidate time $b$. To grow the magnitude of their transverse momentum, additional forces must supply a redshift-weighted absolute radial impulse. Uniform axial fields, including the magnetic guide, preserve the transverse comoving magnitude. The extra electric and magnetic fields produce radial acceleration of magnitude at most $\sqrt{2}Q|\Psi|$ in mass units. The factor becomes one for a common one-dimensional perpendicular wavevector in the neutral sector because the relevant electric and magnetic deflections are orthogonal. The general result in the table uses $\sqrt{2}$.

At time $t$, the specified particles are a subset of those then present. Weighting the field by their positions can increase its RMS above the spatial RMS. Its RMS with respect to the final number measure is at most $\sqrt{C D(t)/D(T)}\|\Psi(t)\|_2$. The triangle inequality for RMS values integrated over time gives the following upper bound on their radial impulse, normalized using all particles present at $T$:

$$
J_{\rm max}(b,T)=\frac{\sqrt{2}Q\sqrt C}{a(T)}
\int_b^T\frac{dt}{a(t)}\sqrt{\frac{D(t)}{D(T)}}
\left[a_0^2\|\Psi(0)\|_2+\frac{Q\sqrt C}{K}\int_{t_0}^t\frac{D(s)}{a(s)}ds\right].
$$

This estimate allows the field to be strongest where there are many particles. It uses no independence of the particle positions and fields.

## Allowing for collisions in denser regions

Let $h_1(t)$ be the upper bound on the binary hazard computed with the mean opposite-sign target density and the maximum of the prescribed collision kernel over all relative momenta. Even moving through overdense regions, a particle sees hazard at most $Ch_1(t)$ under the local density cap and the binary Boltzmann law. Thus the guaranteed unhit share from births after $b$ is at least

$$
f(b,C)=\frac{1}{D(T)}\int_b^T\exp\left[-C\int_s^T h_1(t)dt\right]\,dD(s).
$$

A transverse thermal CDF $F_T$ within discrepancy $\epsilon$ permits the lowest unhit particles to have even smaller momenta than the exact thermal quantiles. Grant every scattered particle unlimited redistribution. With initial transverse momentum bounded by $p_*=\sqrt{1.05^2-1}$, a necessary RMS impulse is still

$$
J_{\rm min}(b,T)=\left[\int_0^{[f(b,C)-\epsilon]_+}
   ([F_T^{-1}(u)-p_*]_+)^2 du\right]^{1/2}.
$$

All momenta and impulses here use units of $mc$. If $J_{\rm min}>J_{\rm max}$ for any cohort cut, the assumed cap cannot reach the target. The code scans native birth knots and increases $C$ until this contradiction disappears. For fixed $C$ the same calculation gives a minimum initial RMS field if a seed is allowed. For instance, at the 1% marker under $C=1$ on the refined histories the necessary initial strengths are at least 0.00127 and 0.00090 in wall units. These permissive initial-energy floors grant subsequent particle currents their maximal effect.

The numerical upper integrals use upper endpoint masses, lower endpoint radii and the entire cell containing the requested marker. The lower number fractions discard all births in that last cell. Quantile costs use lower rectangles with two grid resolutions. Direct integration at each limiting witness provides a separate cross-check. The estimates assume the measured monotonic expansion and birth histories and the adopted interpolation of the recorded target temperature. [Implementation](clumping.py) · [Single-harmonic current bound](../causal/CURRENT_CAP.md)
