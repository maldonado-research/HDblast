# Does the emitted radiation control the expansion?

After the scalar-controlled interaction turns off, massless daughters should redshift like radiation. Gravity may still be governed by other degrees of freedom. We test both statements on the same four trajectories at common proper times $s=2.055$ through $2.10$. The latest field-clock switch-off in the baseline histories occurs at $s=2.050090$. We use the simulated positive strength $y^2=0.01$ to report the daughter abundance. This is a test value. [Protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7001073449)

On the fine quiet shell with the static calibration and pair cutoff 120, the hatted daughter energy density at readout is 0.11930. The total modeled matter density is 97.54677. Their ratio is **0.1223%**; the denominator includes the shell matter and omits the vacuum tension and gravitational terms. Across the tested arms, meshes, two cutoffs and two calibrations the fraction ranges from 0.1170% to 0.1231%. These densities use the study’s finite band and subtraction prescription. There is no measured temperature associated with this calculation. [Numerical results](RADIATION.json)

At leading order write $\rho_D=y^2\rho_{D,1}+O(y^4)$. The zeroth-order daughter density vanishes. After switch-off $a_0^4\rho_{D,1}$ should therefore be constant, with $a_0$ the zero-coupling scale factor. The maximum fractional drift on 81 points of the late interval is below $7.0\times10^{-11}$ over all tested cases. This checks the propagation of the modeled particles on the archived background.

For a spatially flat, expanding four-dimensional solution of the ordinary Friedmann equation containing only separately conserved radiation, $Q=a^2H/H_0$ is constant. In the fine quiet positive-coupling run $Q$ falls from 0.2785304 to 0.2697979, a **3.135% decrease**. For any constant $C$, the smaller achievable maximum of the two endpoint fractional errors $|C/Q_i-1|$ is $|Q_1-Q_2|/(Q_1+Q_2)$. This follows by choosing the harmonic mean of the two positive measurements so the two errors have equal magnitude and opposite sign. Even the best radiation-only normalization misses at least one of the two expansion rates by **1.593%** at its measured scale factor. This comparison applies to the stated interval, geometry and matter description. Its zero-coupling counterpart gives nearly the same discrepancy. [Other four-dimensional sources](../comparator/CANONICAL.md)

## Direct density term in the junction equation

In the native constraint the square involving the matter density is $r_b^2(\sigma+\rho_\mathrm{native})^2/36$. Holding the tension, existing density, scalar field, velocity and Weyl contribution fixed, the derivative of this square from the daughter density alone is

$$
D_1=\frac{r_b^2}{18}(\sigma+\rho_{\mathrm{native},0})\frac{b}{r_b}\,\widehat\rho_{D,1}.
$$

The hatted density converts to native units by $b/r_b$. In the fine quiet static model at cutoff 120, $D_1/h_0^2=+7.2241\times10^{-4}$, where $h_0$ is the baseline expansion rate in $H_0$ units at the same readout. The coupled derivative of $h^2/h_0^2$ is $-3.2297\times10^{-4}$. Their signs differ. The density square is one identifiable term in a constraint that also contains parent, scalar and gravitational response. Assigning it a fraction of the complete gravitational effect would require a distinct, dynamically consistent counterfactual. [Independent expansion and flux check](../audit/AUDIT.md)

The calculation in [radiation.py](radiation.py) reads each run on its own proper clock, reconstructs the scale factor with the measured expansion rate, and compares the signed-strength central difference to the archived zero-coupling history. The normalized scale factor is the same convention used by the native solver.
