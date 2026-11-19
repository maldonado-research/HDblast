# Does the five-dimensional solution support the measured change in expansion?

There are two checks on the expanding shell. The field in the interior supplies an estimate of work that can be compared to the matter work on the boundary. The scalar field, pressure, density and junction conditions predict how the expansion should change during the run. Both comparisons use the full recorded trajectories and a common proper-time readout $s=2.10$.

## A direct check of expansion from the trace equation

Write $h=H/H_0$, $v=d\phi_b/ds$, $r_b=1/H_0$ and $a=\exp(A_b-A_{\rm ref})$. The symbols $\rho,p,J$ below include the native conversion $b/r_b$ from hatted matter units. In conformal coordinates, eliminate $A_{ZZ}$ from the $A_{TT}$ evolution equation using the Hamiltonian constraint. Inserting the three junction gradients then gives

$$
\frac{1}{a^2}\frac{d(a^2h)}{ds}
= r_b^2\left[
\frac{(\sigma+\rho)(2\sigma-\rho-3p)}{36}
+\frac{U(\phi_b)}3-\frac{(\sigma_{,\phi}+J)^2}{24}\right]
-\frac{v^2}{6}.
$$

The check integrates the right-hand side from the common fork and divides by the final $a^2$, separately on each run. The initial term is $a_i^2h_i/a_f^2$. The current, stress, tension/potential and scalar kinetic terms are recorded separately. The observed endpoint $h_f$ enters only when comparing with the prediction. The check uses the measured scalar and scale-factor histories, so it tests a necessary gravitational evolution identity on the returned trajectory. It supplies an independent check on computing $h$ from the simulated metric velocities. [Analysis code](audit.py) · [Complete results](AUDIT.json)

On the fine quiet wall, with static calibration and comoving cutoff 120, the change in $h$ per $y^2$ is about $-2.85\times10^{-5}$. Of the integrated terms, the change in the explicit scalar-current term is about $-2.60\times10^{-5}$ and the term linear in the stress trace changes by about $-7.74\times10^{-6}$. The changing field, tension, volume and normalization also enter the sum. Treating either contribution as a complete energy source would omit those terms.

**Direct contribution of an isotropic massless daughter.** Its stress perturbation has $(\delta\rho,\delta p)=(x,x/3)$. Thus it cancels from the term proportional to $\sigma(\rho-3p)$. The separate first-order contribution from the quadratic stress term is exactly $-r_b^2(\rho+p)x/12$, where $\rho,p$ denote the unperturbed matter. We integrated that term on the quiet histories. At cutoff 120 it contributes about $-1.60\times10^{-9}$ per $y^2$ on the fine mesh, about $0.0056\%$ of the measured signed static slope. This percentage describes one term of this trace equation. The radiation density still enters the junctions and the reconstructed Friedmann equation. Responses of the parent, scalar current, field and metric include their interaction with the daughter sector. The trace equation also requires an initial expansion rate.

## Bulk scalar flux and boundary matter work

The bulk estimate differentiates quadratic and quartic polynomials fitted to the nearest three or five strictly interior values of $\phi(T,Z)$. Their derivatives are taken at the fixed wall $Z=0$. These fits omit the boundary value and the imposed boundary derivative. For the reflected pair of bulks with normal $n=e^{-B}\partial_Z$ toward the shell, the inferred current of scalar energy into the whole wall is

$$F_{\rm bulk}=-2\int a^3 e^{-B_b}\phi_T\phi_Z\,dT.$$

Subtracting the varying tension rate inside the integrand leaves the estimate of matter work. Independently, the shell gives $W_m=\int a^3 J\phi_T\,dT$. We compared both estimates with the native scalar-work accumulator and with the change in comoving energy plus pressure work. Each signed-strength history is first integrated to the same proper time. The centered derivative is computed afterward. Positive $W_m$ is work into the matter system in the specified orientation. Pressure work and tension work are distinguished throughout. The spatial check uses a different estimator on fields produced by the same boundary-value problem. [Earlier interface test](https://github.com/maldonado-research/HDblast/blob/a25c97eed48d64a0868b6120d7902158a320ce6a/research/HDBLAST_SOURCE_WORK_20261118/interface/INTERFACE.md) · [Integrated signed results](AUDIT.json)

For the quiet fine static case at cutoff 120 the change in matter work per $y^2$ is about $-1.14\times10^{-5}$ in native units per normalized comoving volume. The quartic interior measurement also resolves this negative change. Across the fine-grid primary strength pairs at both cutoffs and both original calibrations, its largest absolute discrepancy from the current integral is $1.12\times10^{-10}$ in the native flux units. The largest fine-mesh difference between the integrated trace and the measured expansion slope across the same primary pairs is $5.7\times10^{-10}$. It is the variation of the net scalar work on all the shell matter, which already receives work in the zero-coupling run.

## Why the stage at which a constraint is recorded matters

Let $C_H$ denote the signed Hamiltonian residual as returned by the native code, with zero as its continuum target. The exact semi-discrete trace uses actual discrete normal gradients and includes $-r_b^2e^{-2B_b}C_H/6$. The stored bulk rows sample the state at the beginning of each Runge–Kutta step. Their integrals are available in AUDIT.json under `stencil_correction`, `hamiltonian_correction` and `remainder`. A sampled correction can differ substantially from that contributing to the update.

An independent quiet zero-coupling replay recorded the actual four stages per time step. At every stage it checked the measured expansion acceleration against the sum of the junction prediction, the boundary gradient discrepancy and the raw signed Hamiltonian contribution. The fine-grid maximum algebra error was $5.2\times10^{-14}$ or less in $d(a^2h)/dT$. Up to the final full native step ($s\simeq2.10533$), the raw Hamiltonian contribution integrated using all four stages was $6.39\times10^{-9}$. Sampling only the first stage and integrating by the trapezoid rule instead gave $2.16\times10^{-5}$. Integrating the actual acceleration over the four stages reproduced the observed change in $a^2h$ to $7.2\times10^{-12}$. The full-step sum, the pointwise constraint and the continuum trace check answer different numerical questions. The final full-step proper times in STAGES.json differ slightly between grids. The paired results in AUDIT.json all use $s=2.10$. [Instrumented replay](stages.py) · [Two-grid stage diagnostics](STAGES.json)

## Interpretation and replay

These are finite-grid, finite-cutoff checks. The signed auxiliary strengths approximate derivatives at zero coupling. To this perturbative order the emission amplitude is proportional to $y$ on the original trajectory. A shift of that trajectory proportional to $y^2$ changes the amplitude at order $y^3$ and the emitted energy at order $y^4$. Detecting a leading change in geometry is a different task from calculating that higher-order change in yield. Inferring a pre-Big-Bang origin from an observed universe would require a calibrated physical stress and a discriminator against alternative four-dimensional histories.

With the native trajectories regenerated in bulk/runs, invoke the following from the study directory. The second command instruments two zero-coupling histories directly and can reuse the same trusted native geometry cache as the bulk replay.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB audit/audit.py bulk/runs --output /tmp/AUDIT.json
OPENBLAS_NUM_THREADS=1 python -IB audit/stages.py --bulk bulk --output /tmp/STAGES.json
~~~
