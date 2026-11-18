# Signed energy at the shell

The initial-area normalization is common to each pulse and quiet pair. Set $V=\exp[3(A_b(T)-A_b(0))]$, and count positive work into the shell. In the identity conformal chart the measured two-copy energy flux rate per initial comoving volume, multiplied by $\kappa_5^2$, is

$$P_{\rm bulk}=-2 V e^{-B_b}\phi_T\phi_Z.$$

On this preproduction branch the scalar junction is $\phi_Z=-e^{B_b}\sigma_\phi/2$, giving $P_\sigma=V\sigma_\phi\phi_T$. Its time integral can be checked independently by integration by parts. With $\Delta\sigma(T)=\sigma(\phi_b(T))-\sigma(\phi_b(0))$,

$$\int_0^T P_\sigma\,dt=V(T)\Delta\sigma(T)-\int_0^T 3 A_T V\Delta\sigma\,dt.$$

Subtracting the initial tension well before differencing avoids the cancellation of two nearly equal extensive terms. The code calculates its difference using the exact cubic tension polynomial. The bulk-side diagnostic fits a quadratic to three strictly interior samples and a quartic to five. Neither uses the wall value or the imposed Neumann derivative. Both include the measured wall time derivative and geometry from the same numerical evolution. All time integrals use the composite trapezoid. [Calculator](interface.py)

At the prespecified conformal time $T=0.15$ the primary positive pulse's full quartic bulk-side integral on the finest grid is −0.00034275510. The junction integral is −0.00034275510 and the pressure-corrected endpoint result is −0.00034275567. The quiet solver has a small numerical drift, with a signed integral of +0.00000113236 at that cut. The magnitude of the pulse's net integral is only 6.33% of the integral of the absolute rate, 0.00541853. The cumulative history in the figure subtracts the quiet curve at each common conformal time. The shell field crosses neither the production threshold nor a new matter state during the experiment. [Full four-grid report](INTERFACE.json)

For a comparison at equal elapsed shell proper time, the normalized clock is $s=\int_0^T e^{B_b(t)}dt/\rho_b$. Integrate each pulse up to the time at which its $s$ equals the quiet value at $T=0.15$, then subtract the quiet work at its own cut. Each common-clock cut used below falls before the direct metric-control characteristic. The code also reports paired differences at the positive pulse's own clock for all stored windows where the quiet history covers that clock.

| Scalar amplitude | Excess work, positive pulse | Excess work, negative pulse | Half-difference | Half-sum |
|---:|---:|---:|---:|---:|
| 0.001 | −0.00034394146 | +0.00033896915 | −0.00034145530 | −0.00000248615 |
| 0.0005 | −0.00017134746 | +0.00017010438 | −0.00017072592 | −0.00000062154 |

Doubling the amplitude multiplies the odd part by 2.00002 and the even part by 4.00001 on the finest grid. The next coarser ratios are 2.00002 and 3.99998. These finite comparisons are consistent with linear and quadratic leading terms over the tested amplitudes. The quiet scalar has a nonzero normal gradient. Its product with the perturbed wall velocity already gives a term linear in the pulse amplitude, so signed work can reverse under amplitude reversal. The half-sum includes the response of the geometry and tension and does not isolate a flux of newly created particles.

Endpoint closure at one discretization does not set a continuum error bound. For instance the positive-pulse net flux at $T=0.15$ changes by about $3.0\times10^{-6}$ between the two finest grids. The time derivative, wall gradient and bulk-side and junction routes all describe the same numerical solution. [Convergence and independent characteristic test](evolution/EVOLUTION.md)
