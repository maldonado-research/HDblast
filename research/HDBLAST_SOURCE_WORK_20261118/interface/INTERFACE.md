# Which way does scalar energy cross the shell?

Over the rolling-shell interval studied, the **net scalar energy flux across the two copies of the bulk points out of the shell**. Over the same interval the changing scalar mass does positive work on the modeled shell matter, including the damped field and its massless reservoir. A much larger decline of the tension along the trajectory accounts for the difference in the integrated junction balance. Resolving the small matter term requires subtracting the tension and resolving rapid oscillations in the boundary gradient. [Flux and signed work](INTERFACE.json)

## Geometry, normalization and independent measurement

The native conformal chart has transverse metric $e^{2B}(-dT^2+dZ^2)$, planar scale factor $e^A$, and the shell sits at the right endpoint $Z=0$. The chart transformation used for the static reference preserves this conformal metric form. All postfork records in this check use full fields. On one retained copy the future-directed wall tangent is $u=e^{-B}\partial_T$ and the normal toward the wall is $n=e^{-B}\partial_Z$. The scalar energy current along this normal is $-T_{ab}u^a n^b=-e^{-2B}\phi_T\phi_Z$. The potential term drops out because $u\cdot n=0$. The two reflected copies supply the factor of two. With $a=e^{A_b-A_{*,0}}$ normalized by the original quiet crossing, physical wall time $d\tau=e^{B_b}dT$, and $s=H_0\tau$, the signed cumulative scalar flux into the wall per normalized comoving volume in the solver's $\kappa_5^2$ units is

$$
F=-2\int a^3 e^{-B_b}\phi_T\phi_Z\,dT .
$$

Here the friction coefficient is $Y=0$. The scalar junction prescribes $\phi_Z=-e^{B_b}(\sigma_{,\phi}+J)/2$, with $J=(b/r_b)\hat j$, $b=6.49\times10^{-6}$, and $r_b=7.83828554$. The junction predicts $F=W_\sigma+W_m$, where $W_\sigma=\int a^3\sigma_{,\phi}\,d\phi$ and $W_m=(b/r_b)\int a^3\hat j\,d\phi$. The tension contribution includes its pressure work when the comoving volume grows. The matter source work equals its comoving energy increase plus its pressure work. Decay transfers energy between the massive and radiation components. [B1 equations and signs](https://github.com/maldonado-research/HDblast/blob/b5be2cb3c044457e20a9251658fd7e36f5e65864/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py) · [Earlier pulse flux convention](https://doi.org/10.5281/zenodo.24354840)

The bulk-side measurement extrapolates $\phi_Z$ from three and five **strictly interior scalar field values** by differentiating their quadratic and quartic Lagrange interpolants at zero. It reads each actual $Z_i$ from the native grid. Neither fit uses the shell field value, a ghost value or the imposed Neumann derivative. The shell $B_b,A_b,\phi_T$ in the common prefactor come from the same evolution. To avoid subtracting large integrated quantities, the code subtracts the tension rate point by point before integrating the small inferred matter rate. The independent matter estimate uses the change in the scalar-work accumulator advanced in the native RK4 system. The entire evolution uses the junction law, so agreement checks the spatially resolved numerical solution and energy normalization. It provides no independent physical validation of this boundary model.

## Measured work from the fork to the common readout

All arms fork at $T=1.11$ and end this comparison at the same proper time $s=2.10$. The fine-grid fork is at $s=1.30371715$. The quiet arm reaches the final proper cut at $T=1.58416457$. This is the 24-node, $y=1$ reciprocal model with the specified instantaneous subtraction, not the stable free field evaluated on the historical prescribed trajectory. Its initial preparation is at $T=0.70$. The integrals here begin at the fork to isolate the subsequent shared-clock comparison.

For quiet, the fine-grid quartic estimate gives $F=-0.9532254634$. The tension term is $W_\sigma=-0.9533554074$. The interior residual gives $W_m=+0.0001299439$ and the native scalar-work increment gives $+0.0001299477$. Of this latter work, $+0.0001275099$ goes into the change in comoving matter energy and $+0.00000243786$ into pressure work. Matter therefore receives about **0.01363% of the magnitude of the tension work** in these units. All four arms have negative net flux and positive total matter work over this interval. A positive number in the matter-work column alone would give the wrong orientation for the total bulk flux.

| Pulse relative to quiet | Fine-grid change in matter work, from the interior flux | Fine-grid change from the native matter source | Fractional discrepancy of that change | Change in total flux, fine grid |
|---|---:|---:|---:|---:|
| Near positive | $-2.48654\times10^{-7}$ | $-2.48482\times10^{-7}$ | 0.0692% | $+1.213\times10^{-4}$ |
| Near negative | $+1.25725\times10^{-6}$ | $+1.25739\times10^{-6}$ | 0.0108% | $+1.621\times10^{-4}$ |
| Distant compensation, positive front | $+7.67601\times10^{-8}$ | $+7.69207\times10^{-8}$ | 0.209% | $+2.098\times10^{-4}$ |

Positive total-flux differences mean less outgoing flux than in quiet. They include the tension response. The near positive front reduces the matter work even as the magnitude of the total outward flux falls. Two-grid quartic estimates give all three signs, with relative discrepancies of 0.420%, 0.0635% and 1.89% on the coarse mesh. The independent three-point gradient gives all three signs at both resolutions, but overestimates the smallest coarse-grid excess by roughly a factor of four. The two grids do not certify a continuum error bound.

## Why recording cadence matters

The boundary evolution imposes a derivative on the boundary stencil. Strictly interior extrapolation becomes noisy as the pulse crosses the last few grid points. The instantaneous fine-grid discrepancy in the pulse's flux rate can reach $0.00218$ in these units. We saved every RK4 time step, $\Delta T=0.00025$ on the fine grid and 0.0005 on the coarse grid. The two fine-grid quadratures of each signed quartic excess then agree to within $5.1\times10^{-13}$. If every fourth fine-grid record is retained and Simpson weights are used, the apparent near-positive excess is $+1.47\times10^{-6}$ and its sign is wrong. Alternating quadrature weights applied at the sparser cadence magnify oscillations. The quartic estimate with full cadence resolves a small integrated transfer. Instantaneous residuals and coarsely sampled outputs do not have the same accuracy. [Primitive samples for eight runs](SAMPLES.zip)

The recorded $T,s,A_b,\phi_b,h,\hat j$, integrated scalar work, massive energy and bath energy agree with the original coupled dataset at each of its samples. The script verifies these values before computing the work. Both available meshes are simulations with the same cutoff, effective bath, subtraction and geometry. Energy transfer after the fork is conditional on that construction. The globally prepared pulse energy, any prehistory of the shell, the full quantum stress and the origin of our Universe remain separate questions.

Run the standard-library replay from this directory. To regenerate the inside-wall histories in a temporary directory use the second command with NumPy and SciPy. The source repository must include the four historical B1 modules, the matched-data builder and each mesh’s historical quiet summary. Its four historical evolution modules and initial-data builder are hash-checked. The recording code augments the diagnostic row and samples each native time step. Each regeneration uses a new temporary working directory.

~~~sh
python -I -B interface.py
python -I -B run_interface.py /path/to/source/repo /tmp/interface-samples.zip --jobs 2
python -I -B interface.py /tmp/interface-samples.zip --walls ../coupled/WALLS.zip
~~~
