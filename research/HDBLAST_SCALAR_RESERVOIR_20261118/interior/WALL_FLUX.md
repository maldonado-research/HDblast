# Measuring the scalar energy flux at the material wall

In the selected five-dimensional B1 simulation, the scalar field carries energy **from the wall into the bulk** when integrated over $1\leq s=H_0\tau\leq4$. The particle-plus-radiation system on the wall gains energy during this interval. These signs are compatible because the wall's field-dependent tension releases a much larger amount. The fine-grid scalar gradient measured from five strictly interior nodes gives a total wall energy input of $-592.510291$. Subtracting the integrated tension work of $-592.625808$ gives $+0.115517$ for matter. The independently tabulated matter-side source gives $+0.115503$. These numbers apply to the specified simulation and clock.

## What the gradient measures

The solver evolves a metric with radial and temporal line element $e^{2B}(-dT^2+dz^2)$ and a fixed wall at $z=0$, approached from $z<0$. Let $a=\exp(A_b-A_b(0))$, $r_b=H_0^{-1}$, $v_h=r_b e^{-B_b}\partial_T\phi_b$ and $w=e^{-B_b}\partial_z\phi_b$. The outward unit normal to the simulated domain points toward increasing $z$. For a canonical scalar, the two reflected copies supply the following signed scalar energy flux into the wall per unit $s$ and reference comoving volume,
$$
 f_{\rm bulk}=-2r_b a^3 w v_h.
$$
The rate uses the same $\kappa_5^2/H_0$ surface-density scale as the matter account. The simulated scalar boundary condition requires $w=-(\sigma_{,\phi}+J)/2$. It predicts $f_{\rm junction}=r_ba^3(\sigma_{,\phi}+J)v_h$, with the two summands representing tension and matter. Scalar potential energy makes no separate contribution to this mixed normal-time stress component. [Bulk evolution and wall law](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py) · [Matter equations](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/shell_matter.py)

At ordinary output times the rerun saves the actual radial coordinates and evolved fields near the wall. For a degree-$p$ check we fit $p+1$ scalar samples on the physical interior, differentiate that polynomial at $z=0$, and multiply by the local inverse lapse to obtain $w$. This fit excludes the wall scalar value, ghost cells and the prescribed junction gradient. The replay also implements fits that include the wall scalar value. The underlying solver enforces its usual boundary condition throughout the evolution. The comparison checks how accurately the evolved field realizes that condition in a derivative constructed independently after the run.

## Sensitivity to radial resolution

| Integrated wall input | Wall spacing 0.001 | Wall spacing 0.0005 |
| --- | ---: | ---: |
| Interior linear fit | −601.375894 | −597.054272 |
| Interior quadratic fit | −592.361189 | −592.522521 |
| Interior cubic fit | −592.298712 | −592.509227 |
| Interior quartic fit | −592.307192 | −592.510291 |
| Junction total | −592.307252 | −592.510305 |
| Junction tension contribution | −592.422718 | −592.625808 |
| Junction matter contribution | +0.115466 | +0.115503 |

Every fit gives the same total direction. Resolving the smaller positive matter work requires higher accuracy. For example, subtracting the fine-grid tension from the linear-fit flux would give $-4.428464$. In the quartic fit, the signed discrepancy is 0.0522% and 0.0124% of the integrated matter work on the coarse and fine grids. The corresponding integrals of the *absolute* rate discrepancies are 0.359% and 0.0597% of the absolute matter-rate integral. In the first eighth of the window the matter source is zero while the fine-grid flux integral differs from the tension work by $1.73\times10^{-5}$. Such residuals cannot establish energy transfer at an arbitrary time. Two spacings and the original output cadence do not determine a continuum error bound. [Eight intervals and both stencil families](WALL_FLUX.json)

## Reproduction

The cell has $\phi_*=0.5$, $G=100$, $y=1$, $b=6.49\times10^{-6}$ and detuning $0.1$. The wrapper authenticates four published solver modules, builds the initial state afresh, and inserts a read-only output operation. It retains the original numerical parameters while stopping at coordinate $T=2.43$, just past the comparison window. Each grid supplied 244 output times. At each original output through the first sample at or beyond $s=4$, eleven archived wall columns, including both clocks, three background/field values, expansion, scalar velocity, radiation, particle density, matter force and tension derivative, were reproduced exactly. For integration the extraction retains all 153 interior sample times, adds the two interval endpoints by linear interpolation of each weighted rate, and integrates the resulting 155-point series by trapezoids. [Trajectory checks and source hashes](WALL_RERUN.json) · [Interior-gradient extraction](extract_wall_flux.py)

Replay the public CSVs with standard Python. For full regeneration, install NumPy and SciPy, use an HDblast checkout containing the four authenticated files and choose fresh output paths. The historical time-series archives are optional for rerunning the field evolution. The resulting flux concerns this simulated interface. Inferring a cosmological heat source or a later temperature history requires additional physics and evidence.

~~~sh
python -I -B wall_flux.py > /tmp/hdblast-wall-flux.json
cmp WALL_FLUX.json /tmp/hdblast-wall-flux.json

python -B run_instrumented.py /path/to/HDblast /tmp/b1-coarse coarse
python -I -B extract_wall_flux.py /tmp/b1-coarse coarse /tmp/b1-wall-coarse.csv
cmp wall_flux_coarse.csv /tmp/b1-wall-coarse.csv
python -B run_instrumented.py /path/to/HDblast /tmp/b1-fine fine
python -I -B extract_wall_flux.py /tmp/b1-fine fine /tmp/b1-wall-fine.csv
cmp wall_flux_fine.csv /tmp/b1-wall-fine.csv
~~~
