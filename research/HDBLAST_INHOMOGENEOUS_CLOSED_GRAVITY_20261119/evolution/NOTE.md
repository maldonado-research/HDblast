# Coupled local Einstein--scalar calculation

`evolve.py` evolves the constrained data on the **whole closed three-sphere**, with Gaussian proper time `x`, unit lapse and zero shift. It imports `../initial/constraint.py`; its `--previous` path must name the preparation package containing `nonlinear/nonlinear.py` and `prior/`. Reproduce all results with:

```sh
python evolution/evolve.py --previous previous --output evolution/results.json
```

The potential is the same fitted potential as before, including the same consistently removed tiny stationary-point slope. No quadratic truncation of the scalar force is used. Initial data are `q=b+A U_2(cos chi)`, zero scalar velocity and zero extrinsic curvature, with the conformal factor obtained from the nonlinear Hamiltonian constraint. The standard family has `b=0`; a second family chooses `b` to make the **proper-volume** initial mean zero.

## Formulation and regularity

The metric is

`ds²=-dx²+exp(2L)dchi²+exp(2M)sin²(chi)dOmega²`,

with `L=M+(1-y²)D`, `y=cos chi`. Thus the two radial/angular scale factors coincide smoothly at each pole; their difference vanishes quadratically. The code analytically factors apparent divisions by `sin² chi` and evolves smooth, reflection-even Chebyshev-U series. This retains the reflection symmetry of the chosen zonal n=2 initial harmonic without changing the S³ topology or volume measure. It subtracts the exact closed de Sitter background from metric variables, making the stationary control exact.

With `alpha=Ldot`, `beta=Mdot`, `theta=alpha+2 beta`, the equations are

- `alpha_dot=-Ric_rad-theta alpha+exp(-2L)q_chi²+3u`,
- `beta_dot=-Ric_ang-theta beta+3u`,
- `q_dot=Pi`, `Pi_dot=Delta_g q-theta Pi-3u_q`.

The independently monitored constraints are `R3+theta²-alpha²-2beta²-Pi²-|grad q|²-6u=0` and the radial momentum constraint. A separate trace Raychaudhuri residual is also stored. In particular, the observable volume radius is `aV=(V/2pi²)^(1/3)`, with

`aV_dot/aV=<theta>/3`,

`aV_ddot/aV=<u-Pi²/3-2(alpha-beta)²/9>+2 Var(theta)/9`.

Brackets always use the evolving proper volume. The script records the proper mean **change** from its actual initial mean; a zero round-sphere mean is not silently identified with a zero proper mean.

## Results at equal proper time x=0.12

The table averages the `+A` and `-A` mean changes. This removes the leading signed contamination from the initial proper-volume weighting.

| Branch phi_m | A | Coupled even mean change | Fixed-geometry even mean change |
|---|---:|---:|---:|
| +0.5 | 0.001 | +1.57458025e-5 | +1.57472391e-5 |
| +0.5 | 0.0005 | +3.93644577e-6 | +3.93653555e-6 |
| -0.5 | 0.001 | -1.09253914e-6 | -1.09260939e-6 |
| -0.5 | 0.0005 | -2.73142417e-7 | -2.73146808e-7 |

The even change divided by A² converges to approximately `+15.7458` and `-1.09257`. Gravitational coupling changes these fixed-geometry comparisons by less than 0.01% at this time. The proper-mean-centered A=0.001 pair gives `+1.57464147e-5` and `-1.09257285e-6`; its initial proper means are below `1.3e-19`. Thus the induced mean survives removal of the small initial proper-mean offset.

For the standard A=0.001 pairs, the even acceleration difference `aV_ddot/aV-h²` is `-7.17882167e-4` and `-2.55547480e-4`; the even volume expansion-rate difference is `-2.41100304e-5` and `-1.02384395e-5`. The scalar spreads are approximately `0.00279919` and `0.00211209`, much larger than their means. A larger volume than the reference background at x=0.12 reflects the different constraint initial slice; its expansion rate and acceleration are lower. These computations support a nonzero induced mean, and do not establish homogenization or faster volume expansion.

## Controlled extension and numerical checks

The optional extended runs stop at the earlier study's peak-field caps, 0.04 (+ branch) and 0.08 (- branch). All four signed amplitudes reach those caps:

| Branch | +0.001 | -0.001 | +0.0005 | -0.0005 |
|---|---:|---:|---:|---:|
| +0.5 stop x | 0.229040 | 0.240206 | 0.279427 | 0.290380 |
| -0.5 stop x | 0.352728 | 0.349075 | 0.415247 | 0.411599 |

Endpoint times differ, so their means must not be combined into equal-time signed averages. `controlled_extensions` also records all runs at common branch times 0.2290400646 and 0.3490750466, plus each run's fixed-geometry comparison at its exact same proper time. The main `at_probe` remains x=0.12.

Primary degree 32 / 128 quadrature-node runs are compared with degrees 16 and 24 at x=0.12 and degree 24 at the caps. Dense proper-volume diagnostics and constraints use independent higher quadrature. Halving the maximum time step from 0.006 to 0.003 and tightening relative tolerance by ten also verifies the main and extended runs.

Across all committed runs:

- Max Hamiltonian residual divided by `6h²`: `8.6e-13`; max physical momentum residual divided by `h²`: `1.1e-12`.
- Largest spatial-refinement mean difference: `7.5e-15`; largest temporal-refinement mean difference: `5.1e-16`.
- Largest refined cap-time difference: `2.0e-13`.
- An independently integrated homogeneous closed-FRW scalar solution agrees in mean, volume and H below `7e-16`.
- The zero control has exactly zero scalar and metric deviations.
- The largest final metric logarithmic deviation is `8.35e-4`; the smallest spatial-volume Jacobian ratio to the reference chart is `0.997499`. The geodesic chart remains regular over this interval.

At the negative-branch cap, individual local kinetic/potential shifts reach about 15.3% of the hilltop potential, while proper-volume averaged individual terms remain below 1.7%. These are full nonlinear evolutions within the inherited local potential model. The extension is a local-time result through x=0.41525, not a statement about subsequent global evolution or a homogeneous attractor. Numerical audit values and pass/fail bounds are stored under `verification` in `results.json`.
