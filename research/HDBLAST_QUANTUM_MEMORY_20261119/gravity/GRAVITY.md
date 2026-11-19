# Testing four high-frequency modes in the five-dimensional evolution

The test compares three finite quantum baths with the same parent and the same incoming initial data. The low-only model keeps the original 64-node band $q<\omega<60$. The full reference adds a separate 64-node frequency quadrature on $60<\omega<120$. The compressed model replaces only that upper band by four positive nodes per parent momentum. Four nodes were selected before inspecting their gravitational errors. The earlier static compensation is applied to each model's complete spectral measure so that all three have zero residual frozen restoring-force shift. The compressed upper band also reproduces the reference static and kinetic spectral moments. [Registered study](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7015073318)

The matter input includes the causal parent variation, daughter energy, interaction, and compensating potential and its associated stress and scalar current. The vacuum fluctuations and the driven retarded response come from the same quantum bath. Each signed run evolves the five-dimensional Einstein–scalar variables and the original quantum parent at every Runge–Kutta stage, as in the earlier native study. The explicit first-order matter source is evaluated on that arm's own archived zero-coupling history. [Native bulk solver and source normalization](https://doi.org/10.5281/zenodo.24360096)

## What the subtraction measures

For model $m$, let $h=H/H_0$. We evaluate the centered slope $S_m(s)=[h_m(+\epsilon,s)-h_m(-\epsilon,s)]/(2\epsilon)$ in $\lambda=y^2$ at the same proper clock $s$ on each trajectory. The high band's exact increment is $I=S_{64}-S_{\rm low}$. The four-node prediction is $I_4=S_4-S_{\rm low}$ and its error is $E=I_4-I=S_4-S_{64}$. Additivity follows at first order on the common unperturbed geometry. A comparison of $S_4$ and $S_{64}$ as a fraction of the entire response can hide a larger error in $I$. The high increment includes the response of the parent, the full boundary conditions and the bulk to the upper band. Isolated daughter radiation contributes only one part.

The simulations use fine quiet, coarse quiet and fine near-minus initial data. They comprise eighteen signed evolutions at step $\epsilon=0.01$ and six additional fine-quiet evolutions at $0.005$. The recorder stores every native time step. For each run we interpolate the native observables cubically in its own advancing proper time, with a linear-interpolation endpoint check. Along $1.35\le s\le2.10$ we report $\|E\|_2/\|I\|_2$, using $ds$ and 1,201 common sampling times. The endpoint relative error is $|E(2.10)|/|I(2.10)|$. The denominator retains its own physical scale even where the positive and negative contributions of the band nearly cancel. [Geometry histories, energy, sources and numerical controls](GRAVITY.json)

## Measured gravitational increment

Adding the resolved high band increases the fine quiet expansion coefficient from $-2.90069\times10^{-5}$ to $-2.85329\times10^{-5}$ at the primary strength step. The total interaction still reduces expansion relative to the common zero-coupling trajectory. The upper band offsets part of that reduction, giving $+4.73957\times10^{-7}$ in $\partial_\lambda(H/H_0)$, or about $+0.0268$ ppm in the expansion rate at illustrative $y=0.1$. Four positive nodes recover the direction on each tested history.

| Grid and incoming history | Signed step | Reference increment $I$, in $10^{-7}$ | Four-node increment $I_4$, in $10^{-7}$ | Endpoint $\lvert E/I\rvert$ | History $\|E\|_2/\|I\|_2$ |
|---|---:|---:|---:|---:|---:|
| Fine quiet | 0.010 | +4.73957 | +4.73762 | 0.0410% | 0.0139% |
| Fine quiet | 0.005 | +4.74676 | +4.75249 | 0.1206% | 0.0234% |
| Coarse quiet | 0.010 | +4.84740 | +4.85225 | 0.1001% | 0.0106% |
| Fine near-minus | 0.010 | +4.64791 | +4.67325 | 0.5452% | 0.2018% |

On the held-out near-minus pulse, the four-node expansion excess is $2.5341\times10^{-9}$, and integrating the wall acceleration gives $2.5731\times10^{-9}$. The high band contributes 0.9596445 units of outgoing energy in the reference pulse. The four nodes give 0.9688022, an overestimate of 0.9543%. The corresponding quiet high-band energy is 0.9326764 and the four-node energy is lower by 0.00241%. These are coefficients in $y^2H_0^4$ per normalized comoving volume. The measured total-energy or total-expansion errors would be smaller because the low band dominates both totals. [Endpoint coefficients for all sources and geometry](GRAVITY.json)

The absolute quiet high-band increment changes by about 2.3% between the two native meshes at step 0.01. On a fixed fine grid it changes by $7.20\times10^{-10}$ when the strength step is halved. The full split-band slopes on all three histories agree with the previously published, unsplit direct-frequency calculation to within $6.5\times10^{-11}$ on this implementation. The tested geometry supports the leading response and the comparison of the finite models. A continuum error estimate would require further spatial refinement.

## Resolving the small quiet-shell discrepancy

On fine quiet, subtracting the four-node and full metric slopes gives $-1.94\times10^{-10}$ at step 0.01 and $+5.72\times10^{-10}$ at step 0.005. Those two subtractions alone leave the residual's sign unresolved. An integrated wall trace supplies a separate estimate using the live scalar velocity, potential, tension, density, pressure and signed scalar current in $d(a^2h)/ds$. This gives $+5.11\times10^{-10}$ and $+5.03\times10^{-10}$ respectively. The reconstructed trace has its own junction, constraint and time-discretization errors. [Earlier trace derivation and test](https://doi.org/10.5281/zenodo.24360096)

After seeing the initial subtraction, we registered a further numerical check. At first order the derivative of a solution with respect to an additive source is linear. We formed the difference of the complete four-node and reference source tables, including all time derivatives, and applied that single signed difference to four new native evolutions of the same quiet initial data. The auxiliary strengths are $\pm1$ and $\pm0.5$, applied only to the small difference. They are numerical probes of the tangent equations. Matching each signed pair at its own proper time cancels the zero-coupling solution. [Exploratory control protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7015211027)

| Auxiliary step on the difference source | Endpoint slope from metric | Endpoint slope from integrated trace | Error relative to reference high increment | Relative history error |
|---|---:|---:|---:|---:|
| 1.0 | $+5.229\times10^{-10}$ | $+5.201\times10^{-10}$ | 0.1103% | 0.00933% |
| 0.5 | $+5.165\times10^{-10}$ | $+5.199\times10^{-10}$ | 0.1090% | 0.00933% |

Linear interpolation of the auxiliary endpoint instead of cubic changes the half-step metric slope by $4.1\times10^{-13}$. The small isolated positive response and the independently integrated trace put the fine-grid quiet discrepancy at about **0.11% of the high band's endpoint effect** in this finite model. The extrapolated final digits would depend on roundoff, the difference scheme and spatial resolution. Larger errors on the near-minus history limit the extent to which that quiet result transfers to another incoming perturbation.

## Reproduction

Pass the corresponding nine generated static source tables and the pinned upstream repository to the native test. They can be produced using the source script with the three indicated cases and high orders 0, 4 and 64. The native solver authenticates its upstream files and checks the exact preturn equality of all signed trajectories on each fixed grid and arm. Each native history records the SHA-256 of the source table it actually consumed. A new floating-point environment can produce numerically equivalent source tables with different hashes.

~~~sh
python -IB source/source.py --native /path/to/repo/research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native --cases fine_quiet coarse_quiet fine_minus --high-orders 0 4 64 --output /tmp/hdblast-memory-tables --jobs 2
python -IB gravity/gravity.py --repo /path/to/repo --tables /tmp/hdblast-memory-tables --output /tmp/GRAVITY.json --jobs 2
~~~

The optional `--work` directory retains the 28 full-cadence time series and exact input tables. A trusted, compatible geometry cache can be provided with `--cache`. The experiment uses comoving regulators, a finite Gaussian parent and the specified interaction clock. The result describes the tested histories and momentum band.
