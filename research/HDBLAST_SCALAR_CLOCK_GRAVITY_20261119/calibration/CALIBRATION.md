# Can identical outgoing pair energy accompany opposite expansion responses?

For the quiet shell, the direction of the leading expansion change depends on how the parent field responds virtually to the fermion pairs. We tested a third finite Hamiltonian on both archived spatial grids. A calculation using only the first two Hamiltonians predicts the sign for the third. This was an exploratory follow-up after inspecting the first two slopes. The third coefficient and the test were recorded before reading its simulated result. [Test design](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7001212633)

## The family and its prediction

Write $\lambda=y^2$ for the physical nonnegative coupling strength. In each retained parent mode, the static completion added in the preceding energy study is

$$
H_{\mathrm{ct},q}=\frac{\kappa}{2}\lambda a F(\phi_b)^2 S_q\chi_q^2,
\qquad S_q=\sum_j g_{qj}^2/\omega_{qj}^2.
$$

The sum and all numerical results use the finite comoving cutoff $\Lambda=120$. Here $\kappa$ multiplies a specified compensation of the change in the oscillator’s restoring force. It is separate from the native constraint-damping coefficient. The bare, statically matched and third Hamiltonians have $\kappa=0,1,3$ respectively. Completing the square in the bath coordinates produces the static shift $-\lambda a F^2 S_q\chi_q^2/2$. Thus $\kappa=3$ overcompensates this shift in the finite Hamiltonian. It gives a model with a different parent susceptibility. Treating the three as physically interchangeable would require a common matching condition for that susceptibility and the corresponding stress and scalar current.

At first order in $\lambda$ the explicit canonical parent response, counterterm work, pressure and scalar current are affine in $\kappa$. The unperturbed background and parent covariance are common, since every added term is zero at $\lambda=0$. Linearizing the full bulk evolution about that same solution therefore gives the prediction

$$
S(\kappa)= S(0)+\kappa\,[S(1)-S(0)],
\qquad S=\left.\frac{\partial (H/H_0)}{\partial\lambda}\right|_{\lambda=0,s=2.10}.
$$

The daughter-energy coefficient at order $\lambda$ is the same across the family. A change in the order-$\lambda$ parent response enters the outgoing energy only at higher order. The strongest test here is whether evolving the complete geometry at the third calibration recovers the predicted positive sign, after both reference evolutions gave a negative sign.

## Two-grid test

The third calculation builds $\delta\rho$, $\delta p$, $\delta j$ and their coordinate-conversion inputs from the bare and static **source histories**. It then reruns the coupled five-dimensional solver at signed auxiliary strengths $\pm0.01$ and $\pm0.005$. The outgoing radiation sources agree pointwise between the bare and static inputs. Expansion is read independently from the metric velocity at matched proper time on each resulting trajectory. Feeding a linear combination of two *final expansion rates* back into a summary would bypass this test. The table gives the centered differences with step 0.005. Reported pair energies are coefficients of $y^2H_0^4$ per comoving volume.

| Grid | $S(0)$ | $S(1)$ | Direct $S(3)$ | $\kappa$ at zero response, inferred from $0,1$ | Pair energy, all three |
|---|---:|---:|---:|---:|---:|
| Coarse | $-4.48044\times10^{-5}$ | $-2.85309\times10^{-5}$ | $+4.01591\times10^{-6}$ | 2.75321 | 22.51758 |
| Fine | $-4.49152\times10^{-5}$ | $-2.85332\times10^{-5}$ | $+4.23216\times10^{-6}$ | 2.74173 | 22.51057 |

The sign reversal survives halving the strength step on both grids. The largest absolute difference between the held-out direct slope and its affine prediction is $1.20\times10^{-9}$ over these four centered tests. The extra $O(\lambda^2)$ Richardson extrapolation of the centered slopes gives a fine-grid difference of $1.98\times10^{-9}$, illustrating the precision available when differencing the metric. The daughter-energy coefficients across the three calibrations differ numerically by less than $2.2\times10^{-10}$ within each grid and strength pair. The spread in $S(3)$ between the grids is about 5% of the fine-grid value. The control establishes a finite-model sign difference. Determining which susceptibility a real shell would have requires independent physical input. The three couplings yield the same leading quiet-wall radiation while predicting opposite signs for the change in its expansion relative to the uncoupled solution. [Numerical controls](CALIBRATION.json) · [Code and reconstruction](calibration.py)

With the source tables and the reference bulk/runs available, run this command from the study root. A trusted geometry cache built by the bulk replay can be supplied with `--cache`. The optional `--intermediates` directory retains the third model’s source tables and eight additional histories.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB calibration/calibration.py --bulk bulk --tables source/tables --runs bulk/runs --output /tmp/CALIBRATION.json
~~~
