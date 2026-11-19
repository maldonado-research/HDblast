# Five-dimensional expansion after the second pulse

At proper time $s=2.16$, the second scalar-clock pulse recovers 37.5% of the first pulse’s decrease in the quiet shell’s expansion-rate coefficient. It recovers 34.3% on the shell perturbed by a negative incoming bulk pulse. The combined response on both shells is positive at the midpoint of the late pulse and negative after switch-off.

The entries below are coefficients of $y^2$ in $h=H/H_0$, evaluated at the same proper time on the two signed evolutions. They are multiplied by $10^6$ for readability. The incremental column comes from its own pair of native solves with the source difference.

| Wall | Proper time | First only | Both pulses | Increment from second pulse |
|---|---:|---:|---:|---:|
| Fine quiet | 2.10 | −28.533 | +0.525 | +29.058 |
| Fine quiet | 2.16 | −29.085 | −18.182 | +10.903 |
| Fine near-negative | 2.10 | −29.049 | +0.110 | +29.159 |
| Fine near-negative | 2.16 | −29.616 | −19.446 | +10.169 |
| Coarse quiet | 2.10 | | | +29.020 |
| Coarse quiet | 2.16 | | | +10.892 |

At illustrative $y=0.1$, the quiet wall has $h_0=0.167667$ at the final readout. Its late pulse adds $6.50\times10^{-7}$ to the fractional expansion rate, and the combined interaction contributes $-1.08\times10^{-6}$. The near-negative values are $+6.06\times10^{-7}$ and $-1.16\times10^{-6}$. The late coefficients of the log scale factor at the same time are $+1.634\times10^{-6}$ and $+1.612\times10^{-6}$. The corresponding changes in comoving daughter energy per $y^2$ are $+13.64794$ and $+13.68378$ in native energy units.

## Native evolution and controls

The finite pair band has $q<\omega<120$. The source calculation retains parent–daughter phase information across the two disjoint scalar-clock gates and uses the same static contact prescription throughout. For unit positive gate amplitudes the late direction is the both-pulses source minus the first-only source. It includes interference. The difference retains all first-variation energy density, pressure, scalar current and their native time derivatives. These enter the scalar and both metric junction conditions at every Runge–Kutta stage, with bulk geometry and the baseline parent state evolving. A first-order explicit source tabulated on the reference trajectory suffices at order $y^2$, since varying its arguments contributes at order $y^4$. The negative signed step is a numerical derivative parameter.

The finer grid uses native spatial spacings 0.0005 and 0.002, with both doubled on the coarse grid. Each run starts from its archived postfork state, uses native fourth-order time stepping, and stops after $s=2.162$. The last archived proper time exceeds 2.167 on all three walls. The previous driver and its vendor files are authenticated before each worker runs. The single driver edit extends the terminal time from 2.105 to 2.162. Wall geometry, evolution equations and diagnostics come from the same pinned implementation. Endpoint derivatives use cubic interpolation in each signed run’s own proper time. JSON also records first-order interpolation, integrated trace results, shell constraints and the incremental expansion history.

- The primary step is $\pm 0.01$. Halving it changes the final quiet and near-negative incremental coefficients by $+4.09\times10^{-10}$ and $-3.83\times10^{-10}$. Subtracting the two independently evolved total responses agrees with the direct late direction to $1.1\times10^{-11}$ on quiet and $6.4\times10^{-10}$ on near-negative.
- An independent diagnostic integrates $d(a^2h)/ds$ from the native wall stress, scalar velocity, tension and potential. Its primary incremental coefficients are $1.090306\times10^{-5}$ and $1.016906\times10^{-5}$ at $s=2.16$. Each differs from its direct metric readout by less than $8\times10^{-11}$. The trace uses the same evolved geometry and supplies an integrated equation check.
- The coarse–fine quiet differences in the increment are 0.130% at $s=2.10$ and 0.104% at $s=2.16$. Two spatial resolutions support a consistency check. Continuum extrapolation would require additional resolutions.
- The second gate starts at $s=2.055000$ on fine quiet, $2.054314$ on fine near-negative and $2.054894$ on coarse quiet. Every explicit increment component and its tabulated slope vanish through onset. Positive and negative native trajectories are bit-identical on all 1,801, 1,799 and 901 pre-gate rows respectively for both metric and matter diagnostics. Equal-sign first-only and both-pulse trajectories also coincide there on the two fine walls.

## Reproduction

From the study directory, generate the unit-amplitude first, second and both source tables for fine quiet, fine near-negative and coarse quiet. Pass their common directory to the native solver. It contains the per-arm event times in `GATES.json` and the `first`, `second` and `both` subdirectories.

```sh
python -IB gravity/gravity.py --repo /path/to/upstream --tables /path/to/tables \
  --jobs 2 --output /tmp/repump-gravity.json --compare gravity/GRAVITY.json
```

The upstream repository snapshot is `44370269720415893dcd8d5bc802697e629d7655`. The replay uses the earlier quantum-memory gravity driver and scalar-clock live solver, authenticated vendored sources, native fork and reference histories. It runs 18 native evolutions. Supplying `--work` retains the signed native histories and source coefficients. The committed JSON records input hashes and readouts. The comparison checks numerical observables for newly generated tables. Each native result is checked against the hash of the source table it consumed.
