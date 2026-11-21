# Elastic encounter histories

The radiation calculations follow the actual accepted collisions of the dedicated-interaction shell. We replayed the particle-current calculation at wall times 6.67 and 7.05, recording the positive and negative daughters immediately before and after each accepted elastic collision. The recorder consumes no random numbers. All eight replays recover the original event totals, the fractions of particles ever hit, and the final gas energies to within $4.3\times10^{-15}$ in relative terms. [Compact results](EVENTS.json) · [Source](replay.py)

| Wall proper time | Pairs | Steps | Accepted collisions | Median pair mass / daughter mass | Median pair Lorentz factor in wall frame | Fraction within one running width |
|---:|---:|---:|---:|---:|---:|---:|
| 6.67 | 16,384 | 4,096 | 51,015 | 8,893 | 1.86 | 51.0% |
| 7.05 | 16,384 | 8,192 | 154,876 | 9,295 | 1.92 | 46.5% |

The collision kernel has a resonance mass $M/m=8300$ and a denominator $(s-M^2)^2+W(s)^2$ in daughter-mass units, with $W(s)=0.30(s+2)\sqrt{1-4/s}$. The last column counts encounters satisfying $|s-M^2|\le W(s)$. At 7.05, 42.3% of the fine sample falls above this band. Those encounters account for 49.9% of the sum of incoming pair energies over events. A particle can collide repeatedly, so this sum measures the weight of different encounters and cannot be interpreted as an independent supply of energy.

At 7.05 the 1st and 99th percentiles of the pair mass are 4,727 and 60,084 daughter masses. The spread matters for a radiation calculation that depends on the invariant collision energy. Only 0.074% of events in the fine sample revisit an identical pair of simulation tracks. Individual tracks often participate in several events. Each endpoint starts from the specified birth law, with a sampling weight fixed by its own final population. Track identifiers apply within one file.

## What the generated files contain

The per-run files are transient NPZ arrays named by the exact label in the parent kinetic study. The incident and scattered particles have equal unit rest mass. Each event has the following entries.

| Arrays | Meaning |
|---|---|
| `t`, `a` | Accepted Poisson event time in wall proper units $L$, and radius $a/L$ used for the collision. The split particle solver evaluates momenta and radius after the wall step's force update, at the right-edge radius. |
| `i`, `j` | Zero-based global track IDs. The incoming positive track is even and the negative partner is odd. Adjacent tracks $2k,2k+1$ share a birth. |
| `p_in`, `p_out` | Spatial momenta in the orthonormal wall frame, divided by daughter mass, stored in positive-then-negative order. Shape is event count by 2 by 3. |
| `E_in`, `E_out` | Both time components in the same order, computed as $E/m=\sqrt{1+|\mathbf p/m|^2}$. |
| `s`, `mu` | The invariant squared pair mass in daughter-mass units and the center-of-momentum scattering cosine from the incoming positive particle to the outgoing positive particle. |

The per-track arrays `birth` and `p_final` give every track's original birth time and final wall-frame momentum, including unhit tracks. The arrays `wall_t`, `wall_E` and `wall_B` give the parent simulation's physical field traces at the wall integration nodes. Metadata include the endpoint proper time `time`, radius `endpoint_a`, pair count `pairs`, representative weight `weight` $=D/(2N_\mathrm{pairs})$, and the conformal endpoint `gas` $=\sum w a_f\gamma_f$. The parent values of `action`, `charge_to_mass`, `daughter_mass_L`, `channel_fraction` and the gauge coefficient `K` are included for normalization.

We reconstruct the invariant mass independently from each pair of incoming and outgoing on-shell four-vectors. We also recover the reported scattering cosine from the invariant momentum transfer using $\cos\theta=1+2t_{\rm Mandelstam}/(s-4)$. The largest discrepancy in $s$ is under $3.4\times10^{-13}$ in relative units. The largest cosine difference is below $1.5\times10^{-13}$. Individual event four-momentum balance is checked in the wall frame, and event times must follow the births of both participating tracks. Subsequent photon estimators use this fixed elastic ensemble. Generating emitted photons with recoil would require a new collision evolution.

~~~sh
python -I -B events/replay.py --axial "$axial" --polarization "$polarization" \
    --scattering "$scattering" --history "$extended/dedicated_half.npz" \
    --output "$results/EVENTS.json" --artifacts "$results/events" --jobs 2
~~~

The command uses the pinned wall-current, polarization and scattering study directories. Generate the extended input using the wall-current study's `fields/extract.py` with `--cases dedicated_half --end 7.0505`. The test reuses the parent solver, birth sample, random seeds and collision kernel. The archive tracks summary statistics and the generating source. The detailed arrays can be regenerated for independent analyses.
