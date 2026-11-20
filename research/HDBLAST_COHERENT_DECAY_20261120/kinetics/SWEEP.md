# Can new fermion states become available quickly enough?

A time-integrated count of available momentum states can conceal a brief supply shortage. We resolve those intervals in two large-mass, cycle-averaged trajectories. This calculation uses continuum phase-space density, a homogeneous force, a sharp two-body birth shell and collisionless transport at injection. It tests the supply of newly swept states, summed over emission directions. The underlying fluid accepts the prescribed source even if these conditions predict blocking. [Finite magnetic-state counts](MANY.md) · [Parameters of the averaged histories](../envelope/ENVELOPE.json)

## Local rate and electric reversals

In comoving momentum coordinates the birth sphere has radius $k=a\mu\sqrt{R^2-1}$ and a center translating along one axis with $\dot b=\mu Q aX$. Write $c=|\dot k|$ and $d=|\dot b|$. A surface element at polar angle $u=\cos\vartheta$ sweeps volume at a rate proportional to $|\dot k+\dot b u|$. The favorable upper capacity, allowing a vacant state on either side of the surface, is
$$
\dot N_{\rm sweep}=12 k^2 F(c,d),\qquad
F(c,d)=\tfrac12\int_{-1}^{1}|c+du|du
=\begin{cases}c,&d\le c,\\(d^2+c^2)/(2d),&d>c.\end{cases}
$$
The factor 12 counts both spins of all six sign/flavor populations on the closed sphere. Equivalently the Weyl density per spin and species is $2\pi^2/(2\pi)^3$ per unit comoving momentum volume. The requested rate is $2\pi^2 S_0\dot D_c/\mu$. The action $S_0$ and the rest-energy injection $\dot D_c$ follow the same conventions as the finite-state calculation.

When the electric field reverses, the instantaneous sweep can be much slower than its average. Each native solver step uses the actual change in cumulative injected rest, the difference of the birth-shell radii, and the linearly interpolated signed electric acceleration. For an aggregate rate $P$, the local deficit is $[P-\dot N_{\rm sweep}]_+$. The code integrates this positive part analytically over the substep, including crossings of zero, $|d|=c$, and $|d|=P/(12k^2)+\sqrt{[P/(12k^2)]^2-c^2}$ when a deficit exists. The slowly changing radius and injection rate use their step-centered values. Two closed-form electric ramps validate this integration. The comparison with halved time step is an independent accuracy check for the sharp peaks. [Integrator](sweep.hpp) · [Replay](sweep.py)

## Daughter mass $mL=4\times10^{10}$

The action is $S_0=3.609e+48$, the earlier joint-screen example. At this action the thermal drag target requires 5.77% of its single-channel ceiling at the tightest reference state. The table covers the central 99.998% of total production, from cumulative fraction 0.00001 to 0.99999. On the half-step grid this spans proper wall times 6.252192 to 15.912045.

| Resolution | Peak demand / capacity, interpolating the field | Peak at native step endpoints only | Percentage of total production locally unfilled | Largest demand / capacity in full fixed bins of widths 0.0001, 0.001 and 0.01 |
|---|---:|---:|---:|---|
| native | 285.1 | 95.59 | 0.105164% | 39.25, 5.027, 0.4495 |
| half step | 285.1 | 147.8 | 0.105163% | 49.11, 5.133, 0.4508 |

The largest interpolated ratio on the half-step grid is 285.06, at proper time 6.618352 and cumulative fraction 0.351139. At the smallest action consistent with the reference hot-drag rate, the largest ratio is 16.455. Across this delivery window the conditional instantaneous action ceiling is therefore $S_0=1.266e+46$. Meeting the hot drag at that limiting action would require **1.65e+03%** of the thermal elastic unitarity ceiling. This exceeds the available local phase-space flux at the example action. The half-step trajectory has 5 signed electric reversals inside the tested interval. For the nested 1%–99% window (6.263501 to 10.122976), the maximal ratio is 285.06 and the estimated unmet share is 0.105163% of all production.

## Daughter mass $mL=4\times10^{12}$

The action is $S_0=1.042e+54$, five times the reference hot-drag floor, chosen before this sweep. At this action the thermal drag target requires 20.00% of its single-channel ceiling at the tightest reference state. The table covers the central 99.998% of total production, from cumulative fraction 0.00001 to 0.99999. On the half-step grid this spans proper wall times 6.252188 to 15.914101.

| Resolution | Peak demand / capacity, interpolating the field | Peak at native step endpoints only | Percentage of total production locally unfilled | Largest demand / capacity in full fixed bins of widths 0.0001, 0.001 and 0.01 |
|---|---:|---:|---:|---|
| native | 0.8228 | 0.3417 | 0% | 0.1111, 0.01449, 0.001297 |
| half step | 0.8228 | 0.5207 | 0% | 0.1372, 0.01476, 0.0013 |

The largest interpolated ratio on the half-step grid is 0.82277, at proper time 6.618348 and cumulative fraction 0.351140. At the smallest action consistent with the reference hot-drag rate, the largest ratio is 0.16455. Across this delivery window the conditional instantaneous action ceiling is therefore $S_0=1.266e+54$. Meeting the hot drag at that limiting action would require **16.5%** of the thermal elastic unitarity ceiling. The specified action passes this angle-integrated local sweep screen at both temporal resolutions. The half-step trajectory has 5 signed electric reversals inside the tested interval. For the nested 1%–99% window (6.263497 to 10.122974), the maximal ratio is 0.82277 and the estimated unmet share is 0% of all production.

The two discarded extreme tails together contain 0.002% of the final comoving particle number. The first and last partial clock bins are omitted from fixed-bin maxima. Each complete bin reports the ratio of integrated demand to integrated capacity; the positive instantaneous deficit counts over-capacity requests even when spare capacity exists elsewhere in the same bin. The reported partial-wave fraction uses the hot small-drift s-wave benchmark and its numerical action floor. The rate screen is angle-integrated. A favorable sum can still hide crowded directions, angular selection rules or previously occupied states. Collisions, other final states and a finite production bandwidth change the instantaneous test.

The high-mass continuum calculation is motivated by the many occupied magnetic branches, with an independent finite-state check at $4\times10^{10}$. It assumes instantaneous, sharp quasiparticle energies in the prescribed electromagnetic background and separates production from later collisions. It provides a conditional necessary test of this kinetic interpretation. A quantum kinetic treatment would have to track the actual occupations, spectral width, accessible scattering channels and field response.

Run both cases with their native and halved time steps from the capsule root.

```sh
python -I -B kinetics/sweep.py --envelope envelope/ENVELOPE.json --output /tmp/SWEEP.json
```

[All windows, grid checks and normalization choices](SWEEP.json)
