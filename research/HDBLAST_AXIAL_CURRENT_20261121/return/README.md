# How much radiation from the shell comes back?

A change in the shell's electric current changes the electric field it radiates into the five-dimensional bulk. Curvature can scatter some of that radiation back onto the shell. We followed the difference field in the same finite radial region as the source simulation. In the three colliding late runs the axial field still reverses direction when the particles supply their own current. [Numerical results](RETURN.json)

The outer tube fixes the incoming wave emitted there. Its **perturbation** is zero. The tube therefore transmits the extra outward radiation without directly reflecting it. Scattering distributed throughout the bulk can return a signal from arbitrarily near the wall. On the sampled geometry from just before first birth to wall time 6.67, crossing the whole region one way takes between 0.0020000 and 0.0020036 model time units. A bound valid for any time profile of the three selected electric coefficients limits each inward characteristic to **0.0912% of the peak change in the corresponding wall electric field** through time 6.67. It gives 0.0816% on the shorter window ending at 6.45 and 0.1063% on the longer window ending at 7.05. The bound applies to the fixed geometry and absorbing tube in this study.

## What changes in the particle calculation?

At time 6.45 the field along the incident pulse still points in its original direction. At 6.67 it points the other way in all three colliding samples. The fixed-inflow column uses the same particle calculation with the original incident wave. The return column allows radiation emitted by the altered current to scatter off the curved bulk and reach the wall again. The times are in model units. [Paired runs](RETURN.json)

| Simulated pairs | Time steps | Final axial field, fixed inflow | Final axial field, with return | First reversal with return | Mean encounters per particle with return |
|---:|---:|---:|---:|---:|---:|
| 4,096 | 2,048 | +0.016164 | +0.016143 | 6.481747 | 3.107 |
| 8,192 | 2,048 | +0.016936 | +0.016864 | 6.481911 | 3.175 |
| 16,384 | 4,096 | +0.016460 | +0.016482 | 6.482531 | 3.121 |

Adding return shifted the first reversal by less than 0.00000064 within each of these paired runs. The variation among the three samples is larger. In the largest run, 82.27% of the particles had scattered by 6.67. Its 20-bin distance from the prescribed two-stream target changed from 47.61% to 47.96%. The shifts across the three samples have both signs and reach 0.59 percentage points. In the 4,096-pair collisionless control, the field had already reversed twice by 6.67. Elastic scattering affects the timing here much more than the return wave. [Complete diagnostics](RETURN.json)

## Extending the run through time 7.05

At 7.05 the reference pulse has delivered 97.86% of its net work. The coupled return runs use 8,192 pairs at two time resolutions and 16,384 pairs at the finer time resolution. All three cross zero a second time near 6.718. This second crossing also occurs with fixed incident characteristics. The resulting particle populations are still far from sharing pressure equally among the three axes. [Extended radial and kinetic results](EXTENDED.json)

| Simulated pairs | Time steps | Second reversal with return | Mean encounters per particle | Common thermal target distance | Largest principal pressure divided by energy |
|---:|---:|---:|---:|---:|---:|
| 8,192 | 4,096 | 6.717683 | 9.376 | 30.40% | 0.764 |
| 8,192 | 8,192 | 6.717688 | 9.610 | 29.89% | 0.767 |
| 16,384 | 8,192 | 6.717533 | 10.248 | 28.67% | 0.750 |

The 20-bin distances in this table compare all particles to a single thermal gas with the same measured energy and momentum. In each run the final self-consistency residual in the incident waveform is at most $5.1\times10^{-20}$ in model field units. Doubling the radial mesh resolution or halving the radial time step on these final wall traces changes the returning characteristic by at most $2.4\times10^{-12}$. We also repeated the 16,384-pair simulation with the doubled radial mesh feeding the particle currents. It gave the same 167,902 accepted collisions and a final axial field agreeing within $2\times10^{-16}$. Its final incoming-wave residual fell below $7\times10^{-19}$.

The late electric field is sensitive to the particle realization. With 16,384 pairs, allowing the wave to return raised the final axial field from $-1.43\times10^{-5}$ to $+0.00668$. In the last 0.01 time units its root-mean-square magnitude rose from $2.14\times10^{-5}$ to 0.00657. The two 8,192-pair return runs instead gave magnitudes of $4.17\times10^{-5}$ and $1.70\times10^{-5}$ on this same window. A small bound on the returning wave does not bound the ultimate change of a nonlinear particle trajectory with random discrete collisions. A prediction for the late field amplitude requires more particle realizations and convergence work. [Paired endpoint and field results](EXTENDED.json)

## The return problem

The inherited radial Maxwell equations use the round three-sphere's lowest curl magnitude, two. The three homogeneous components obey the same linear equations on the prescribed geometry. The axial reference is the original fluid solution. The two additional reference fields are the propagated source waves before adding the electric current of the new particle calculation. For each direction use the change $X$ from that reference as boundary data. Define $I$ and $O$ as changes in the inward and outward bulk characteristics, $U$ as the tangential gauge-potential coefficient, $E=I+O$ and $Z=I-O$. The proper radial coordinate runs from the wall at zero to the outer tube at depth $d=0.002$. At the inner edge $O=X-I$ and at the outer edge $I=0$. All three difference fields start at zero.

Write $N$ for lapse, $r$ for areal radius and $u=r_t/N$. Set $G=N r_x/r$, $H=N u/r$ and $V=4N/r^3$. The equations used for propagation are

$$
\begin{split}
I_t-N I_x &= \tfrac12[(2N_x+3G-3H)I-(G+H)O - VU],\\
O_t+N O_x &= \tfrac12[(G-H)I-(2N_x+3G+3H)O - VU],\\
U_t &= N r(I+O).
\end{split}
$$

We recover the metric on each radial node from the constrained parent solver. On the postbirth window the source background uses six radial elements with four Gauss–Lobatto nodes per element. We interpolate the metric between recorded times and reuse its spatial cubic inside each element. The Maxwell solver uses the inherited upwind radial fluxes and fourth-order Runge–Kutta time stepping. A spatial control subdivides each element and resolves the same interpolated geometry. A second control halves the time step. [Metric and native control extraction](geometry.py) · [Independent radial solver](collar.cpp)

For given $X(t)$ this yields a retarded operator $I=\mathcal R[X]$. The wall simulation adds this wave to the original incident characteristic and recomputes its three currents, particles and fields to produce a new $X(t)$. We iterate this procedure with the same particle tags and random number streams. In a finite realization a small forcing change can alter an accepted collision, so we check the actual waveform residual at the final iteration. The solver allows the particle current and all three homogeneous electromagnetic components to react. The radial metric, lapse, production history and imposed outer incident source are prescribed. This test leaves their gravitational response to the changed matter distribution open. [Coupled iteration](coupled.py)

## Why the bound covers any incident waveform in this family

Let $Y_I,Y_O$ be the suprema of $|I|,|O|$ over space and the time interval of length $T$. Each characteristic leaves this finite interval through a spatial boundary or its initial surface in at most $\tau=d/N_{\min}$. The potential satisfies $|U|\le T\max(Nr)(Y_I+Y_O)$. Take absolute values and integrate the two characteristic equations along these paths. If $M$ contains the resulting upper bounds on the transfer coefficients, the pair of suprema satisfies the componentwise inequality

$$
\binom{Y_I}{Y_O}\le M\binom{Y_I}{Y_O} + \binom{0}{\|X\|_\infty}.
$$

The inner boundary contributes a unit coefficient to $M_{OI}$. All entries are nonnegative. The reported $M$ has spectral radius below one, so $(1-M)^{-1}$ is nonnegative. Its element in the first row and second column bounds $Y_I/\|X\|_\infty$. The calculation bounds the coefficient $|G+H|=N|r_x+u|/r$ directly and bounds the other coefficients by their separate suprema. Each recorded spatial cubic is converted to the Bernstein basis, whose coefficients enclose all values within the element. Linear interpolation in time preserves these enclosures. The margins guard ordinary floating-point rounding. The continuum statement concerns these interpolated coefficients and this Maxwell initial/boundary problem. It does not bound the error in the background Einstein solver. [Bound construction and outputs](response.py)

## Energy and independent numerical checks

For the difference wave alone the radial energy density and outward flux are

$$
q=r^3(I^2+O^2)+2U^2/r,\qquad
f=Nr^3(O^2-I^2).
$$

The time-dependent metric exchanges energy with this wave. Substitution in the equations above gives

$$
q_t+f_x=r^3N_x EZ +\tfrac12Nu r^2(Z^2-E^2)-2Nu U^2/r^2.
$$

We accumulate both boundary fluxes, this volume work, and the nonnegative loss from upwind jumps. The archive reports the residual of this balance on every run. This ledger is quadratic in the difference field and excludes interference with the nonzero background. The wall calculation separately keeps the energy transactions for the full physical fields.

The independent native solver evolves the characteristic equations on the background solver's own internal Runge–Kutta stages. We drive its inner edge with a unit smooth turn-on and oscillation, starting from zero, and compare the returned characteristic against the separate solver. Its peak is about $7.64\times10^{-4}$. The largest absolute discrepancy over the interval is below $1.4\times10^{-7}$, including the interpolation of that test waveform. Halving the return time step and doubling the radial elements preserve the result. Two independently unforced components remain zero.

For reproduction, the main study generates the dedicated field history and particle trajectories. The two commands below can also be run after that stage.

~~~sh
python -I -B return/geometry.py --polarization /path/to/HDBLAST_POLARIZATION_20261121 \
    --output-dir /tmp/axial-geometry
python -I -B return/coupled.py --polarization /path/to/HDBLAST_POLARIZATION_20261121 \
    --scattering /path/to/HDBLAST_SCATTERING_20261121 --history /path/to/dedicated_half.npz \
    --live-output /path/to/live --geometry /tmp/axial-geometry/dedicated.npz \
    --output /tmp/RETURN.json --jobs 2
~~~

For the extended study first produce the 7.05 field history and the matching particle runs with the main study's `live/follow.py`. Repeat the geometry command with `--end 7.05` and a separate output directory. Repeat the coupled command with `--extended`, the extended history and geometry, the directory containing `CROSSINGS.json` and `cross-snapshots`, and `--output /tmp/EXTENDED.json`. Both invocations accept `--snapshots DIRECTORY` to keep the final particles and fields for additional diagnostics.
