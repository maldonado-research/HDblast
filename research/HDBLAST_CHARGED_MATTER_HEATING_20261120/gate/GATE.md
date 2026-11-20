# A local, delayed onset of electric response

Can the charged gas heat before contraction resumes if the material has no preset switch-on time? In this test the material responds to the local scalar and electromagnetic fields on the wall. On the standard 48-element early mesh, internal charged-gas energy reaches 91.51% of added wall energy. The gas crosses 90% at radius 9.6188 and remains above it as the radius grows by another 1.485% to 9.7616. The wall then starts shrinking. The prescribed constitutive law and the perfect-fluid assumptions remain central to this result.

## Matter and local trigger

Use the same closed wall, source pulse, scalar potential, positive Maxwell kinetic function and six opposite-charge streams as the coupled charged-gas study. The neutral precursor rest density at the initial radius 2 is 0.0001. Its rest energy gravitates from the start. No physical carrier species, binding energy or ionization cross section has been specified. The terminal charge-to-mass parameter is $Q=10^8$ in the model's projected electromagnetic units and the drag rate is $\nu=6$ per unit of wall proper time. The local electric and magnetic amplitudes measured in the wall frame are $X$ and $B$.

The constitutive rule has a memory. At the first accepted time step where the scalar has descended its potential, $\psi\ge0.38$, and $B^2/(B^2+X^2)<\epsilon$ with nonzero denominator, a mobility variable latches on. Its age $u$ then grows at one unit per unit of wall proper time. The force and the Maxwell current both use $Q_{\rm eff}=Q S(u/w)$, where $S$ is the parent solver's smooth monotonic step with $S(0)=0$, $S(1)=1$ and flat endpoints. The primary inputs are $\epsilon=0.01$ and $w=0.5$. Once active, the age continues to grow even if the magnetic share rises. A local observer determines the onset from field measurements and the material's subsequent response from elapsed proper time. The scalar cut, field threshold, irreversible latch and mobility law are material assumptions.

Time-dependent mobility represents the strength of electrical contact. Before contact the initially cold precursor is neutral to the force. During the ramp there is no carrier-number injection, mass jump, binding-energy term or impulsive transfer. The equal positive and negative streams have zero net charge density and their homogeneous spatial current is divergence free. The bulk-to-wall Maxwell evolution subtracts $jX$ from the sheet and adds the identical work $jX$ to matter, using the same $Q_{\rm eff}$ at every Runge–Kutta stage. This gives a homogeneous, energy-accounted transport law. Microscopic charge conservation for conversion between neutral and ionized species, the energy needed to free charges and a gauge-invariant action for that conversion would need additional fields. The run leaves open whether an actual material could realize the chosen response.

## Heating and sensitivity

The 48-element wall bounces at $\tau=0.86697$, completes its scalar descent at 4.51324 and first meets the local mobility condition at 10.02608. At contact the scalar is 0.75988, the magnetic fraction is 0.00999984, and the gas has exactly zero electric work and temperature. At 10.82926 the internal-energy fraction of the full added wall density passes 90%. It peaks at 0.91512 as the wall turns around at 12.42350. Internal energy here is $n_{\rm lab}m(h-1-\theta)$ for the common streams, where $h$ is rest-frame specific enthalpy and $\theta=T/m$. This counts random energy in the stream rest frames, using wall-frame particle number. It excludes particle rest mass and ordered drift. The background wall tension is excluded from the denominator. The gas closure assumes local equilibrium within each stream and drag between streams. No photon spectrum is evolved.

| Electrical condition and ramp | Latch time | Peak internal fraction in expansion | Expanding e-folds above 90% |
|---|---:|---:|---:|
| Magnetic fraction below 1%, width 0.50, 48 elements | 10.02608 | 0.91512 | 0.01474 |
| Same, 64 elements | 10.02603 | 0.91513 | 0.01474 |
| Below 1%, width 0.25, 24 elements | 10.02675 | 0.91259 | 0.01675 |
| Below 1%, width 1.00, 24 elements | 10.02675 | 0.91657 | 0.01083 |
| Below 5%, width 0.50, 24 elements | 9.44088 | 0.89179 | 0 |
| Below 10%, width 0.50, 24 elements | 8.99960 | 0.85591 | 0 |
| Below 0.01%, width 0.50, 24 elements | 10.45469 | 0.91447 | 0.00782 |

Turning off the force or drag gives zero internal heating. The collisionless control still develops a directed current and satisfies $a u_{\rm drift}=\int aQ_{\rm eff} X\,d\tau$ to a normalized residual below $2\times10^{-10}$. An impossible scalar threshold or zero magnetic-fraction cut keeps the neutral rest energy present without electrical contact. With no source pulse the wall contracts through the inner stopping radius with zero charge work or temperature.

Lowering the scalar cut to 0.004 opens the same 1% magnetic gate at 0.01327, while the wall is still in the early loading phase. Two resolutions follow this counterexample through $\tau=0.12$, where $T/m$ is 0.03848 and the accumulated electric work is about $5.37\times10^{-5}$. These short runs establish that the scalar cut matters for the onset and do not determine a late outcome. Because the early Maxwell coefficient is small, this check also needs a time step resolving the cold plasma frequency. The patch adds the cap $\Delta\tau\,Q_{\rm eff}\sqrt{\rho_{\rm rest}/K}\le0.1$ when the carriers are active.

At 48 elements the peak gas work error is $1.6\times10^{-13}$, the comoving sheet error is $9.4\times10^{-10}$, the total wall work error is $4.6\times10^{-10}$ and the Friedmann residual is $1.3\times10^{-10}$. Across the collar the enclosed mass remains above 0.0965 and the lapse above 0.9977. Refining from 48 to 64 early elements changes peak internal fraction by $8.7\times10^{-6}$. A separate test changes the late radial mesh from 6 to 10 elements and the time-step factor from 0.08 to 0.05 on the same 24-element early mesh. The peak fraction changes by $1.2\times10^{-6}$. The latch is detected on accepted steps and the smooth ramp then evolves continuously. The growth interval, and sensitivity to field threshold, remain properties of one selected material and initial state.

## Reproduce

Run `python -I -B gate/gate.py --output /tmp/GATE.json --jobs 2` from the study root. Python, a C++17 compiler named `g++` and the standard `patch` utility are required. The program authenticates the exact parent solver, applies a small source diff in a temporary directory and recalculates the stored cases and balance checks. [Parent fluid equations](../coupled/COUPLED.md) · [Implementation diff](mobility.diff) · [Runner](gate.py) · [Results](GATE.json)
