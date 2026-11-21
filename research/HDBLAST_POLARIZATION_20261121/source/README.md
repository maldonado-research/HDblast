# Can the shell supply the added fields and still produce daughters?

The emitting tube sends the same total energy in two trials. In the first, both added electric patterns arrive from an early pulse. In the second, the tube holds back about one part in 1,600 of its energy for later. The first trial produces charged particles. In the second, the scalar field barely enters its first preparation interval and the decay channel stays closed through wall time 8. All times and energies below use the parent model's dimensionless units. [Numerical results](SOURCE.json)

## Inputs and scope

The numerical model starts from the dedicated-interaction, $mL=4\times10^{12}$ benchmark of the parent study. The original axial mode has curl magnitude 2. Two further angular patterns on the round three-sphere have unit mean square and vanish in their pairwise inner products with each other and the axial pattern. At equal curl magnitude they occupy different real quadratures. The new pulses have coefficients $\epsilon/\sqrt{2}$ each and the original pulse has coefficient $\sqrt{1-\epsilon^2}$. Thus $\epsilon=0.08$ diverts **0.64%** of the angularly integrated early source power. This prescription allows the local intensity to vary over the emitting tube. The origin of these particular pulse shapes remains an input.

Each mode propagates through the bulk and the moving boundary. The simulation updates the mean scalar, radius, expansion, bulk mass, radial constraint fields and original matter variables. The new modes enter the mean stress and the scalar's equation of motion at every integration stage. A separate diagnostic switches these terms off while propagating the new waves and reducing the original pulse. Both of those diagnostic runs missed the opening through time 6.36. [Implementation](backreaction.py) · [Angularly averaged stress](closure.hpp)

The calculation projects the geometry and matter onto a round shell at each time. Angular fluctuations in the scalar and metric, shear and transverse daughter currents are omitted here. An angular fluctuation linear in the added amplitude can contribute to the mean at quadratic order. Consequently, the leading physical correction to the five-dimensional mean could differ from this calculation. The conditional daughter law after opening uses the original axial current. The pre-opening results precede that current. [Parent geometry and source construction](../README.md)

## The scalar's energy budget

Write $t$ for wall proper time, $a$ for wall radius, $H=\dot a/a$, $\psi$ for the homogeneous scalar and $K(\psi)$ for its electromagnetic coupling. Orthonormality removes the cross terms in the spatial average. The gauge contribution to $\ddot\psi$ is

$$
\frac{K'(\psi)}{2}\sum_{j=0}^2(E_j^2-B_j^2).
$$

Work taken from the scalar by the two extra patterns is therefore

$$
W_{\rm extra}(t)=\frac12\int_0^t a^3\,\dot K
       \sum_{j=1}^2(B_j^2-E_j^2)\,dt'.
$$

Positive $W_{\rm extra}$ is energy transferred to the extra fields. The parent oscillator, the axial field and the scalar's pressure work also change between trials. The simulation accounts for their energy transfers separately. A value of $W_{\rm extra}$ larger than the incident tube energy is possible because the evolving scalar can supply energy to the fields.

For each added mode the independent field ledger uses the radiation scaling $a^4$. Here $i_j$ and $o_j$ are the incoming and returning wave characteristics at the wall.

$$
\frac12 a^4 K(E_j^2+B_j^2)
=2\int_0^t a^4(i_j^2-o_j^2)\,dt'
 +\frac12\int_0^t a^4\dot K(B_j^2-E_j^2)\,dt'.
$$

The physical source energy reported here is $2\int r_{\rm tube}^3\sum_j i_{j,\rm tube}^2\,d\tau_{\rm tube}$ per unit solid angle. Returning energy and the wall's energy ledger are separate diagnostics. For the matched delayed runs the simulation recalibrates the later pulse against the evolving tube radius and lapse. [Calibrated waveforms](../fields/waveforms.hpp)

## Results

| Added early patterns | Emitted energy share | First daughter, wall time | Scalar comoving energy at first birth | Scalar work given to added modes by first birth |
|---|---:|---:|---:|---:|
| None | 0 | 6.25161 | 11.5273 | 0 |
| 2 and 2 | 0.64% | 6.25161 | 11.5273 | 0.38731 |
| 5 and 7 | 0.16% | 6.25599 | 10.4677 | 1.92109 |
| 5 and 7 | 0.64% | 6.26983 | 7.44247 | 7.55612 |
| 5 and 7 | 2.56% | 6.38781 | 4.33994 | 25.98388 |
| 5 and 7 | 4.00% | None through 8 | — | — |

First birth is recorded at the last step with zero accepted daughter rest before that rest becomes positive. The same-mode rotation recovers the original first-birth time and scalar energy to roundoff, as all three equally shaped, lowest-curl waves have the same radial gain before any daughter current appears. At 0.64% in modes 5 and 7, halving the time step gives first birth at 6.26985 and scalar work 7.55589. Doubling the radial cell counts gives 6.26848 and 7.57820. The scalar's own comoving energy at first birth shifts to 7.03888 on the doubled radial grid. Its mesh dependence is appreciably larger than the relative changes in the integrated work.

At 4% the scalar reverses inside the last change of gauge. Its native, half-step and doubled-radial peak values are 0.527966, 0.527966 and 0.527939 against an opening value of 0.530000. The smaller 2.56% diversion still produced daughters. These samples apply to the specified modes over the simulated interval. They leave the intervening amplitudes open.

For the timing test, each of two independent mode-7 patterns initially receives amplitude 0.025. Emitting both early produces daughters at 6.25585. Moving the second pattern to a compact tube pulse of duration 0.008 centered at tube proper time 6.2100 gives it 0.0624952% of the matched total energy. Both trials deliver about 0.97803744. The three matched delayed runs reproduce the target to relative error below $2\times10^{-12}$. They all remain below $\psi=0.36034$ through wall time 8. Preparation spans 0.36 to 0.37. Final release starts at 0.5265. Immediately before late emission the scalar is already at 0.359882. Withholding that same part of the early pulse and then switching the late emitter off entirely gives a final scalar of 0.36032464933, within $4\times10^{-10}$ of the matched delayed run. For these pulses the timing sensitivity comes from reducing the early input. Delivering the energy later scarcely changes the endpoint.

For each run, the archive records the time at which the scalar is largest, its direction of travel in the final release region, shell and scalar energy balances, the emitting tube mass error, and a per-mode field ledger. The zero-input control remains exactly zero. The inherited radial mass check has a larger numerical error than the local field ledgers. At 0.64% in modes 5 and 7 it improves from $5.46\times10^{-5}$ to $6.44\times10^{-6}$ on the doubled radial grid. [Results for every pulse and mesh](SOURCE.json)

## Reproduction

The capsule's replay supplies the pinned parent model and the field calibration table. To run this component separately, use Python 3 with NumPy and a C++17 compiler.

~~~sh
python source/backreaction.py --parent /path/to/replayed-parent \
  --profiles /path/to/fields/FIELDS.json --output-dir /tmp/hdblast-mean --jobs 2
~~~

The command writes the checked scalars to SOURCE.json. For runs reaching the parent study's fixed reference electrical-work marker, it also writes the birth histories and scaled modal fields in source_traces. The optional --cases argument accepts case names from the JSON for a shorter run.
