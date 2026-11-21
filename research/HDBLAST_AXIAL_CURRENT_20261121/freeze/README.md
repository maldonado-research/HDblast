# Which particles get the collisions?

A box experiment isolates one feature of the proposed interaction. Take momenta from a live shell simulation, keep the volume fixed, switch off particle births and driving fields, and let the same charged population scatter. The registered vector resonance strongly favors some pair energies. The control gives all pair energies the same total cross section. Both branches keep opposite-charge selection, the relativistic relative-flux factor, the same energy-dependent angular law in the pair center-of-momentum frame, and exact elastic two-body kinematics. The two branches start with identical particles and random seeds. The cross-section normalization in the control sets its clock.

We stop at identical numbers of **accepted** encounters. Each event counts for both participants. An average of eight still allows a particle to have zero. The endpoints below use 8,192 positive-negative pairs sampled from each snapshot. The high-energy tenth is selected once using the initial shell-frame particle energies. The pressure column gives the largest principal pressure divided by energy density in the gas rest frame. A relativistic, directionally uniform gas approaches one third on each axis.

| Box start (shell time) | Total cross section | Count gap to a common thermal gas | Largest pressure / energy | Encounters in initial highest-energy tenth | Particles hit |
|---|---|---:|---:|---:|---:|
| 6.27855 | Vector resonance | 22.39% | 0.647 | 4.66 | 95.62% |
| 6.27855 | Constant | 1.63% | 0.337 | 8.03 | 99.97% |
| 6.67 | Vector resonance | 26.40% | 0.731 | 2.32 | 91.43% |
| 6.67 | Constant | 2.35% | 0.345 | 7.92 | 99.94% |
| 7.05 | Vector resonance | 20.49% | 0.640 | 2.49 | 85.58% |
| 7.05 | Constant | 3.66% | 0.344 | 7.98 | 99.97% |

The gap is the total-variation distance between particle counts in the inherited 20 energy-and-direction bins and a single Maxwell–Jüttner distribution matched to the combined number, energy and momentum. For instance, a gap of 20.49% means at least that fraction of the sample probability would have to move between bins to match the specified target. Differences within bins are invisible. The control is a diagnostic construction. A small gap on this grid leaves finer departures from equilibrium open.

These are additional collisions after each snapshot. The two branches use different elapsed times to reach each count. In the normalized box clock, eight encounters take about 855 units for the resonance and 16 for the constant cross section at the late start. The comparison isolates how cross-section selection affects energy sharing at fixed encounter count. A separate calculation follows the shell on its physical clock. After sixteen additional encounters at the late start, the count gaps on the finer sample are 15.08% and 1.33%, respectively.

A second subsample with 4,096 pairs and another collision seed gives the same ordering at every start. The replay checks the exact collision count and total energy and momentum at each saved checkpoint. Across twelve runs the largest reported relative energy error in a single event is below $9 \times 10^{-16}$. The JSON includes both sample sizes, normalized elapsed time, coverage, tagged counts, energy-only gaps, and pressure eigenvalues at mean counts 0, 1, 2, 4, 8 and 16.

Run the study's top-level replay to reconstruct the two radiation-history inputs and generate the live snapshots before running this test. The calculation and complete results are in [compare.py](compare.py), [box.cpp](box.cpp), and [FREEZE.json](FREEZE.json).
