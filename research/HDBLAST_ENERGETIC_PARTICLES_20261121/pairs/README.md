# Instantaneous collision participation of the highest-energy tenth

This supplemental measurement holds each set of live dedicated-pole momenta fixed and sums over every allowed pair of one positive and one negative particle in the representative species. For each pair, sum its two contributions to the total number of collision participants per unit time. Include one contribution for each member of the pair that belongs to the highest-energy tenth. The resulting share divided by that tenth's actual particle fraction is its per-particle collision rate relative to the full population. Both comparison kernels use precisely the same input. Equal particle weights, number density and the constant cross section's arbitrary normalization cancel in the participation ratios. The energy share includes rest energy.

| Wall time | Gas energy in top tenth | Share of instantaneous participations, vector pole | Same share, constant cross section | Vector-pole rate per high-energy tag / overall mean | Constant-cross-section rate per high-energy tag / overall mean |
|---|---:|---:|---:|---:|---:|
| 6.67 | 31.12% | 0.934% | 8.704% | 0.0934 | 0.8704 |
| 7.05 | 36.37% | 1.192% | 10.010% | 0.1192 | 1.0009 |

The finest endpoint of each row contains 16,384 particles of each charge and 268,435,456 eligible pairs. At 7.05 the vector-pole rate per high-energy particle was 0.1122–0.1192 times the population mean across four particle and time-step choices. Under a constant cross section the range was 0.9988–1.0017. At 6.67 the ranges were 0.0884–0.0972 and 0.8704–0.8829, respectively. [Eight measurements](RESULTS.json) · [Ready-to-run calculation](measure.py)

[Source snapshots and replay at PR 200](https://github.com/maldonado-research/HDblast/tree/32b30bd651578e269ccd921a734bd772c1564ba2/research/HDBLAST_AXIAL_CURRENT_20261121/live) · [Registered vector-pole kernel](https://github.com/maldonado-research/HDblast/blob/32b30bd651578e269ccd921a734bd772c1564ba2/research/HDBLAST_SCATTERING_20261121/rate/README.md)

The second kernel takes a constant total cross section times the Møller flux $v_M=\sqrt{g^2-1}/(E_1E_2)$, where $g=E_1E_2-p_1\cdot p_2$. To avoid subtracting large momenta for nearby beams, the code uses $g-1=2\sinh^2[(\eta_1-\eta_2)/2]+|p_1||p_2||n_1-n_2|^2/2$ with radial rapidities $\eta=\operatorname{asinh}|p|$ and unit directions $n$. A 50-by-50 subset cross-checked against direct, extended-precision four-momentum arithmetic gave relative errors below $10^{-15}$ for both full-pair kernel means and absolute errors below $10^{-15}$ for both participant shares.

The ratios measure the instantaneous interaction preference of the sampled state. They neither count the individual particles' earlier encounters nor evolve the state with a substituted force or collision rule. The run-to-run range records finite numerical realizations. No sampling of possible partners enters each endpoint sum.

From the repository root, replace OUT with the output directory from the primary study replay and run

```sh
python -I -B research/HDBLAST_ENERGETIC_PARTICLES_20261121/pairs/measure.py OUT/live/snapshots/t6.67000*particle-pole.npz \
  OUT/live/cross-snapshots/t7.05000*particle-pole.npz --output RESULTS.json
```
