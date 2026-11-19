# Spectator trajectories through actual first minimum passage

We propagate the inhomogeneous daughter modes continuously from the initial slice to just before and just after the actual first minimum passage on the $\phi_m=-0.5$ branch. Both seed signs are included, with the physical kinetic weight $B(\Theta)=1/f(\Theta)$ and the same-geometry $Z=1$ control. This extends the earlier [finite-band propagation](../dynamic/DYNAMIC.md), whose endpoint vacuum and field-Hamiltonian conventions are retained exactly.

The source imports the frozen dynamic routines without modifying them. From the study root:

```sh
OPENBLAS_NUM_THREADS=1 python -IB passage/passage.py --gravity gravity --output /tmp/PASSAGE.json
```

[passage.py](passage.py) · [PASSAGE.json](PASSAGE.json)

## Independent crossing and common trajectory

The gravity solution is freshly evolved until the first zero of the canonical minimum offset minus the spatial maximum of the parent displacement in the roll direction. The event grid includes the poles and equator. An independent dense spatial grid verifies the extremum at the crossing. We then continue that same solution for another 0.005 proper-time units. The symmetry-fixed first worldline is the pole for the positive seed and the equator for the negative seed.

| Seed | Independently computed $x_*$ | Before | After | Signed field distance before | Signed field distance after |
|---:|---:|---:|---:|---:|---:|
| $+0.001$ | 0.488395411319240 | 0.483395411319240 | 0.493395411319240 | +0.0109398 | −0.0108663 |
| $-0.001$ | 0.589134244832117 | 0.584134244832117 | 0.594134244832117 | +0.0107757 | −0.0106930 |

Positive distance means that no point has yet reached the minimum in the roll direction; negative distance means that part of the slice has passed it. The independent crossing times agree with the refined gravity event calculation. The daughter modes are propagated from $x=0$ once per kinetic choice and sampled at all three times $x_*-.005,x_*,x_*+.005$. We do not restart the quantum state at the crossing.

## Low-mode trajectories

The primary radial band is $n=2,4,\ldots,36$, for one fixed $\ell=2,m$ block. The following are occupations under the instantaneous original-field Hamiltonian convention. Each $N_j=\sum_i|\beta_{ji}|^2$ includes every initial mode in this finite band.

| Seed | Kinetic choice | $N_2$, before → after | $N_4$, before → after | $N_8$, before → after |
|---:|---|---:|---:|---:|
| $+0.001$ | $B=1/f$ | 0.00697954 → 0.00730092 | 0.00465449 → 0.00485116 | 0.00109930 → 0.00117418 |
| $+0.001$ | $Z=1$ | 0.00846597 → 0.00899518 | 0.00370906 → 0.00380907 | 0.000564515 → 0.000591517 |
| $-0.001$ | $B=1/f$ | 0.0185682 → 0.0194586 | 0.00734193 → 0.00842688 | 0.00817452 → 0.0102037 |
| $-0.001$ | $Z=1$ | 0.0138873 → 0.0143680 | 0.00424577 → 0.00429173 | 0.00108822 → 0.00129871 |

The physical-weight differences $N(x_*+.005)-N(x_*-.005)$ are $(0.000321372,0.000196673,0.0000748856)$ for the positive seed and $(0.000890333,0.00108495,0.00202916)$ for the negative seed, for $n=(2,4,8)$. The JSON also records the values at the crossing and the off-diagonal contribution to each low-mode occupation.

The full off-diagonal $\beta$ sum is explicitly finite-band:
$\sum_{i,j\in\{2,4,\ldots,36\},i\ne j}|\beta_{ji}|^2$, the 306 off-diagonal entries of an 18-by-18 matrix, without angular multiplicity or averaging. It changes from 0.00127059 to 0.00148080 for the positive physical-weight seed and from 0.0468293 to 0.0565224 for the negative seed. Increasing the band changes this aggregate, so the convergence claim below concerns the listed low-mode occupations, not an all-mode sum.

## Controls

The primary gravity degree 64 and 160 projection nodes are compared with degree 96 and 256 nodes. A separate run tightens the gravity time tolerance and halves its maximum step. The daughter time integration is independently refined, and its quadrature changes from 192 to 272 points. The event-time or event-distance refinement error is below $8.10\times10^{-14}$; the maximum sampled normalized Hamiltonian or physical momentum constraint is $3.52\times10^{-10}$. All sampled parent fields remain more than 0.2307 canonical field units inside the fitted domain.

Canonical commutators and Bogoliubov identities hold to $2.83\times10^{-12}$. Across the three reported low modes, every time, quadrature and gravity refinement changes occupations by less than $1.53\times10^{-14}$ in absolute units. Daughter bands 18,24,36,48 are compared. The largest low-mode difference from band 24 to 36 is $1.05\times10^{-10}$, and from 36 to 48 is $3.41\times10^{-16}$. The coarser band 18 gives a largest difference $6.62\times10^{-8}$.

These are mode trajectories over a finite time stretch that straddles a physical first-passage event. The change also includes expansion, the evolving inhomogeneous geometry and kinetic weight, and the changing endpoint Hamiltonian basis. It does not isolate excess power caused by passage through the minimum, establish a coherent parent oscillator, define unique finite-time particles, or supply a UV-summed radiation energy or backreaction.
