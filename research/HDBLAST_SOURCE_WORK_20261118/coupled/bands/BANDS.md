# Extending the gravitating momentum band

The original preparation follows one real scalar on $0\le k\le8$. We repeated all four coupled, decaying cases with upper edges of 12 and 14. These wider-band calculations were chosen after the nominal results were available. Each uses the same preparation at $T=0.70$, couplings, pulse construction and two bulk meshes. The wider domains use 32 Gauss–Legendre nodes. The modes, pressure, reciprocal scalar force and daughters of the additional momenta all contribute during the evolution. Each different domain defines a different finite-band matter model. [Nominal experiment](../COUPLED.md) · [Wider-band readouts](BAND_RESULTS.json)

| Fine-grid band at $s=2.10$ | Quiet $\rho_R/H_0^4$ | Near positive change | Near negative change | Distant positive change | Quiet daughter term / $H^2$ |
|---|---:|---:|---:|---:|---:|
| $0\le k\le8$ | 37.1576 | −0.1714% | +1.0167% | +0.0994% | 0.22512% |
| $0\le k\le12$ | 56.0699 | −0.2554% | +1.4051% | +0.1137% | 0.33922% |
| $0\le k\le14$ | 59.7506 | −0.2826% | +1.5272% | +0.1168% | 0.36142% |

Counting and coupling the wider band materially raises the daughter density and fraction. The sign of each contrast persists. The largest native daughter-radiation term divided by $H^2$ over all 16 wider-band runs and the recorded times through $s=2.10$ is about **0.367%**, on the coarse near-negative run at the readout. This is a model-specific expansion contribution for these finite preparations. The three edges do not extrapolate to the continuum, establish a thermal bath or limit other species.

We checked the physical frequency of the outermost member of each band along every history. On each saved segment, with $A_b$ and $\phi_b$ interpolated linearly in proper time, the squared frequency $E_K^2=K^2\exp[-2(A_b-A_{*,0})]+10^4(\phi_b-0.5)^2$ is convex. Its maximum therefore occurs at an endpoint, including the interpolated cut. Across both grids the interpolant is bounded above by $51.888H_0$ for $K=12$ and $53.085H_0$ for $K=14$, versus $M_5=53.6107H_0$. The latter is about 99.0% of the chosen gravitational scale, so omitted short-distance physics is particularly relevant. The bound applies to the prescribed interpolation and says nothing about rapid unsampled bulk dynamics. The native geometric decomposition and the instantaneous subtraction used for stress and force are identical to the nominal calculation. [Samples and independent frequency check](BAND_SAMPLES.zip) · [Verification](verify_bands.py)

Run the compact check with standard Python. Regenerating the native experiments uses the same NumPy/SciPy environment as the parent study. The adapter changes the quadrature endpoint and its frequency diagnostic. The original native driver and quantum sector supply the rest of the equations.

~~~sh
python -I -B verify_bands.py
python -I -B run_band.py /path/to/repository /tmp/b1-wider --edge 12
python -I -B run_band.py /path/to/repository /tmp/b1-wider --edge 14
python -I -B band_results.py /tmp/b1-wider /path/to/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118
~~~
