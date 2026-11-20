# The field seen by charged test matter

The neutral background is the primary finite-emitter Maxwell wall of the preceding study. The wall has an initial radius of two curvature lengths and has no added matter. The dense file records the incident pulse, bounce, scalar exit, electromagnetic era and the beginning of recollapse. A signed wall potential and its conformal-time derivative can be recovered as $A=a^2 B/2$ and $A'=a^2 X$. The electric amplitude in one of the three equally excited Abelian directions is $X/\sqrt{3}$. Changing that factor would redefine the particle charge used by the calculations.

The standalone generator compiles the original C++17 solver with only the reporting schedule replaced. It reruns the 64-element primary calculation on both sides of the wall. It extracts a sample every 0.00002 units of wall time through 0.08, then every 0.0005, plus the actual final step. The 29,558 samples include the metric, potential, electric and magnetic fields, and source quantities. The integrated bulk and wall continue to use their original timesteps. The two reported electromagnetic landmarks agree with the preceding release within $2\times10^{-7}$ absolute. The study verifier checks all 17 columns.

~~~sh
python -I -B background/generate.py --output-dir /tmp/wall-trace
~~~

[Historical source at the pinned revision](https://github.com/maldonado-research/HDblast/blob/aee9b3ed9c1f25257df2fcecbd8e47ceb2c7656a/research/HDBLAST_CURVED_MAXWELL_SCALAR_WALL_20261120/loading/forward.cpp) · [Previous release](https://doi.org/10.5281/zenodo.24387565)
