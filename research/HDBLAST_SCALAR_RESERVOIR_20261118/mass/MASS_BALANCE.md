# What changes in the bulk next to the shell?

A rerun of the B1 particle-production simulation gives two ways to follow the geometry inside the bulk. The metric gives a mass parameter on each planar symmetry orbit. The scalar stress gives its time and space derivatives through Einstein’s equations. In the fine run, the mass difference between the wall and the interior surface at $z=-0.04$ changes from **+2,266.21** at shell time $s=1$ to **−24,466.56** at $s=4$. Integrating the scalar stress reproduces this change along two independent paths through the same computed solution. The sign reversal concerns a quasi-local geometric quantity that includes the scalar potential. Identifying an energy source for a hot universe also requires a physical coupling, a specified history, and an independent thermal prediction. [Numerical results](MASS_BALANCE.json) · [Generalized Misner–Sharp mass and unified first law](https://arxiv.org/abs/0709.1199)

## Normalization and field identity

Take one of the two smooth copies. Its metric is $g=h_{ab}dx^a dx^b+R^2d\mathbf y^2$ with $R=e^A$, $h=e^{2B}(-dT^2+dz^2)$, and three planar orbit coordinates $\mathbf y$. Use the coordinates and normalization of the archived solution. For this fixed reference orbit metric the curvature parameter is $k=0$. Define

$$
\mu=R^2[k-(DR)^2]=R^4e^{-2B}(A_T^2-A_z^2).
$$

Here $D$ is the covariant derivative of $h$. In the usual five-dimensional Einstein normalization, the associated mass for a fiducial orbit-coordinate volume $V_3$ is $3V_3\mu/(2\kappa_5^2)$. All numbers in this note report $\mu$ with the native $A$, including its original additive constant. They refer to one smooth copy. A rescaling of planar orbit coordinates changes the convention for $R$ and the numerical normalization. The shell's comoving energy units used in the matter ledger are different.

For three symmetry directions the orbit-space components of Einstein’s tensor satisfy

$$
G_{ab}=-\frac3R D_aD_bR+\frac3R h_{ab}D^2R
+\frac3{R^2}h_{ab}\big[(DR)^2-k\big].
$$

Subtract the two-dimensional mixed trace times the identity and contract with $2R^3D_bR/3$. The second derivative of $R$ gives $-R^2D_a(DR)^2$, and the algebraic term gives $2R[k-(DR)^2]D_aR$. In the B1 normalization $G_{ab}=T_{ab}$, where $T$ is the canonical scalar stress multiplied by $\kappa_5^2$. Consequently

$$
D_a\mu=\frac23 R^3 (T_a{}^b-\delta_a{}^b T_c{}^c)D_bR .
$$

The trace here runs over $T,z$ only. The reported mass retains all of the model's scalar potential $U(\phi)$. As a normalization check, planar vacuum with Schwarzschild–AdS function $f(R)=R^2/\ell^2-m$ and $U=-6/\ell^2$ gives $\mu=m-R^4/\ell^2$. The right side is $-(4/\ell^2)R^3D_aR$, exactly its derivative. Using a cosmological-constant-subtracted convention would change both the mass and its balance. [Einstein-limit normalization in Maeda and Nozawa](https://arxiv.org/abs/0709.1199)

## Independent metric and stress calculations

Set $K=e^{-2B}(\phi_T^2+\phi_z^2)/2$ and $P=e^{-2B}\phi_T\phi_z$. The orbit-space scalar trace is $-2U$, and the two components are

$$
\partial_T\mu=\frac23R^4[(U-K)A_T+PA_z],\qquad
\partial_z\mu=\frac23R^4[(U+K)A_z-PA_T].
$$

The left sides of the numerical tests use differences of the geometric $\mu$. The right sides use potential, kinetic and mixed scalar terms and metric first derivatives. Both use one shared numerical geometry. The component involving $P$ at a fixed radius is only one contribution to the time variation. For instance, on the fine-grid wall the three integrated terms are respectively −41,102.35, −86,525.10 and +71,031.78. The signed sum is −56,595.67. Interpreting this number as energy delivered to radiation would ignore the other boundaries, the potential, and the normalization of the mass. [Complete component integrals](MASS_BALANCE.json)

At the two shell proper times $s=H_0\tau$ we use the corresponding constant-$T$ slices throughout the collar. The native conformal clock runs approximately from $T=0.8873$ to $2.4109$. The radial surfaces are fixed **coordinate** tubes with native nodes at $z=0,-0.02,-0.04$ within rounding. Their areal radii evolve. The temporal integrals use $dT$. The scalar spatial derivatives and $A_z$ at the wall are extrapolated from the five strictly interior field nodes using the derivative of the interpolating polynomial. The boundary Neumann values and ghost cells do not enter this derivative. Interior grid points use five local values. Values at the time cuts are interpolated linearly for each derived column using the sampled wall clock. We retain all native interior time and radial nodes between the endpoints. Spatial and temporal integrals use the composite trapezoid.

| Fixed surface | Coarse geometric change | Coarse integrated stress | Fine geometric change | Fine integrated stress |
| --- | ---: | ---: | ---: | ---: |
| $z=0$ | −56,569.52 | −56,570.31 | −56,594.87 | −56,595.67 |
| $z=-0.02$ | −40,567.23 | −40,567.99 | −40,587.49 | −40,588.08 |
| $z=-0.04$ | −29,845.56 | −29,846.15 | −29,862.10 | −29,862.57 |

For the fine-grid $z=-0.04$ collar, the change in the geometric wall-minus-inner contrast is **−26,732.769**. Subtracting the two temporal stress integrals gives **−26,733.107**. Subtracting the initial radial stress integral from the final one gives **−26,732.772**. The radial comparison also tests each endpoint directly.

| Fine-grid collar | Initial geometric difference | Initial radial integral | Final geometric difference | Final radial integral |
| --- | ---: | ---: | ---: | ---: |
| $[-0.04,0]$ | +2,266.212 | +2,266.198 | −24,466.557 | −24,466.574 |
| $[-0.02,0]$ | +2,358.219 | +2,358.171 | −13,649.158 | −13,649.145 |

For five-point derivatives the temporal discrepancies on the three individual surfaces are at most 20 parts per million on each grid, measured relative to the absolute geometric change. Radial endpoint discrepancies are at most 84 parts per million on the coarse grid and 21 on the fine grid. Repeating the extraction with three-point derivatives gives maximum temporal discrepancies of 0.108% and 0.0276%, and maximum radial discrepancies of 0.255% and 0.0606%, respectively. This sensitivity matters most at the extrapolated boundary. The report also gives the effect of thinning the integration nodes by a factor of two. These are numerical comparisons on a finite model domain. The sign of this quasi-local mass has no general positivity guarantee in the prescribed potential and planar convention.

## Reproduce the two paths

Run the following from this directory with Python 3.10 or later. Only the public tables and Python's standard library are needed for the balance.

```sh
python -I -B mass_balance.py > reproduced.json
cmp reproduced.json MASS_BALANCE.json
```

The time tables contain the geometric $\mu$ and three local scalar terms at each surface. The radial tables include the same geometric quantity and three radial scalar terms on every native mesh point from $z=-0.04$ to the wall at both cuts. The replay checks common endpoint values across the two tables. [Metric and stress tables](data) · [Standard-library replay](mass_balance.py)

To regenerate the measurements, supply an HDblast checkout containing the B1 producers authenticated by the observer script. The runner uses NumPy and SciPy. It configures the original numerical source and parameters, adds passive field recording at the existing diagnostic cadence, and stops after the final cut. Use fresh output directories for each grid. The example below also extracts the four public tables. A different platform can change floating-point digits.

```sh
python -I -B ../interior/run_instrumented.py /path/to/HDblast work-coarse coarse
python -I -B ../interior/run_instrumented.py /path/to/HDblast work-fine fine
python -I -B derive_mass.py work-coarse/run/collar regenerated coarse
python -I -B derive_mass.py work-fine/run/collar regenerated fine
python -I -B mass_balance.py regenerated
```

Pass `--order 3` to each extraction with a separate output directory to reproduce the lower-order differentiation check. [Numerical observer and source pins](../interior/run_instrumented.py) · [Table derivation](derive_mass.py) · [Original B1 solver](https://github.com/maldonado-research/HDblast/blob/a011dee2b490faa56ae81d3165b4b5fa5bbfe636/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py)

[Primary source for the mass convention and cosmological-constant treatment](MASS_REFERENCE.md)
