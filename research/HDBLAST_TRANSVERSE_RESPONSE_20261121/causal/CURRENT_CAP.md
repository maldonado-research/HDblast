# Can one field pattern turn the injected particles quickly enough?

Give one initially empty electromagnetic Fourier mode the largest current the charged particles can supply. Let the particles bunch as strongly as they like and move at any speed up to the speed of light. On the two prescribed histories even this allowance is too small to produce the archived early thermal spread with the original homogeneous bound on collision rates. This test covers nonlinear trajectories in the selected mode. It requires the recorded expansion and total birth density and excludes an imposed incoming transverse wave.

| At 1% of whole-run absolute electrical work, finer wall time step | Dedicated pole | Shared pole |
|---|---:|---:|
| Maximum single-mode RMS broadening allowance over particles already born | $311.33\,mc$ | $546.76\,mc$ |
| Minimum required by the thermal transverse distribution under the preceding bounded-hazard model | $541.65\,mc$ | $674.10\,mc$ |
| Minimum fraction requiring a collision or another deflection mechanism for that thermal distribution, given the single-mode force cap alone | 96.16% | 79.03% |
| Largest distinct fraction the inherited homogeneous, any-target collision bound allows to have been hit | 66.24% | 43.52% |
| Minimum gap in the cumulative transverse-momentum distribution under both assumptions | 29.92 points | 35.51 points |
| Minimum number of Fourier pairs needed to get within **one percentage point** at all transverse momenta at the interpolated marker | 9 | 4 |

The minimum counts permit every mode to attain its own maximum current and all forces to add in the same direction along every trajectory. Meeting those counts does not establish a sufficient field. The cap on electromagnetic broadening and the minimum required hit fraction in the third row need only the global rest inventory. Combining them with a maximum binary-collision count uses an additional assumption about the local target density. Particle clumping can change the binary rate. [Machine-readable values, both wall step sizes](CURRENT_CAP.json) · [Density-dependent multimode test](../global/clumping.py) · [Earlier thermal-impulse derivation](https://github.com/maldonado-research/HDblast/blob/0d73635c128efb988f92a591bf07ddb19ee731a3/research/HDBLAST_KINETIC_TEST_20261121/ballistic/KICK.md)

## Maxwell response to one real Fourier pair

Use a periodic local wall patch with the usual linear constitutive relation of constant $K=2.54$. The free Maxwell equations propagate comoving fields $a^2(\mathbf E,\mathbf B)$ in conformal time $d\eta=dt/a$ by free rotations for transverse polarizations. The longitudinal electric polarization accumulates its current. Its homogeneous external axial drive is already included in the ballistic background. The extra field is initially zero and receives only the included particle currents. The total comoving daughter rest weight is $D$, in the same effective current normalization as the parent wall, $|J|\le QD/a^3$ for the spatially averaged absolute current. Taking this full normalization for each mode is generous. A real Fourier harmonic includes both complex coefficients $\mathbf k$ and $-\mathbf k$. Our convention for either coefficient is a spatial average with $\exp(-i\mathbf k\cdot\mathbf x)$.

Let $\mathbf n$ point along a tag's momentum perpendicular to its field axis $\hat z$, and $V$ denote its axial velocity. Only the component of force along $\mathbf n$ can increase the magnitude of transverse momentum. The axial guide can rotate the transverse momentum at arbitrary speed. To bound radial broadening, consider one source element with velocity $\mathbf u$, and set $c=\mathbf n\cdot\hat k$, $d=\hat z\cdot\hat k$ and $\delta=k(\eta-\eta')$. Write $\mathbf n_L=c\hat k$, $\mathbf n_T=\mathbf n-\mathbf n_L$, and let $\phi$ denote the spatial phase from source to tag. Apart from signs and the source normalization, the contribution of this real conjugate pair to the radial Lorentz force is

$$
2\mathbf u\cdot(\mathbf A\cos\phi+\mathbf B\sin\phi),\qquad
\mathbf A=\mathbf n_L+\cos\delta\,\mathbf n_T,\qquad
\mathbf B=V\sin\delta(c\hat z-d\mathbf n).
$$

Here the first term comes from the electric field, including its longitudinal component, and the second from the magnetic field. The transverse part of the tag velocity is parallel to $\mathbf n$ and drops out of the radial magnetic force. The operator taking $(\cos\phi,\sin\phi)$ to the vector in parentheses has norm at most one. To see this, form its two-column Gram matrix. Its first diagonal entry is at most one with deficit $\sin^2\delta(1-c^2)$. The second is at most one since $c^2+d^2\le1$. The largest norm occurs at $|V|=1$ since reducing $|V|$ contracts the second column. At $|V|=1$ the determinant of the identity minus the Gram matrix is

$$
\sin^2\delta\,(1-c^2-d^2)(1-\sin^2\delta\,c^2)\ge0.
$$

This gives a pointwise radial bound even with a charge-density perturbation, arbitrary wave orientation and particles bunching at different phases. At zero transverse momentum one can use any transverse unit vector in the same bound. Summing the real particle current elements and integrating their forcing gives the radial-force allowance per unit charge

$$
F_{\mathrm{cap}}(s)=\frac{2Q}{K a(s)^2}\int_{t_0}^{s}\frac{D(r)}{a(r)}dr.
$$

For total particle orbits the currents entering this equation may vary nonlinearly. The bound uses the fixed global particle inventory and the speed limit at each time. It concerns the radial projection of the Lorentz force. The full vector impulse can also rotate particles inside a region of small transverse momentum.

## Cohort and collision tests

For a tag born at $t_b$, the integrated radial-broadening allowance at time $T$ in mass units is $J_b=Q\int_{t_b}^{T}a(s)F_{\mathrm{cap}}(s)ds/a(T)$. An unhit tag obeys $p_\perp/m\le b+J_b$, $b=\sqrt{1.05^2-1}$. For $M$ nonzero Fourier pairs the triangle inequality grants $b+M J_b$. This holds regardless of whether harmonics occupy multiple axes or have phase-correlated currents. Each Fourier pair is counted once; no independent-random-phase assumption is used.

If the probability of remaining unhit is at least $f_0$ and the transverse distribution must be thermal with CDF $F$, the smallest RMS radial broadening compatible with that target grants all hit tags free redistribution, and assigns the lowest-momentum fraction $f_0$ to the unhit tags. Its square is at least $\int_0^{f_0}(F^{-1}(u)-b)_+^2du$. This is the same quantile integral evaluated in the previous study and applies directly to the broadening allowance. The table compares that lower bound with $\sqrt{\int J_b^2dD/D(T)}$.

For an independent CDF test, the fraction of all births with cap below a selected momentum $p$ is $G(p)=\Pr(b+M J_b\le p)$. At least $G(p)-F(p)$ of all particles must have collided or received another kick if their final transverse distribution equals $F$. Maximizing over cuts gives the third row with $M=1$, independently of any proposed rate of collisions. When the homogeneous binary hazard bound used in the parent study applies, replace each birth weight in $G$ by its lower probability of remaining unhit. The resulting function is a compulsory lower bound on the actual CDF. The script records the largest excess over the thermal CDF, and the first mode count for which the lower bound permits a stated tolerance. Exact-CDF counts at zero tolerance are especially sensitive to the latest cells of continuous injection. The quoted one-percentage-point counts agree at both wall resolutions when the interpolated marker is used. An additional test treats the marker as lying anywhere in its final recorded time cell. It grants the full cell’s field growth and collision exposure, excludes all births in that cell from the CDF floor, and bounds all normalizations by their less restrictive endpoints. The required counts are then at least **eight and four** on the finer grids, and **eight and three** on the original grids. In all four cases the single-mode gap remains greater than 29 percentage points.

In each time cell the source is bounded above by the largest endpoint $D$ divided by the smallest endpoint $a$. The free-force integral uses that envelope and grants every birth a whole interval from the earliest permitted time, which makes the result conservative for monotone interpolated histories. The survival test additionally maximizes the target density and gives each cohort the entire cell's collision exposure. Piecewise-linear trapezoidal checks give RMS allowances of $311.02$ and $546.31$ on the finer histories. These are slightly below the upper enclosures in the table. The real-mode Gram identity is also checked on 4,900 synthetic source orientations and phases. [Calculation](current_cap.py)

Run from the study directory with `python -I -B causal/current_cap.py --history history --output /tmp/CURRENT_CAP.json`.
