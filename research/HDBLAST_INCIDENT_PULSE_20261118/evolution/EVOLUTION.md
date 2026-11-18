# A compact bulk pulse produces a time-resolved response at the shell

We propagated matched pulse and quiet initial data in the five-dimensional B1 Einstein–scalar equations. The main scalar pulse and the changes to the initial geometry required by its constraints occupy a finite region. On the finest grid the positive pulse raises the shell field by 0.00104929 above its quiet comparison, with its peak at conformal time 0.12019. A pulse of the opposite sign lowers the field by 0.00104889. Neither approaches the B1 production crossing at 0.5. [Wall measurements](wall_histories.zip) · [Numerical comparisons](REPORT.json)

## What arrived and when

Take the identity conformal chart with the wall at $Z=0$ and increasing $Z$ directed toward it. The unperturbed configuration has $A=B_0(Z)+T$, $B=B_0(Z)$ and $\phi=\phi_0(Z)$. Both histories start with the same tuned tension, `dc=0`, no particles, and no material radiation. The primary scalar profile has height $\eta=0.001$ and support $[-0.14,-0.10]$. Its initial scalar time derivative is minus the derivative of its perturbation in $Z$. Solving both Einstein constraints supplies the accompanying initial metric. Three changes in $B_T$ behind the pulse make the exterior match the quiet solution to numerical precision at the designated cut. Their rightmost edge is $Z=-0.152$. [Numerical pulse coefficients and matching cuts](profiles.json)

The scalar and metric disturbances associated with the incident pulse first extend to $Z=-0.10$. Its leading radial null characteristic reaches the wall at $T=0.10$. Before $T=0.152$ no signal from the added metric-rate controls can reach the wall in this chart. The diagnostic for that statement is a finite-patch run with exactly the same scalar pulse and wall-side geometry, omitting the three controls. The two fine-grid shell field values differ by at most $1.14\times10^{-15}$ through $T=0.15$. At $T=0.18$ the difference is about $5.44\times10^{-5}$. The omission describes a separate finite-patch comparator with differing outer geometry. It does not supply a matched exterior.

On the finest position grid the positive response stays numerically indistinguishable from quiet through $T=0.095$. At the leading characteristic itself, the sampled difference is $5.72\times10^{-8}$, down from $1.30\times10^{-5}$ at the coarsest resolution. At $T=0.105$ the fine-grid difference is $+0.00029456$. At the peak the largest actual shell scalar is $0.00166608$. The figures below give the maximum magnitude of the *difference from a paired quiet wall*. They use common conformal time. The proper shell clock in units of the initial Hubble time is $s=\int e^{B_b}dT/r_b$. For the fine positive pulse its peak node has $s\simeq0.1201823$. The initial geometry, incoming wave, and choice of time interval concern this local experiment. [Archived physical model](https://github.com/maldonado-research/HDblast/tree/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source)

| Spatial step | Positive peak response | Conformal peak time | Quiet scalar drift, maximum | Maximum normalized Hamiltonian residual, pulse | Maximum normalized momentum residual, pulse |
|---:|---:|---:|---:|---:|---:|
| 0.001 | 0.00105138 | 0.12100 | 3.80e-5 | 2.37e-4 | 1.68e-3 |
| 0.0005 | 0.00104969 | 0.12000 | 9.59e-6 | 9.24e-5 | 8.56e-4 |
| 0.00025 | 0.00104928 | 0.120125 | 2.41e-6 | 3.55e-5 | 3.65e-4 |
| 0.000125 | 0.00104929 | 0.1201875 | 6.03e-7 | 9.81e-6 | 1.11e-4 |

The maximum differences between the positive response curves on successive common time nodes are $5.13\times10^{-5}$, $2.21\times10^{-5}$ and $7.20\times10^{-6}$. The small change in the height of the sampled peak understates the discrepancy in the full time history. Residuals and grid comparisons measure numerical sensitivity. A rigorous time-continuous error enclosure remains open.

For signed pulses let $r_+$ and $r_-$ denote their changes from quiet. On $0.10\le T\le0.15$ the ratio of the supremum of $|(r_++r_-)/2|$ to the supremum of $|(r_+-r_-)/2|$ is $0.000784$ on the fine grid. The independently matched half-amplitude runs give $0.000392$. Doubling either half-amplitude response reproduces its full-amplitude response over this window to about 0.0392% of the corresponding full peak in maximum norm. These values characterize the nearly linear response of the scalar field on this interval. The instantaneous flux and integrated tension response involve other quantities. [Half-amplitude initial profiles](half_profiles.json)

## Additional matched data

A second example keeps $B_T=0$ initially. Four scalar bumps supply the profile, with their amplitudes chosen to match the geometry at both cuts. Its finest paired response peaks at $+0.00009581$ at $T=0.125$. The leading bump first has causal access to the wall at $T=0.015$. Here the spatial changes in the initial metric and in $A_T$ still follow from the constraints.

A third example preserves the initial areal radius and its expansion and uses a mass-matching condition. We evolved its computed near-root $\eta=0.00202094949502$ with a dedicated quiet reference. The wall radius reconstructed from the quiet first integral is 7.838285538074238. The value read from the earlier simulation summary is 7.838285538073243. The generator uses the reconstructed value in both runs of this pair. The remaining mass mismatch at the computed root is about $8.75\times10^{-17}$. The two finest grids give peak shell changes of $0.002164816$ and $0.002164832$, both at $T=0.1205$. The fine grid reaches an actual shell scalar of $0.00278162$. These are finite numerical evolutions of the computed near-root. [Mass-matched propagation wrapper](evolve_matched_mass.py)

## Equations and numerical domain

We use the registered sextic bulk potential $U(\phi)$ and the registered scalar-dependent tension $\sigma(\phi)$. Write $V=e^{2B}U$ and $K=A_T^2-A_Z^2$. In each smooth interior the equations are

$$
\begin{aligned}
 A_{TT}-A_{ZZ}&=-3K+2V/3,\\
 B_{TT}-B_{ZZ}&=3K-(\phi_T^2-\phi_Z^2)/2-V/3,\\
 \phi_{TT}-\phi_{ZZ}&=-3(A_T\phi_T-A_Z\phi_Z)-e^{2B}U_\phi .
\end{aligned}
$$

Before the first production event the three wall conditions give $A_Z=B_Z=e^B\sigma/6$ and $\phi_Z=-e^B\sigma_\phi/2$. The numerical branch starts below $\phi_*=0.5$ and checks the field at every time step. A crossing would require evolving the original B1 matter system and its corresponding boundary terms. All the reported runs stay below the crossing, so the B1 cohort-creation channel is inactive over the tested interval. [Original equations, constraints and junction code](https://github.com/maldonado-research/HDblast/blob/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py)

The position scheme uses centered second differences in $Z$, centered positions in $T$ and three local corrections to time-center the velocity-dependent force. The initial step uses the initial velocity and acceleration. At the wall, a cubic Hermite formula uses the prescribed derivative and the nearest two interior values for the second derivative. We evolve the full fields in each run. The quiet drift in the table therefore tests the discretization on a nonconstant static spatial profile.

The initial patch is $[-0.46,0]$. Each step discards the leftmost point, so the update at every retained point has all its numerical ancestors in the initial patch. The time step is one half the spatial step. Through $T=0.18$ the shell and fixed near-wall constraint monitor remain inside the retained domain. The scheme requires no artificial far-end condition. It propagates the unmatched initial tails without smoothing them. The same rule applies to the finite-patch comparators. At every time sample, both original constraints are evaluated on interior nodes at $-0.05\le Z<0$, excluding the wall. Each is divided pointwise by the sum of the absolute terms in the original constraint and the maximum over those nodes is recorded.

A separate implementation transports $\partial_T u\pm\partial_Z u$ along the two null characteristics, for each field $u=A,B,\phi$. It integrates the semilinear force using the two endpoint values along the ray and sets the incoming characteristic at the wall from the outgoing characteristic and the actual Neumann data. This uses the same initial fields and a different discretization of the evolution and boundary. The finest positive peak is $0.0010492626$ at $T=0.120125$. The maximum discrepancy from the position scheme across their common time nodes is $3.10\times10^{-6}$. In this characteristic scheme the two spatially separated compensation cases coincide through $T=0.15$ at all three tested spacings. Neither method supplies an independent physical realization. [Characteristic program](check_characteristics.py) · [Cross-check results](CHARACTERISTICS.json)

## Reproduction and table schema

The ZIP contains 46 CSV files and ten associated metadata files. The seven cases in `runs/g1`, `g2`, `g4` and `g8` use spatial steps of $0.001$, $0.0005$, $0.00025$ and $0.000125$. `delayed` and `delayed_negative` are the primary matched pair, each compared with the same `quiet`. `uncompensated` and `uncompensated_negative` are the corresponding finite-patch cases without the compensating metric controls. `scalar_only` is the four-bump matched case. `previous_pulse` is the earlier one-bump example whose metric failed to match the same quiet exterior. The two cases in each `mass_runs` directory use their own quiet normalization. The five cases in each `half_runs` directory include quiet, two matched signs and their two uncompensated equivalents, at the two finest spacings.

Each CSV gives the coordinate time, $A_b,B_b,\phi_b,A_{T,b},B_{T,b},\phi_{T,b}$, six scalar samples strictly inside the bulk at $Z=-j\Delta Z$ for $j=1,\ldots,6$, and the two normalized constraint diagnostics. Output cadence is $\Delta T=\Delta Z/2$. Every CSV has a SHA-256 listed in its directory's metadata. The compact report checks these data and can be reproduced with the Python standard library.

~~~sh
python -I -B report_evolution.py wall_histories.zip > replay.json
cmp REPORT.json replay.json
~~~

For field regeneration install NumPy and SciPy and supply the [pinned B1 summary](https://raw.githubusercontent.com/maldonado-research/HDblast/07317e899a7b8ecb3367add9d6e5a27d6e777aea/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/runs/main/main_ps0.5_G100_y1_bcut_dzf5e-4_summary.json). The source digest is enforced by the solver. In the following examples save the summary as `B1_summary.json` and replace the mass module location with the supplied `mass_match.py`. Repeat using the other published spacings for the full comparison.

~~~sh
python -I -B evolve_patch.py B1_summary.json profiles.json regenerated/runs/g8 --spacing 0.000125
python -I -B evolve_patch.py B1_summary.json half_profiles.json regenerated/half_runs/g8 --spacing 0.000125
python -I -B evolve_matched_mass.py B1_summary.json /path/to/mass_match.py regenerated/mass_runs/g8 --spacing 0.000125
python -I -B check_characteristics.py B1_summary.json profiles.json regenerated/runs > characteristic-replay.json
python -I -B report_evolution.py regenerated/runs > complete-replay.json
~~~

The characteristic comparison needs the three finest primary directories. The combined report needs all ten run directories. Full numerical regeneration on other platforms may differ in the last floating-point digits. [Primary evolution](evolve_patch.py) · [Standard-library replay](report_evolution.py) · [Packed histories](wall_histories.zip)
