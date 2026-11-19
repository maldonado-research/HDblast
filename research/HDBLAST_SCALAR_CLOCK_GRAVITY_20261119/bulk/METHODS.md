# The shell and bulk response

The observable is the fractional expansion response to a weak Yukawa coupling at a fixed proper time on each wall. Write $\lambda=y^2$, $s=H_0\tau$, $h=H/H_0$, and $a=\exp(A_b-A_{\rm ref})$. The reported derivative is

$$
D_h(s)=\frac{h_{+\epsilon}(s)-h_{-\epsilon}(s)}{2\epsilon\,h_0(s)}.
$$

The subscripts indicate separate five-dimensional integrations. The reference $h_0$ belongs to the same spatial grid and pulse arm. Thus $D_h$ equals $\partial_\lambda\log h$ to the accuracy of the difference. Negative $\lambda$ is a signed numerical probe. The leading prediction at physical $y$ is $y^2D_h$. Source terms of higher microscopic order are outside the calculation.

## Dynamical equations and clock

Each integrator begins from the archived Gaussian parent and six bulk Cauchy arrays at $T=1.11$. The quiet arm uses those data as stored. The three disturbances use the archived, constraint-matched near-positive, near-negative and distant-compensator preparation. All eight initial conditions have the same native parent prescription and a vanishing added source. The two native meshes contain 4,291 and 8,581 spatial points. The code retains the full spatial chart with nominal length 16, native spatial stencils, background, tension, potential, constraint damping, fixed left nodes and fourth-order Runge–Kutta evolution. The parent source normalization is $b/r_b$, with $b=6.49\times10^{-6}$, $r_b=7.838285538\ldots$ and $H_0=1/r_b$. [Native background and Gaussian forks](https://doi.org/10.5281/zenodo.24357388)

At every stage the native Gaussian shell reevaluates its freely evolving covariances and its instantaneous subtraction on the current wall, then adds the correlated parent correction, daughter stress, interaction and chosen local completion tabulated on the appropriate zero-coupling history. For a native state $U$ and right-hand side $R_0$, the order-$\lambda$ equation is $\dot U_1=D R_0(U_0)U_1+R_1[U_0]$. Holding $R_1$ at the corresponding baseline coordinate time in the signed solves implements its inhomogeneous term. The covariance, geometry, homogeneous response and constraint-damping terms evolve at live stage values. A physical finite-$y$ model would also require higher-order stress and its response to the evolving geometry. [Finite-band energy and scalar-clock formulation](https://doi.org/10.5281/zenodo.24358059)

For each grid and pulse, the baseline map $s_0(T)$ uses cubic Hermite interpolation of the native archived pairs $(T,s)$ and the lapse $ds/dT=e^{B_b}/r_b$. A second cubic Hermite interpolant supplies $(\hat\rho_2,\hat p_2,\hat j_2,\hat\rho_{d,2})$ and its proper-time derivative. Its conformal derivative is the derivative of that same interpolant times $ds_0/dT$. In the live source totals, the scale factor and field dependence of the native free parent and vacuum subtraction use their own current time derivatives. The signed sources enter the native junctions,

$$
A_Z=\frac{e^B}{6}(\sigma+\rho),\qquad
B_Z=\frac{e^B}{6}(\sigma-2\rho-3p),\qquad
\phi_Z=-\frac{e^B}{2}(\sigma_{,\phi}+J).
$$

Here $\rho,p,J$ include the $b/r_b$ conversion. The code computes all three time-differentiated junction values. The native momentum constraint uses $\partial_T A_Z$ in its velocity stencil. The spatial evolution consumes all three junctions. The scalar and pressure work accumulators use the total hatted sources at every stage. Radiation and remaining stress are exported separately as $\hat\rho_d$ and $\hat\rho-\hat\rho_d$. The field Wy_over_H02 is the native Weyl term reconstructed from the live constraint identity.

Each solution integrates its own proper clock by the native trapezoid rule on the RK time steps, $\Delta T=0.0005$ or $0.00025$. This clock is distinct from the baseline source lookup. Readouts use cubic interpolation in that solution's own proper clock, including perturbed times for both signs. Linear interpolation provides an additional sensitivity check. The common endpoint is $s=2.10$. Every native RK record is saved until the first step at or after $s=2.105$.

## Released data

- READOUTS.csv gives 32 central differences at $\epsilon=0.01$ and 10 additional differences at $\epsilon=0.005$. The smaller differences cover both grids for all four quiet prescriptions and the distant pulse at static cutoff 120. A row gives the fractional expansion response, logarithmic scale response, planar-area and comoving-volume responses, signed differences of each observable and the actual forward finite-difference shift. The letter d in a column name means a derivative per unit $\lambda$.
- HISTORIES.csv covers nine common proper times. DIAGNOSTICS.csv includes control and perturbed near-wall constraint maxima, a preturn test, the two interpolation estimates, and the integrated energy-work residual. Maxima of normalized constraints are absolute numerical diagnostics. Interpreting the response of a signed constraint requires signed arrays at matched locations and times.
- BASELINE_CHECK.json gives the maximum absolute difference from the retained native archive at every overlapping archived sample through the simulation window. SOURCE_GRID_TEST.json records the quiet static cutoff-120 result before and after inserting native background knots into the separate source-integration mesh. Floating-point subtraction of nearby native solutions limits the final digits of the centered slopes.
- RUNS.zip contains the 92 individual full-cadence histories. Its format stores byte-shuffled arrays and reversible bitwise XOR differences from each corresponding zero run. The decoded NumPy arrays preserve every bit, including all strictly interior scalar samples for the independent flux check. The vendor/SHA256SUMS file identifies the immutable input files. The upstream solver modules and fork bytes come from the frozen HDblast source. The references ZIP contains eight original CSV members from the earlier native archive with shortened directory names. The two scale references come from the corresponding earlier rolling quiet summaries.

Run these commands from this directory with Python, NumPy and SciPy. The first group recovers the archived histories and recomputes all readout tables. The last command repeats every native integration in a separate output directory with the canonical source tables.

~~~sh
python -IB pack_runs.py unpack --runs /tmp/hdblast-live-decoded
python -IB summarize.py --runs /tmp/hdblast-live-decoded --output /tmp/hdblast-live-readouts
python -IB verify.py --runs /tmp/hdblast-live-decoded
python -IB ensemble.py --output /tmp/hdblast-live-recomputed --jobs 2
~~~

A single fine-grid signed benchmark uses `python -IB live.py fine quiet 0.01 120 static --output /tmp/hdblast-quiet`. A clean first run solves the two static vacuum shell shooting problems before evolving the postfork fields. The cache resides in the system temporary directory. Repeated runs reuse it. Changing the numerical libraries can change the last digits. The stored source hashes identify the exact driving tables used to construct the release.

## Measured expansion and controls

All 32 primary calculations give a negative expansion response at $s=2.10$. The unperturbed fine quiet wall has $h=0.17668876$ there, and continues expanding throughout the experiment. At the same proper time its scale factor decreases and its shell scalar increases relative to the zero-coupling trajectory. The entries below are fine-quiet coefficients per unit $y^2$. The relative shifts at $y=0.1$ are one hundredth of the corresponding coefficients at leading order. [All geometrical and matter readouts](READOUTS.csv)

| Finite prescription | Pair cutoff | $\partial_\lambda\log h$ | $\partial_\lambda\log a$ | $\partial_\lambda\phi_b$ |
|---|---:|---:|---:|---:|
| Bare | 60 | −0.00018639 | −0.000013880 | +0.000011352 |
| Bare | 120 | −0.00025421 | −0.000035634 | +0.000021087 |
| Static | 60 | −0.00016417 | −0.000006464 | +0.000008145 |
| Static | 120 | −0.00016149 | −0.000004752 | +0.000007815 |

For the static cutoff-120 fine wall, the near-positive pulse reduces the magnitude of the fractional expansion response relative to quiet by $6.86\times10^{-7}$ per unit $y^2$ and the near-negative arm increases its magnitude by $2.87\times10^{-6}$. These two directions persist for both meshes and all four prescriptions. The far arm gives a smaller negative contrast with static subtraction, $-1.67\times10^{-8}$ at step 0.01 on the fine grid. This changes to $-1.17\times10^{-8}$ at step 0.005 and to $+1.30\times10^{-6}$ at cutoff 120 with the bare prescription and step 0.01. Signs and magnitudes here refer to differences of fractional responses. [Matched-time contrasts](READOUTS.csv)

On the fine quiet wall, the static prescriptions have a short period of enhanced expansion during the open gate. Sampling every 0.001 in proper time gives a peak fractional coefficient of $+3.08\times10^{-5}$ near $s=1.682$ at cutoff 120 and $+1.25\times10^{-5}$ near $s=1.685$ at cutoff 60. The sign at the final readout is common to all arms. [Matched-time histories](HISTORIES.csv)

The eight zero-coupling controls reproduce the original samples to within $1.5\times10^{-13}$ in absolute units over the compared fields, sources, work and constraint norms. The preturn arrays checked below $s=1.345$ coincide exactly. For the ten individually repeated slopes, halving the signed step changes $\partial_\lambda h$ by at most $3.15\times10^{-5}$ relatively. Linear and cubic matched-time interpolants produce slopes within $9.0\times10^{-5}$ of one another relatively. On the static cutoff-120 quiet case, adding native knots to the source integration changes the fine- and coarse-grid fractional coefficients by $1.21\times10^{-9}$ and $-5.19\times10^{-10}$, respectively. A separate run generated the static background in a fresh numerical cache and recovered the coarse zero-coupling arrays bit for bit. [Native archival controls](BASELINE_CHECK.json) · [Step and constraint diagnostics](DIAGNOSTICS.csv) · [Source discretization comparison](SOURCE_GRID_TEST.json)

The integrated Ward balance also reveals the order at which the prescribed perturbation ceases to define a closed finite-strength system. On the fine quiet static cutoff-120 arm, $\Delta(a^3\hat\rho)+W_p-W_\phi$ at $s=2.10$ is $9.79\times10^{-7}$ for both signs of step 0.01. It drops to about $2.43\times10^{-7}$ at both signs of 0.005. The zero-coupling residual is $-2.06\times10^{-9}$. This scaling is consistent with an even, higher-order residual in the chosen completion of the tangent equation. The independent near-wall constraint and scalar-flux analyses test the native fields directly. [Full-step native samples](RUNS.zip)
