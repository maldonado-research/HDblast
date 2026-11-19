# A separable state reverses all three excess pair-count contrasts

For the archived shell histories, the same smooth phase correlation between opposite parent momenta reverses the signs of the plus, minus and far pulse contrasts. It works for both switch windows, both wall grids and both pair-energy ranges. In units with the archived dimensionless comoving momentum, one such preparation is

$$
n_q=1,\qquad c_q=\langle a_{\mathbf q}a_{-\mathbf q}\rangle
=r\exp[i(q^2/20-2.43)],\qquad 0<q<8.
$$

The parent modes have real positive field amplitude at preparation. The unoccupied modes outside this band are in their vacuum. The phase is in radians. We chose the quadratic after inspecting the coarse, wide-window unscreened far response on the baseline quadrature, rounded its coefficients, and kept that expression fixed when doubling the numerical resolution.

## What reverses

Each contrast is the pulse-minus-quiet difference in the count *in excess of the matched vacuum on each wall*, per mean initial parent in a comoving band, divided by the squared coupling. All cases use the same leading-order switched Yukawa response. The screened count integrates pair conformal energy from the kinematic threshold $q$ to 60. The full count integrates the resolved spectrum above $q$. Raw counts include the switched vacuum contribution, and proper densities also depend on the final volume. Their contrasts require their respective normalizations. [Observable and oscillator conventions](../coherence/COHERENCE.md)

For the fine wall and wide switch, the full excess counts give the following signed differences at $r=1$.

| Pulse minus quiet | Incoherent $c=0$ | Same smooth $c(q)$ in all arms |
|---|---:|---:|
| Plus | $-2.65712\times10^{-4}$ | $+3.94444\times10^{-5}$ |
| Minus | $+1.73418\times10^{-3}$ | $-1.90387\times10^{-4}$ |
| Far | $+2.05821\times10^{-4}$ | $-1.14912\times10^{-5}$ |

A classical mixture of product coherent states prepares these correlations. For each unordered opposite-momentum pair, draw independent circular complex Gaussian variables $A$ and $B$ of unit mean square. Assign coherent amplitudes $A$ and $c_q A^*+\sqrt{1-|c_q|^2}\,B$. Averaging these product states gives each mode mean occupation 1 and anomalous covariance $c_q$. The distribution is homogeneous and isotropic, and independent between distinct pairs. At $r=0.98$ the covariance is strictly within the separable disk while the smallest reversal is still 3.08% of that row's original contrast magnitude. The walls prescribe the evolution. Equal initial occupation fixes the free energy expectation in the instantaneous diagonal basis on a given wall. The initial pressure and scalar force can depend on the preparation.

## How much coherence is required?

Index by $i$ the 24 combinations of pulse, grid, window and spectral range. With $D=\int K(|F_+|^2+|F_-|^2)d\omega$, $Z=\int KF_+F_-d\omega$, let $d_i=\sum_q w_q\Delta D_{iq}$ denote the incoherent contrast. The differences are pulse minus quiet. The Gauss weights include $3q^2/8^3$. Define $V_{iq}=-2w_q\Delta Z_{iq}/d_i$. A preparation $z_q=c_q$ changes the contrast to

$$
d_i[1-m_i(z)],\qquad m_i(z)=\mathrm{Re}\sum_q V_{iq} z_q.
$$

All signs reverse precisely when $\min_i m_i>1$. Compact convex minimax gives the sharp unit-disk bound on the chosen quadrature,

$$
t=\max_{|z_q|\le1}\min_i m_i(z)
 =\min_{p_i\ge0,\;\sum_i p_i=1}\sum_q\left|\sum_i p_iV_{iq}\right|.
$$

For any simplex vector $p$, its right-hand expression is an upper bound. The aligned state $z_q=\overline{(pV)_q}/|(pV)_q|$ supplies a lower bound by direct evaluation on every row. A zero component can be set to zero. The stored simplex weights and state attain the same value to numerical precision in each group. The joint 24-row bound uses unit weight on the coarse, wide-window unscreened far contrast. If $z$ must be a single complex constant, the same argument applies after summing each row over momentum.

| Constraint on a common preparation | Baseline | Refined |
|---|---:|---:|
| Best attainable minimum $m_i$ with $|c_q|\le1$ | 1.05571238 | 1.05571236 |
| Minimum $m_i$ of the fixed quadratic phase with $r=1$ | 1.05182111 | 1.05182108 |
| Best attainable minimum $m_i$ with constant complex $c$ | 0.71979132 | 0.71979122 |
| Upper bound with constant phase and arbitrary $0\le |c_q|\le1$ | 0.725258 | 0.725166 |

The sharp refined threshold on a uniform ceiling $|c_q|\le r$ follows by linear scaling. Strict simultaneous reversal is attainable for $r>0.94722771$ and impossible below the threshold on this discrete momentum grid. Equality cancels at least one contrast. The rounded quadratic requires $r>0.95073204$. A constant phase fails at separable amplitudes even when allowed to vary the amplitude with momentum. For one fine, narrow-window screened minus row, the largest possible factor is only 0.72517. For that row the best amplitude at fixed phase enables exactly the momentum nodes with positive projected response. All distinct positive subsets needed for the optimum can be enumerated between the angles at which a projection changes sign.

The result establishes preparation dependence in a specified shell calculation. It leaves the physical mechanism that could generate such a phase open. The numerical bounds certify the displayed finite grids. The two resolutions test their stability. The refinement does not supply a continuum error bound.

## Reproduce

From the parent research package, with the repository root in the shell variable `repo_root`, run the following commands with NumPy and SciPy installed.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B phase_feasibility/feasibility.py "$repo_root" > /tmp/hdblast-phase-baseline.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B phase_feasibility/feasibility.py "$repo_root" --nodes 48 --samples 2048 --padding 32 > /tmp/hdblast-phase-refined.json
```

The baseline uses 24 Gauss momentum nodes, 1024 switch time samples and Fourier padding 16. The refined run doubles each. Both runs use the coarse and fine input walls, and both use the independent 96-node finite-interval energy integral for the screen. The largest change over the 24 fixed-phase margin factors is $3.55\times10^{-7}$. Each output includes every signed contrast, state, dual weight, and subgroup result. [Baseline](CERTIFICATES.json) · [Refinement](REFINED.json) · [Calculation](feasibility.py)
