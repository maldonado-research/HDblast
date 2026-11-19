# Direct Fourier check and squeezed-parent counterexample

On the fine stored walls, a single isotropic squeezed state reverses the three wide-window pulse contrasts in the number and energy of outgoing pairs from the parent band. The change occurs at fixed initial occupation and fixed mean free energy. The tested phase was chosen after examining the response to exhibit a counterexample. The inputs do not predict that this state was prepared.

## Independent evolution and integration

For the pulse, quiet wall, and their common preparation segment, `transfer.py` reconstructs the conformal clock from the published proper-time wall data. It advances the canonical pair $(\chi,\Pi)$ with the exact constant-coefficient oscillator transfer at each conformal midpoint. This propagator preserves phase-space area and has second-order global accuracy. Each mode starts with real positive $\chi_i=(2a_i^3E_i)^{-1/2}$ and $\Pi_i=-i a_i^3E_i\chi_i$. All phase choices below refer to this common initial basis.

For switched mode $g=f a\chi$, we evaluate $F_\pm=\int g(\eta)e^{\pm i\omega\eta}d\eta$ directly, using Gauss–Legendre panels in $\omega$ and composite Simpson quadrature in $q$. The spin-summed Yukawa vertex weights the massless two-body phase space by $W=(\omega^2-q^2)/(8\pi^2)$ for $\omega\ge q$. An extra factor $\omega/a_f$ gives the pair energy. The full numerical band ends at $\omega=768$. For individual quiet/near-positive modes the extra energy integrated from 384 to 768 is below $7.8\times10^{-12}$ of the total on the finest run. This is a tail diagnostic. For the number integral the separate time-domain Parseval and spacelike correction check has maximum absolute discrepancy below $1.7\times10^{-8}$.

Across the point momenta $q=0,2,4,6,8$, the uniform and tapered occupations, both windows and both frequency cuts, the largest relative discrepancy in inclusive number or energy against the refined proper-time DOP853 and Fourier-grid calculation is $4.57\times10^{-8}$ on the two checked fine walls. The largest discrepancy in the near-positive pulse contrast is $1.31\times10^{-8}$ percentage point. [Transfer results](TRANSFER.json)

## Two squeezed-state experiments

Let $n=1$ for $q\le8$, with $n=0$ outside. The two-mode squeezed vacuum has $|\langle a_{\mathbf q}a_{-\mathbf q}\rangle|=\sqrt{n(n+1)}$. The common phase profile is a function of $|\mathbf q|$ and preserves translation and rotation invariance. The free instantaneous energy depends only on $n$. Other components of stress and the later response can depend on phase.

For each total pair momentum $+\mathbf q$, the physical squared amplitude in the occupied band is

$$
S_\theta=\left|\sqrt{2}\,F_-^* + e^{i\theta(q)} F_+\right|^2 .
$$

The matched-vacuum weight is $V=|F_-|^2$. In the number-diagonal reference state the excess weight is $I=|F_+|^2+|F_-|^2$. The script integrates $S_\theta$ and $V$ separately on each arm for a constant $\theta=-0.5$ and for $\theta(q)=-2.2+3(q/8)^2$ radians. It checks the covariance expansion pointwise. Physical number or energy density includes the arm-specific factor $n_c/a_f^3$. Here “absolute” denotes the switched daughter yield before the matched-vacuum subtraction, summed over the modeled band of pair total momenta $q\le8$. Contributions by parent-vacuum modes outside this band would enter an unrestricted absolute observable.

The table gives fine, wide-window contrasts at $s_f=2.10$, with the quiet number-diagonal **excess** in the same observable as denominator for every entry. Entries are $100$ times that ratio. Each row uses the same state on all four arms.

| Observable and state | Near positive | Near negative | Distant positive |
|---|---:|---:|---:|
| Excess pairs per input parent, number diagonal | −0.089192 | +0.582117 | +0.069088 |
| Excess pairs per input parent, constant phase | +0.005310 | −0.025092 | −0.001707 |
| Absolute band pairs per input parent, quadratic phase | +0.017869 | −0.104621 | −0.010874 |
| Absolute band energy per final physical volume, quadratic phase | +0.023214 | −0.110176 | −0.005929 |

The quadratic phase reverses the corresponding number-diagonal contrasts in both windows and with and without the $\omega\le60$ screen. The signs hold for the absolute band yields and matched-vacuum excesses in number and energy, measured per input parent or per final physical volume. The smallest contrast magnitude in these checks is 0.005899 percentage point. Doubling the momentum resolution from 64 to 128 intervals and the conformal evolution resolution from 16,384 to 32,768 steps shifts the quadratic-phase contrasts by at most $1.06\times10^{-7}$ percentage point. The two phase experiments use only the fine wall data.

These are differences of experiments on prescribed geometries at order $y^2$. The same initial mean energy does not impose the same initial scalar force, pressure or later backreaction. [Squeezed-state results](SQUEEZE.json)

## Sharp phase bounds

After frequency integration define $C(q)=2\sqrt{2}\int W F_+F_-\,d\omega$ and include $\omega/a_f$ in $W$ for energy. At fixed $n$, a contrast is $D+\Re\langle e^{i\theta(q)}\Delta C(q)\rangle_q$, where $D$ is its phase-independent part and the angle brackets denote the normalized $q^2dq$ average. For an absolute observable $D$ includes the matched-vacuum difference. For a final physical density both the fixed part and the cross term include the corresponding arm-specific dilution. A single constant phase has exact extrema $D\pm|\langle\Delta C\rangle_q|$. Allowing arbitrary isotropic phases gives $D\pm\langle|\Delta C|\rangle_q$, attained separately for each observable by phase alignment or antialignment at each $q$. The analytic statement applies to the continuous integrals. The reported endpoints are quadrature estimates.

For the absolute wide-window pair count, the independent sharp intervals in table units are $[-0.274298,0.018749]$ for near positive, $[-0.108256,1.724448]$ for near negative, and $[-0.012636,0.194284]$ for distant positive. The quadratic phase in the table is a single explicit state inside all three intervals. [Phase ranges and refinements](SQUEEZE.json)

## Reproduction

Use a checkout containing the two digest-pinned historical inputs named in `transfer.py`. The scripts require NumPy and SciPy with `scipy.integrate.cumulative_simpson`. From this directory run

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B transfer.py /path/to/HDblast > transfer-replay.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B squeeze.py /path/to/HDblast > squeeze-replay.json
~~~

`TRANSFER.json` includes four numerical refinements. `SQUEEZE.json` includes three, both proper-time windows, both frequency cuts and all four fine wall trajectories. The files retain absolute, vacuum-subtracted and volume-diluted coefficients in separate named fields.
