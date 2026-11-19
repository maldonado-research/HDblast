# How much high-frequency energy can any parent preparation suppress?

On the prescribed quiet wall, the archived preparation emits 0.932676 units of proper high-band pair energy per comoving volume. **At the same initial free-parent energy, the smallest attainable yield is 0.919004, a 1.466% reduction.** Without that energy budget, a squeezed preparation can reduce the quiet yield to 0.00000327169. Its initial parent energy is 83% higher than the archived preparation. The distinction matters when attributing a change in emission to the quantum state. Both comparisons use prescribed backgrounds. [Unrestricted test](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7015238836) · [Equal-budget test](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7015331720) · [Numerical results](FLOOR.json) · [Calculation](floor.py)

## Variational bound

For one parent momentum, use dimensionless canonical coordinates $r=(\sqrt I\chi_0,\Pi_0/\sqrt I)$, where $I=a_0^3\sqrt{(q/a_0)^2+m_0^2}$. The fork vacuum has covariance $V=\mathbf{1}/2$. Let $v(\omega)$ contain the complex drive amplitudes for these two initial coordinates. The positive energy weight $d\mu$ includes the parent momentum quadrature weight and $\omega(\omega^2-q^2)d\omega/(8\pi^2)$. For the interval $60<\omega<120$, define

$$
H_{ij}=\int d\mu\,v_i^* v_j,\qquad G=\operatorname{Re}H,\qquad c=-\operatorname{Im}H_{12}.
$$

The conformal ordered energy contribution of a centered parent mode is $E_q(V)=\operatorname{tr}(GV)+c$. This is the ordered two-point result in the preceding quantum calculation. A nonzero mean contributes a further nonnegative quadratic form. The uncertainty relation requires positive $V$ and $\det V\ge 1/4$. Taking the trace and determinant of $G^{1/2}VG^{1/2}$ gives

$$
E_q\ge 2\sqrt{\det G\det V}+c\ge\sqrt{\det G}+c,\qquad
V_* = \frac{\sqrt{\det G}}{2}G^{-1}.
$$

The covariance $V_*$ attains equality and is realized by a pure squeezed Gaussian. The bound uses only second moments, so the same minimum applies over all admissible parent states for these linear probes. Each native momentum can realize its own optimum, and the sum gives the total achievable floor on the prescribed wall. The code divides the conformal energy by $a(2.10)$ for its proper readout. [Ordered spectrum and the two canonical amplitudes](../quantum/SPECTRUM.md)

## Why the continuum floor is strictly positive

Since $\det G=c^2+\det H$, strict positivity follows if the complex amplitude Gram matrix $H$ is positive definite. If $H$ had a null vector, some fixed complex linear combination of $v_1(\omega),v_2(\omega)$ would vanish throughout the open high-frequency band. Both amplitudes are Fourier transforms of compactly supported conformal-time drives $aFh,aFk$, up to constant canonical normalizations. Their transforms are entire. The identity theorem then makes that combination zero at every frequency, and uniqueness of the Fourier transform makes the same combination of the original drives zero wherever $F$ is positive. The free fundamental modes $h,k$ have a nonzero Wronskian throughout the open gate, which precludes proportionality. This argument proves a positive floor for each mode on these continuous positive-weight bands. Its magnitude depends on the wall and the gate.

For a stable numerical evaluation, the verifier writes $H=M^\dagger M$, where each frequency row of the two-column amplitude matrix is multiplied by the square root of its quadrature weight. It evaluates $\det H$ from the product of the two singular values squared. Whenever $c<0$ it uses

$$
\sqrt{\det G}+c
=\frac{\det H}{\sqrt{\det G}-c}.
$$

On the quiet wall, the smaller singular value of $M$ is roughly $0.000595$–$0.000655$ of the larger one. Their near proportionality accounts for the small quiet floor. The corresponding ratios on the three pulse histories span $0.00363$–$0.00751$. The real canonical Gram condition numbers range from 5.74 to 34.6 across the four arms. A direct sum of squared ordered amplitudes and an independent two-parameter numerical minimization agree with every mode floor to $3\times10^{-11}$ relative in the reference run. Joint refinement from three to five Gaussian nodes per native time interval and from 48 to 96 frequency nodes changes each floor by less than $6\times10^{-9}$ relative.

## Four fixed-wall tests

All entries in the middle columns are leading proper pair-energy coefficients per unit comoving volume, with $y^2$ factored out. The last column compares the common quiet-optimal preparation with the preparation optimized separately for the curve in that row.

| Wall | Retained preparation | Wall-specific floor | Quiet-optimal preparation | Excess above wall-specific floor |
|---|---:|---:|---:|---:|
| Quiet | 0.932676 | 0.00000327169 | 0.00000327169 | 0% |
| Plus | 0.929598 | 0.000225644 | 0.000226653 | 0.447% |
| Minus | 0.959644 | 0.000317842 | 0.000342848 | 7.87% |
| Far | 0.936804 | 0.000204574 | 0.000212162 | 3.71% |

The quiet optimum drops the quiet emission to about 3.51 millionths of the retained value. It is a deliberately chosen preparation that uses knowledge of the later waveform. The calculation gives a kinematic minimum. It provides no physical selection mechanism for that state.

At the fork the free parent has proper Hamiltonian $\sum_q w_q E_{0q}\,\operatorname{tr}V_q/2$, where $E_{0q}=\sqrt{(q/a_0)^2+m_0^2}$. Its vacuum, retained, and quiet-optimal energies are 94.91007, 94.93201, and 174.13829 per comoving volume. The quiet-optimal state thus carries **79.22822 above vacuum** and **79.20628 above the retained preparation**. The fork energies are zeroth order in $y$. The table gives coefficients at order $y^2$, and those numbers cannot be equated to the cost of preparing the state. These comparisons prescribe the same background throughout each state variation. A calculation allowing that extra initial parent energy to gravitate would require new constraint-consistent initial data and a new evolution.

## With the same free-parent energy budget

The large initial energy cost motivated a second test after the unrestricted result was known. It constrains the total fork energy to the archived value $E_b=94.93200739$, which is only $0.02193965$ above the vacuum on the native momentum grid. Each momentum can use a different amount of that common budget. For $k_q=w_q E_{0q}/2$ the constraint is $\sum_q k_q\operatorname{tr}V_q\le E_b$. A nonnegative multiplier $\lambda$ shifts each matrix to $G_q+\lambda k_q\mathbf{1}$. The previous trace minimization gives

$$
V_q(\lambda)=\frac{\sqrt{\det(G_q+\lambda k_q\mathbf{1})}}{2}
                 (G_q+\lambda k_q\mathbf{1})^{-1},\qquad
D(\lambda)=\sum_q\left[\sqrt{\det(G_q+\lambda k_q\mathbf{1})}+c_q\right]-\lambda E_b .
$$

$D(\lambda)$ is a lower bound on any feasible conformal emitted energy. The initial energy of $V_q(\lambda)$ decreases as the multiplier increases. This follows from concavity of the minimum of affine functions of $\lambda$, whose derivative is the initial energy excess over the budget. The calculation chooses one common multiplier that reaches the budget. The feasible covariance then attains the lower bound, proving global optimality for this constraint. The code computes the direct ordered energy by summing positive frequency powers. After converting to proper energy, each direct answer agrees with its dual bound within $2.4\times10^{-15}$ in the reference run.

| Wall | Archived emission | Equal-budget minimum | Reduction from archived | Same quiet-optimal covariance |
|---|---:|---:|---:|---:|
| Quiet | 0.932676411 | 0.919004494 | 1.466% | 0.919004494 |
| Plus | 0.929598266 | 0.915995801 | 1.463% | 0.915995812 |
| Minus | 0.959644470 | 0.945689283 | 1.454% | 0.945689393 |
| Far | 0.936803624 | 0.923133876 | 1.459% | 0.923133926 |

The four registered histories share the same fork. Using the quiet-optimal equal-budget preparation unchanged on the other arms gives relative excesses of $1.16\times10^{-8}$, $1.16\times10^{-7}$ and $5.42\times10^{-8}$ above their own optima. Within this comparison, reshaping the initial covariance at fixed total energy changes the minimum yield only modestly. Large suppression in the unrestricted calculation uses a different preparation-energy regime. The constraint fixes the total initial free-parent energy on the retained finite momentum grid. Equality of that one total does not fix the full stress tensor or propagate a modified cosmology.

## Replay

From the study root, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB floor/floor.py --repo ../.. --output /tmp/FLOOR.json
~~~

The code reuses the Fourier amplitudes in the sibling quantum calculation. The JSON includes all 24 optimal covariances in the dimensionless canonical coordinates above, their weighted per-mode proper floors, conditioning checks, energy costs and the common-preparation control.
