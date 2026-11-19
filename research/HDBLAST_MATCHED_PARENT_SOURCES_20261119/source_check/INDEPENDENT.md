# Matched initial sources on the fine wall

With one parent per mode up to $q=8$, matching the initial pressure and scalar force allows a change of sign in the excess pair count for each near pulse separately. The distant pulse retains its original sign for every allowed coherence on the fine wide-window calculation. Each pulse also retains its original sign in the absolute count from the occupied band. These statements concern the number per input parent, at order $y^2$, on the four prescribed fine walls.

| Maximum factor opposing the diagonal contrast | Near positive | Near negative | Distant positive |
|---|---:|---:|---:|
| Excess pairs, separable cap | 0.793389 | 0.751599 | 0.704395 |
| Excess pairs, quantum cap | 1.122021 | 1.062921 | 0.996165 |
| Excess pairs, purely imaginary coherence, separable cap | 0.704277 | 0.673260 | 0.644197 |
| Excess pairs, purely imaginary coherence, quantum cap | 0.995998 | 0.952133 | 0.911032 |
| Absolute band pairs, separable cap | 0.553819 | 0.541419 | 0.535822 |
| Absolute band pairs, quantum cap | 0.783219 | 0.765682 | 0.757767 |

The reversal factor is the negative ratio of the anomalous contribution to the signed number-diagonal contrast in the same quantity. A factor above one reverses its sign. The absolute band count includes the matched-vacuum contribution on each arm from modes with $q\le8$. Its denominator is the pulse-minus-quiet absolute number-diagonal count in that band.

## Support calculation

Write $E_i(q)=\sqrt{m_i^2+q^2/a_i^2}$, $w(q)=3q^2/8^3$, and $h(q)=(m_i/E_i,(q/8)^2m_i/E_i)$. We fix the initial mean occupation $n(q)=1$ in the band. Matching the initial scalar force and pressure is equivalent to $\int w h\,\Re c\,dq=0$ in the common real-positive instantaneous-mode convention. The mean free energy is fixed by the occupation. Subsequent source evolution may differ between preparations.

Let $D$ be the signed incoherent excess contrast and $K(q)=-2\Delta\int W F_+F_-\,d\omega/D$. The reversal factor is $\Re\int w Kc\,dq$. For any real two-vector $\lambda$, the source condition gives the upper bound $R\int w|K-\lambda\cdot h|\,dq$, where $R=1$ is the separable cap and $R=\sqrt{2}$ is the quantum positivity cap. We minimize over the two multipliers, align a trial covariance to the complex conjugate of the residual, remove any residual moments by a two-dimensional real projection, and rescale it into the allowed disk. This constructs an attainable lower value on the numerical quadrature as well. With a purely imaginary covariance the source moments vanish modewise, and the support reduces to $R\int w|\Im K|\,dq$. These are separate optimizations for each pulse. The minimizing multipliers and both numerical certificates are retained. [Full numerical results](SOURCE_TRANSFER.json)

We use the frozen canonical-transfer solver and positive-square Fourier integrals from the preceding moving-parent study, the wide proper-time window, and the numerical spectrum through $\omega=768$. We refine the conformal propagator and the response grid in $q$, spline each response and integrate with a composite Gauss rule, splitting the imaginary-only integrand at its zeros. At 128 momentum intervals the maximal quantum dual–witness gap is $1.8\times10^{-11}$ and the two normalized moment residuals are below $1.7\times10^{-16}$. Going from 64 to 128 intervals changes all six displayed bound factors by less than $6.9\times10^{-8}$. These are numerical bounds on reconstructed response curves for the stated histories.

From this directory run the following against an HDblast checkout that contains the released moving-parent study and its historical wall data. The driver verifies both released Python helper digests before loading either helper.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B source_transfer.py /path/to/HDblast > replay.json
~~~
