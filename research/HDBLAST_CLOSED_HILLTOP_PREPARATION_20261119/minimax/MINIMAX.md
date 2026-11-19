# How long can a quantum state keep the hilltop smooth?

Take a real canonical scalar on the stationary closed background. The field thresholds here are **absolute** amplitudes in reduced Planck units. The spatial root-mean-square (RMS) is the square root of the quantum expectation of the spatial mean of the squared inhomogeneous field. The experiment asks whether one state can keep that number below a specified tolerance for every time between the bounce and a later observation. A state can have a small spatial variance at a single time by squeezing the field.

Use real harmonics with unit spatial RMS on the unit three-sphere. Let $q_{n j}$ be the dimensionless amplitude, with $j=1,\ldots,(n+1)^2$. The nine components at $n=2$ are the first multiplet used in the constrained rolling calculation. With $x=H_{\rm top}\tau$ and $a_b=1/h_*$, the bounce commutator for each component is

$$
[q(0),q'(0)]=i\kappa,\qquad
\kappa=\frac{(H_{\rm top}/M_{\rm Pl})^2}{2\pi^2a_b^3}.
$$

Thus the propagator $q(x)=A_n(x)q(0)+B_n(x)q'(0)$ gives $[q(0),q(T)]=i\kappa B_n(T)$. Every component considered in the following all-time result is unstable throughout the interval. It obeys $q''+3h_*\tanh(h_*x)q'-F_n q=0$ with $F_n=\mu^2-n(n+2)h_*^2/\cosh^2(h_*x)>0$.

## An exact full-interval optimum

For a nonzero waiting time $T$ and a mode set $I$ with unit weights, the smallest possible peak spatial RMS satisfies

$$
\inf_\rho\;\max_{0\le x\le T}
\left[\sum_{n\in I,j}\langle q_{nj}(x)^2\rangle_\rho\right]^{1/2}
=
\left[\frac{\kappa}{2}\sum_{n\in I}(n+1)^2 B_n(T)\right]^{1/2}.
$$

The infimum ranges over states of these real-line oscillators with finite second moments. A filter multiplies each term by its squared weight.

1. Robertson's inequality gives $\Delta q(0)\,\Delta q(T)\ge s_n=\kappa B_n(T)/2$ for each component. The larger of the two endpoint sums of variances is at least their average, which is at least $\sum_{n,j}s_n$. This argument allows entanglement, inter-mode correlations and mixed states. Nonzero means also contribute to RMS.
2. For a fixed component, write $q(x)=u(x)q(0)+w(x)q(T)$. The two endpoint functions have respective boundary values $(1,0)$ and $(0,1)$. The elementary maximum principle gives $u,w\ge0$ and $u+w\le1$. A positive interior maximum of the latter above 1 would require a nonpositive second derivative precisely where the equation requires a positive one. A negative interior minimum is excluded by the reverse argument. The strictly positive solution $B_n(x)$ for $x>0$ gives uniqueness of this boundary interpolation.
3. A zero-mean Gaussian with equal endpoint variances $s_n$ and zero symmetrized endpoint covariance attains the uncertainty product. It exists since $q(0),q(T)$ are a scaled canonical pair. At intermediate times its variance is $(u^2+w^2)s_n\le s_n$. A product of these Gaussians gives equality in the formula. The worst time really is an endpoint for this optimizing state. The state can depend on the requested final time.

This is a result in the linear scalar sector. For a theory with extra quantum-gravity constraints the kinematical lower bound still holds wherever the indicated field operators and their commutator are represented. Attainment on the full physical state space would require checking its constraints. For the rolling four-dimensional background the independent reduced symplectic form gives an endpoint obstruction directly; the stationary proof of all-time attainment is specific to the equation above. [Constrained rolling calculation](../rolling/ROLLING.md)

## Numerical consequences

The table uses tolerance $0.01M_{\rm Pl}$. The required time is the earlier marginal classical launch's time to reach a displacement of $10^{-5}M_{\rm Pl}$, measured separately for each potential. It is only part of the time needed to reach radiation equality. The two latest end times are solutions of the equation above for either nine or 280 components (sum of degeneracies for $2\le n\le8$).

| Minimum | Detuning | Required time | Latest with $n=2$ | Latest with $2\le n\le8$ |
|---|---:|---:|---:|---:|
| $+0.5$ | $10^{-3}$ | 3.382 | 0.365 | 0.126 |
| $-0.5$ | $10^{-3}$ | 3.028 | 0.431 | 0.140 |
| $+0.5$ | $10^{-10}$ | 8.755 | 1.580 | 1.406 |
| $-0.5$ | $10^{-10}$ | 8.401 | 1.944 | 1.753 |
| $+0.5$ | $10^{-20}$ | 16.430 | 3.352 | 3.186 |
| $-0.5$ | $10^{-20}$ | 16.077 | 4.143 | 3.957 |

Smaller detuning reduces the bounce commutator. It also increases the required residence time. Across these cases the gain in smooth time falls short. Values obtained by propagating large formal linear amplitudes to the end of the wait are useful for checking the equations. The physical inference uses the much earlier **crossing of a small tolerance**. The study does not predict the amplitude of fully nonlinear ripples from a linear extrapolation.

## Initial kinetic-energy budget

Impose a second condition on the nine $n=2$ components. Let their summed mean kinetic energy be at most a fraction $\eta$ of $3M_{\rm Pl}^2H_{\rm top}^2$ at the bounce, so $\sum_j\langle q_j'(0)^2\rangle/6\le\eta$. Rotation or permutation averaging makes the nine diagonal covariance blocks equal while preserving summed budgets. Conversely a product of identical Gaussians realizes any admissible one-component covariance. Write $D=d^2/9$, $A=A_2(T)$, $B=B_2(T)$ and $s=\kappa B/2$. A necessary condition for both endpoints to satisfy the spatial RMS budget is $s\le D$. The least required initial velocity variance per component is

$$
V_{\rm min}=\begin{cases}
\kappa^2A^2/(4D), & (s/D)^2\le 1-A^{-2},\\
[D(1+A^2)-2A\sqrt{D^2-s^2}]/B^2, & \text{otherwise}.
\end{cases}
$$

For endpoint variances $X,Y\le D$, the maximal allowed positive covariance is $Z=\sqrt{XY-s^2}$. The resulting velocity variance is $(Y+A^2X-2AZ)/B^2$. At fixed $Y$ the optimum is $X=Y/A^2+s^2/Y$ when feasible, otherwise $X=D$ on the active upper boundary. The optimal $Y$ is $D$ in both cases since $A\ge1$. The endpoint construction also controls the entire interval by the positive interpolation above. Hence $9V_{\rm min}/6\le\eta$ is both necessary and sufficient in this linear nine-mode experiment. At $T=0$ the minimum is $\kappa^2/(4D)$ per component.

For $d=0.01$, imposing $\eta=0.01$ leaves the table's $n=2$ horizons attainable. Tightening the budget to $10^{-4}$ at the largest detuning shortens them to 0.237 and 0.291. The calculation bounds the positive initial kinetic contribution only. It places no upper bound on later stress, excludes higher modes from this budget and does not assume a renormalized total vacuum energy.

The code verifies the momentum optimization against numerical minimization, tests the endpoint interpolation by an independent backwards integration, and checks the maximum principle throughout 401 sampled times per run in addition to the proof. [Machine results for three field tolerances and two budgets](MINIMAX.json) · [Calculation](minimax.py)
