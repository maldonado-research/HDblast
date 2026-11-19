# Finite velocity at the closed homogeneous bounce

A classical initial velocity shifts the field value that can linger near the hilltop. At detuning 0.001, choosing a compensating displacement from the linear equations misses the required interval even at a velocity of 0.001. The center curves in opposite directions for the two calibrated potentials. The one-sided width changes much less than the center across the velocities sampled here. [Results](MANIFOLD.json) · [Radiation event checks](EVENTS.json)

## What is held fixed

Write $q=\Theta-\Theta_*$, $x=H_{\rm top}t$, $p=dq/dx$, $h=H/H_{\rm top}$, $l=\log a$ and $U=V_E/V_{\rm top}$. The bounce has zero radiation, unit positive curvature, $p(0)=v$ and $h(0)=0$. On this surface the Einstein constraint fixes $a_b=[U(q_0)+v^2/6]^{-1/2}$. Before the harmonic decay gate the scalar and scale factor obey

$$
q'=p,\quad p'=-3hp-3U_q,\quad h'=U-p^2/3-h^2,\quad l'=h .
$$

These equations propagate the Friedmann constraint through $h=0$. Each tested launch has positive $h'(0)$. The maximum velocity magnitude $\sqrt{0.06}$ gives initial scalar kinetic energy equal to 1% of the reference hilltop energy. Here it specifies classical examples. A bound on a quantum mean kinetic energy permits a velocity distribution with tails beyond the sampled range.

## The stable curve

For fixed $v$ a boundary-value solve picks $s(v)=q(0)$ such that the scalar approaches the stationary solution. Define $h_*=\sqrt{U(0)}$, $\mu^2=-3U''(0)$ and $\lambda_-=[-3h_*-\sqrt{9h_*^2+4\mu^2}]/2$. At $x=5$ the numerical boundary condition is $p=\lambda_-q$. Repeating at $x=6$ at tighter tolerance changes the selected centers by at most $5.4\times10^{-15}$ in four sampled cases per branch. This accuracy is sufficient to display the bend. Resolving the width requires a calculation in relative coordinates.

At the origin set $q=vS+v^2 T+O(v^3)$. The geometry changes at second order in $v$. Define $L y=y''+3h_*\tanh(h_*x)y'-\mu^2y$. Then $LS=0$ and $LT=-3U'''(0)S^2/2$ with $S'(0)=1$, $T'(0)=0$ and decaying endpoints. This gives $s(v)=-Bv+Cv^2+O(v^3)$ with $C=T(0)$. In the interpolated potential the stationary point is within $2\times10^{-15}$ of a spline knot. Adjacent cubics give $U'''$ values differing by 0.00108 or less. The finite-difference coefficients at $|v|=0.001$ agree with the origin calculations to 9 parts per million or better. The numbers below refer to the specified interpolated potential.

| Potential minimum | $B$ | $C$ | One-sided zero-velocity distance $d_0$ | $\sqrt{d_0/\lvert C\rvert}$ |
|---|---:|---:|---:|---:|
| $+0.5$ | 0.0692117492 | −0.008259871 | $6.410\times10^{-25}$ | $8.81\times10^{-12}$ |
| $-0.5$ | 0.0837801575 | +0.001549225 | $1.254\times10^{-19}$ | $9.00\times10^{-9}$ |

The last column estimates the speed at which the quadratic shift reaches the archived zero-velocity displacement. At $v=0.001$ the linear prediction misses the center by $8.26\times10^{-9}$ and $1.55\times10^{-9}$, respectively. Even the quadratic approximation misses by about $1.8\times10^{-12}$ and $1.2\times10^{-12}$. This distinction matters when assigning a precision to a putative launch quadrature.

## Width relative to the exact center

Let $\epsilon=q_0-s(v)$ be a signed displacement at fixed velocity on the bounce surface. The previous study found a threshold on the side rolling toward each calibrated minimum. That side has $\epsilon<0$ for the positive minimum and $\epsilon>0$ for the negative minimum. The exactly centered solution approaches the top indefinitely and produces no oscillating parent in this model. Outcomes on the other side require their own potential and endpoint tests.

We evolved the deviations from the stable solution, divided by $\epsilon$, from the bounce to the same expansion radius as the critical reference history. The initial normalized displacement is one and the normalized initial momentum and Hubble-rate perturbations are zero. Differences of potential and force are evaluated by integrating their derivatives between the two fields. The difference in initial log-radius uses the Friedmann constraint and a logarithm accurate for small arguments. This represents $\epsilon$ independently even when both unshifted field values would round to the same number.

We adjusted $|\epsilon|=d(v)$ until the outgoing field at the archived matching radius equaled the value obtained by propagating the zero-velocity critical seed through the same nonlinear equations. The matching field magnitude is about $10^{-5}$. At all sampled velocities the matching momentum agrees to better than $2\times10^{-10}$ relative. The Friedmann residual along the unperturbed central trajectory is below $9\times10^{-12}$. Direct integrations at a resolvable field separation of $2\times10^{-7}$ check the normalized deviations to relative accuracy $1.3\times10^{-7}$ or better at $x=0.55$.

| Potential minimum | $v$ | $s(v)$ | $d(v)/d_0$ |
|---|---:|---:|---:|
| $+0.5$ | −0.244949 | +0.01648275 | 0.881678 |
| $+0.5$ | +0.244949 | −0.01747742 | 1.112196 |
| $-0.5$ | −0.244949 | +0.02063332 | 1.006762 |
| $-0.5$ | +0.244949 | −0.02044595 | 0.971028 |

As a downstream check at both listed speeds on each branch we actually continued two launches per case through the first minimum, the two-cycle gate, and the first minimum-crossing at least 24 units later. From there the published averaged decay equations determine the first turnaround and any earlier equality. At $0.94d(v)$ radiation equals matter first. At $1.06d(v)$ turnaround comes first. For the positive minimum the radiation/matter ratios at the first maximum are 1.01569 and 0.98548. For the negative minimum they are 1.01948 and 0.98204. Similar results at the two velocity signs follow from matching the outgoing histories to the same field and radius. The tabulation samples velocities and brackets the event locally. It supplies no measure on initial conditions. [Velocity sample and numerical checks](MANIFOLD.json) · [Eight event outcomes](EVENTS.json)

## Direct launch below machine precision

For $v=-0.001$ we also shot in the absolute field coordinate using 50-decimal arithmetic. A Taylor integrator evolves the coupled field and geometry with a local expansion of the two neighboring spline cells. We move the decaying endpoint from $x=1.3$ to $2.25$ and then $2.85$ and solve for the center at each stage. We form the displaced field by adding the signed threshold separation in decimal arithmetic and integrate it directly to the target time. Repeating with 20 and 24 steps per unit time, time-series orders 28 and 32, and local field polynomials through orders 10 and 12 tests discretization and truncation. The refined runs reproduce the outgoing field to within $3.5\times10^{-10}$ relative and the target log-radius to within $1.1\times10^{-12}$. These digits test a fixed numerical potential. The uncertainty in the physical hilltop position is much larger than the formal width. [Direct shoots and refinement](PRECISION.json)

## Reproduce

From the study directory with NumPy and SciPy installed, run

~~~sh
python -I -B manifold/manifold.py --prior prior --output /tmp/manifold.json
python -I -B manifold/precision.py --prior prior --manifold /tmp/manifold.json --output /tmp/precision.json
python -I -B manifold/events.py --prior prior --manifold /tmp/manifold.json --output /tmp/events.json
~~~

The extraction module rebuilds the cubic interpolant from the bundled archived field and Planck-factor data. The main output records SHA-256 hashes of the four prior inputs. The event check additionally uses the supplied prior hilltop and curvature solvers.
