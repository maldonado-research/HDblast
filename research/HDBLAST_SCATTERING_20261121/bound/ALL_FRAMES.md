# Can any single scattering frame make elastic turns sufficient?

Even if scatterers could move in **any direction at any subluminal speed**, simply turning the particles while preserving each particle's energy in that one frame leaves a measurable discrepancy from the target hot stream. In every frame, the simulated and thermal fractions below some energy cutoff differ by at least **five percentage points**. This holds separately for positive and negative daughters, so each charge could choose its own fixed frame. It holds for both interaction histories and both source costs, in the 4,096-pair and 16,384-pair samples. The cutoff may vary with the frame and charge. We give a conservative bound on the largest cumulative count difference. The program ends each branch of the search once it proves 5%. [Results for all sixteen populations](ALL_FRAMES.json)

The comparison uses the prescribed moving thermal reference and counts all particles at the common wall time of the original snapshot. A fixed-frame magnetic deflection or other elastic turn preserves each such energy count. A binary collision changes individual energies in a particular frame whenever its center of momentum moves relative to that frame. Different pairs can have different centers of momentum. The calculation bounds what a single fixed frame can explain in these recorded samples. It allows scattering during production or energy exchange among particles to have different outcomes.

## A finite check of infinitely many boosts

Work in the target stream's rest frame. For any candidate scattering-center velocity introduce a vector $\mathbf b$ with $|\mathbf b|<1$, with sign chosen so the particle energy divided by $\theta\Gamma_{\rm rel}$ is

$$
Y_{\mathbf b}=\frac{\gamma'+\mathbf b\cdot\mathbf p'}{\theta}.
$$

The same positive scale applies to the simulated and thermal populations. We enclose the closed three-dimensional unit ball inside the cube $[-1,1]^3$ and repeatedly bisect all three coordinates. A cube wholly outside the ball can be discarded. In a cube centered at $\mathbf c$ with halfwidth $h$, the energy in each particle tag lies in

$$
Y_{\mathbf c}\;\pm\;h\,(|p'_x|+|p'_y|+|p'_z|)/\theta .
$$

This gives upper and lower cumulative counts valid everywhere in the cube. If the upper simulated count lies at least 0.05 below the lower thermal count at some threshold, or the lower simulated count at least 0.05 above the upper thermal count, the entire cube passes. Otherwise we subdivide. All sixteen checks ended at depths of five to seven, with between 2,196 and 4,287 passing cubes per population. The algorithm visits all intersecting cubes. It includes the limiting directions arbitrarily close to light speed.

## Thermal bounds in a velocity cube

Put $\epsilon=1/\theta$, $q=|\mathbf p'|/\theta$, and let $\mathbf n$ be the unit momentum direction. The exact scaled thermal energy is $\sqrt{q^2+\epsilon^2}+q\mathbf b\cdot\mathbf n$. Relative to $q(1+\mathbf b\cdot\mathbf n)$ it lies higher by at most $\epsilon$. The normalized radial distribution and wall-time angular weights together differ in total variation from their massless limits by at most $(e^\epsilon-1)/2+\epsilon/8$, as derived for the [axial calculation](README.md). In the massless limit $q$ has unit-scale gamma density of shape three and the directional distribution on the sphere is $(1+\beta n_z)/(4\pi)$, where $\beta$ is the reference drift speed.

For $C=|\mathbf c|>0$ let $m=\mathbf n\cdot\mathbf c/C$. Its density after integrating azimuth is $(1+H m)/2$, with $H=\beta c_z/C$. For any velocity in the cube,

$$
1+C m-\sqrt{3}h\;\le\;1+\mathbf b\cdot\mathbf n\;\le\;1+C m+\sqrt{3}h .
$$

These limits give two radial CDF bounds that hold without selecting a particular direction inside the velocity cube. For either limit $k=L+C m$, take the integral of the massless probability above a positive threshold $x$. For $k>0$ its value at that angle is $e^{-x/k}[1+x/k+(x/k)^2/2]$. For $k\le 0$ the above-threshold probability is zero. An antiderivative in $k$ with the angular weight folded in is

$$
A(k)=e^{-x/k}\left[(C-H L)(k+x/2)+\frac H2 k(k+x)\right],\quad k>0.
$$

Set $A(k)=0$ otherwise. The angularly averaged survival probability is $[A(L+C)-A(L-C)]/(2C^2)$. For $C=0$ the angular multiplier is constant. The solver applies the nonzero-mass shift and the total-variation allowance to these upper and lower CDFs. It widens the floating-point bounds by a further $10^{-7}$. Actual checks select one of 96 thresholds spaced by 0.125 and publish the most restrictive cube, its sign of discrepancy, and its threshold.

## Independent controls and replication

An independent angular quadrature checks the closed-form primitive at 72 inputs covering small speeds, either dipole sign, a multiplier crossing zero, and the light-speed limit. An independent quadrature of the **massive** Maxwell–Jüttner density checks the reported discrepancy at physical points in each limiting cube, including its physical corners. The proof covers the deposited sample momenta and the specified thermal reference. Agreement between original and refined samples addresses the finite sampling empirically. The result is conditional on instantaneous turns after the same recorded state and the selected marker in each simulation.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B bound/all_frames.py --snapshots inputs --output output/ALL_FRAMES.json
~~~

[Adaptive boost cover and closed-form thermal bounds](all_frames.cpp) · [Independent quadratures and snapshot checks](all_frames.py)
