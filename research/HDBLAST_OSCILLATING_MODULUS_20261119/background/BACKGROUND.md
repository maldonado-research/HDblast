# Does the scalar reach a lasting oscillatory stage before the shell turns around?

With the archived hilltop displacement and positive spatial curvature, both expanding branches reach a maximum size before the scalar completes two nearly harmonic cycles. The first maximum occurs at $a=1.20851$ for the minimum at $\phi_m=+0.5$ and $a=1.30115$ for $\phi_m=-0.5$. The positive-minimum branch subsequently contracts, expands briefly again, and reaches a larger second maximum. With zero spatial curvature the scalar on both branches loses energy to expansion, becomes trapped near its minimum and oscillates for many cycles. [Registered comparison](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7018476842) · [Primary histories](BACKGROUND.json) · [Independent checks](CHECKS.json)

## Inputs and equations

This is the four-dimensional Einstein-frame shell-modulus effective theory. The archived modified Jordan-frame detuning is $t(1-\phi/\phi_m)^2(1+b\phi)$ with $b=c_*+2/\phi_m$ and $t=10^{-3}$. The archived flow tables give $\phi(\Theta)$ and $f(\Theta)$ for the canonical scalar $\Theta$ measured in reduced Planck units. The potential $V_E$ is the detuning divided by $f^2$. We reconstruct it and its derivative from the same twice-continuously differentiable cubic interpolants of the original nodal values. The squared zero at the minimum is retained exactly as a function of $\phi$. The archived nearest-node masses are $m/H_{\rm top}=7.68360$ and $16.54254$. The reconstructed local curvatures at the exact zero give $7.68136$ and $16.55313$. The weak width uses the archived mass and coupling, and the tests of harmonicity use the reconstructed curvature. [Archived potential implementation](../archive/blast_common.py) · [Original benchmarks](../archive/BLAST_REHEATING_RESULTS.json)

Let $x=H_{\rm top}t_E$ denote Einstein-frame proper time in hilltop Hubble units. Write $u=V_E/V_{\rm top}$, $R=\rho_r/V_{\rm top}$, $v=d\Theta/dx$, $h=H_E/H_{\rm top}$ and $V_{\rm top}=3M_{\rm Pl}^2H_{\rm top}^2$. For $\gamma=\Gamma/H_{\rm top}$ the homogeneous equations are

$$
\Theta'=v,\qquad
v'=-(3h+\gamma)v-3u_{,\Theta},\qquad
R'=-4hR+\gamma v^2/3,\qquad
(\ln a)'=h,\qquad
h'=k/a^2-v^2/2-2R.
$$

The Friedmann constraint is $h^2=v^2/6+u+R-k/a^2$. The scalar energy decreases through the decay term at the same rate that the radiation sector gains energy. This follows from differentiating $S=v^2/6+u$. All solutions start with $R=0$. The width describes a stipulated transfer to one real massless minimally coupled scalar species near a harmonic minimum. At this detuning its values are $\gamma=3.15294\times10^{-5}$ and $5.55356\times10^{-5}$. Radiation here means massless-particle energy. Its thermal distribution and interactions with ordinary matter require a separate calculation.

The archive puts $H_{\rm top}/M_{\rm Pl}=0.00881344$ and displaces $\Theta$ by $H_{\rm top}/(2\pi M_{\rm Pl})$ from zero toward the specified minimum. The flat case starts at $a=1$ with the archive's growing-mode velocity. The closed case starts at rest with $k=1$, $h=0$ and $a=u^{-1/2}$. This initial zero of $h$ is a local minimum of the scale factor. We evolve the signed expansion rate and locate its subsequent positive-to-negative zero directly.

## When to turn on the decay

The constant weak width applies to a small-amplitude oscillator. In the primary calculation the width turns on only after two successive qualifying cycles. Each cycle extends between consecutive crossings of the minimum in the original direction. Over a cycle we require $|h|/m_*\le0.01$ and $|u/u_{\rm quad}-1|\le0.10$ at both scalar turning points and across a 129-point phase grid, with $m_*=m/H_{\rm top}$ and $u_{\rm quad}=m_*^2(\Theta-\Theta_m)^2/6$. The width is then held constant. Turning it on from launch and setting it to zero for the whole evolution supply separate controls. The early always-on calculation extends an effective friction law into a regime where its microscopic rate has yet to be established.

| Flat launch | First crossing $x$ | First certifiably trapped crossing $x$ | Width activation $x$ | Matched crossing $x$ | $h/m_*$ at match | Cycle-averaged $w$ at match |
|---|---:|---:|---:|---:|---:|---:|
| $\phi_m=+0.5$ | 0.52072 | 2.26944 | 10.49061 | 51.39980 | 0.00167830 | −0.00002894 |
| $\phi_m=-0.5$ | 0.53418 | 1.72561 | 4.78273 | 50.72337 | 0.00079005 | −0.00000826 |

The trapping test requires the scalar to cross its minimum with energy below the potential at both $\phi_m-0.1$ and $\phi_m+0.1$. As long as the universe continues to expand with passive damping, its monotone scalar energy then forbids an exit through either boundary. For the average $w$, we divide the time integral of scalar pressure by the time integral of scalar energy over a full cycle. The matched state is taken at the end of the first full cycle beginning after $x=50$. Its measured radiation fractions are about 0.0906% and 0.166%, respectively. These late states, including the radiation already produced, serve as inputs to the averaged continuation.

## The closed shell and the limits of a dust estimate

| Closed launch | First maximum $x$ | First maximum $a$ | Scalar crossings before maximum | Radiation fraction there with width from launch | Archive dust estimate for maximum $a$ |
|---|---:|---:|---:|---:|---:|
| $\phi_m=+0.5$ | 1.26089 | 1.20851 | 2 | $2.5594\times10^{-5}$ | 1.24898 |
| $\phi_m=-0.5$ | 1.61629 | 1.30115 | 5 | $5.8105\times10^{-5}$ | 2.06058 |

The first two columns show the delayed-width solutions. Their radiation is zero by the prescribed switch rule. Activating friction at launch shifts $x$ at the first maximum by less than $2.4\times10^{-5}$. The positive-minimum run reaches that maximum before even one complete cycle. The negative-minimum run completes two cycles, but the latter departs from a quadratic potential by as much as 47% and its integrated $w$ is about −0.091. Replacing those early trajectories by dust loses the pressure response measured directly here.

A secondary continuation followed the two control equations beyond the first maximum. We stopped on the first occurrence of $|\phi|=0.9$, $a=a_{\rm initial}/2$, or 12 more units of $x$. For the positive minimum with zero width, $h$ returns to positive at $x=1.39189$, $a=1.20570$, then changes back to negative at $x=1.60354$, $a=1.21544$. When $\phi$ reaches +0.9 at $x=3.07531$, $a=0.91028$ and the shell is contracting again. For the negative minimum we found no intervening re-expansion before $\phi$ reached +0.9 at $x=3.44157$, $a=0.75944$. The always-on control gives the same sequences of sign changes. This finite field boundary marks the end of this follow-up integration. It provides no determination of the ultimate fate of the contracting branch.

## Verification

The primary calculation integrates the acceleration of the expansion rate independently of the Friedmann constraint. The largest absolute constraint residual is below $7.2\times10^{-11}$ across the twelve main runs. The maximum residual in the integrated total comoving energy and pressure-work identity is below $1.5\times10^{-8}$. The radiation-injection integral agrees with $a^4R$ to $8.5\times10^{-15}$ in the reference run.

A second solver uses the positive Friedmann root directly on the flat branch and evolves $a''=a(u-v^2/3-R)$ on the closed branch. It uses tighter integration tolerances and smaller steps. It also verifies the count, times and sizes of the post-maximum extrema. The undamped flat first-crossing times agree with the saved archive within $8.4\times10^{-5}$ in $x$. The two undamped closed large-displacement controls ($|\Theta_0|=0.08$) recover the archive's maximum scale factors within $4.7\times10^{-7}$. The archive used sampled-time readouts, whereas the present calculation finds zeros by event location. [Numerical differences and continuation boundaries](CHECKS.json)

From the study root, reproduce the two files with

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB background/background.py --archive archive --output /tmp/BACKGROUND.json --checks-output /tmp/BACKGROUND_CHECKS.json
~~~
