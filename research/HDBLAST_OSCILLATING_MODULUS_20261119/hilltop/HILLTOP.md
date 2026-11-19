# How close to the hilltop must a closed bounce start?

In this four-dimensional effective model, the closed universe **can** expand far enough for decay radiation to match the oscillating scalar before the first turnaround. Its launch must start extraordinarily close to the unstable maximum of the potential. At detuning $t=10^{-3}$ the critical initial displacement is about **$4.57\times10^{-22}$** times the archive's $H_{\rm top}/(2\pi M_{\rm Pl})$ displacement on the $+0.5$ branch. On the $-0.5$ branch it is **$8.94\times10^{-17}$** times that reference. On either branch the original referenced displacement leads to recollapse before the decay can supply half the energy. Here radiation equality means $\rho_r=\rho_\sigma$. It is an energy comparison and provides no temperature measurement. [Measured threshold and controls](HILLTOP.json) · [Archived-launch evolution](../background/BACKGROUND.json) · [Curvature criterion](../curvature/CURVATURE.md)

## The time-symmetric closed launch

This is an exploratory extension of the evolution of the two archived launches. Set $K=1$, $x=H_{\rm top}\tau$, $a=H_{\rm top}a_{\rm physical}$ and $u=V_E/V_{\rm top}$. Initialize zero scalar velocity, zero expansion rate and zero radiation at the bounce. Friedmann's constraint fixes $a_i=u(\Theta_i)^{-1/2}$. Move the initial scalar only, on the side of the maximum that rolls into the chosen minimum. The perturbative width is the calibrated value for **one real massless minimal scalar**.

The code finds the maximum $\Theta_*$ of the reconstructed spline potential by solving $u'(\Theta_*)=0$. Its coordinate is within $2\times10^{-15}$ of zero for these input tables. Many of the displacements in this test are far below the precision of the original interpolant's absolute location. We represent them as a **logarithmic offset from that extremum** and work in the linearized equations until the displacement is numerically resolved. This gives a well-defined classical boundary-value experiment on the specified potential. Changes to that potential, quantum fluctuations or a different initial state may change its physical interpretation.

At the top $h_*^2=u(\Theta_*)$ and the closed background is $a(x)=\cosh(h_*x)/h_*$. Define $R(x)$ as the linear growth of the displacement with $R(0)=1$, $R'(0)=0$. It obeys
$$
R''+3h_*\tanh(h_*x) R'-\mu^2 R=0,\qquad \mu^2=-3u''(\Theta_*)>0.
$$
We evolve $\ln R$ and $R'/R$. A handoff at time $x_m$ with displacement magnitude $d_m$ then corresponds to $\ln|\delta_i|=\ln d_m-\ln R(x_m)$. The primary calculation uses $d_m=10^{-5}$. At handoff we take $a(x_m)$ from the closed solution, restore the full nonlinear potential and velocity, and take the positive expansion rate from Friedmann's constraint. For the nonlinear expanding history the square root is evaluated from the constraint at every integration stage. In a larger-displacement overlap where a direct, time-symmetric nonlinear bounce is numerically resolvable, matching the linear evolution at $2\times10^{-4}$ gives velocity differences below 0.038% and scale-factor log differences below $1.4\times10^{-5}$ across the two branches.

## Decay and a measured handoff

The scalar and spatial curvature evolve without a decay term until two consecutive oscillations pass the earlier harmonicity checks. These require at most 10% fractional error of the quadratic potential and $|H|/M\le0.01$ throughout each cycle. The measured rest-frame width then damps the scalar; the energy loss feeds free massless radiation. We continue through several more complete oscillations. At the first minimum crossing toward the original rolling direction at least $24/H_{\rm top}$ after the **first** minimum crossing, we pass the measured scale, scalar density and radiation density to the independently solved oscillation-averaged equations. Positive spatial curvature stays at its original value throughout:
$$
\mathcal M=a^3 e_\sigma,\qquad g=(\Gamma/H_{\rm top})\mathcal M,\qquad b_0=a/\mathcal M,\qquad r_0=e_r/e_\sigma ,
$$
where $e$ denotes density divided by $V_{\rm top}$. We solve for the initial displacement where radiation equality coincides with the averaged first maximum. Direct conformal-time envelope integration checks the angular-constraint solver. We varied the matching age to 12 and 36 and also matched at the next scalar turning point. The late evolution is oscillation averaged. These tests measure sensitivity to that interface; they do not give a rigorous bound for the exact scalar evolution through the eventual turnaround.

| Minimum | Detuning $t$ | $\log_{10}[|\delta_{\rm crit}|/(H_{\rm top}/2\pi M_{\rm Pl})]$ | $\ln a$ at first minimum |
|---|---:|---:|---:|
| $+0.5$ | $10^{-3}$ | −21.340 | 3.567 |
| $-0.5$ | $10^{-3}$ | −16.049 | 3.330 |
| $+0.5$ | $10^{-10}$ | −48.116 | 8.940 |
| $-0.5$ | $10^{-10}$ | −36.975 | 8.702 |
| $+0.5$ | $10^{-20}$ | −86.364 | 16.615 |
| $-0.5$ | $10^{-20}$ | −66.866 | 16.378 |

Local brackets around each quoted boundary give opposite first events. At $t=10^{-3}$ on the positive branch, the offsets with logarithms −21.330 and −21.350 give radiation-to-scalar ratios of 0.99415 and 1.00589 at the first turnaround. The negative-branch brackets are −16.040 and −16.057, with the same ratios to four significant figures. The larger initial offset rolls off earlier and reaches maximum expansion first. The smaller offset spends longer near the maximum and reaches radiation equality while still expanding. Varying $d_m$ between $10^{-6}$ and $10^{-4}$ changes the $t=10^{-3}$ log boundary by less than 0.00013 on the positive branch and 0.000017 on the negative branch. From the measured local slope of the crossover, changing the late handoff ages and phases shifts the estimated logarithmic boundary by less than 0.0006. The data also include the signs of all these comparisons on each side of the boundary. These numerical differences are specific to the matched model. Corrections to the effective potential and the use of the archived partial width are separate limitations.

An independent integration tracks the tiny offset logarithmically from the time-symmetric bounce through a resolved displacement of 0.03. At the two primary $t=10^{-3}$ thresholds it changes the scale-factor logarithm there by at most $2.22\times10^{-6}$. Its change to $\ln(g/g_{\rm critical})$ at the late reference crossing is below $6.7\times10^{-6}$. An analytic hypergeometric solution of the closed linear equation checks the growth logs directly. [Independent launch checks](CHECK.md)

For example, at $t=10^{-3}$ on the $+0.5$ critical history the decay gate opens at $x=14.254$. At the primary envelope match $x=28.985$ the actual closed solution has $a=337.76$, $e_\sigma=0.00068685$ and $e_r/e_\sigma=0.00036802$. Its $g=0.834415$ differs from the vanishing-start limit 0.835082 because both measured starting ratios enter the envelope.

## Why weaker decay requires a much finer displacement

At fixed $d_m$ and late departure the solution grows with exponent $\lambda=(-3h_*+\sqrt{9h_*^2+4\mu^2})/2$, while $a$ grows as $\exp(h_*x)$. Consequently $a$ at exit scales as $|\delta_i|^{-h_*/\lambda}$. In the weak-width limit, at fixed oscillation phase the density after exit approaches a common value, the retained curvature and radiation fractions become small, and the critical $g$ approaches its early-matter limit. This gives
$$
|\delta_{\rm crit}|\propto(\Gamma/H_{\rm top})^{\lambda/(3h_*)},\qquad
\frac{|\delta_{\rm crit}|}{H_{\rm top}/(2\pi M_{\rm Pl})}\propto t^{\lambda/(3h_*)-1/2}.
$$
In this family $\Gamma/H_{\rm top}\propto t$ and $H_{\rm top}/M_{\rm Pl}\propto\sqrt t$. The predicted exponents of the second relation are 3.82476 and 2.98909. Independent nonlinear solutions separated by ten decades in detuning, $10^{-10}$ and $10^{-20}$, give 3.82476 and 2.98909. This identifies the classical source of the fine tuning. The reference $H_{\rm top}/(2\pi M_{\rm Pl})$ is a comparison scale chosen in the archive. Assigning a quantum probability to an offset this small would require a state, a fluctuation treatment for the tachyonic hilltop and a measure over initial data.

## Reproduction

From the study root, regenerate both branches and all three detunings:
~~~sh
OPENBLAS_NUM_THREADS=1 python -IB hilltop/hilltop.py --archive archive --output /tmp/HILLTOP.json
~~~
The calculation uses the archived background tables, the nonlinear potential reconstruction and the two expanding dust-radiation solvers. The complete output records both sides of every threshold, the activation time, the measured scale and densities, all six handoff choices and the finite-displacement overlap.
