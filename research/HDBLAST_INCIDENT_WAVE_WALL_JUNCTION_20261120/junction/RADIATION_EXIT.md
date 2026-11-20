# A finite scalar plateau followed by radiation on the expanding wall

A finite pulse deposits energy into a canonical scalar on an initially empty, tension-tuned closed wall. In one dimensionless example the wall first contracts, bounces, expands while more than 90% of its material energy is in a positive component with pressure one-third of its density, then recollapses. The expanding interval above 90% lasts 0.320 e-folds, a 37.7% increase in radius. The plateau shape and transfer law are chosen inputs. This example provides neither a temperature nor an observationally calibrated radiation history.

## What changes after the bounce

Use the symmetric five-dimensional AdS–Vaidya junction and the finite source and absorption law in the scalar-capture calculation. Initially $a=2$, $H=0$, $\mu=4$, and $\psi=\dot\psi=e_r=0$. Wall tension is separate from material energy. The baseline pulse has integrated force $0.485$ and ends at proper time $\tau=0.01$. The scalar climbs a smooth step from $\psi=0.004$ to $0.012$, reaching a plateau of height $U=0.111$. At the bounce $\tau=0.794582$, $a=1.899038$, $\psi=0.117629$. The wall is already moving outward at $a=2.05$ when the late solver takes over. Both descriptions have exactly the same equations up to this point. [Scalar capture and absorption law](SCALAR_CAPTURE.md)

Give the plateau a finite width. Define $S(x)=0$ for $x\leq0$, $S(x)=1$ for $x\geq1$, and in between use the smooth step $S(x)=1/[1+\exp(1/x-1/(1-x))]$. The potential becomes the original uphill step minus $U S(x)$, with $x=(\psi-0.36)/0.02$. Set $\Gamma(\psi)=4 S(x)$ and keep it on after the downhill section. The force and the friction switch on smoothly. Forcing from the source has ended and the bulk mass on each side stays constant. The late equations with $\pi=\dot\psi$ are

$$
\dot\pi=-3H\pi-V_\psi-\Gamma\pi,\qquad
\dot e_r=\Gamma\pi^2-4H e_r,\qquad p_r=e_r/3 .
$$

Thus $d(a^4e_r)/d\tau=a^4\Gamma\pi^2\geq0$. The scalar loses exactly the energy that the radiation gains internally. Total material density and pressure are $e=\pi^2/2+V+e_r$ and $p=\pi^2/2-V+e_r/3$, with $e+p=\pi^2+4e_r/3\ge0$. They obey $\dot e+3H(e+p)=0$ after the pulse. The total stress sources the junction equations

$$
H^2=2e+e^2-\frac{1}{a^2}+\frac{\mu}{a^4},\qquad
\frac{\ddot a}{a}=(1+e)(1-2e-3p)-1-\frac{\mu}{a^4}.
$$

$\Gamma$ is an effective transfer coefficient. The perfect-fluid closure specifies the pressure, with homogeneity and isotropy imposed. Particle species, interactions, thermalization, the radiation spectrum, light element abundances and a present-day scale do not follow from this closure.

## How long does the radiation phase last?

| Baseline event | Wall time $\tau$ | Radius $a$ | $H$ | Wall radiation fraction $e_r/e$ |
|---|---:|---:|---:|---:|
| Finite pulse ends | 0.010000 | 1.998979 | −0.121428 | 0 |
| Bounce | 0.794582 | 1.899038 | 0 | 0 |
| Downhill begins | 3.293949 | 3.124644 | 0.368327 | 0 |
| Potential reaches zero | 3.461178 | 3.326322 | 0.364495 | 0.14725 |
| Radiation crosses 90% | 3.716587 | 3.601809 | 0.271443 | 0.90000 |
| Maximum radius | 7.080160 | 4.959739 | 0 | $>0.999999$ |
| Returns inward to integration endpoint | 11.485806 | 2.050000 | −1.267841 | $>0.999999$ |

The e-fold count is $\ln(4.959739/3.601809)=0.319917$. Here 90% refers to material energy on the wall. At the 90% crossing, the linear radiation term $2e_r$ supplies 86.04% of the sum $2e+e^2+\mu/a^4$, while the magnitude of the negative curvature term $1/a^2$ already equals 51.13% of that sum. At maximum radius curvature cancels the sum. These values explain why a radiation-rich wall can recollapse so soon. In this model the continued inward integration ends at $a=2.05$.

## Sensitivity at fixed potential and fixed transfer law

| Integrated scalar force | $\Gamma$ after exit | Limiting $\psi$ on an indefinitely extended plateau | Expanding e-folds above 90% | Subsequent outcome |
|---:|---:|---:|---:|---|
| 0 | 4 | — | 0 | No capture. Contracts to the inner source. |
| 0.483 | 4 | 0.356433 | 0 | Remains below the exit. Integrated to $a=50$. |
| 0.484 | 4 | 0.371543 | 0.591130 | Maximum radius 8.864075, followed by contraction. |
| 0.485 | 4 | 0.386645 | 0.319917 | Maximum radius 4.959739, followed by contraction. |
| 0.486 | 4 | 0.401782 | 0.204510 | Maximum radius 3.812094, followed by contraction. |
| 0.485 | 0 | 0.386645 | 0 | Bounces, converts the potential into scalar motion and recollapses. |

There is an independent test of whether the field can reach the exit. On the outward plateau $\Pi=a^3\pi$ is constant. With $z=1/a$, $e(z)=U+\Pi^2 z^6/2$ and $F(z)=2e+e^2-z^2+\mu z^4$, the limiting field for a permanently flat plateau would be

$$
\psi_\infty=\psi_j+\Pi\int_0^{1/a_j}\frac{z^2\,dz}{\sqrt{F(z)}} .
$$

The subscript $j$ denotes the outward handover at $a_j=2.05$. In these controls $H_j>0$ and $a^3 d(H^2)/da=2-4\mu/a^2-6(1+U)\Pi^2/a^4-3\Pi^4/a^{10}$ is already positive there and increases. Positive $\Pi$ then moves the field monotonically toward $\psi_\infty$. If $\psi_\infty\leq0.36$, the finite exit is never reached. If it exceeds 0.36, the original and extended plateau solutions coincide up to the first crossing. A root in the tested force bracket $[0.4832,0.4833]$ places this threshold at $0.483235835$. This is an exit threshold for the specified pulse family, distinct from the earlier condition for a bounce. The sensitivity follows from the finite distance a freely rolling scalar travels in an expanding universe.

## Reproduction and numerical checks

The calculation uses the temporal Israel equation through both reversals and checks the unsquared spatial junction and Friedmann constraint. In the baseline late integration the peak absolute Friedmann error is $1.8\times10^{-12}$, the unsquared error is below $9\times10^{-13}$, the comoving total energy balance closes within $2.6\times10^{-11}$, and the comoving radiation-transfer identity closes within $6\times10^{-13}$. Halving the maximum step and tightening the tolerance changes the radiation e-fold count by less than $6\times10^{-11}$. Radiation is nonnegative on the solver's accepted states. Dense interpolation at the exact smooth switch can return a negative value below $5\times10^{-26}$ in magnitude.

After $\psi=0.38$, the potential vanishes and $\Gamma$ becomes constant. The computed scalar momentum obeys $\pi=\pi_d(a_d/a)^3\exp[-\Gamma(\tau-\tau_d)]$ within $3\times10^{-16}$ in the baseline, where subscript $d$ denotes the end of the descent. For an independently evaluated radiation-only turning radius insert $R=a^4e_r$ at the maximum into $a^2-(2R+\mu)-(R/a^2)^2=0$. The result agrees with the measured radius to $2\times10^{-11}$. At the emitter $r=1.6$, $\mu\leq4$ gives $f=1+r^2-\mu/r^2\geq1.9975$, and $f$ increases outward for $\mu\ge0$. The source and the simulated wall stay in the untrapped region. The interior below the emitter is excised; the emitter itself is prescribed. [Calculation](radiation_exit.py) · [Machine-readable checkpoints](JUNCTION.json)
