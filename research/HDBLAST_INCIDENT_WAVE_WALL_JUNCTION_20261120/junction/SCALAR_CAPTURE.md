# Incident radiation can trigger a bounce if the wall stores energy in a scalar potential

A closed wall initially at its maximum radius contracts under its own gravity. In this example, an incident pulse deposits energy in a scalar field on the wall. The field moves up a potential, its pressure becomes negative, and the contracting wall reaches a minimum radius and expands. Removing the pulse or removing the potential erases the reversal. This is a conditional solution of five-dimensional Einstein junction equations with a phenomenological absorption law.

## What supplies the energy?

Use the [two-interior retarded Vaidya geometry and normalization](JUNCTION.md). Both sides have the same retained AdS interior with curvature radius $\ell=1$, closed three-sphere slices, tuned tension, and an emitting tube at $r_s=1.6$. The initial conditions are $a=2$, $H=0$, mass parameter $\mu=4$ on each side, scalar $\psi=0$ and velocity $\pi=0$. The initial scalar energy vanishes. The initial acceleration per radius is $-1/4$. This is a pre-existing wall of finite radius.

Choose a canonical homogeneous scalar with $e=\pi^2/2+V$, $p=\pi^2/2-V$, where densities are divided by the wall tension. Its potential is zero up to $\psi=0.004$, rises smoothly and reaches a positive constant $U$ at $\psi=0.012$. On the ramp use $s=(\psi-0.004)/0.008$ and $V/U=[1+\exp(1/s-1/(1-s))]^{-1}$. This potential is a time-independent choice in the model. The positive potential energy reached in the run has to be supplied by the evolving system.

In wall proper time $\tau$, set $j=Jb(\tau/d)/(d\int_0^1 b\,dx)$ with $d=0.01$ and $b(x)=\exp[4-1/(x(1-x))]$ on $(0,1)$, zero elsewhere. The phenomenological interaction directs the incident energy into ordered scalar motion with $\dot\psi=\pi$, $\dot\pi=-3H\pi-V'+j$ and absorbed power per proper volume $Q=\pi j$. The full pulse finishes while $\psi<0.004$. During the drive $(a^3\pi)^{\cdot}=a^3j\ge0$, so $Q\ge0$ exactly. After the pulse $j=Q=0$ and the scalar moves according to its free canonical equation. Its null energy condition $e+p=\pi^2\ge0$ holds throughout.

Writing $B=1+e+H$, the ordinary two-interior branch gives

$$\dot a=aH,\qquad \dot u=\frac{1}{aB},\qquad \dot\mu=-2a^4 B Q,$$

$$H^2=e(2+e)+\frac{\mu}{a^4}-\frac{1}{a^2},\qquad
\frac{\ddot a}{a}=(1+e)(1-2e-3p)-1-\frac{\mu}{a^4}-Q.$$

The model takes a compact proper-time interaction and reconstructs the outgoing pulse $\mu(u)$ via the monotone retarded clock. The tube stress $\kappa_5^2T_{uu}=-3\mu_u/(2r_s^3)$ is nonnegative. The transfer obeys $\dot e+3H(e+p)=Q$ with the entire source term supplied by the decreasing bulk mass. There is no second energy source in this accounting. A microscopic action generating this specific partition of the incident radiation remains open.

For $J=0.485$, $U=0.111$, the radiating tube stays outside the apparent horizon, $\min f(r_s)=1.9975$, and the positive mass remaining on each side is 0.24294027. Each emitting side uses 93.9265% of its initial mass parameter in a retarded interval 0.00500797. The proper energy emitted from the two tubes is 1.14726117 and the proper energy incident on the wall is 0.94084519, in units of $\lambda\,{\rm Vol}(S^3_{\rm unit})\,\ell^3$. Propagation and motion shift the locally measured energy. The first law is $\int a^3Q\,d\tau=\Delta(a^3e)+\int p\,d(a^3)$. At the pulse end the internal energy is 0.94205975 and the contraction work in this convention is $-0.00121456$. After the drive, pressure work continues while the flux vanishes.

| Event | $\tau$ | $a$ | $H$ | $e$ | $p$ |
| --- | ---: | ---: | ---: | ---: | ---: |
| Initial wall | 0 | 2 | 0 | 0 | 0 |
| Incoming pulse ends | 0.010000 | 1.998979 | −0.121428 | 0.117938 | +0.117938 |
| Field reaches the plateau | 0.044705 | 1.990252 | −0.125129 | 0.119209 | −0.102791 |
| Contraction reverses | 0.794582 | 1.899038 | 0 | 0.121878 | −0.100122 |
| Expansion reaches 2.05 | 1.752238 | 2.050000 | +0.159514 | 0.117874 | −0.104126 |

The pulse initially makes contraction faster. Reversal comes later, after energy has moved from the kinetic term into the potential and the pressure has changed sign.

## An independent turning-point test

On the plateau, $a^3\pi$ and $\mu$ are constant. Put $C=(a^3\pi)^2/2$. The Friedmann constraint then gives

$$\dot a^2=F(a)=U(2+U)a^2-1+\frac{\mu}{a^2}+\frac{2(1+U)C}{a^4}+\frac{C^2}{a^{10}}.$$

For $U,\mu>0$, $F$ is strictly convex at positive $a$ and diverges at both ends. A negative minimum gives two simple roots. A contracting solution that enters the plateau on the allowed large-radius branch reaches the larger root and has $\ddot a=F'(a)/2>0$ there. To observe this event with the source tube present the root must lie above $r_s$. For the primary run $C=0.51019709$ and $F$ has minimum $-0.13716453$ at $a=1.523972$. The larger root is $1.89903817129$, with $\ddot a=+0.316588$. The smaller mathematical root 1.26048344 lies below the tube. The larger root agrees with the time integration to $5\times10^{-12}$.

| $J$ | $U$ | Outcome before reaching $r_s$ |
| ---: | ---: | --- |
| 0 | 0.111 | Scalar remains at rest; contraction to the tube |
| 0.485 | 0 | Energy enters the scalar as kinetic energy; contraction to the tube |
| 0.485 | 0.099 | Reaches plateau, $F_{\min}=+0.017554$; contraction to the tube |
| 0.485 | 0.104 | Minimum radius 1.846382 |
| 0.480 | 0.111 | Minimum radius 1.906914 |
| 0.485 | 0.111 | Minimum radius 1.899038 |
| 0.490 | 0.111 | Minimum radius 1.889876 |
| 0.485 | 0.118 | Minimum radius 1.921560 |
| 0.485 | 0.121 | Field turns back on the ramp; contraction to the tube |

For fixed $J=0.485$ we also solve $F_{\min}(U)=0$ directly after evolving the field onto the plateau. The bracket $0.099<U<0.104$ contains a lower critical height $U_c=0.1012164656$, with a double root at $a_c=1.74137075$ above the source. On the larger-radius contracting branch at this critical point, $F$ vanishes quadratically. The integral of $da/\sqrt{F}$ diverges logarithmically, so the ideal solution approaches $a_c$ asymptotically. The finite-time reversal occurs on the side with two simple roots. This is a local threshold for the stated pulse and ramp in the part of parameter space that reaches the plateau. The scalar was evolved at two tolerances and each one-dimensional root search used absolute tolerance $10^{-12}$ in $U$. The threshold record includes their differences.

The two outer failures at fixed $J$ have different causes. The lower plateau lacks enough negative-pressure support. The higher one requires more climbing energy than the incident pulse and contraction provide. The strict minimum and simple outer root in the central case persist under sufficiently small parameter changes. These examples identify a finite range of successful configurations for this source and potential profile.

The calculation verifies the unsquared spatial junction condition, proper-time normalization, Friedmann constraint, integrated scalar energy balance and emitter fluence. On the primary trajectory all these residuals are below $5\times10^{-12}$ in the units used here. Tightening the numerical tolerance shifts the bounce radius by less than $6\times10^{-13}$. The scalar is a chosen wall degree of freedom and the absorbed energy is initially coherent. The construction gives a finite-radius classical bounce for this matter content. Establishing a hot radiation era, the origin of the emitter, the microscopic scalar interaction and a link to our universe requires additional physics. [Calculation](scalar_capture.py) · [Numerical record](JUNCTION.json)
