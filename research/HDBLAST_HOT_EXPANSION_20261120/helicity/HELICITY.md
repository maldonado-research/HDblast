# Response to inhomogeneous electromagnetic seeds

Both successful homogeneous preparations amplify some initially small electromagnetic waves. A brief scalar interaction with the electromagnetic dual field produces a pronounced dependence on the handedness of the wave. The parity-even preparation amplifies low harmonics in both polarizations. The calculation propagates infinitesimal classical disturbances on two prescribed wall histories. The absolute energy in disturbances depends on their initial state and on the subsequent response of the wall and bulk. [Scalar–dual-field transfer](HELICITY.json) · [Parity-even transfer](EVEN.json)

## Mode equation and measured gain

On the round unit three-sphere, transverse vector harmonics have curl eigenvalues $hk$ for integer $k\ge2$ and helicity $h=\pm1$. For one member of the Abelian triplet, a mode of amplitude $A$ obeys

$$
 (K A')' + (k^2 K + h k\theta') A=0,\qquad d\eta=d\tau/a.
$$

Here $K$ is the time-dependent electromagnetic coefficient, $\theta$ multiplies the parity-odd interaction, and primes denote conformal-time derivatives. The sign labeled negative uses the stated orientation and the measured positive $\theta$ jump. Reversing the sign of the jump on this *same* prescribed history interchanges the helicities. The spectator equation omits charged matter currents, electromagnetic exchange with the bulk, and perturbations of the metric and scalar. Its endpoint data and temporal interpolation come from the coupled wall solver. [Scalar and field history](../axion/HANDOFF.json) · [Parity-even wall model](../exit/growth_forward.cpp)

At each constant-coefficient endpoint use $z=(\sqrt{kK}\,A,\,K A'/\sqrt{kK})$. The comoving energy per normalized real oscillator is $k\,\lVert z\rVert^2/2$. If $z_f=Tz_i$, the largest energy gain across all initial oscillator phases is $G=\sigma_{\max}^2(T)$. The exact equations preserve $\det T=1$. The phase-averaged energy gain of a seed of fixed magnitude is $(G+1/G)/2$, so large gains also affect the phase average. A constant-$K$, constant-$\theta$ wave has $G=1$. Amplification requires an initial perturbation or a specified quantum state to assign an energy.

## Scalar–dual-field preparation

The pulse lasts 0.04661 proper wall time, or 0.00580047 conformal time, while $K$ rises from 0.04 to 6000.04 and $\theta$ rises from 0 to 28000. The control sets $\theta$ to zero on the *same measured* scale factor and $K$ history. This isolates the parity-odd term. In the negative-helicity equation $k^2 K-k\theta'$ becomes negative in part of the pulse for affected modes. The maximum sampled $|\theta'|/K$ is about 39,624. Finite-time evolution determines the actual gain.

| Unit-sphere harmonic $k$ | Negative helicity $G$ | $\theta=0$ control $G$ | Positive helicity $G$ |
| ---: | ---: | ---: | ---: |
| 2 | $3.56\times10^6$ | $1.50\times10^5$ | $3.28\times10^6$ |
| 128 | $3.90\times10^7$ | $1.47\times10^5$ | $1.50\times10^5$ |
| 512 | $9.82\times10^9$ | $1.09\times10^5$ | 22.9 |
| 4096 | $3.08\times10^{19}$ | 9.64 | 2.18 |
| 8192 | $1.63\times10^{22}$ | 1.13 | 1.53 |
| 32768 | $5.48\times10^4$ | 1.006 | 1.027 |
| 65536 | 1.0010 | 1.0002 | 1.0012 |

At $k=8192$ the physical wave number is $k/(aL)$, about $10^3/L$ during the short pulse, where $L$ is the length unit of the wall model. Applying the model at that scale requires an effective description valid there. The table samples nineteen discrete wave numbers. It establishes amplification at the listed scales while leaving the full spectrum, high-frequency tail, and feedback of the fluctuations unresolved.

There is also a check using the actual initial homogeneous magnetic and electric field components. At $k=2$, positive helicity, they give final comoving energy **3,274,991** times the initial value in the spectator evolution. The coupled wall gives **3,277,075**. The final electric shares are **95.63349%** and **95.63343%**, respectively. At larger $k$ the amplified disturbance would draw energy from the prescribed background. A computation with backreaction is needed as soon as its energy is significant.

The solver reconstructs the smooth scalar step from 1,025 measured field-level crossing rows, using the proper field velocity for cubic interpolation. It propagates $(A, K A')$ with 16,384 implicit-midpoint matrix steps. Each Cayley step preserves the Wronskian algebraically. The final 70-digit matrix products have determinant errors below $1.5\times10^{-47}$. Halving the step count changes the largest energy gains by at most 0.027% across the table and the other sampled harmonics. Halving the source sampling changes them by at most 0.021%. A separate wall replay with its Courant factor halved at the same spatial resolution changes them by at most 0.028%. Static and instantaneous shear cases have analytic checks. [Transfer calculation](helicity.py)

## Parity-even preparation and delayed release

In the second candidate the scalar raises $K$ from 0.04 to 3400.04 and later lowers it to 26.04. Propagation starts at proper time 4.0 with constant $K$ before the first preparation. It ends at 6.67, after the second step. The constant-$K$ interval between the two numerical windows spans 0.32155 conformal time. Relative phase across this interval changes the amplification of each harmonic. Both helicities have the same equation. A frozen-$K$ oscillator on the same geometry rotates without changing energy. [Complete finite-band table](EVEN.json)

| $k$ | Preparation alone $G$ | Release alone $G$ | Both in sequence $G$ |
| ---: | ---: | ---: | ---: |
| 2 | 84,948 | 130.57 | $4.88\times10^6$ |
| 4 | 84,788 | 130.56 | $1.09\times10^7$ |
| 128 | 4,549 | 123.47 | 187,410 |
| 256 | 7.48 | 104.27 | 179.13 |
| 512 | 1.18 | 52.25 | 44.44 |
| 1024 | 1.018 | 5.541 | 5.452 |
| 4096 | 1.00001 | 1.028 | 1.028 |
| 8192 | 1.00000 | 1.00016 | 1.00016 |

The single-step columns optimize their initial phases separately. The combined column propagates the same disturbance through both steps. It includes interference across the plateau. The $k=4$ row gives the largest combined gain among the nineteen sampled harmonics. The endpoint-normalized exterior waves do not require an absolute field normalization to report these ratios.

The calculation replays the original 48-element primary through the second step using its built-in dense sampling. It advances the spectator with exact constant-coefficient rotations over short layers using the midpoint $K$ and scale factor. Halving the 16,384 layers per window changes the reported gains by at most 0.0025%. Halving the source sample count changes them by at most 0.0054%. Across all three transfers, the determinant error stays below $1.7\times10^{-10}$. The small high-frequency excesses require tighter absolute precision for use in a particle count. [Two-step calculation](even.py)

From the study root, run the following commands with Python 3 and a C++17 compiler.

~~~sh
python -I -B helicity/helicity.py --output /tmp/HELICITY.json
python -I -B helicity/even.py --output /tmp/EVEN.json
~~~
