# How much magnetic energy can become heat before another doubling?

The closed wall in the preceding calculation prepared a largely magnetic field. The field had to become electric before it could do work on the newly charged particles. Both the field rotation and the growth of the wall depend on conformal time. Comparing their rates gives a limit that applies even if the electric energy could be extracted instantly. [Numerical output](PHASE.json) · [Archived coupled solution](../prior/DECAY.json)

## A bound on electromagnetic energy release

On a unit three-sphere the invariant magnetic mode obeys $b'=-2y$, where $b=a^2B$, $y=-a^2X$ and a prime denotes differentiation with respect to $d\eta=d\tau/a$. With constant electromagnetic coefficient $K$, zero current and no bulk driving, the other Maxwell equation gives $y'=2b$. Thus the vector $(b,y)$ rotates by $2\Delta\eta$ and its comoving energy $E_{\rm EM}=K(b^2+y^2)/2$ is constant. A purely magnetic field needs $\Delta\eta=\pi/4$ to become fully electric.

Let $R_0^2=b_0^2+y_0^2$, choose the sign of the field so that $b_0\ge0$, and let $s=y_0^2/R_0^2$ be the initial electric energy fraction. Even with a particle current, as long as $K$ is constant and $E_{\rm EM}(\eta)\le E_{\rm EM,0}$ throughout the interval, Faraday's law gives $|b'|\le2\sqrt{R_0^2-b^2}$. Starting from $\alpha_0=\arcsin\sqrt{s}$, the earliest possible reduction in $b$ therefore obeys $b/R_0\ge\cos(\alpha_0+2\Delta\eta)$ up to the first quarter turn. The residual field energy includes at least this magnetic component. The net fractional decrease $f$ in the initial sheet energy satisfies

$$
 f\le \sin^2\!\left[\min\!\left(\frac{\pi}{2},\alpha_0+2\Delta\eta\right)\right].
$$

This upper bound uses only the magnetic evolution and the ceiling on field energy. A freely rotating field with favorable initial sign can attain it in the idealized limit of complete, instantaneous extraction of its electric component at the measurement time. Finite currents, delayed charge creation and reversible exchange can give smaller net transfer. A reduction in the sheet energy measures energy made available to other degrees of freedom. An internal heat fraction of the whole wall also depends on the other energy sources and sinks.

## Conformal time in a closed radiation era

Take an expanding closed wall with constant bulk mass $m_b$, zero late vacuum energy and conserved total radiation $E=a^4e>0$. This includes any sheet energy transferred to another component with radiation pressure. The brane Friedmann equation gives

$$
 (a')^2=2E+m_b-a^2+E^2/a^4.
$$

Write $A$ for the first turnaround and $u=a/A$, $\epsilon=E^2/A^6$. Since $A^2=2E+m_b+E^2/A^4$, substitution $u=\sin\theta$ gives the exact advancing phase between two radii,

$$
 \Phi(u_0,u_1)=2\int_{\arcsin u_0}^{\arcsin u_1}
 \frac{\sin^2\theta\,d\theta}{\sqrt{\sin^4\theta+\epsilon(1+\sin^2\theta)}}
 \le2(\arcsin u_1-\arcsin u_0).
$$

For a subsequent growth factor $G$ during expansion, the conversion must occur by $u_1=1/G$. If $u_0>1/G$, the wall has already passed that radius. Otherwise a target net release $f$ requires an initial electric share of at least

$$
 s_{\min}=\sin^2\!\left[\max\!\left(0,\arcsin\sqrt{f}-\Phi(u_0,1/G)\right)\right].
$$

For a pure magnetic start, reaching $f$ leaves at most $1/\sin[(\arcsin\sqrt f)/2]$ subsequent growth even as $u_0$ and the high-density correction approach zero. A complete quarter turn leaves at most a factor of **1.4142**, and 90% release leaves **1.7102**. A subsequent doubling permits at most **75%** release from such a magnetic start. Even in the most favorable limiting geometry a 90% release before a doubling needs an electric seed of at least **4.0192%**. A larger starting radius or positive brane correction reduces the available phase.

## What happened on the previous trajectory

The archived wall completed scalar exit at radius **4.92072** and turned at **9.76159**. Half its maximum radius was **4.88080**, already behind it at exit. At the first fully post-exit sample, proper time 4.52, only **0.06140%** of the sheet was electric. Integrating $2\,d\tau/a$ from that sample to the decay threshold predicts a field angle of **1.53742** radians, compared with **1.53712** from the measured fields. Doubling the trace spacing changes the phase integral by $8.3\times10^{-6}$ radians. Through the entire pre-threshold interval the angle discrepancy stays below $0.00030$ radians. The full evolution includes exchange with the bulk.

| Event in the coupled solution | Radius | Remaining growth factor to turnaround |
| --- | ---: | ---: |
| Scalar exit | 4.92072 | 1.98377 |
| Electric share first reaches 90% after exit, interpolated | 9.08156 | 1.07488 |
| Decay channel opens, 99.887% of the sheet electric | 9.51647 | 1.02576 |
| Internal particle energy reaches 90% of added wall energy | 9.68672 | 1.00773 |

From scalar exit to maximum size, the geometric phase budget is 1.94565 radians. Only 0.43027 radians remain after the decay threshold. The actual charge current strongly changes the field after that threshold.

For a separate all-radiation continuation from the first post-exit sample, assign the entire measured density $e=0.0991781$ at $a=4.93445$ to radiation, hold $m_b=0.350410$ fixed and turn off bulk exchange. The new turnaround is **10.87178**. The preparation occurs at $u_0=0.453876$, leaving only **0.102449 radians** before half of this new maximum radius. A 90% sheet release by that point requires an electric seed of **83.0593%** in the favorable phase. With the measured seed, the ceiling is **1.6100%** net release. Two quadratures, one in $a$ and one in $\theta$, agree within $2\times10^{-15}$ radians. This continuation assumes a particular stress history. Direct heating by the scalar, an electrically prepared sheet before the radiation phase, later field energization or an additional source of expansion changes the applicable test. A rolling scalar coupled to the electromagnetic dual field can prepare electric energy during exit. Its work and stress must be included when assessing those initial conditions.

## Checking the two electrically prepared cases

Both new preparations emerge from scalar exit earlier in the expansion and with much larger electric shares. At their respective first archived post-exit samples, the actual radii are 28.77% of turnaround in the parity-even case and 27.73% in the parity-odd case. For each, the table holds that sample's total density and bulk mass fixed and repeats the separate all-radiation continuation. The exact seed bound applies to sheet energy measured from the specified sample.

| Primary case | Proper time | Matched maximum radius | Available phase before half the matched maximum | Minimum electric seed for 90% release | Measured electric share |
| --- | ---: | ---: | ---: | ---: | ---: |
| Parity-even preparation and release | 6.66 | 36.82085 | 0.47566 | 48.7984% | 99.9993% |
| Transient parity-odd interaction | 5.62 | 31.02257 | 0.50990 | 45.3815% | 95.7370% |

Each prepared state clears the necessary seed fraction for the matched comparison. The actual uninterrupted growth is determined by the fully coupled trajectories. [Parity-even study](../exit/EXIT.md) · [Parity-odd study](../axion/AXION.md)

Run `python -I -B phase.py --output /tmp/phase.json` from this directory to reproduce the numbers with the Python standard library.
