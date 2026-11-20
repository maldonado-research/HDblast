# How much expansion is still available when the particles first get hot?

Heating must finish while the closed wall still has time to expand. The coupled field equations give a necessary test that uses the wall radius and expansion rate at the beginning of each hot interval. Each preparation has its own measured starting density, bulk mass and field phase. All quantities below use the dimensionless normalization of the archived solver. The internal energy fraction refers to energy above the wall tension. [Coupled equations](../prior/decay_forward.cpp) · [Previous trajectory](../prior/DECAY.json)

## An expansion budget for the open wall

Write $a$ for the radius, $H=\dot a/a$, $e,p$ for total added wall energy density and pressure, and $\mu$ for the bulk mass parameter at the wall. The proper-time Friedmann constraint gives

$$
C\equiv a^2[1+(aH)^2]=a^4(2e+e^2)+\mu .
$$

Let $I,O$ denote the two bulk wave amplitudes at the wall, $X=I+O$ the tangential electric field and $B$ the tangential magnetic mode in the solver's convention. The bulk energy flux and normal pressure are $F=2(I^2-O^2)$ and $P_n=I^2+O^2-B^2/2$. The wall and bulk mass equations are

$$
\dot e+3H(e+p)=F,\qquad
\dot\mu=-4a^4[H P_n+(1+e)F/2].
$$

On an expanding segment their flux terms cancel in the derivative of $C$, giving the exact identity

$$
\frac{dC}{d\ln a}=2a^4\left[e-3p-e^2-3ep+B^2-X^2-(I-O)^2\right].
$$

Thus a pulse can change $e$ and $\mu$ substantially while their flux contributions to this combination cancel. Bulk normal stress still affects the evolution. The field term can increase $C$ only when magnetic energy exceeds electric energy.

The three opposite-stream pairs give isotropic stress in total. For lab-frame daughter rest density $D$, drift Lorentz factor $\gamma$, rest-frame pressure per rest energy $\Theta$, and specific enthalpy $h=5\Theta/2+\sqrt{1+9\Theta^2/4}$, the gas trace is $e_g-3p_g=D(h-4\Theta)/\gamma=D/[\gamma(h-\Theta)]$. The second equality follows from $(h-\Theta)(h-4\Theta)=1$. The parent contributes its dust density and the canonical scalar contributes $4V-\dot\psi^2$. The isotropic Maxwell field and radiation have zero trace. With $R=M/(2m)$, $y$ the decayed number fraction, and $D_*=a_{\rm init}^3 D_{\rm parent,init}/R_{\rm init}$ the total comoving daughter rest energy after all parents decay,

$$
e-3p=4V-\dot\psi^2 + D_*a^{-3}\left[R(1-y)+\frac{y}{\gamma(h-\Theta)}\right]
\le 4V+D_*\bar R/a^3,\qquad \bar R=\max(1,\sup R).
$$

Consider an uninterrupted expanding interval with $V=0$, nonnegative $\mu$, and positive $K\ge K_{\min}$. Suppose an internal energy fraction of at least $f$ is to persist. All remaining wall energies are nonnegative, and the sheet density $e_s=K(B^2+X^2)/2$ then obeys $e_s\le(1-f)e$. Here $p\ge0$ and $2a^4e\le C$. For the electric share $s_E=X^2/(B^2+X^2)$ one obtains

$$
\frac{dC}{d\ln a}\le 2D_*\bar R a+
       \frac{2(1-f)}{K}(1-2s_E)_+\, C
       \le 2D_*\bar R a+\delta C,\qquad
\delta=\frac{2(1-f)}{K_{\min}} .
$$

The field term is zero when both fields vanish. If $G=a/a_h\ge1$ and $a_h,H_h$ are the values when the interval begins, integration and $C\ge a^2$ give a necessary condition

$$
G^{2-\delta}\le 1+(a_hH_h)^2+
  \frac{2D_*\bar R}{a_h}\frac{G^{1-\delta}-1}{1-\delta}.
$$

At $\delta=1$ the last factor is $\ln G$. For a dilute parent and large positive $K$, doubling requires $a_hH_h$ close to $\sqrt{3}$. The criterion permits bulk exchange throughout. Its assumptions must continue to hold up to the proposed endpoint. Passing it gives an upper bound on growth. Demonstrating a doubling still requires a trajectory with uninterrupted heating and expansion.

## A separate finite-time check on producing the heat

Writing $E_{\rm int}$ for the density the solver calls intrinsic random energy, the charged-gas decomposition gives $E_{\rm int}\le e_g-D$. Between two expanding samples $u,s$, the charged-gas first law then gives

$$
a_s^3 E_{{\rm int},s}\le
a_u^3(e_g-D)_u+D_*\int_u^s(R-1)\,dy
+\int_u^s a^3 JX\,d\tau
-3\int_u^s H a^3p_g\,d\tau .
$$

With $Q$ the solver's charge-to-mass coefficient and $|v|\le1$, $|a^3J|\le D_*|Q|y$. The electric-work term is bounded above by $D_*|Q|\int y|X|d\tau$. Even a strong electric field has zero current capacity while $y=0$. A birth channel that opens when the field is large still needs enough time, charge density and actual work to reach the desired fraction.

On the previous primary trajectory the doubling condition closes at interpolated proper time **4.7691**, radius **5.3995**. Only **1.460%** of the sheet energy is electric there. The daughter fraction is exactly zero. The first saved point after intrinsic heating crosses 90% occurs at time 11.28 and radius 9.6884. At this point $aH=0.12854$, compared with the required **1.73199** for doubling. The upper bound on subsequent growth is **1.00823**, and integrating to turnaround gives **1.00755** from the same saved point. The ODE locates the original continuous threshold slightly earlier, at time 11.26695.

Between the first saved post-exit state and the first saved hot state the field does **4.5015** units of comoving work. Injected birth kinetic energy is **0.0000340** and pressure work costs **0.04669**. The resulting upper bound on internal energy using actual work is **4.45485**, compared with **4.40583** measured. The phase-space and thermalization assumptions behind the charged fluid remain separate microphysical questions. [Reproducible results](CRITERION.json)

## Comparison of the two new preparations

The first saved hot state provides a common starting point for the bound and measured growth. These rows use each branch's own $K$, density, radius and expansion rate. The independent identity for $e-3p$ agrees with total stress to better than $5\times10^{-17}$ in all three traces.

| Preparation | First saved hot time | $a_hH_h$ | Minimum $a_hH_h$ for doubling | Maximum factor allowed by the bound | Measured factor to turnaround |
| --- | ---: | ---: | ---: | ---: | ---: |
| Previous magnetic release | 11.28 | 0.1285 | 1.7320 | 1.0082 | 1.0076 |
| Two-stage, parity-even release | 6.80 | 3.2911 | 1.7259 | 3.4561 | 3.3232 |
| Rolling scalar coupled to the dual electromagnetic field | 6.54 | 2.5879 | 1.7320 | 2.7745 | 2.6782 |

The two-stage integrator resolves 90% slightly before its next saved row, at radius **10.75233** and $aH=3.31038$. From that event the upper factor is **3.47476** and the observed factor is **3.33937**. From its first saved hot row to turnaround, $C$ drops by **91.976**. The resolved matter trace contributes approximately **−59.186**, the quadratic density and pressure corrections **−32.732**, and $B^2-X^2$ contributes **−0.034**. The small remaining radial-stress integral is sensitive to sampling. The primary improvement is preparing the hot state with a much larger expansion rate times radius. [Two-stage run](../exit/EXIT.json) · [Dual-field run](../axion/GROWTH.json)

The gas has radiation-like pressure in this fluid closure. Throughout the sampled hot interval, $|3p_g/e_g-1|$ stays below $3.4\times10^{-9}$ in both new preparations. The rolling scalar can keep the *total* wall pressure above $e/3$. As a supplementary diagnostic, we additionally require $|3p/e-1|$ to stay within the tolerances below and report the longest run of consecutive qualifying saved outputs. These pressure thresholds were selected after the runs. The original success criterion uses the internal energy fraction alone.

| Additional total-pressure tolerance | Previous growth factor | Two-stage growth factor | Dual-field growth factor |
| --- | ---: | ---: | ---: |
| 1% | No qualifying sample | 1.043 | No qualifying sample |
| 5% | No qualifying sample | 2.353 | 1.658 |
| 10% | No qualifying sample | 3.323 | 2.378 |

Both endpoints of these supplementary intervals are saved physical states. The calculation stops at the last saved expanding state before the continuous turnaround, with no interpolation of the pressure crossing. Intervening sub-sample excursions are unresolved in this check. The strongest tolerance permits only a short interval even on the two-stage branch.

## Phase preparation and scalar work

The tangential magnetic field obeys Faraday's law during the handoff. A local threshold on an electric invariant can change which carriers participate, but a kinematically closed decay still requires an energy source. If a parent mass is made to depend directly on $z=F_{\mu\nu}F^{\mu\nu}/4$, the term $-nM(z)$ in the action also changes the electromagnetic susceptibility from $K$ to $K+nM_{,z}$ and contributes to stress. Preserving the previous Maxwell equation would omit that response.

One optional active preparation uses a scalar-dependent coefficient $\vartheta(\psi)$ multiplying $F\widetilde F$. In the sign convention $\dot{(KX)}|_\vartheta=-\dot\vartheta B$, backreaction requires $\ddot\psi|_\vartheta=+\vartheta_{,\psi}BX$. Their opposite works cancel in total wall energy. Write $b=a^2B$, $x=a^2X$ and $d\eta=d\tau/a$. The surface equations can be written

$$
b'=2x,\qquad Z=Kx+\vartheta b,\qquad
Z'=-2Kb+2\vartheta x+a^3[2(I-O)-J],
$$

where the prime differentiates with respect to $\eta$. For an impulsive change over which $b,Z$ remain nearly fixed and initial $Z/b$ is negligible, the sheet energy approaches $B^2(K+\vartheta^2/K)/2$. One can prescribe $K+\vartheta^2/K=K_{\rm target}(\psi)$ to match a scalar loading curve while moving a chosen portion to electric energy. This is a constitutive design. Actual scalar work, field rotation, retained bulk waves and particle backreaction require solving the full system. The curvature identity above applies to a metric-independent $\vartheta F\widetilde F$ density with the same kinetic stresses. The two-stage parity-even preparation sets the initial phase through its time-varying $K$.

During the actual dual-field transition, the archived 1,025-point handoff spans conformal time **0.0058005**. The comoving magnetic field changes by **−2.327%** and the displacement change is **−0.284%** of the final $\vartheta b$. Holding the starting $Z/b$ fixed in $X/B=(Z/b-\vartheta)/K$ predicts a final ratio of **−4.66664** and electric share of **95.6097%**. The coupled solution gives **−4.67988** and **95.6334%**. Those values describe the sheet on completing the scalar transition. The charged gas gives a separate, later test of energy transfer. [Dense handoff](../axion/HANDOFF.json)

Run `python -I -B criterion.py --output /tmp/criterion.json` from this directory. The script uses the Python standard library and reads the three archived primary traces and the dense dual-field handoff. The event diagnostics use saved states for energetic calculations and interpolate the geometric deadline. The output includes the effect of doubling the quadrature spacing and the measured Friedmann and matter-trace residuals. [Analysis code](criterion.py)
