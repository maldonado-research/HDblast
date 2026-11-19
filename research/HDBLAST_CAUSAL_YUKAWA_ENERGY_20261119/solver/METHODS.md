# Causal energy and work on the retained quiet wall

Each number below is the coefficient of $y^2$ in an energy per unit comoving volume, in units of $H_0^4$. A negative parent entry means that coupling the parent to the daughters lowers its energy relative to evolving the same initial parent covariance on the same wall with $y=0$. The native wall and scalar mass history are external inputs throughout this calculation.

## Finite spectral model and initial state

The clock is proper time $s=H_0\tau$. The native trajectories specify $a=\exp(A_b-A_{\rm ref})$, signed mass $m=100(\phi_b-1/2)$, and their proper derivatives. We use cubic Hermite interpolation of $\log a$ and $m$. The quiet retained wall is used for each preparation. At its fork the parent has native $X=\langle\chi^2\rangle$, $C=\langle\{\chi,\Pi\}/2\rangle$ and $P=\langle\Pi^2\rangle$. The erase control shares $X,P$ and has $C=0$. The matched fork vacuum has $X=1/(2I)$, $P=I/2$, $C=0$, where $I=a^3\sqrt{(q/a)^2+m^2}$ at the fork. The daughters start in their conformal vacuum with zero parent–daughter covariance. Every result uses all 24 native Gauss–Legendre nodes of the measure $q^2dq/(2\pi^2)$ on $0<q<8$.

For each $q$, take a comoving spectral cutoff $\Lambda$ and Gaussian quadrature on $q<\omega<\Lambda$. The massless Dirac vacuum spectrum is $K_q(\omega)=(\omega^2-q^2)/(8\pi^2)$. At order $y^2$, a collection of harmonic coordinates with $g_j^2=2\omega_j K_q(\omega_j)\Delta\omega_j$ realizes the required vacuum two-point function. For $A=a^{-3}$, $B=aq^2+a^3m^2$, the proper-time generator is

$$
H_\chi=\tfrac12(A\Pi^2+B\chi^2),\qquad
H_D=\frac{1}{2a}\sum_j(p_j^2+\omega_j^2x_j^2),\qquad
H_I=y f\chi\sum_jg_jx_j .
$$

The switch is $f(s)=\sin^4[\pi(s-s_l)/(s_r-s_l)]$ on its support and zero outside. The wide interval is $[1.35,2.05]$, the narrow one $[1.45,1.95]$, and the readout is at $s=2.10$. The two finite Hamiltonians are labeled by $k=0$ (bare) and $k=1$ (static). The latter adds

$$
H_{\rm ct}=\tfrac12 k y^2 a f^2 L_q \chi^2,\qquad L_q=\sum_jg_j^2/\omega_j^2 .
$$

This term compensates the static decrease of the parent stiffness obtained by completing the square in the bath coordinates at fixed $a,f$. A varying mass, switch and scale factor can still produce regulator dependence.

## Initial-value equations

Write $X=X_0+y^2X_2+\cdots$ and likewise for $C,P$. The free parent obeys $X_0'=2AC_0, C_0'=AP_0-BX_0, P_0'=-2BC_0$. For each bath mode, $Q_x,Q_p,R_x,R_p$ denote the leading symmetrized covariances of $(\chi x,\chi p,\Pi x,\Pi p)$ divided by $y g_j$. They obey

$$
\begin{aligned}
 Q_x'&=A R_x+Q_p/a, &
 Q_p'&=A R_p-\omega^2Q_x/a-fX_0,\\
 R_x'&=-B Q_x+R_p/a-f/(2\omega), &
 R_p'&=-B Q_p-\omega^2 R_x/a-fC_0 .
\end{aligned}
$$

For $S(z)=\sum_jg_j^2 z_j$, $V=k a f^2 L_q$ and bath conformal energy correction $D$, the causal response equations are

$$
\begin{aligned}
X_2'&=2AC_2,&
C_2'&=AP_2-BX_2-f S(Q_x)-V X_0,\\
P_2'&=-2BC_2-2f S(R_x)-2V C_0,&
D'&=-f S(Q_p).
\end{aligned}
$$

All mixed and second-order initial data vanish at the fork. The implementation evolves the response to $V$ additively, which lets both calibrations share the same free and first-order trajectories. This also explains why the two calibrations give identical daughter energy to this order.

## First law and internal transfer

At any time the four energy coefficients per $q$ are $E_\chi=(A P_2+B X_2)/2$, $E_D=D/a$, $E_I=f S(Q_x)$ and $E_{\rm ct}=V X_0/2$. Write $h=(\log a)'$. Positive work is supplied to these degrees of freedom. Integrating the following partial derivatives of the Hamiltonian gives the accumulated external work independently of the final energies.

$$
\begin{aligned}
\dot W_m &= a^3m m' X_2,\\
\dot W_{a,\chi}&=\tfrac12h[(a q^2+3a^3m^2)X_2-3A P_2],&
\dot W_{a,D}&=-hD/a,&
\dot W_{a,\rm ct}&=hV X_0/2,\\
\dot W_{f,I}&=f' S(Q_x),&
\dot W_{f,\rm ct}&=k a f f' L_q X_0 .
\end{aligned}
$$

The integrated internal currents have rates $J'_{I\to\chi}=-f A S(R_x)$, $J'_{I\to D}=-f S(Q_p)/a$ and $J'_{\rm ct\to\chi}=-V A C_0$. For example $E_\chi=W_m+W_{a,\chi}+J_{I\to\chi}+J_{\rm ct\to\chi}$, $E_D=W_{a,D}+J_{I\to D}$, and $E_I=W_{f,I}-J_{I\to\chi}-J_{I\to D}$. The counterterm obeys $E_{\rm ct}=W_{a,\rm ct}+W_{f,\rm ct}-J_{\rm ct\to\chi}$. Summing gives the total first law. During the switch these free-sector energies include reversible dressing. The cut follows switch-off and gives free outgoing daughter energy, redshifted to that cut.

## Numerical results

The following table uses the fine wall, wide switch and retained parent. Figures at cutoffs 120, 240 and 480 use the most refined settings. Both $W_f$ and $W_a$ include their explicit counterterm derivatives whenever it is present.

| Comoving cutoff | Daughter energy | Bare parent energy | Bare switch work | Static parent energy | Static switch work |
|---:|---:|---:|---:|---:|---:|
| 40 | 11.969184 | −25.821482 | −4.150218 | −25.220288 | −0.706465 |
| 60 | 21.577893 | −39.667496 | −6.651434 | −38.223323 | +1.508294 |
| 80 | 22.492854 | −41.698896 | −12.164219 | −39.064084 | +2.643269 |
| 120 | 22.510569 | −45.035161 | −30.325898 | −38.982968 | +3.542724 |
| 240 | 22.5105747285 | −63.383168 | −131.942526 | −38.821125 | +5.098852 |
| 480 | 22.5105747299 | −137.339069 | −542.943464 | −38.659863 | +7.126116 |

At cutoff 120 the bare mass and expansion works are $+6.955908$ and $+0.845399$, while static calibration gives $−17.299784$ and $−2.715338$. At switch peak $s=1.70$, the interaction energy is $−35.552103$, and the static counterterm contributes $+22.300116$. Both interaction terms vanish at switch-off. At cutoff 480 the static mass and expansion works are $−19.158191$ and $−4.117214$. Increasing the cutoff from 240 to 480 changes the final daughter energy by $1.43\times10^{-9}$ and the calibrated switch work by $+2.027265$. This scan demonstrates residual sensitivity to the chosen finite spectral model after static subtraction.

At cutoff 120 the matched fork-vacuum daughter energy is 23.100116 and the retained value is 22.510569. Subtracting the fork vacuum gives the following signed excess budget, normalized in the same units.

| Fine retained minus fork vacuum | Parent energy | Daughter energy | Switch work | Mass work | Expansion work |
|---|---:|---:|---:|---:|---:|
| Bare | +1.342396 | −0.589547 | +0.412449 | +0.345361 | −0.004961 |
| Static | +1.020595 | −0.589547 | −0.062831 | +0.437866 | +0.056013 |

## Verification and reproduction

A piecewise uniform proper-time mesh places nodes exactly at the fork, switch-on, peak, switch-off and cut. Fourth-order Runge–Kutta advances covariances, the bath energy and each work/current integral. We evaluate background values at every Runge–Kutta stage. The default frequency and time resolutions are 64 and 1,024. The identifiers containing `N` and `T` give explicit quadrature orders and time step counts. At cutoff 120, separate refinements to 96 modes and 2,048 steps distinguish integration errors. For the three states and both schemes the largest endpoint component change is $2.23\times10^{-5}$ on doubling the time steps and $4.39\times10^{-12}$ on increasing frequency order alone. At cutoffs 240 and 480 we increase both resolutions and obtain largest endpoint changes $5.27\times10^{-6}$ and $1.30\times10^{-6}$.

The total first-law residual is evaluated at 65 uniform step indices and at each of the five named checkpoints. Its maximum absolute value across all saved cases is $9.89\times10^{-8}$. Across all four component laws the maximum is $1.05\times10^{-7}$. These test energy accounting. The finer-step comparison separately quantifies solution accuracy. At cutoff 60 the earlier switched-emission calculation gave retained daughter energy 21.5778929119, compared with 21.5778928947 from the present hierarchy. The companion independent check calculates outgoing energy from free parent modes and frequency amplitudes.

From this directory run the following with Python 3.10+ and NumPy/SciPy. Replace the native directory argument with the location of the prior capsule's `FORK.json` and `WALLS.zip`.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B budget.py ../../HDBLAST_COVARIANCE_FEEDBACK_20261119/native --output BUDGETS.json
python -I -B check.py BUDGETS.json --output CHECKS.json
~~~

`BUDGETS.json` records the input SHA-256 digests, energies, external work, internal transfers, residuals and scale factor at every named checkpoint. The fine, wide, cutoff-60 run additionally saves intermediate histories for all three states and both calibrations. `CHECKS.json` reports componentwise integration differences. [Causal solver](budget.py) · [Budgets](BUDGETS.json) · [Identity and refinement checks](CHECKS.json)
