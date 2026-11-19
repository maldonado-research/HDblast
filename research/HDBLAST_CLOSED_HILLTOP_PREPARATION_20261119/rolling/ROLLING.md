# Scalar and metric perturbations on a rolling closed bounce

For the two previously tuned bounce trajectories at each of three detunings, the stationary closed-de Sitter calculation accurately predicts the early field response of every scalar harmonic tested, $n=2,\ldots,16$. At a homogeneous displacement of $10^{-5}M_{\rm Pl}$ from the hilltop, the rolling and stationary field-transfer rows differ by at most 0.00399% for the positive minimum and 0.000530% for the negative minimum. At $10^{-2}M_{\rm Pl}$ the differences reach 3.93% and 0.554%. These comparisons use equal elapsed proper time, starting at the bounce. They apply to the response to infinitesimal initial data. [Machine-readable rows](ROLLING.json)

## Field and metric

Use $x=H_{\rm ref}t_E$, $a=H_{\rm ref}a_{\rm physical}$, $u=V/(3M_{\rm Pl}^2 H_{\rm ref}^2)$, $h=(\log a)'$, and $v=\Theta'$, with the canonical field expressed in reduced Planck units. Here $H_{\rm ref}$ is the archive hilltop scale for the indicated detuning. The exact reconstructed stationary scale obeys $h_s^2=u(\Theta_s)\simeq1$. The signed squared mass in these units is $3u_{,\Theta\Theta}$, equal to $-207.2548854$ or $-140.9655115$ at the two stationary points. For the rolling trajectory the background obeys $v'=-3hv-3u_{,\Theta}$, $h'=a^{-2}-v^2/2$ and $h^2=u+v^2/6-a^{-2}$. Radiation and the oscillator decay width are inactive throughout this early interval.

Choose the Newtonian slicing with metric $H_{\rm ref}^2ds^2=-(1+2\Phi Y)dx^2+a^2(1-2\Phi Y)d\Omega_3^2$. The canonical scalar has zero scalar anisotropic stress, and $n\ge2$ admits equal Newtonian potentials. The harmonic obeys $\Delta Y=-k_n^2Y$, $k_n^2=n(n+2)$. Use any real unit-RMS member of its multiplet, $\int_{S^3}Y^2d\Omega_3/(2\pi^2)=1$. Write the field perturbation as $q(x)Y$. In this slicing $q$ and $\Phi$ coincide with the shear-free gauge-invariant field perturbation and Bardeen potential. The Einstein momentum equation and perturbed Klein–Gordon equation give

$$
\Phi'=-h\Phi+\tfrac12vq,\qquad
q''+3hq'+(k_n^2/a^2+3u_{,\Theta\Theta})q
=4v\Phi'-6u_{,\Theta}\Phi .
$$

The time-time Einstein equation and the momentum constraint imply

$$
C=\left[2(k_n^2-3)/a^2-v^2\right]\Phi+v q'-v'q=0 .
$$

Differentiating gives $C'=-3hC$, so $a^3 C$ is conserved in exact arithmetic. At the time-symmetric bounce, $v=h=0$, $a_0=u(\Theta_0)^{-1/2}$, $\Phi'_0=0$ and $\Phi_0=v'_0 a_0^2q_0/[2(k_n^2-3)]$. These relations permit an arbitrary infinitesimal initial pair $(q_0,q'_0)$. The evolution stays regular through the bounce and through a possible zero of $2(k_n^2-3)/a^2-v^2$. The same calculation on the stationary background has $v=u_{,\Theta}=0$. The constraint then fixes $\Phi=0$ and the field evolves autonomously with $a_s=\cosh(h_s x)/h_s$. The exceptional homogeneous and $n=1$ sectors require their own constraint analysis.

## Numerical comparison

Each recorded row $T_q=(q^{(1)},q^{(2)})$ maps the initial pair $(q_0,q'_0)$ to the mode at a saved endpoint, and analogous rows give $q'$ and $\Phi$. The reported fractional change is $\|T_{q,\rm roll}-T_{q,s}\|_2/\|T_{q,s}\|_2$ in the specified classical coordinates. A specially correlated initial pair can cancel a row, so the row comparison is an operator-scale difference. The figures below are maxima across $n=2,\ldots,16$ and the three detunings. Each line uses the event time on its own rolling trajectory.

| Homogeneous displacement in $M_{\rm Pl}$ | Field change, positive minimum | Field change, negative minimum | Velocity change, positive minimum | Velocity change, negative minimum |
|---:|---:|---:|---:|---:|
| $10^{-5}$ | 0.003984% | 0.000530% | 0.007982% | 0.001064% |
| $10^{-4}$ | 0.03981% | 0.005277% | 0.07970% | 0.01060% |
| $10^{-3}$ | 0.3974% | 0.05292% | 0.7941% | 0.1064% |
| $10^{-2}$ | 3.925% | 0.5538% | 7.724% | 1.136% |

A control evolves a test field on the same rolling background with the metric source in its field equation set to zero. Restoring the metric source changes the field transfer by at most $3.0\times10^{-11}$ at $10^{-5}M_{\rm Pl}$ and $3.0\times10^{-5}$ at $10^{-2}M_{\rm Pl}$, across both minima. Most of the finite change in the table already appears in the field equation with rolling scalar mass and background geometry.

At detuning $10^{-3}$ and displacement $10^{-5}M_{\rm Pl}$, the $n=2$ field rows are $(1.17021\times10^{19},8.25864\times10^{17})$ and $(5.61461\times10^{13},4.84131\times10^{12})$, for positive and negative minima respectively. A quantum prediction for a rolling background also requires a state on the scalar-gravity reduced phase space. The transfer rows specify the dependence on initial data. The canonical area below gives a state-independent two-endpoint bound. Neither a classical row norm nor a growing linear root by itself determines nonlinear evolution.

## What the canonical bounce implies for two-time quantum fluctuations

For a unit-RMS real harmonic define $F=2\pi^2(M_{\rm Pl}/H_{\rm ref})^2$. Set $W(A,B)=A_1B_2-A_2B_1$ for two classical perturbations, and take the symplectic orientation $dq\wedge dp$. In Newtonian gauge the metric shift vanishes and the lapse perturbation is $+\Phi Y$. The scalar's canonical momentum contains $\sqrt{\gamma}/N$. Its linear variation for this mode is $Fa^3(q'-4v\Phi)$ after harmonic integration. For the spatial conformal variable $b=\ln a-\Phi Y$ the gravitational momentum density is $\pi_b=2\gamma_{ij}\pi^{ij}=-2M_{\rm Pl}^2\sqrt{\gamma}\,K$, with extrinsic curvature $K_{ij}=(2N)^{-1}\partial_t\gamma_{ij}$. Its linear variation after integration is $6Fa^3(\Phi'+4h\Phi)$. Pairing this with $\delta b=-\Phi$ and adding the matter variation gives

$$
\frac{\Omega}{F}=a^3\left[W(q,q')-4v W(q,\Phi)-6W(\Phi,\Phi')\right]
=a^3\left[W(q,q')-vW(q,\Phi)\right].
$$

The last equality uses the momentum constraint. The perturbation equations give $(\partial_x+3h)W(q,q')=(-4hv-6u_{,\Theta})W(q,\Phi)$. Wedging the time-time constraint with $\Phi$ gives $vW(q',\Phi)=v'W(q,\Phi)$. Consequently $(\partial_x+3h)[vW(q,\Phi)]=(2v'+2hv)W(q,\Phi)$, which equals the same right-hand side by the background field equation. This proves conservation, including at $v=0$ and when $2(k_n^2-3)/a^2-v^2=0$.

A direct two-variable reduction supplies another check. Write $P=q'-v\Phi$, $\alpha=2(k_n^2-3)/a^2$ and $\omega=k_n^2/a^2+3u_{,\Theta\Theta}$. For $n\ge2$ and finite $a$, $\alpha>0$ and the constraint fixes $\Phi=(v'q-vP)/\alpha$. The remaining equations become $q'=P+v\Phi$ and $P'=-3hP-(\omega-3v^2/2)q+v'\Phi$. Substitution of the constraint gives a matrix of trace $-3h$, again conserving $a^3 W(q,P)$.

At the time-symmetric bounce the conserved area becomes $\Omega=Fa_0^3W(q,q')$. With $\hbar=1$ the linearized canonical field thus has $[q_0,q'_0]=i\kappa$, $\kappa=(H_{\rm ref}/M_{\rm Pl})^2/(2\pi^2a_0^3)$. For $q(X)=A(X)q_0+B(X)q'_0$ it follows that $[q_0,q(X)]=i\kappa B(X)$. Robertson's inequality bounds the product of the two endpoint standard deviations from below by $\kappa|B|/2$ for each real harmonic. For the full $n=2$ multiplet let $S(x)=\sum_{i=1}^9{\rm Var}\,q_i(x)$. This equals the contribution to the spatially averaged connected field variance under the unit-RMS convention. Cauchy–Schwarz on the nine pairs gives

$$
\max\{\sqrt{S(0)},\sqrt{S(X)}\}\ \ge\ \sqrt{9\kappa |B_2(X)|/2}.
$$

The corresponding one-real-harmonic floor is one third of this multiplet floor. On the unconstrained collection of linear oscillators, a product Gaussian saturates the two-endpoint inequality by equalizing endpoint variances and setting each symmetrized endpoint covariance to zero. This statement needs no prescription for a cosmological vacuum. Extension to a physical state of full compact gravity additionally requires the gravitational constraints at higher order, including any integrated constraints associated with background symmetries. The rolling calculation establishes a necessary linear two-endpoint bound. It leaves the intervening-time optimum, nonlinear evolution, and absolute preparation probabilities open.

| Minimum | Detuning | $x$ when the two-endpoint $n=2$ RMS floor reaches $0.01M_{\rm Pl}$ | Homogeneous displacement there in $M_{\rm Pl}$ | $n=2$ endpoint floor at prior match in $M_{\rm Pl}$ |
|---|---:|---:|---:|---:|
| $+0.5$ | $10^{-3}$ | 0.364848 | $5.66\times10^{-23}$ | $3.82417\times10^6$ |
| $+0.5$ | $10^{-10}$ | 1.580413 | $3.54\times10^{-46}$ | $1.66076\times10^{18}$ |
| $+0.5$ | $10^{-20}$ | 3.352454 | $2.05\times10^{-79}$ | $6.98413\times10^{34}$ |
| $-0.5$ | $10^{-3}$ | 0.431479 | $9.45\times10^{-18}$ | $9.25908\times10^3$ |
| $-0.5$ | $10^{-10}$ | 1.944437 | $4.32\times10^{-35}$ | $4.77883\times10^{12}$ |
| $-0.5$ | $10^{-20}$ | 4.142873 | $5.63\times10^{-60}$ | $1.33281\times10^{25}$ |

Each value uses the constrained rolling $B_2$. At these early crossings the rolling and independently evaluated stationary crossing times agree within $1.4\times10^{-14}$ in $x$, inside the integration accuracy. Samples 1% later yield floors of $0.01024$–$0.01243M_{\rm Pl}$. Large formal late-time floors describe the continued linear equations and do not give a nonlinear fluctuation amplitude. The initial radius $a_0=u(\Theta_0)^{-1/2}$ is unity at displayed precision for all six, and $\kappa\simeq3.93515\times10^{-6}(t/10^{-3})$ for detuning $t$.

A local numerical test initializes two independent balanced perturbations spanning the time-time constraint plane at nonzero background velocity. It transports them over a short rolling interval for $n=2,8$ on every background. The uncoupled field expression $a^3W(q,q')$ changes by up to 2.10 times the initial full area. The full scalar–metric area stays constant to $1.08\times10^{-10}$ across all twelve checks. Nine of these tests cross $2(k_n^2-3)/a^2-v^2=0$. Starting the checks locally avoids catastrophic cancellation between two columns that grew almost parallel over a long unstable interval.

## Compensating initial velocity

Four further integrations test the two endpoint fluctuation bound when the background starts on the compensated stable trajectory with $v_i=\pm\sqrt{0.06}$ and $t=10^{-3}$. The input is the signed recentered bounce field supplied by the homogeneous calculation. This trajectory approaches the hilltop during the test. Its background displacement at the fluctuation crossing has magnitude $7.86\times10^{-5}$–$1.05\times10^{-4}M_{\rm Pl}$. [Homogeneous centers](../manifold/MANIFOLD.json) · [Numerical data](VELOCITY.json)

Now $a_i=(u_i+v_i^2/6)^{-1/2}$, $D_i=\alpha_i-v_i^2$ and $\Phi_i=(v'_iq_i-v_iq'_i)/D_i$, with $\alpha_i=2(k_n^2-3)/a_i^2$. For these four regular initial charts, the canonical area at launch is $Fa_i^3(\alpha_i/D_i)W(q,q')$. The commutator for $q_i,q'_i$ thus has coefficient $\kappa_{\rm eff}=(1-v_i^2/\alpha_i)/(Fa_i^3)$. Its extra factor is 0.9940004–0.9940007 for $n=2$. A second propagation in $q,P=q'-v\Phi$ uses the canonical coefficient $1/(Fa_i^3)$ directly. The transfer rows computed in the two coordinate systems agree after the initial-coordinate transformation to a relative difference below $1.70\times10^{-13}$. The full symplectic area drifts by less than $9.71\times10^{-12}$.

| Minimum | $v_i$ | $x$ when the two-endpoint multiplet RMS floor reaches $0.01M_{\rm Pl}$ | First minimum crossing, bracket pair | Two-cycle gate, bracket pair |
|---|---:|---:|---:|---:|
| $+0.5$ | $-\sqrt{0.06}$ | 0.360735 | 4.28020–4.28946 | 14.25025–14.25951 |
| $+0.5$ | $+\sqrt{0.06}$ | 0.369206 | 4.28028–4.28954 | 14.25034–14.25959 |
| $-0.5$ | $-\sqrt{0.06}$ | 0.432323 | 4.03034–4.04182 | 8.27904–8.29052 |
| $-0.5$ | $+\sqrt{0.06}$ | 0.430692 | 4.03033–4.04181 | 8.27903–8.29050 |

The threshold distance from the specified center is below a double-precision increment at that field value. For this early calculation the center itself is integrated directly. Repeating with an initial offset of $+10^{-9}$ and $-10^{-9}M_{\rm Pl}$ changes the crossing times by at most $2.43\times10^{-8}$ in $x$. The calibrated threshold distances are $5.65\times10^{-25}$–$1.27\times10^{-19}M_{\rm Pl}$. The later event clocks in the table come independently from the two relative launches at 0.94 and 1.06 times the threshold distance. Only these four early compensated examples were tested. [Paired homogeneous events](../manifold/EVENTS.json)

## Checks and reproduction

The second evolution uses the Einstein spatial trace

$$
\Phi''+4h\Phi'+(2h'+3h^2-a^{-2})\Phi
=\tfrac12(vq'-v^2\Phi-3u_{,\Theta}q)
$$

together with the field equation and the same constrained bounce data. These two evolutions share the homogeneous solution. Each run tracks the normalized time-time constraint in the momentum evolution and the normalized momentum constraint in the trace evolution. For a residual row these diagnostics divide its Euclidean norm by the norm of the sum of the absolute componentwise terms. The homogeneous Friedmann residual uses the original dimensionless energy normalization.

For stability the code evolves the logarithm of the homogeneous offset magnitude and its log derivative. Within $2\times10^{-5}M_{\rm Pl}$ it uses the prior root-centered polynomial to evaluate $u_{,\Theta}/(\Theta-\Theta_s)$. It divides all transfer rows by the same analytic factor $\exp(\lambda_+x)$ and divides the metric rows by the signed background displacement as well. Here $\lambda_+=(\sqrt{9h_s^2-12u_{,\Theta\Theta}(\Theta_s)}-3h_s)/2$. This scaling resolves the smallest starting metric response. Every rolling seed has nonzero homogeneous displacement. The stationary comparison evolves the decoupled field with zero metric perturbation. All persisted transfer rows restore their original units.

Across all six runs, the maximum scaled time-time residual was $8.42\times10^{-10}$, the maximum scaled momentum residual of the spatial-trace evolution was $7.75\times10^{-12}$ and the maximum Friedmann residual was $1.91\times10^{-13}$. The two metric systems agreed in every endpoint field, velocity and metric transfer row within $4.86\times10^{-13}$ in fractional row norm. Repeating with `--rtol 3e-12 --max-step .017` changed any endpoint time by at most $3.68\times10^{-12}$ in $x$. The rolling $0.01M_{\rm Pl}$ crossing time changed by at most $5.78\times10^{-15}$. Rolling field, velocity and metric transfer rows changed by at most $1.93\times10^{-10}$, $1.72\times10^{-9}$ and $3.03\times10^{-10}$ respectively.

Run from the study root with Python, NumPy and SciPy. The prior capsule supplies hilltop/hilltop.py, hilltop/check.py, hilltop/HILLTOP.json, background/background.py, curvature/curvature.py and the two pinned archive files. DOP853 uses relative tolerance $10^{-11}$, absolute tolerance $10^{-15}$ for rescaled variables, and maximum step 0.035 in $x$.

```sh
OPENBLAS_NUM_THREADS=1 python -IB rolling/rolling.py --prior prior --output /tmp/ROLLING.json
OPENBLAS_NUM_THREADS=1 python -IB rolling/velocity.py --prior prior --manifold manifold/MANIFOLD.json --output /tmp/VELOCITY.json
```
