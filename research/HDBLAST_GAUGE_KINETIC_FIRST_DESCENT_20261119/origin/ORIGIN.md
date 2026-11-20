# Does the five-dimensional wall provide a photon coupling?

The archived wall action fixes its four-dimensional gravitational coupling and the field that moves the wall. It specifies no bulk Maxwell action or scalar-dependent gauge kinetic function. A constant-coupling Maxwell field confined to the wall has zero direct response to its changing Weyl factor. Adding a bulk gauge field introduces a different question. Does the putative photon have a finite norm in the extra dimension, and, if so, does its normalization change as the wall moves? The answers depend on the geometry and on the added gauge action. [Archived action and field tables](../prior/prior/gravity/previous/prior/archive/blast_common.py) · [Original dynamics report](../prior/origin/BLAST_DYNAMICS_REPORT.md)

## Reduction and normalization

Write the added action as $-\int d^5X\sqrt{-g_5}\,G(\phi)F_{MN}F^{MN}/(4g_5^2)$, where $g_5$ is a constant gauge coupling. Take the metric to be $e^{2A}g_{\mu\nu}dx^\mu dx^\nu+e^{2B}dy^2$. For a static leaf assume $A_y=0$, a real, even gauge profile $A_\mu(x,y)=a_\mu(x)\psi(y)$ and a source-free Neumann condition at the wall. The radial equation is

$$-\partial_y\big(G e^{2A-B}\partial_y\psi\big)=m^2G e^B\psi.$$

A constant profile solves this equation at zero mass. Normalize the gauge coordinate once so that this profile equals one at the wall. If it has finite norm and the zero mode may be retained, the reduced action is $-\int\sqrt{-g_4}\,I^2 F_{\mu\nu}F^{\mu\nu}/4$, where

$$ I^2=\frac{1}{g_5^2}\sum_{\rm copies}\int G(\phi)e^B\,dy .$$

The factors $e^{4A}$ from the volume and $e^{-4A}$ from the two inverse four-metrics cancel. For $G=1$ on a finite interval, the coefficient measures proper length. A finite brane-localized kinetic term adds its own coefficient. A fixed proper length and fixed microscopic coupling give a constant coefficient. The Einstein–Hilbert coefficient instead integrates $e^{2A+B}$. In four dimensions the same cancellation removes a universal Weyl dependence from a minimally coupled Maxwell field living directly on the brane. These results require the specified field action and boundary conditions. A coupling to charged species, an anomaly, or a varying brane kinetic term would require additional physical input. [Compact warped example, Eqs. 4.6–4.9](https://arxiv.org/abs/0705.3844) · [Scalar-weighted thick-brane example](https://arxiv.org/abs/1402.6480)

## The two archived geometries

**Static detuned shell.** The registered chart is $ds^2=dv^2+\rho(v)^2ds_{\rm dS_4,unit}^2$ between a horizon at $v=0$ and a fixed shell at $v_b=8.22820155 L_0$. The even constant vector has finite two-copy radial norm $g_5^2 I^2=2v_b=16.45640310 L_0$ on this retained patch. Gravity instead gives $M_4^2/M_5^3=2\int_0^{v_b}(\rho/\rho_b)^2dv=2.06773278 L_0$, which we recomputed from the regular cone initialization and Einstein–scalar equations. A single static length supplies no derivative with respect to a rolling modulus. Normalizability on the chosen patch does not by itself specify continuation of all field configurations across the Lorentzian horizon. [Registered solution](../prior/origin/REGISTERED_SHELL_FLOAT_SOLUTION.json) · [Five-dimensional equations](../prior/origin/hdblast_background.py)

The static minimal Maxwell radial equation reads $-(\rho^2\psi')'=\lambda\psi$. Set $dz=dv/\rho$ and $J=\sqrt\rho$. The canonical function $J\psi$ has potential $J_{zz}/J=\rho_v^2/4+\rho\rho_{vv}/2$. Near the regular cone $\rho=v+O(v^3)$ and $z\to-\infty$, so the potential tends to $1/4$ exponentially. The continuum starts at physical vector mass $H_{\rm shell}/2$, with $H_{\rm shell}=1/\rho_b$. The tensor continuum of this patch starts at $3H_{\rm shell}/2$. This calculation locates the vector continuum threshold, leaving the possible presence of additional discrete vector modes below it open. [Archived tensor calculation](../prior/origin/TENSOR_SECTOR_AND_CALIBRATION.md)

**Tuned leaves underlying the rolling effective theory.** On the regular branch $\phi(y)=-\tanh y$, the retained proper coordinate runs from $y_b$ to $+\infty$. Its warp exponent obeys $A_y=-W/3$, with $W=1-\phi+\phi^3/3$. Gravity localizes because $2\int_{y_b}^\infty e^{2(A-A_b)}dy$ converges. The program checks that this integral agrees with the archived $f(\Theta)$ within $10^{-10}$ at launch and at both first crossings. The minimal bulk photon instead has norm $2\int_{y_b}^\infty dy/g_5^2=\infty$. Thus the massless constant profile cannot supply a normalized four-dimensional photon on these leaves. In particular one cannot substitute either $f(\Theta)$ or the proper length of the separately registered detuned shell for the gauge norm on the rolling leaves. [Leaf geometry](../prior/prior/gravity/previous/prior/archive/blast_common.py)

## A photon confined to the wall

An independent gauge field living on the shell can have a local action
$S_{\rm wall}=-\tfrac14\int d^4x\sqrt{-h}\,J(\phi_b)^2 F_{\mu\nu}F^{\mu\nu}_{(h)}$, where $h$ is the induced metric and $J^2>0$. It needs no bulk normalization integral. A Weyl transformation from $h$ to the Einstein-frame metric cancels between the measure and the raised indices. The remaining coefficient is $I(\Theta)=J[\phi_b(\Theta)]$, where the positive square root fixes the sign convention. One can therefore specify a positive brane coefficient with constant values before and after a smooth transition in field space, as in the switched test. The wall action leaves its contrast, sign and shape as inputs. Varying it adds $J^2(F_{\mu\alpha}F_\nu{}^\alpha-h_{\mu\nu}F^2/4)$ to the surface stress and $-(\partial_{\phi_b}J^2)F^2/4$ to the surface action derivative per proper four-volume. A coupled solution must include both contributions in the metric and scalar matching conditions. [Switched test protocol](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7025474763)

The pure bulk constant-profile projection on the two tuned leaves is more restrictive. Since $d\phi/dy=-(1-\phi^2)$ on the tanh branch, a scalar weight independent of shell location gives
$$ I^2(\phi_b)=\frac{2}{g_5^2}\int_{-1}^{\phi_b}\frac{G(u)}{1-u^2}\,du,\qquad
\frac{dI^2}{d\phi_b}=\frac{2G(\phi_b)}{g_5^2(1-\phi_b^2)}\geq0 .$$
This assumes convergence at the throat and no additional brane kinetic term. A strictly positive bulk weight produces an increasing coefficient on the open branch. Enforcing an exactly flat interval of this bulk coefficient would require the local Maxwell kinetic weight to vanish there. Positive contrast on both oppositely directed rolls also cannot come from one such monotone bulk coefficient. A brane function can have constant plateaus at positive kinetic normalization and can vary in either direction. The smooth bulk example below and the switched brane test thus describe different interactions.

## A specified bulk extension

As an explicit alternative assume $G(\phi)=(1-\phi^2)^2=W_{,\phi}^2$. The archived wall field supplies the function's argument. The choice to multiply $F^2$ by this nonnegative function is a new microscopic assumption. On the open tanh branch it is positive, decays into the throat, and gives a normalizable constant vector. Writing $\phi_b=\phi(y_b)$ produces the following analytic coefficient per copy,

$$ P(\phi_b)=\int_{y_b}^{\infty}\operatorname{sech}^4y\,dy
=\frac23+\phi_b-\frac{\phi_b^3}{3},\qquad
\frac{I(\Theta)}{I(\Theta_i)}=\sqrt{\frac{P[\phi_b(\Theta)]}{P[\phi_b(\Theta_i)]}} .$$

Direct radial quadrature checks this expression. The canonical-field slope is $(\log I)_{,\Theta}=(1-\phi_b^2)\phi_{b,\Theta}/(2P)$. We evaluated it on the unseeded homogeneous trajectories of the two modified-potential branches, using their recorded displaced launches. At either bounce $d\log I/dx=0$ because the scalar starts at rest. Both endpoints below are *moving* first crossings. The variable $x=H_0t_E$ is the archived Einstein proper time. The normalization in both branches is $H_0/M_{\rm Pl}=0.00881344$. [Canonical field calibration](../prior/prior/gravity/previous/prior/background/background.py) · [Rolling trajectory](../prior/homogeneous/results.json)

| First minimum $\phi_b$ | $\Theta_*$ | $\log(I_*/I_i)$ | $(\log I)_{,\Theta}$ at crossing | $d\log I/dx$ at crossing |
|---|---:|---:|---:|---:|
| $+0.5$ | −0.385805 | +0.260028 | −0.354042 | +0.794669 |
| $-0.5$ | +0.296939 | −0.579973 | −3.325056 | −7.535520 |

The localized projection has a limitation even for a free bulk field. With $dz=e^{-A}dy$, its canonical radial function is $J\psi$, $J=(G e^A)^{1/2}$. At the throat $e^A\sim e^{-k y}$ with $k=5/9$, and $G\sim 16e^{-4y}$. Thus $J\sim z^{-s}$ with $s=(1+4/k)/2=4.1$, and $J_{zz}/J\sim 20.91/z^2$. The bulk continuum starts at zero mass even though the constant mode has finite norm. There is no finite mass gap that automatically justifies neglecting continuum excitation by a moving wall. The moving normal boundary condition and overlaps with bulk modes require a separate five-dimensional evolution. Also, $G\to0$ in the deep throat would require interaction control in any completion containing charged matter. Changing $G$ or adding an independent brane kinetic term changes the response.

## Energy of one prepared wave at a moving crossing

We tested what this particular zero-mode projection would predict if it were a valid four-dimensional action. In a round universe the Hamiltonian of one real transverse harmonic with $L\ge2$ is $[p^2/I^2+L^2 I^2q^2]/(2a)$. Prepare the positive-frequency ground state of this original Hamiltonian at the time-symmetric bounce and compare with the state obtained by adding one excitation in this same harmonic. The energy difference coincides with that of the normalized complex classical solution. It is positive and finite, with no vacuum subtraction ambiguity or sum over unprepared modes. The phase-averaged classical excitation obeys the same ratio. In the adiabatic reference the excitation retains its action and has energy $E_{\rm ad}=E_i a_i/a$.

For $d\eta=dx/a$ and $s_x=\partial_x\ln I$, write the instantaneous ladder variable as $b=e^{-iL\eta}(\alpha b_i+\beta b_i^\dagger)$. The evolution equations are $\dot\alpha=s_x e^{2iL\eta}\beta^*$ and $\dot\beta=s_x e^{2iL\eta}\alpha^*$ in $x$. They imply $|\alpha|^2-|\beta|^2=1$, and the measured ratio is $E_{\mathrm{add}}/E_{\rm ad}=1+2|\beta|^2$. A second solver integrates the original canonical variables and evaluates their physical stress energy directly. This definition permits an energy measurement while the coupling is changing. It specifies no asymptotic photon count.

| First minimum | $E_{\rm ad}/E_i$ | $E_{\rm add}/E_{\rm ad}$ at $L=2$ | $E_{\rm add}/E_i$ at $L=2$ | $E_{\rm add}/E_i$ at $L=16$ |
|---|---:|---:|---:|---:|
| $+0.5$ | 0.900254 | 1.122932 | 1.010924 | 0.902259 |
| $-0.5$ | 0.880029 | 1.695341 | 1.491950 | 0.998147 |

The lowest allowed harmonic gains about 1.09% or 49.19% of its initial energy on the two prescribed rolling backgrounds. A sufficiently high-frequency excitation instead approximately follows ordinary redshift. For any such number-state addition $\operatorname{arsinh}|\beta|\le\int |d\log I|$. On these monotone legs this makes the all-frequency upper bounds on $E_{\rm add}/E_{\rm ad}$ equal to $\cosh(2|\Delta\log I|)$, or 1.138305 and 1.751631. This bound is for the specified unsqueezed incident covariance. Coherent waves of specified phase can show an interference term.

Since $d\log I/d\eta$ at the endpoint is nonzero and vanishes at launch, integration by parts gives $L^2(E_{\rm add}/E_{\rm ad}-1)\to (\partial_\eta\log I)_*^2/2$. The high-frequency endpoint term depends on the measurement prescription. Summing instantaneous zero-point occupations over all frequencies here would produce a quadratic ultraviolet sensitivity. The measured single-excitation difference needs no such sum. These finite-time projected results do not determine the absolute quantum stress, real photon abundance, thermalization, or the effect of leakage into the five-dimensional continuum. [Numerical results and controls](ORIGIN.json) · [Replay](origin.py)
