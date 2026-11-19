# A full-frequency check of the finite oscillation train

A four-cycle source produces 1.37500014 times the energy estimated from the flat long-duration width, already at leading order in the drive amplitude. The corresponding eight- and sixteen-cycle ratios are 1.09040179 and 1.02239118. The following time-domain calculation evaluates the continuum energy over all frequencies and supplies an independent reference for the discrete-sphere comparison. [Numbers](SPECTRAL.json) · [Replay](spectral.py)

## Normalization on the closed sphere

Use Einstein proper time and let $A=F^{-1/2}/F_m^{-1/2}$. The constant rescaling sets the daughter normalization outside the externally driven interval. On a static sphere of physical radius $R$, each real normalized harmonic has action $\frac12\int dt\,A^2(\dot\chi^2-\omega^2\chi^2)$. The scalar eigenfrequencies are $\omega_n=\sqrt{n(n+2)}/R$ with $(n+1)^2$ real oscillators at level $n$. The homogeneous degree of freedom has conserved $A^2\dot\chi_0$ and the Hamiltonian of a free particle in the stationary exteriors. It needs its own state. The outgoing excitation energy and number of individual particles in the nonzero modes are

$$
 \rho=\frac1{2\pi^2R^3}\sum_{n\ge1}(n+1)^2\omega_n|\beta_n|^2,\qquad
 n_X=\frac1{2\pi^2R^3}\sum_{n\ge1}(n+1)^2|\beta_n|^2 .
$$

The pair count is half the particle count. These quantities subtract the static vacuum contribution. The factor $(n+1)^2\,dn/(2\pi^2R^3)$ tends to $k^2\,dk/(2\pi^2)$ as $R$ grows at fixed physical momentum.

Set the parent mass $m=1$, use $t=m t_{\rm physical}$ and write $k=\omega/m$, $T=2\pi N$, and $A=\exp[\epsilon f(t)]$ with $f=\sin^4(\pi t/T)\cos t$ supported on $[0,T]$. Energy density and number density below are in units of $m^4$ and $m^3$. The tangent coupling is $\epsilon=\alpha\Phi/M_P$, with dimensional scalar amplitude $\Phi$ and $\alpha=-\tfrac12(d\ln F/d\Theta)_m$. The canonical oscillator $u=A\chi$ obeys $u''+[k^2-A''/A]u=0$. In terms of $u=\alpha_k e^{-ikt}+\beta_k e^{ikt}$ and $A\chi'=-ik(\alpha_k e^{-ikt}-\beta_k e^{ikt})$, its exact equations are

$$
 \alpha_k'=W e^{2ikt}\beta_k,\quad
 \beta_k'=W e^{-2ikt}\alpha_k,\quad W=(\log A)'.
$$

They preserve $|\alpha_k|^2-|\beta_k|^2=1$. Starting from the incident vacuum gives $\beta_{1,k}=2ik\epsilon\widetilde f(2k)$ with $\widetilde f(\nu)=\int f(t)e^{-i\nu t}\,dt$. The exact $\beta$ is odd in $\epsilon$ for this exponential drive. The time-ordered series bounds $|\beta-\beta_1|\le\sinh B-B$, $B=\int|W|\,dt$, at every $k$. An occupation computed at first order in the amplitude means $|\beta_1|^2$, of order $\epsilon^2$.

## Continuum energy from an autocorrelation

Set $g=f''$, extended by zero, and $C(h)=\int g(t)g(t+h)\,dt$. Parseval gives $n_{X,1}=\epsilon^2 C(0)/(16\pi)$. The energy is $\epsilon^2/(32\pi^2)$ times $\int_0^\infty \nu |\widetilde g(\nu)|^2\,d\nu$. Using $\int_{-\infty}^\infty [1-\cos(\nu h)]\,dh/h^2=\pi |\nu|$, the energy becomes

$$
 \rho_1=\frac{\epsilon^2}{16\pi^2}\int_0^\infty\frac{C(0)-C(h)}{h^2}\,dh
        =-\frac{\epsilon^2}{16\pi^2}\int_0^T\frac{C'(h)}h\,dh .
$$

For integer $N$, write $g(t)=\sum_a b_a\cos(at)$ on the pulse, where $a$ runs over the nonzero distinct frequencies $|1+j/N|$. The coefficients are $-a^2$ times the corresponding entries of $(1,-4,6,-4,1)/16$, merging repeated frequencies. They sum to zero. The diagonal autocorrelation integral of two unit cosines at $a$ is $(T-h)\cos(ah)/2-\sin(ah)/(2a)$. For two distinct frequencies it is $[b\sin(bh)-a\sin(ah)]/(a^2-b^2)$. Set $H_a=\operatorname{Ci}(aT)-\ln(aT)$ and $S_a=\operatorname{Si}(aT)$. Integrating the differentiated expressions gives the full energy as $\epsilon^2/(16\pi^2)$ times

$$
 \sum_a b_a^2\left(H_a+\frac{aT}{2}S_a\right)
 +2\sum_{a<b}b_a b_b\frac{a^2 H_a-b^2 H_b}{a^2-b^2}.
$$

The arbitrary common constant in $H_a$ cancels. The replay also evaluates the original time integrals directly for $N=4$. Both versions give $4.69970750649\times10^{-8}$ at $\epsilon=0.001$.

Once the continuum resonance is narrow, the parent rest-frame width $\Gamma=\alpha^2m^3/(32\pi M_P^2)$ gives $\rho_{\rm width}=\Gamma (m^2\Phi^2/2)\int w^2\,dt_{\rm physical}$, where $w=\sin^4(\pi t/T)$. In $m^4$ units this is $\epsilon^2(35T/128)/(64\pi)$. The same result follows by approximating $\nu^5$ by one on the positive carrier line, where $\widetilde f(\nu)\simeq\widetilde w(\nu-1)/2$, and integrating its square by Parseval. This limit requires both a narrow continuum line and sufficiently weak excitation as the train length grows. On a fixed sphere, a long train resolves individual levels. A fixed off-resonant mode has $\beta_1=O(\epsilon T^{-4})$, whereas an exactly resonant mode at $k=1/2$ has $|\beta_1|=3|\epsilon|T/16$ for the durations tested. The order in which the radius and the train duration increase therefore matters.

## Ultraviolet control

For $b=1/N$, the Fourier transform of $w$ is $(1-e^{-isT})3b^4/[2i\,s(s^2-b^2)(s^2-4b^2)]$, with removable poles. The two carrier shifts give $\widetilde f$. For $\nu\ge\nu_c=2k_c>1+2b$,

$$
 |\widetilde f(\nu)|\le B_c\nu^{-5},\qquad
 B_c=\frac{3b^4}{(1-1/\nu_c)^5[1-b^2/(\nu_c-1)^2][1-4b^2/(\nu_c-1)^2]}.
$$

The omitted continuum Born energy above $k_c$ is at most $\epsilon^2B_c^2/(128\pi^2\nu_c^4)$. At $k_c=4$, even the four-cycle bound is $1.03\times10^{-16}$ in $m^4$ units at the stated drive amplitude. The replay includes the analogous discrete bound obtained by bounding $\sum_{j\ge j_0}j^{-5}$ with $j_0^{-5}+j_0^{-4}/4$ for $j=n+1$. Endpoint differentiability matters. A compact $\sin^2$ envelope times this same cosine would instead give jumps in $f''$, a $k^{-2}$ amplitude envelope and logarithmically divergent continuum Born energy.

For the exact oscillator, set $V=A''/A$ and $\Omega=\sqrt{k^2-V}$ at sufficiently large $k$. The instantaneous-frequency Bogoliubov equations couple the two components with $\Omega'/2\Omega=O(k^{-2})$. Both $V$ and $V'$ vanish at the exterior endpoints. Integrating the outgoing component once by parts against $e^{-2i\int\Omega}$ gives $|\beta|=O(k^{-3})$, because $V''$ is piecewise integrable and the Bogoliubov components are bounded by $\exp\int|\Omega'/2\Omega|$. This suffices for absolute convergence of the exact energy in flat and closed space. The preceding explicit tail constants apply at Born order.

The prescribed source supplies energy for the whole train. A first passage in the fitted inhomogeneous geometry supplies only the local field history. Computing daughter production there requires a specified quantum state, propagation through the changing geometry and a particle or stress-energy measurement convention. A later radiation fraction would also require a backreaction calculation. Collisionless production alone supplies no thermalization time or temperature. [Previously calibrated flat width](../modulus/widths/WIDTHS.json) · [Finite round-sphere computation](../spectrum/SPECTRUM.json)
