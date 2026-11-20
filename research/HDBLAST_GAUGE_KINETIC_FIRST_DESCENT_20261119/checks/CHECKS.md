# How much energy can a smooth change in the photon coupling produce?

A smooth change in the electromagnetic kinetic coefficient can populate transverse waves even on a round closed spacetime. An exactly solvable example with $\Delta\ln I=0.62832$ gives $2.94896$ outgoing photons in expectation on the entire three-sphere. Its physical energy density above the outgoing vacuum is $0.50803213/a^4$ in units $\hbar=c=1$. A second example ends after a finite conformal time and admits a rigorous upper bound on the energy in all omitted harmonics. Both histories are specified test functions. They provide independent checks on the propagation and normalization used for the rolling backgrounds.

## Source, spectrum and interchange of electric and magnetic fields

In conformal time $\eta$ the scale factor cancels from the four-dimensional action $-\frac14\int\sqrt{-g}\,I^2 F_{\mu\nu}F^{\mu\nu}$. Each real transverse harmonic on the unit three-sphere therefore has action $I^2(A'^2-k^2A^2)/2$, canonical momentum $P=I^2 A'$ and Hamiltonian
$\mathcal H_k=E_k+B_k$ with $E_k=|P/I|^2/2$ and $B_k=k^2|IA|^2/2$. Here $k=2,3,\ldots$ and each shell has $2(k^2-1)$ real oscillators including both polarizations. No transverse zero mode occurs. The energy per oscillator divided by $a$ is its proper energy. From Hamilton's equations, with $u=(\ln I)'$,

$$
\mathcal H'_k=\left.\partial_\eta\mathcal H_k\right|_{A,P}=2u(B_k-E_k).
$$

Thus $\dot\rho_\gamma+4H\rho_\gamma=2\dot{\ln I}(\rho_B-\rho_E)$, with the dot in proper time. It holds directly for classical waves, for each finite oscillator sum, and for differences of sufficiently regular states in the same theory. A dynamical scalar driving $I$ receives the opposite exchange. Expansion accounts separately for the $4H\rho_\gamma$ dilution. The work can change sign during a pulse even when the final excess from the in-vacuum is positive.

Define $X=IA$, $Y=P/I$ and $\alpha e^{-ik\eta}=(\sqrt{k}X+iY/\sqrt{k})/\sqrt{2}$, $\beta e^{ik\eta}=(\sqrt{k}X-iY/\sqrt{k})/\sqrt{2}$. The positive incoming vacuum mode starts with $\alpha=1,\ \beta=0$ in this convention. The evolution is

$$
\alpha'=u e^{2ik\eta}\beta,\qquad
\beta'=u e^{-2ik\eta}\alpha,\qquad |\alpha|^2-|\beta|^2=1 .
$$

Replacing $I$ by a constant times $I^{-1}$ sends $u\to-u$, leaves $\alpha$ fixed and reverses $\beta$. The electric and magnetic contributions interchange and the final occupation $N_k=|\beta_k|^2$ agrees for the two profiles. Also $d|\beta|/d\eta\le |u|\sqrt{1+|\beta|^2}$ wherever differentiable. Integration gives the nonperturbative bound $N_k\le\sinh^2(\int |d\ln I|)$ on each mode, including paths that reverse direction. A monotone switch uses simply $|\Delta\ln I|$. This uniform bound alone cannot control the full high-frequency sum.

With constant initial and final coefficients the in and out positive frequencies are fixed on the conformal cylinder. The late physical energy above the out vacuum satisfies

$$
a^4\rho_{\rm excess}=\frac{1}{\pi^2}\sum_{k=2}^\infty k(k^2-1)N_k .
$$

The vacuum $k/2$ per oscillator cancels in this comparison. An initial state diagonal in the incoming number basis with mean occupation $m_k$ has outgoing occupation $m_k+(2m_k+1)N_k$. Neither the vacuum nor a squeezed out state is a thermal bath.

## A soluble smooth switch

Consider a conformal cylinder covering the whole real $\eta$ line, with $I(\eta)=\exp[2g\arctan e^{\eta/\tau}]$, $\tau>0$. It tends to $1$ and $e^{\pi g}$ at past and future infinity. In this benchmark $u=g\,\operatorname{sech}(\eta/\tau)/\tau$. It gives an ideal smooth scattering problem with asymptotically constant couplings. A pulse with exact finite plateaus is treated below.

Eliminate $\beta$ and use $y=[1+\tanh(\eta/\tau)]/2$ and $\kappa=k\tau$. The regular incoming solution satisfies

$$
y(1-y)\alpha_{yy}+\big[\tfrac12-i\kappa-y\big]\alpha_y-g^2\alpha=0,\qquad
\alpha={}_2F_1(ig,-ig;\tfrac12-i\kappa;y).
$$

The other solution at $y=0$ would supply an incoming negative-frequency component. Gauss's summation at $y=1$ gives $\alpha_{\rm out}=\Gamma(c)^2/[\Gamma(c-ig)\Gamma(c+ig)]$, $c=1/2-i\kappa$. Using $|\Gamma(1/2+ib)|^2=\pi/\cosh(\pi b)$ then gives the exact result

$$
N_k=\frac{\sinh^2(\pi g)}{\cosh^2(\pi k\tau)}.
$$

At $g=0.2$, $\tau=0.18$, the lowest harmonic has $N_2=0.15361429$. Keeping shells through $k=96$ excludes less than $5.66\times10^{-43}$ in $a^4\rho$. This last bound uses $\operatorname{sech}^2 z\le4e^{-2z}$ and a geometric upper bound on consecutive terms $k^3e^{-2\pi k\tau}$. Independent integration of the original $(A,P)$ equations and instantaneous work from $-35\tau$ to $+35\tau$ checks the formula for six harmonics through $k=16$. The largest absolute discrepancies in occupation and work balance are $1.7\times10^{-15}$ and $3.8\times10^{-15}$ respectively. The reciprocal coupling gives the same occupations to $1.3\times10^{-15}$.

## Bounding the entire uncomputed tail for a finite switch

The rescaled wave $v=IA$ satisfies $v''+(k^2-V)v=0$ with $V=I''/I$. Assume $I>0$ is constant off an interval and is sufficiently smooth for $V,V'$ to vanish at its endpoints and $V''$ to be integrable without singular boundary terms. Define
$J(\eta)=\int_{-\infty}^\eta V(s)\,ds$, $A_1=\int|V|\,d\eta$, and

$$
R_2=-V''+2V^2+J V'-\tfrac12VJ^2,\qquad B_2=\int |R_2|\,d\eta .
$$

An incoming positive-frequency trial wave, with its constant $1/\sqrt{2k}$ normalization suppressed, is

$$
v_2=e^{-ik\eta}
 \left[1-\frac{J}{2ik}+\frac{-V+J^2/2}{(2ik)^2}\right] .
$$

It matches both incoming data and remains purely positive frequency after the switch. Substitution into the exact equation leaves the residual $e^{-ik\eta}R_2/(2ik)^2$. For the difference $d=v-v_2$, the retarded integral has kernel $\sin[k(\eta-s)]/k$ multiplying $Vd-\mathrm{residual}$. Iterating this Volterra equation bounds $\int |Vd|+\int|\mathrm{residual}|$ by $\exp(A_1/k)\int|\mathrm{residual}|$. Projecting the outgoing integral onto $e^{ik\eta}$ supplies a further factor $1/(2k)$. Consequently

$$
|\beta_k|\le \frac{e^{A_1/k}B_2}{(2k)^3},\qquad
0\le a^4\rho_{k>K}\le \frac{e^{2A_1/(K+1)}B_2^2}{128\pi^2 K^2}.
$$

The second inequality uses $\sum_{k>K}k^{-3}\le1/(2K^2)$ and includes both polarizations. A smooth compact switch admits arbitrarily many further integrations by the recurrence $f_0=1$, $f_{j+1}'=f_j''-Vf_j$ with higher coefficients initially zero. Truncation after $j$ gives $|\beta_k|\le e^{A_1/k}\|f_j''-Vf_j\|_1/(2k)^{j+1}$. Each constant depends on the actual profile and its conformal-time derivatives. A piecewise jump in $I$ instead squeezes every oscillator by the same magnitude at matching and produces divergent outgoing energy if extrapolated to arbitrarily high $k$.

For a finite numerical example take $y=\eta/2\in[0,1]$ and $d\ln I/dy=(630/4)y^4(1-y)^4$, extended by zero outside. This monotone $C^4$ coupling changes $\ln I$ by $1/4$. The program bounds the integrals for $A_1$ and $B_2$ by converting the two piecewise-polynomial integrands to exact rational Bernstein coefficients. After interval subdivision the integral of each absolute polynomial is bounded by the interval width times the average absolute Bernstein coefficient, using positivity and equal integrals of the basis functions. The resulting rational upper bounds are approximately $0.61935238$ and $11.49202623$. Direct propagation through $k=128$ gives $a^4\rho_{\le128}=0.01164346869$ and the remaining tail is below $0.000006443$. The tail inequality is analytic. The displayed partial sum uses numerical mode evolution. The rational bounds pertain only to this specified finite switch.

## What requires renormalization during the switch?

At a constant final $I$, the excess compares two states on the same outgoing geometry. The common state-independent local subtractions cancel, and the finite positive sum measures actual excitation energy relative to that specified reference. The absolute stress of the conformal photon vacuum can contain curvature and topology contributions. During a switch, the local electromagnetic vacuum terms and the scalar operator sourcing $I$ also depend on the spacetime-dependent coupling and its derivatives. A self-consistent absolute quantum stress and source require compatible renormalization of both operators, with local counterterms and their finite choices specified. Truncating the instantaneous oscillator energy to a frequency band does not supply those choices. [Study literature](../literature.md)

From the study root, run

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B checks/checks.py --output /tmp/CHECKS.json
~~~

[Calculation](checks.py) · [Results, exact rational bounds and cutoff checks](CHECKS.json)
