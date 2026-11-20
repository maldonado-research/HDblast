# How much heat can the finite supply of charges produce?

In the principal run, electrical work supplies **56.5102** units of comoving energy by the first hot point. Decay injects **0.0000558** units above daughter rest energy by then. The model's ability to heat depends on whether charges appear while there is an electric field that can do work on them. A second test asks whether a finite observation time and field strength are compatible with interpreting the model's decay as a local vacuum process. [Coupled solution](../coupled/JOINED.json) · [Audit results](EXPOSURE.json)

## Energy, number and a finite-time bound

We use dimensionless wall proper time $\tau$, scale factor $a$, daughter mass $\mu=mL$, parent mass $M$ in the same units, and $R=M/(2\mu)$. Let $q_\chi=\chi-v(\psi)$, $e_\chi=(\dot\chi^2+M^2q_\chi^2)/2$ and $p_\chi=(\dot\chi^2-M^2q_\chi^2)/2$. In the modeled local damping law,

$$
\beta^2=(1-R^{-2})_+,\qquad \Gamma=\kappa R\beta^3,\qquad P_\chi=\Gamma\dot\chi^2.
$$

The dissipated power $P_\chi$ is nonnegative. The scalar receives $U_{,\psi}\dot\psi$ from the driver, spends $3H(e_\chi+p_\chi)$ in expansion, and loses $P_\chi$ to pairs. The gas gains $P_\chi+JX$ and spends $3H(e_d+p_d)$. Here $X$ is the tangential electric field. The moving minimum and varying mass transfer energy internally through the driver's opposite force. Electric work also appears with the opposite sign in the electromagnetic equation.

The model assumes creation in pairs with total energy $M$ per pair, zero mean directed birth momentum, and a fixed rest energy $2\mu$. Thus the comoving daughter rest energy obeys $\dot D=a^3P_\chi/R$. If $N_{\rm pair}$ is the dimensionless integrated pair count in the solver, $D=2\mu N_{\rm pair}$. The physical number over a whole three-sphere adds the factor $2\pi^2S_0$ of the previous normalization. The comoving rest-energy birth rate also equals $-2\mu\,\dot I_\chi|_{\rm damp}$ for instantaneous coherent action $I_\chi=a^3e_\chi/M$. Work by the driver and pressure can separately change $I_\chi$. Equating the dissipated energy divided by $M$ to a pair count is the model's on-shell birth assumption. [Parent normalization](https://github.com/maldonado-research/HDblast/pull/191)

Let $W_{\rm kin}=\int a^3P_\chi(1-1/R)d\tau$ be true birth energy above rest, $W_J=\int a^3JX\,d\tau$ signed electrical work, and $W_p=3\int Ha^3p_d\,d\tau$ pressure work. For the positive Taub–Mathews gas, its intrinsic random energy density $e_{\rm int}$ obeys $a^3e_{\rm int}\le a^3e_d-D$. Starting with empty daughter and parent fields, the gas energy identity and the speed bound on current therefore give

$$
a^3e_{\rm int}+W_p\le W_{\rm kin}+W_J,\qquad
|W_J|\le |Q|\int D|X|\,d\tau .
$$

The second inequality credits every available charge moving at light speed in the most favorable direction. With an upper mass ratio of 1.05, even the most energetic allowed individual birth gives at most 5% of its daughter's rest energy as kinetic energy.

For an observation time $t$, set $A=W_{\rm kin}/\kappa$ and $B=\int_0^t D|X|d\tau/\kappa$ when $\kappa>0$. Substituting the positive number source and exchanging the two integrals gives

$$
A=\int_0^t a_s^3\dot\chi_s^2\beta_s^3(R_s-1)\,ds,\qquad
B=\int_0^t a_s^3\dot\chi_s^2\beta_s^3\left(\int_s^t |X(u)|\,du\right)ds .
$$

A carrier created late receives credit only for the electric exposure still available after its birth. With $\alpha_\Sigma=\kappa/\mu$, a necessary condition for the specified heat is $[a^3e_{\rm int}+W_p]/\mu\le\alpha_\Sigma(A+|Q|B)$. For the measured signed electrical work alone, $|Q|\alpha_\Sigma\ge |W_J|/(\mu B)$. The solver accumulates $A$ and $B$ at its internal steps. These bounds test whether a choice of parameters can reproduce the *specified field, background and condensate histories*. Changing parameters in a new evolution can alter those histories.

| Principal run | First hot point | First hot doubling |
| --- | ---: | ---: |
| Proper time | 6.87263 | 14.35652 |
| Comoving daughter rest energy | 0.001118 | 0.002125 |
| $A$ | $1.39457\times 10^{-6}$ | $2.65315\times10^{-6}$ |
| $B$ | $5.85399\times 10^{-8}$ | $5.91519\times10^{-8}$ |
| Measured $W_J$ | 56.51022 | 56.51028 |
| Speed-limited $W_J$ ceiling at the actual $Q$ | 468.32 | 473.22 |
| Required $\lvert Q\rvert\alpha_\Sigma$ from work | $2.41332\times10^6$ | $2.38836\times10^6$ |

The first hot point uses 12.07% of its generous current ceiling. The run sets $Q=2\times10^8$, $\alpha_\Sigma=0.1$. At the first saved hot state, $\tau=6.88$, the independently evaluated inequality for attaining a 90% internal share gives $|Q|\alpha_\Sigma\ge2.4110\times10^6$ if $\alpha_\Sigma\le1$. For this bound we credited the maximum allowed birth contribution $A$ at fixed history. Using the same saved state for all terms preserves the small difference between gas and internal energies across the threshold.

The exposure bound itself assumes a single coherent parent, fixed daughter rest mass, the prescribed current law and causal speeds. Applying it jointly to a microscopic Yukawa model adds assumptions about its flavors and coupling normalization. For three equal Yukawa couplings, the largest mass shift on the coherent mean and the maximum allowed local gauge coupling bound $|Q|\sqrt{\alpha_\Sigma}$ separately. With unequal couplings, a bound on the largest shift applies to at least one flavor. It constrains a three-equal-fluid interpretation only if all three modeled channels must retain their masses. On the principal trajectory the largest saved $|\chi|/\sqrt K$ between opening and first hot with positive birth power is $1.34495\times10^{-6}$. For a Yukawa interaction $y_i\Phi\bar{\Psi}_i\Psi_i$, $\Phi L=\sqrt{S_0}\,\chi$, so the fractional shift of the positive bare mass is $|y_i|\sqrt{S_0}|\chi|/\mu$. Three equal couplings obey $3y_i^2=4\pi\alpha_\Sigma$. Eliminating $S_0$ gives $|\delta m_i|/m=|Q\chi|\sqrt{\alpha_\Sigma/(K\alpha_g)}$. Even allowing a 5% shift, $\alpha_g\le1$ and $\alpha_\Sigma\le1$, these samples require $|Q|\alpha_\Sigma\le |Q|\sqrt{\alpha_\Sigma}\le37{,}176.1$. The measured electrical work requires at least $2{,}413{,}322.5$, a **64.9-fold incompatibility** on that finite-mass trajectory. A stronger shift, a change in which flavors are active, interactions that protect the mass or a different evolution changes the test. [Finite-mass continuation](../envelope/ENVELOPE.md)

## Threshold resolution and strong electric fields

A perturbative scalar Yukawa width in flat space gives $\Gamma/M=\alpha_\Sigma\beta^3/2$, using the sum over accessible flavors. At the final mass ratio the principal run gives **0.001417** for $\Gamma/M$ and **0.02977** for $\Gamma/(M-2\mu)$. The bare mass changes at $|\dot M|/M^2=0.01083$ at opening. These numbers concern the prescribed masses and local vacuum width. A time-dependent mass, occupations, coherent driving and electric fields can affect actual production.

One way to quantify finite time at a threshold is the accumulated phase $\Xi(t)=\int (M-2\mu)_+d\tau$. A value much less than one gives poor frequency resolution of the excess energy. The crossover to a locally resolved channel depends on the actual time-dependent quantum evolution. The principal run reaches $\Xi=25.11$ radians by its first hot point. The sampled rest ledger assigns about **5.36%** of the rest mass present then to the interval before one radian and **32.34%** to the interval before one cycle. The historical low-mass, fast-damping example reaches only **0.1802 radians** by its first hot point and has $\Gamma/M=28.35$ on the plateau. A weak, sharply resolved vacuum-decay width is outside its domain there. For a fixed homogeneous oscillator, $\Gamma+3H>2M$ is the overdamped regime. Its slowly decaying root approaches $-M^2/(\Gamma+3H)$ for large damping, so increasing a friction coefficient can leave potential energy stored longer.

The same charge that supplies electrical work also couples directly to the field. The force law $d(p/m)/d\tau=QX$ in the absence of drag and expansion fixes $|eE|/m^2=|QX|/\mu$. At scalar exit the principal run has $|X|=0.21172$, an electric share of 99.99984%, and **$|eE|/m^2=1.059\times10^5$**. The flat, constant electric field factor $\exp[-\pi m^2/|eE|]$ would have exponent magnitude $2.97\times10^{-5}$. Scalar exit precedes the decay opening by 0.003838 wall times. For scale comparison, the instantaneous inverse square root of $|eE|$ there is $7.68\times10^{-6}$ wall times. The changing field, duration, curvature, available modes, occupation and backreaction enter a field-induced production calculation. The local field strength leaves the actual Schwinger pair yield undetermined.

## A necessary condition on continued hot expansion

The open-wall curvature budget from the hot-growth study can include a coherent parent with oscillatory pressure. With wall density and pressure $e,p$, bulk mass parameter $m_b$ and incident/outgoing amplitudes $I_b,O_b$, the exact identity is

$$
C=a^2[1+(aH)^2]=a^4(2e+e^2)+m_b,\qquad
\frac{dC}{d\ln a}=2a^4[e-3p-e(e+3p)+B_{\rm mag}^2-X^2-(I_b-O_b)^2] .
$$

For nonnegative scalar potentials, all canonical scalars together contribute at most four times their energy density to $e-3p$. The positive charged-gas closure satisfies $3p_d\ge e_{\rm int}$. The sheet contributes no trace and $B_{\rm mag}^2\le 2 e_{\rm sheet}/K$. If the internal fraction stays above $f\ge 3/4$ during expansion, $K\ge1/2$ and $m_b\ge0$, then $e-3p+B_{\rm mag}^2\le4(1-f)e$, $3p\ge(4f-3)e$ and $2a^4e\le C$. Integrating $dC/d\ln a\le4(1-f)C$ bounds a continuous growth factor $G$ by $G^{4f-2}\le 1+(a_hH_h)^2$, where $h$ denotes the start of the interval. A 90%-hot doubling needs $a_hH_h\ge1.4253$. This general bound allows negative instantaneous scalar pressure and permits bulk exchange. Passing it gives only an upper limit on possible growth. [Open-wall derivation](https://github.com/maldonado-research/HDblast/pull/190)

Run `python -I -B audit/exposure.py --output /tmp/exposure.json` from the study root. The script reports event-level native work and separate heat inequalities at actual saved states. It also checks the gas and scalar energy identities and rest/pair-number normalization on every supplied trace. On the principal trace $D=2\mu N_{\rm pair}$ agrees within $9.5\times10^{-16}$ in comoving units. The native integrator checks and reconstructed saved-state residuals use different time samplings. [Analysis code](exposure.py)
