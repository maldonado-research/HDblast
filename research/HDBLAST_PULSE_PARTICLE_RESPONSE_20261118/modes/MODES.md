# A massive quantum spectator on the simulated shell history

The incident pulse changes a quantum field's Gaussian state even though the classical shell field stays far below the B1 particle-creation crossing. For a comoving momentum of $10 H_0$ at the initial scale factor, the mismatch between the positive-pulse state and its quiet counterpart reaches $2.69082\times10^{-6}$ at elapsed proper time $H_0\tau=0.15$. The same diagnostic was zero to recorded precision at 0.09. The comparison uses a specified free spectator and a common initial vacuum. [Input histories and model](https://doi.org/10.5281/zenodo.24354840) · [Replayed diagnostics](SUMMARY.json)

## Field, clocks and quantum state

Consider one real scalar spectator $\chi$ on the spatially flat four-dimensional induced geometry. Choose the curvature coupling $\xi=1/6$ in its Klein–Gordon action and the B1 mass law $\mu=m_\chi/H_0=100|\phi_b-0.5|$. This is a new test-field calculation on the published classical metrics. It omits a feedback term in the gravitational equations, interactions between spectator quanta, decays and thermalization. Common comoving coordinates have $a(0)=1$. With $s=H_0\tau$, $a=\exp(A_b-A_b(0))$, $ds/dT=\exp(B_b)/\rho_b$ and $d\eta=ds/a$, define $\kappa=k/H_0$. Throughout this note primes denote derivatives with respect to dimensionless conformal time $\eta$.

In dimensionless oscillator units the canonical mode $v=a\chi$ and $p=v'$ satisfy

$$v'=p,\quad p'=-\omega^2v,\qquad \omega^2=\kappa^2+a^2\mu^2,\qquad vp^*-v^*p=i.$$

There is no mass crossing in the input pulse. The differential equation still changes when the time dependence of $a\mu$ changes. All runs begin with the same field and geometry at the shell. The primary seed is the instantaneous vacuum, $v_0=(2\omega_0)^{-1/2}$, $p_0=-i\omega_0v_0$. The alternative WKB preparation includes the amplitude derivative, $p_0=(-i\omega_0-D_0)v_0$ where $D=\omega'/(2\omega)$. At the endpoint we analyze either seed with both the instantaneous basis and this WKB basis using $n_D=|\omega v-i(p+Dv)|^2/(2\omega)$, setting $D=0$ for the instantaneous basis. These finite-time choices concern a finite band; neither is an all-frequency renormalization prescription.

The primary integration band $5\le\kappa\le18$ has sampled physical frequency $\sqrt{(\kappa/a)^2+\mu^2}\le53.0833$ in units of $H_0$ on the positive trajectory through $s=0.15$. For the historical $b=6.49\times10^{-6}$, $M_5/H_0=b^{-1/3}=53.6107$. Using $\omega_{\rm phys}<M_5$ represents an additional conservative screen and the upper band edge is close to this scale. The wider $5\le\kappa\le80$ and individual $\kappa=25,50,100$ are formal sensitivity calculations. They do not obey the same screen. Results use one real field with no multiplicity factor.

## Relative squeezing and the best common analyzer

Compare the two states at the same $s$ and comoving $\kappa$ by identifying their normalized canonical coordinate pairs $(v,p)$. For the pulse and quiet complex positive-frequency solutions $z_p,z_q$, define

$$N_{\rm rel}=|v_qp_p-p_qv_p|^2.$$

This expression is invariant under independent mode rephasings and under a common real canonical transformation because it is the modulus of a two-dimensional determinant. The identification of variables across two different metrics matters. Identifying the original field and its momentum $(\chi,\Pi)=(v/a,ap-a'v)$ instead applies different canonical maps on the two backgrounds. The corresponding finite-band integrals at $s=0.15$ are $0.000245743$ and $0.000245586$, respectively. Each integral uses $d\nu=\kappa^2d\kappa/(2\pi^2)$, the density of modes per initial normalized comoving volume. $N_{\rm rel}$ compares two prepared states; it does not assign all the difference to produced pulse particles.

A sharp common-analyzer bound quantifies this distinction. In either shared phase-space identification the mode covariance is $C=\bigl(\begin{smallmatrix}|v|^2&\mathrm{Re}(v^*p)\\ \mathrm{Re}(v^*p)&|p|^2\end{smallmatrix}\bigr)$ with $\det C=1/4$. Direct expansion gives $\det(C_p+C_q)=1+N_{\rm rel}$. A common centered oscillator analyzer is a positive quadratic form $G$ with $\det G=1$, defining $n_G=(\mathrm{tr}(GC)-1)/2$. For $S=C_p+C_q$ the trace–determinant inequality gives $n_{G,p}+n_{G,q}+1\ge\sqrt{\det S}$. The choice $G=\sqrt{\det S}\,S^{-1}$ attains equality and gives the same count for both states, since $\det C_p=\det C_q$. Consequently

$$\inf_G\max(n_{G,p},n_{G,q})=\frac{\sqrt{1+N_{\rm rel}}-1}{2}.$$

At $\kappa=10$ and $s=0.15$ the positive-pulse pair has a floor of $6.72705\times10^{-7}$. A shared choice cannot make both states vacuum. Each state individually admits a centered pure-state vacuum basis of its own. This bound concerns independent centered analyzers of each real oscillator; a specified detector has its own response functional. The standard-library audit checks the determinant identity and both midpoint counts directly on 504 pairs of stored complex modes. [Covariance audit](audit_modes.py)

## Results at a common proper time

The primary endpoint is $s=0.15$, which occurs at $T=0.15001005$ for the positive pulse, before the direct arrival of its remote compensating controls at $T=0.152$. The quiet run also expands, so its instantaneous occupation is nonzero for the chosen seed.

| Initial scalar amplitude | $N_{\rm rel}$ at $\kappa=10$ | Pulse minus quiet instantaneous $n$ at $\kappa=10$ | $\int N_{\rm rel}\,d\nu$, $5\le\kappa\le18$ | $\Delta\langle\chi^2\rangle_{\rm band}/H_0^2$ |
|---:|---:|---:|---:|---:|
| +0.001 | $2.69082\times10^{-6}$ | $+3.21189\times10^{-5}$ | 0.000245743 | +0.000413647 |
| −0.001 | $2.68508\times10^{-6}$ | $-2.69382\times10^{-5}$ | 0.000245241 | −0.000392500 |
| +0.0005 | $6.72347\times10^{-7}$ | $+1.54118\times10^{-5}$ | 0.0000614044 | +0.000204177 |
| −0.0005 | $6.71629\times10^{-7}$ | $-1.41167\times10^{-5}$ | 0.0000613416 | −0.000198890 |

With an instantaneous initial state at $\kappa=10$, the positive-pulse and quiet instantaneous endpoint occupations are 0.000117675 and 0.0000855560. Changing only the endpoint analyzer to WKB gives a difference of $+1.70428\times10^{-5}$. This dependence is distinct from the state-pair mismatch. Changing the common initial preparation to WKB changes the positive-pulse band-integrated $N_{\rm rel}$ from 0.000245743 to 0.000246095.

The field variance comes directly from $\langle\chi^2\rangle_{\rm band}/H_0^2=\int |v/a|^2\,d\nu$. It requires no endpoint particle labels. The positive pulse raises the band-filtered variance by **0.0690%** relative to its quiet value of 0.59934846. The alternative WKB preparation gives a change of +0.000368632. Doubling the pulse amplitude between the two tested pairs multiplies the sign-even part of integrated mismatch by 3.99999 and the sign-odd variance change by 2.00003. These are finite-amplitude comparisons.

## Stress interpretation and numerical checks

For one conformally coupled scalar on the homogeneous flat background the spatially averaged band energy is $\rho_{\rm band}/H_0^4=\int (|p|^2+\omega^2|v|^2)/(2a^4)\,d\nu$. In the instantaneous basis this equals $\int\omega(n+1/2)/a^4\,d\nu$. On the primary positive pulse the excitation terms are 0.350031 and 0.253664 for pulse and quiet respectively. Their difference is +0.096367. The corresponding change of the finite-band zero-point term is −0.556301, giving a total finite-band difference of **−0.459935** in these units. The sign will depend on the comparison: the negative pulse gives an excitation difference of −0.0808124 and a total band difference of +0.477372. This calculation retains a fixed comoving band on both trajectories. It supplies neither a renormalized all-mode stress tensor nor the source's energy loss. Placing the field variance or the integrated relative squeezing into a Friedmann equation as radiation would require new physical input.

The solver uses cubic Hermite interpolation of $A_b,B_b,\phi_b$ and their saved time derivatives. It evolves the real linear transfer matrix of each oscillator with time-ordered, fourth-order commutator-free exponentials using the two Gauss nodes per coordinate interval. Proper and conformal clocks use the same Gauss rule, clipping at the specified common proper time. The maximum Wronskian error among the sampled modes is $2.26\times10^{-14}$. The finest two bulk grids change the primary integrated relative squeezing by 0.0138%. Switching to a second-order midpoint propagator on the same finest history changes it by 0.000333%. Dropping the remote compensating controls changes the finest integrated result by less than $2\times10^{-17}$ at this cut. A manufactured time-varying oscillator with an independent power-series solution gives error reduction factors near 16 for the fourth-order propagator and near 4 for midpoint. These are consistency checks at finite resolution. [Mode report](MODES.json) · [Integrator benchmark](BENCHMARK.json)

With NumPy installed, pass the preceding release's `evolution/wall_histories.zip` to the solver. The covariance audit itself needs only Python's standard library.

~~~sh
python -I -B spectator.py /path/to/wall_histories.zip > regenerated-modes.json
python -I -B audit_modes.py regenerated-modes.json
python -I -B benchmark.py
~~~
