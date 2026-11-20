# Can the late electric field heat a cold population of charges?

We place a dilute, initially cold, neutral population of pre-existing positive and negative carriers on the archived expanding wall. Elastic deflections can turn the electric-field-driven streams into random motion. In this example the largest random particle energy among eight tested rates amounts to **0.00000526% of the wall electromagnetic energy** at the first expansion maximum. The two oppositely streaming species have zero total momentum even in the collisionless limit. Each species then stays cold in its own moving frame.

## Inputs and omissions

Use the finite-emitter run with $K_{\rm late}=13{,}200$, no extra radiation fluid, and the same wall time, charge-field and tension units as the parent study. The particles enter this test at $t_0=4.6$, the first archived tenth of a model time unit after the scalar completes its descent. The wall has $a_0=5.08648$ there and turns around at $t_f=12.42163$, $a_f=9.74401$. We specify a mass $m=0.1$ and charge-vector norm $|q|=1$ in model units. The two species together have initial rest energy $n_0m=10^{-4}e(t_0)$. Their total proper number density thereafter is $n=n_0(a_0/a)^3$. Charge creation and loss, photon production and chemical reactions are outside this calculation. This density, charge and mass have no measured cosmological normalization.

The three Abelian fields of the input have equal amplitudes and span orthogonal directions. The archived $X$ obeys $e_{\rm electric}=KX^2/2$ for the **sum** of the three. A charge vector of norm $|q|$ couples to one combination with electric force $|q|X/\sqrt3$. Its effective electric and magnetic fields are parallel to one invariant direction on the round three-sphere. For comoving momentum $\mathbf P=a\mathbf p$, the geometric spatial rotation from the round-sphere invariant frame contracts an antisymmetric connection with $p^j p^k$ and vanishes. The magnetic force rotates transverse momentum about the selected direction while preserving $P_\parallel$ and $|\mathbf P|$. Consequently between deflections

$$
dP_\parallel= b(t)\,dt,\qquad b=q\,\dot{\mathcal A}/\sqrt3,
\qquad |\mathbf P_\perp|=\mathrm{constant},
\qquad \mathcal A=a^2 B_{\rm mag}/2,\quad \dot{\mathcal A}=aX .
$$

Here the sign of $q$ selects either member of the neutral pair. The kinetic equations use the signed field from the imported solution. Hubble redshift is included through $\mathbf p=\mathbf P/a$. The directional rotations are Poisson events at a **chosen** constant rate $\nu$ per unit wall time. Each event samples a uniform direction at unchanged local momentum magnitude, as for elastic scattering from an ideal infinitely heavy, comoving background. Such centers can absorb momentum with negligible recoil energy. No microscopic collision rate or thermodynamic heat reservoir has been derived. Deflections do work only by changing the subsequent motion in the electric field.

## An exact momentum diagnostic

Write $u=\langle P_\parallel\rangle$, $S=\langle P^2\rangle$, $D=\langle P_\parallel^2-P_{\perp,1}^2\rangle$, $C=\langle P^2P_\parallel\rangle$ and $F=\langle P^4\rangle$ for the positively charged carriers, with zero initial moments. Azimuthal symmetry makes the two transverse components equivalent. The Markov collision generator gives

$$
\dot u=b-\nu u,\quad \dot S=2bu,\quad
\dot D=2bu-\nu D,\quad
\dot C=b[S+2(S+2D)/3]-\nu C,\quad
\dot F=4bC .
$$

In particular the trace of the momentum covariance obeys $d(S-u^2)/dt=2\nu u^2$. All these **moment** identities hold at any charge-to-mass ratio under the assumed wall-frame rate and elastic rule. For the specified moderately slow carriers, the directed and random kinetic energies per carrier are respectively $u^2/(2ma^2)$ and $(S-u^2)/(2ma^2)$. Deflections preserve the instantaneous total energy. In the nonrelativistic accuracy of those tables, the change in mean total kinetic energy per carrier equals $\int bu/(ma^2)\,dt-\int HS/(ma^2)\,dt$. The code checks this balance, the variance identity and collisionless transport separately. The maximum allowed $p/m$ is 0.03336 on this window. The nonrelativistic single-particle kinetic energy exceeds the exact one by at most 0.028% relative here.

| Rate $\nu$ | Fraction of carriers' kinetic energy in random motion at turnaround | Random energy divided by all electromagnetic energy on the wall |
|---:|---:|---:|
| 0 | 0 | 0 |
| 0.03 | 0.1261 | $1.66\times10^{-8}$ |
| 0.1 | 0.3493 | $3.94\times10^{-8}$ |
| 0.3 | 0.6769 | $5.26\times10^{-8}$ |
| 1 | 0.9149 | $3.11\times10^{-8}$ |
| 3 | 0.9755 | $1.21\times10^{-8}$ |
| 10 | 0.9930 | $3.78\times10^{-9}$ |
| 30 | 0.9977 | $1.27\times10^{-9}$ |

Frequent deflections also suppress the current doing the electrical work. A large *fraction* of random motion among the particles can therefore coexist with very little total energy converted. At $\nu=0.3$ the central fourth-moment ratio is 0.977 relative to an isotropic Maxwell distribution with the same summed variance. Nevertheless, $e^{-\nu(t_f-t_0)}=9.57\%$ experience no deflections and remain a directed monoenergetic component. The parallel-minus-one-transverse central pressure difference is 28.8% of the summed random pressure components. This provides two explicit checks against reading a near-Maxwell fourth moment as equilibrium. Number-changing reactions, interactions that equilibrate particle energies, and the electromagnetic fluctuation spectrum are absent.

## A ceiling on energy capture by pre-existing charges

An elastic deflection preserves $|\mathbf P|$. Each interval of electric acceleration can increase it by at most $|q|\,|d\mathcal A|/\sqrt3$, even relativistically. Thus with $L=\mathrm{TV}_{[t_0,t_f]}(\mathcal A)/(\sqrt3a_f)$ and $r= n_f m/e_{\rm EM}(t_f)$,

$$
\frac{\rho_{\rm kinetic}(t_f)}{e_{\rm EM}(t_f)}
\le r\left[\sqrt{1+ (|q|L/m)^2}-1\right] .
$$

On this background the reconstructed potential decreases from $0.0404722$ to $-0.0158273$, with numerical total variation $0.0562995$. Integrating $|aX|$ independently agrees within $1.3\times10^{-10}$. With $r=2.53218\times10^{-4}$, even the maximum relativistic pickup at $|q|/m=10$ is only $1.4085\times10^{-7}$ of the wall electromagnetic energy. For this monotone field the collisionless stream attains the bound as directed energy. Keeping the same seed rest abundance, reaching even **1%** would require $|q|/m\ge12{,}135$ in these model units. This is a necessary condition on the imported history. It does not justify extending the weak, one-way tracer approximation to large charges. Additional carriers supplied by pair creation or injection change this bound.

At the chosen dilute abundance the total rest plus kinetic energy density of the tracers stays below $2.34\times10^{-4}$ of the surface matter. The two charges add to the electrical current while their momentum densities cancel. A check on the scale of this perturbation integrates the bound $a^2|j|\le a^2|q|n P_{\max}/\sqrt{m^2a^2+P_{\max}^2}$ over the window, with $P_{\max}=|q|\,\mathrm{TV}_{[t_0,t]}(\mathcal A)/\sqrt3$. The resulting impulse is $2.46\times10^{-7}$ or less of the largest projected background displacement $a^2K|X|/\sqrt3$. This is a scale check for the extra term in the integrated Maxwell equation. The calculation imports the geometry, scalar and gauge histories and leaves their dynamical response open.

## Replay

Run `python -I -B kinetic.py --input ../background/trace.csv.gz --output KINETIC.json` from this folder. Python's standard library is sufficient. Piecewise-linear signed potential and piecewise-log-linear radius are used within each input interval. Taking every ninth trace point changes the $\nu=0.3$ random energy by $1.1\times10^{-7}$ relative. An independent seeded sample of 65,536 carriers follows electric acceleration and individual random deflections exactly on those interpolants. It agrees with the moment calculation's random energy to 0.023% and its normalized fourth moment to 0.059%. The calculation records the sample difference separately from the numerical time-step check. [Solver](kinetic.py) · [Full results](KINETIC.json)
