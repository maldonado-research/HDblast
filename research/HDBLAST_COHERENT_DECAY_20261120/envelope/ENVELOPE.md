# A finite-mass test of the charged-particle interpretation

Changing the daughter mass changes the relation between the dimensionless fluid calculation and a microscopic particle model. This calculation asks whether a heavier family can still heat the expanding wall while keeping two proposed couplings small. It keeps the energy taken from the source and returned to the driver in the evolution. The result depends on varying the width and amplitude of the parent source with mass.

## Resolving the source, averaging the free oscillations

Write $\mu=mL$, $M=2\mu R$ in wall units and $I=a^3\rho_\chi/M$. The source changes the minimum of the oscillator while $\psi$ crosses a small prescribed interval. Before and during that interval the numerical evolution resolves $\chi$ and its velocity. All histories begin with the parent and daughter sectors empty. The new family fixes $\mu v_0=0.008$ and $\mu\,\Delta\psi=0.24$. Near the drive, the change of coordinates $u=\mu(\psi-\psi_s)$, $c=\mu\chi$, $s=\mu(\tau-\tau_s)$ keeps $c$, $dc/ds$ and the dimensionless bump shape finite. It also keeps the $\chi$ force on the driver in its fast-time equation. Changing only the masses while holding the source fixed would define a different experiment.

After the minimum has returned to zero, the code matches at the first resolved negative-to-positive zero crossing of parent pressure. It copies the actual energy into $I$. The subsequent proper-time equations use
$$
\dot I=-\Gamma I,\qquad \rho_\chi=MI/a^3,\qquad
\langle p_\chi\rangle=0,\qquad P=\Gamma\rho_\chi .
$$
For a harmonic cycle $\langle\dot\chi^2\rangle=\rho_\chi$ and $\langle\chi^2\rangle=\rho_\chi/M^2$. Hence the driver work from a changing mass is $\dot M\rho_\chi/M$, with the opposite reaction in the driver equation. The Friedmann equation, electric field, gas and dynamical finite bulk use the evolving sources throughout. The code terminates with an error if the drive is re-entered. The approximation assumes slow relative variation on one radian of parent phase. Its output records the largest $|H|/M$, $|\dot M|/M^2$, and $\Gamma/M$ after matching. A full phase-resolved solver at the same spatial and temporal settings supplies the benchmark at $\mu=400$. The high-mass family extrapolates a post-drive approximation checked at that anchor. It does not implement fermion dynamics during the drive.

For a harmonic amplitude $A=\sqrt{2\rho_\chi}/M$, a random phase gives $\langle\chi^2\rangle=A^2/2$. Births in the local law are weighted by $\dot\chi^2$, so the birth-weighted mean is $A^2/4$. The archive includes both the native and averaged integrals against comoving daughter rest production. [Resolved equations](../coupled/JOINED.md) · [Averaged continuation](build.py)

## Joint checks at one action scale

Use the project normalization $S_0=\sigma L^4$ and three equal Dirac flavors sharing a prescribed total width. Then $y_i^2=4\pi\kappa/(3\mu)$, $q=\sqrt3 Q\mu$, $\alpha_g=q^2/(4\pi S_0 K)$ and $|\delta m_i|/m=|y_i|\sqrt{S_0}|\chi|/\mu$. A larger $S_0$ weakens the gauge coupling and increases the mass shift. Equivalently,
$$
\frac{|\delta m_i|}{m}\sqrt{\alpha_g}
=|Q\chi|\sqrt{\frac{\kappa}{\mu K}} .
$$
For fixed trace and equal couplings this product is independent of $S_0$. The screen uses an illustrative gauge limit $\alpha_g\le0.1$ and daughter mass shift at most $0.005m$, one tenth the late threshold margin. These are declared comparison tolerances. Even a small shift can move the exact opening time while $M$ is crossing $2m$.

For the whole-history gauge floor the code uses the analytic lower bound $K\ge0.001/40$ of this positive gauge profile. For the mass ceiling it uses the maximum across accepted time steps of $|v|+\sqrt{(\chi-v)^2+(\dot\chi/M)^2}$ during the resolved source and the oscillation amplitude afterwards. These give a conservative gauge condition and a numerical amplitude check. A separate screen uses $\nu L=18$ and the single elastic s-wave unitarity bound in low-drift, ultrarelativistic hot states. In those states the action must exceed $6\nu\gamma\mu^3\theta^2/(\pi D_{\rm loc})$ if that channel is to produce the fluid drag. The actual cross section can be well below its unitarity ceiling. Each run uses one action to evaluate all three numerical restrictions. The field strength $|eE|/m^2=|QX|/\mu$ is independent of this action. [Finite-history thermal-rate assumptions](../kinetics/KINETICS.md) · [Screen calculation](envelope.py)

## Finite-mass results

The first line uses the charge of the fully resolved reference. The other listed runs have $Q=10^7$ and 24 early and six late radial elements, except where indicated. $G$ is the largest expansion factor continuously above the 90% random-energy threshold. The last column divides the allowed action ceiling for the mass shift by the larger of the gauge and elastic capacity floors. A value above one leaves an interval under these specified tests. A value below one for conservative envelopes alone cannot exclude every microscopic realization.

| $mL$ | $Q$ | $G$ | Largest sampled electric force in mass-squared units | Action ceiling / larger floor |
|---:|---:|---:|---:|---:|
| 400 | $2\times10^8$ | 2.813 | $4.23\times10^6$ | $1.02\times10^{-17}$ |
| 400 | $10^7$ | 2.512 | $2.11\times10^5$ | $4.09\times10^{-15}$ |
| 4,000 | $10^7$ | 2.276 | $2.11\times10^4$ | $2.99\times10^{-12}$ |
| 40,000 | $10^7$ | 2.255 | $2.11\times10^3$ | $3.00\times10^{-10}$ |
| $4\times10^6$ | $10^7$ | 2.253 | 21.15 | $3.00\times10^{-6}$ |
| $4\times10^8$ | $10^7$ | 2.253 | 0.2115 | 0.0300 |
| $4\times10^{10}$ | $10^7$ | 2.253 | 0.002115 | 300.1 |

At $mL=400$ the averaged calculation and phase-resolved 24-by-6 comparison differ by 0.0284% in daughter rest production, 0.0513% in electric work, 0.415% in hot logarithmic radius growth and 0.096% in the birth-weighted $\chi^2/K$ integral. At $4\times10^{10}$ the maximum post-match $|H|/M$, $|\dot M|/M^2$, and $\Gamma/M$ are respectively $6.1\times10^{-12}$, $1.22\times10^{-10}$ and $1.42\times10^{-11}$. At this mass, halving the step gives $G=2.2523$ and using 48 early and eight late elements gives $G=2.2974$. At $4\times10^8$, lowering $Q$ to $3\times10^6$ reduces $G$ to 1.895. The charge still matters even when the parent oscillates rapidly.

With the scaled source the measured peak amplitude is approximately proportional to $1/\mu$ and the macroscopic gas histories at large $\mu$ are similar. At fixed $Q,\kappa$ the gauge floor then scales as $\mu^2$, the mass-shift ceiling as $\mu^5$, and the low-drift s-wave capacity floor as $\mu^3$. This explains the roughly 10,000-fold gain in the ratio between the latter ceiling and floor when the mass rose 100-fold. It relies on the chosen family.

At $mL=4\times10^{10}$ a trial $S_0=3.609\times10^{48}$ gives a conservative whole-history gauge coupling bound $4.24\times10^{-10}$ and a sampled oscillator mass-shift bound $1.21\times10^{-3}$ as a fraction of the daughter mass. The maximum thermal elastic rate at this action in the most demanding tested hot state is 17.3 times the imposed drag. A single elastic channel would therefore need to deliver at least 5.8% of its thermal unitarity ceiling there. Neither small gauge and Yukawa couplings nor room under this ceiling generates that rate. A separate magnetic mode count checks whether enough distinct fermion labels could accept the births across six successive windows. [Available-state bounds and assumptions](../kinetics/MANY.md)

For the high-mass example both the source width in the rolling field ($6\times10^{-12}$) and amplitude ($2\times10^{-13}$) are small by construction. Their microscopic origin, the absolute length scale $L$, the particle species and the actual cross sections remain undetermined. Matching at a pressure crossing reduces the pressure discontinuity. The finite-phase action has corrections of order $|H|/M$ and $|\dot M|/M^2$. Averaging cannot resolve quantum occupations, resonant selection rules or direct electromagnetic pair production.

From the capsule root run `python -I -B envelope/envelope.py --output /tmp/ENVELOPE.json`. [Numerical results](ENVELOPE.json)

## A precision check on the narrowest source

A time-resolved state count at $mL=4\times10^{10}$ found brief moments when too few empty states cross the allowed birth momentum, even though the six cumulative windows passed. Raising the mass by a factor of 100 is a targeted test of that failure. In the large-mass regime a thin resonant shell has area proportional to $\mu^2$ and its momentum sweep speed scales as $\mu$. Requested counts scale as $S_0/\mu$, so the action ceiling from a given instantaneous sweep scales as $\mu^4$. The action needed for hot elastic drag scales as $\mu^3$. The final test holds $Q=10^7$ and sets $S_0$ to **five times** the drag floor. This is a different choice from the geometric mean of the mass and drag bounds in the earlier example. [Time-resolved birth test](../kinetics/SWEEP.md)

At $mL=4\times10^{12}$ the drive spans only $6\times10^{-14}$ in the rolling field coordinate. In a trial using the absolute coordinate, arithmetic rounding changed the logarithmic hot growth from about 0.812 to 0.861. The archived run integrates the local coordinate $z=(\psi-\psi_s)/\Delta\psi$ directly in this interval, with $\dot z=\dot\psi/\Delta\psi$. The force, the driver's opposite reaction, and the wall trajectory remain dynamic. This gives hot growth factors 2.2529, 2.2975 and 2.2524 on the reference, finer spatial grid and half time step. The code also runs the centered formulation at $mL=4\times10^{10}$ as an overlap check. The trial action at $4\times10^{12}$ is $1.042\times10^{54}$. At that action the conservative gauge bound is $1.47\times10^{-11}$, the sampled fractional mass-shift bound is $6.46\times10^{-6}$, and the measured peak electric strength in mass-squared units is $2.12\times10^{-5}$. The physical length $L$ is still free.

An independent integration checks whether the rapidly varying source still transfers the intended energy. In fast time the expansion and smooth exterior forces enter at order $1/\mu$. Setting $g_0=0.008$, $w_0=0.24$, $m_0=1.9$ and writing $b(u)$ for the unit bump on $[0,w_0]$, the leading inner problem obeys
$$
c''=-m_0^2(c-g_0 b),\quad
u''=m_0^2(c-g_0b)g_0 b',\quad
u'^2+c'^2+m_0^2(c-g_0 b)^2=\pi_{\rm in}^2 .
$$
Primes on $u,c$ denote fast time and the prime on $b$ denotes its argument. The pre-drive plateau conserves $a^3\dot\psi$ and fixes $\pi_{\rm in}=0.122878245$ from a separately sampled trajectory. Using $u$ as the independent coordinate eliminates the narrow absolute width. The outgoing neutral energy is $0.0002716758554$. Doubling the inner resolution changes it by $4.1\times10^{-12}$ in relative terms. The high-mass centered runs reproduce it within $2.0\times10^{-6}$ across the three tested grids. This validates source energy transfer even though the wall equations then use a cycle-averaged parent. [Inner source calculation](inner.py)
