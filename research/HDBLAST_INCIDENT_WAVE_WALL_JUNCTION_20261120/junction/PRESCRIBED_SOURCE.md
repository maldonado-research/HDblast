# Can the source specify its pulse in advance?

In the [canonical scalar calculation](SCALAR_CAPTURE.md), the intended scalar force and conservation laws determine the retarded mass profile that supplies its work. For an independent test we specify the source mass directly at emission. Once launched, the prescribed null pulse propagates outward to the wall. We use the same wall potential and compare against a wall with the potential set to zero.

Use the closed symmetric branch with $a=2$, $H=e=\psi=0$ and $\mu=4$ initially, in AdS units. For $s=u/\Delta u$ in $[0,1]$ set $\mu=4[1-\eta(35s^4-84s^5+70s^6-20s^7)]$, and hold the mass fixed before and after. All energies are nonnegative for $0\le\eta\le1$. We assume coherent absorption into the scalar kinetic energy $K$ and select the positive field-velocity branch. This is a property assigned to the wall. The microscopic interaction that would implement it is an open question.

The entire pulse arrives while $V=V'=0$. The retarded clock, local power and matter equations in wall proper time are

$$
\dot u=[a(1+K+H)]^{-1},\quad
Q=-\frac{\mu_u\dot u^2}{2a^3},\quad
\dot K=Q-6HK,\quad \dot\psi=\sqrt{2K}.
$$

Integrating $K$ avoids division by the initially vanishing scalar velocity. We integrate the temporal Israel equation independently of its spatial constraint, then evolve the canonical scalar on the same rising potential as in the other calculation. The stationary emitter at radius 1.6 is outside the initial apparent horizon. The wall remains outside that emitter on the reported bouncing trajectory.

For $\eta=0.94$, $\Delta u=0.005$, $V_{\rm plateau}=0.111$, the wall receives 0.941495 in integrated proper energy, while the two source sides emit a combined proper energy of 1.148051 in the archived units. The mass ends at 0.24 on each side. The scalar reaches the positive plateau, the wall turns at radius 1.899614 and proper time 0.791259, and it passes radius 2.05 while expanding. The $\eta=0.85$ pulse fails to reach the plateau and the wall contracts to the source. At $\eta=0.97$ it reaches the plateau and bounces. For the identical 94% source profile with $V=0$, the wall contracts. This independently specified beam and the force-normalized beam use different profiles and reach slightly different turning points.

We then pass the outgoing solution at radius 2.05 into the independently implemented [finite-plateau continuation](RADIATION_EXIT.md). For the preassigned 94% beam, the radiation-like fluid passes 90% of material energy at proper time 3.7048, radius 3.5914. Expansion continues to radius 4.9336 at time 7.0413 before reversing. The expanding interval above 90% spans 0.3175 e-folds. The model consumes the specified source energy and transfers energy internally according to its chosen damping law. It does not calculate particle populations or a temperature.

[Prescribed-source computation](prescribed_source.py) · [All results](JUNCTION.json)
