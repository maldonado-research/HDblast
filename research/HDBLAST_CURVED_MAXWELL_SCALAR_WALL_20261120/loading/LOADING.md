# Let the Maxwell action determine the energy delivered to a moving wall

Two ordinary anti-de Sitter interiors share a round wall of radius $a$. Each side carries three equally excited, orthogonal Abelian fields. Their total stress is isotropic on the wall. The same potential and fluid transfer switch as in the [preceding study](https://doi.org/10.5281/zenodo.24386416) set the comparison with the earlier scalar-driven bounce. We now derive the scalar force and the stress at the boundary from the gauge action. A separate trial also removes the fluid transfer switch and puts the released energy back into the Maxwell field. [Direct electromagnetic return](RETURN_FIELD.md) · [Forward wave propagation](FORWARD.md)

## Which data are prescribed?

For the calculation in this note, the incoming Maxwell characteristic is specified **at the wall in its proper time**. The departing characteristic, scalar, surface energy, radius and gravitational mass respond dynamically. This choice supplies exact local junction data and by itself leaves the required source at smaller radius unspecified. Reconstructing the surrounding geometry shows that the source must continue to emit weakly after the direct compact wall pulse. The forward experiment in this study instead fixes the pulse at an emitting tube in advance and sets its incoming characteristic to zero after cutoff. [Inverse reconstruction](../identities/COLLAR.md)

Set the AdS length, tuned tension and speed of light to one, with $\kappa_5^2=6$. Use $X=\dot{\mathcal A}/a$, $Y=n\cdot D\mathcal A/a$ for a normal outward from the retained domain, $c=\mathcal A/a$, and $B=2c/a$. The surface electric and magnetic energies sum to $E_s=K(X^2+B^2)/2$ and the Maxwell pressure is $E_s/3$. The equal incoming and outgoing amplitudes per side are $I=(X-Y)/2$ and $O=(X+Y)/2$. Then $q=I^2-O^2$ flows into the wall from one side and $P_n=I^2+O^2-B^2/2$ is the bulk normal pressure. The surface coupling and the canonical scalar equation give

$$
P=KX,\quad \dot c=X-Hc,\quad
\dot P=4I-2X-2HP-4Kc/a^2,\quad
\ddot\psi=\tfrac12 K_{,\psi}(X^2-B^2)-(3H+\Gamma)\dot\psi-V_{,\psi}.
$$

When enabled, the separately postulated fluid satisfies $\dot e_r=\Gamma\dot\psi^2-4He_r$. The total added energy and pressure are $e=\dot\psi^2/2+V+e_r+E_s$ and $p=\dot\psi^2/2-V+(e_r+E_s)/3$. Starting from $a=2$, $H=0$, $\mu=4$, all fields zero and zero wall matter, we integrate

$$
\frac{\ddot a}{a}=(1+e)(1-2e-3p)-1-\frac{\mu}{a^4}-2P_n,\qquad
\dot\mu=-4a^4[H P_n+(1+e)q].
$$

Checks include $H^2=e(2+e)+\mu/a^4-1/a^2$, the unsquared positive-normal junction, the integrated work equation $d(a^3e)+p\,d(a^3)=2a^3q\,d\tau$ and the fluid transfer equation. Positive surface energy and nonnegative surface enthalpy hold for these positive couplings. Initial tension and bulk mass are already present.

## Drives, material and controls

The principal total local incident dose is $D=2\int I^2 d\tau=0.123$. The pulse uses the compact rising shape in [the planar calculation](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7030480769), duration $\delta[3+\log(39)/4+0.0625]$ with $\delta=0.001$. Its scale is $\sqrt{D/\delta}$. Initially $K=\delta$, with logarithmic slope $\sqrt8/(\delta\sqrt{0.12})$ and smooth caps at $K/\delta=1/40,40$. The calibration dose 0.12 stays fixed as the actual dose changes. The earlier potential begins rising at $\psi=0.004$ and reaches 0.111 at 0.012. Its descending edge spans $0.36<\psi<0.38$. In the fluid comparison $\Gamma$ rises from zero to 4 on this edge. The alternative Maxwell return has $\Gamma=0$ and a second positive rise in $K$.

In the fluid comparison, the fully aligned pulse reverses contraction near radius 1.878. The expanding radius grows from 3.059 to 3.717 after the fluid exceeds 90% of added energy, then contracts again. This is 0.195 e-folds. Removing the potential or the gauge interaction prevents the bounce. With damping switched off and the original capped $K$ alone, the independent fluid remains empty. Changes in dose can prevent exit from the plateau during the evolution to the radius cutoff. The base solver ends a case when it reaches radius 1.6 inward, 16 outward, or the specified time limit.

For a coherence control, introduce a second independent triplet odd under interchange of the two sides. This preserves mirror stress while permitting a nonzero odd share of the dose. The odd incoming characteristic is paired with an equal, opposite outgoing amplitude so its boundary potential and electric field vanish. It carries no net energy through the wall and adds $2I_{\rm odd}^2$ to each side's bulk normal pressure. Splitting the dose equally between even and odd sectors reproduces the energy budget of uncorrelated equal-power inputs in this idealized construction. At 99% even dose, the fluid-dominated expanding interval lasts about 0.68 e-folds. At 98% the wall bounces and reaches the radius-16 cutoff without the fluid exit. At 90% even or at equal splitting the wall reaches the inner cutoff during contraction. These are model-specific sensitivities to the incident state.

Two more pulses add a smooth negative lobe with width $40\delta$ or $500\delta$ and zero signed input area, retaining the same total incident dose. With $\Gamma=0$ and the late increment in $K$ set to 2000, neither reproduces the baseline expanding interval dominated by the electromagnetic field. The shorter pulse family reaches the inner cutoff without bouncing. The longer does bounce and crosses the descending edge, but misses that dominance threshold. Coherence, signed area, energy and scalar conversion are different controls.

## A short-pulse check

For a fixed-shape family scale the material width and pulse together, $K_0=\delta$ and $g\propto1/\delta$. In the formal small-$\delta$ limit the radius changes negligibly during the short loading and discharge, while the expansion rate can acquire a finite negative jump. Local magnetic terms and the finite scalar potential ramp do not contribute to the leading impulse. Write $D$ for incoming and $R$ for returned proper dose at the wall in this limit. The positive flux deposits $D-R$ and the normal pressure gives the sum $D+R$. The junction conditions after the impulse are

$$
e_+=D-R,\quad H_+=-(D+R),\quad
\mu_+=4+16[H_+^2-e_+(2+e_+)].
$$

An independent planar action calculation on the rescaled pulse gives $R=0.006119$, $e_+=0.116881$, $H_+=-0.129119$ and $\mu_+=0.307967$. The surface scalar carries the surviving stored energy, so its outgoing velocity is $+\sqrt{2e_+}$. Subsequent free motion across the potential and fluid switch gives 0.19352 expanding e-folds above 90% fluid. The archived finite-width sequence from 0.002 to 0.000125 approaches the limiting turning points. This limit changes a material scale along with pulse duration. At fixed material changing source width is a different physical test.

[Equations and integrated controls](loading.py) · [All results](LOADING.json)
