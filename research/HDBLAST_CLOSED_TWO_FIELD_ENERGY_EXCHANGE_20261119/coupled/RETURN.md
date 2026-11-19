# Exchange through the first polar return

The first passage can give the wrong sign for cumulative direct energy transfer even a short time later. With a displacement-seeded daughter the rolling parent had taken 1.492% of the daughter’s initial energy by first passage. By its return through the same minimum the parent had supplied a net 4.255% of the daughter’s initial energy. The daughter’s total energy nevertheless fell by 17.284% between launch and return as the wave did work on the changing geometry.

We continued the solutions with $A=+0.001$ on the potential with its minimum at $\phi_m=-0.5$. We mark three events on the same polar normal worldline. The parent field passes upward through $q_m$, turns at a zero of its proper-time velocity, then passes back downward through $q_m$. The two crossings mark half a local oscillation. The transfer is integrated over the **entire proper spatial volume** between the times of these events. A local measurement at the pole would require a regional flux calculation. [Initial data, field equations and balance law](COUPLED.md)

## Where the exchange changes sign

Each signed number in the table is direct coupling work into the daughter, divided by its own initial energy and expressed as a percentage. The two phases have closely matched initial integrated energies, approximately $1.667\times10^{-5}$ of the reference energy for $s=0.01$.

| Daughter launch | Seed $s$ | Launch to first crossing | First crossing to turn | Turn to return | Launch to return |
|---|---:|---:|---:|---:|---:|
| Velocity | 0.01 | +0.8511% | +1.0181% | −0.8061% | +1.0631% |
| Velocity | 0.005 | +0.8515% | +1.0177% | −0.8072% | +1.0620% |
| Displacement | 0.01 | −1.4915% | +0.3850% | +5.3615% | +4.2550% |
| Displacement | 0.005 | −1.4915% | +0.3850% | +5.3614% | +4.2549% |

For the displacement launch at $s=0.01$, the daughter does positive pressure work of 21.539% of its initial energy up to the return. For the velocity launch the integrated pressure work is −0.0902% and direct transfer is +1.0631%, leaving the daughter with a net energy gain of 1.1534%. The volume average of the expansion is still positive at both returns, whereas part of the slice has begun contracting. Wave pressure can have either sign. Thus total expansion alone does not fix the sign of the integrated pressure work. The new endpoint does not determine a secular rate, thermalization, or transfer over a complete cycle.

## Event and numerical checks

For $s=0.01$, the velocity launch first crosses at 0.48835603, turns at 0.56406091 and returns at 0.64280676. The displacement launch gives 0.48839417, 0.56409910 and 0.64284508. The unseeded parent returns at 0.64284619. Relative to their matched constant-weight controls, the two matter-weighted waves return earlier by 0.00003609 and 0.0000001073 model-time units, respectively. The output contains the same timing comparisons at the half seed. For constant coupling the direct transfer vanishes at every interval. With zero seed the classical daughter stays zero throughout.

The solver retains the full coupled parent, daughter and Einstein equations with the matched constraint data from the first-passage experiment. Each integration stops at the first downward polar crossing. An angular derivative bound encloses the parent field on the full sphere at every trial stage before evaluation of the fitted potential, and the independent constraint monitors are evaluated at every accepted timestep and every physical event. Separate runs vary angular resolution and the time step. We compute each balance interval directly, compare adjoining intervals with their union and compare independent time quadratures. The nine primary configurations use angular degree 128, and all nine are repeated at degree 160 and with the maximum time step halved. Event times change by less than $6.1\times10^{-13}$ and integrated direct transfer changes by less than $3.4\times10^{-18}$ reference energy units. The primary normalized constraint and Raychaudhuri residuals stay below $1.4\times10^{-10}$. The parent and daughter energy ledgers close to $2.5\times10^{-12}$ and $1.1\times10^{-18}$, respectively. [Numerical data and controls](RETURN.json) · [Replay script](return_pass.py)

Run from the study root with `OPENBLAS_NUM_THREADS=1 python -IB coupled/return_pass.py --prior prior --output coupled/RETURN.json`.
