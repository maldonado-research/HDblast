# Quantum-state differences from continuum retarded work

For two Gaussian parent states propagated on the same prescribed wall, the scalar commutator is independent of preparation. It drops out of the difference of their two-point functions. This makes the change in the quantum pair energy accessible from the causal response of real scalar quadratures, even though the absolute emission of the switched fork vacuum is nonzero.

Let $F=f a\chi_q=U+iV$ denote the switched free basis mode on the retained quiet wall. At the fork take $n=\langle a^\dagger a\rangle$ and $c=\langle aa\rangle$ in the instantaneous physical basis. Write $R$ for the retarded convolution at a fixed invariant subtraction scale, including the two local subtractions derived in [the dispersion calculation](MEMORY.md). With $W_{UU}=-\int U'R[U]d\eta$, $W_{VV}=-\int V'R[V]d\eta$ and $W_{UV}=-\int (U'R[V]+V'R[U])d\eta$, the excess over the matched fork vacuum at one momentum is

$$
\Delta E_{q,\mathrm{conf}}/y^2
=2[(n+\operatorname{Re}c)W_{UU}
 +(n-\operatorname{Re}c)W_{VV} - \operatorname{Im}c\,W_{UV}].
$$

The factors follow from $Fa+F^*a^\dagger=\sqrt{2}(U x - V p)$ and the symmetrized quadrature covariance difference $\Delta S_{xx}=n+\operatorname{Re}c$, $\Delta S_{pp}=n-\operatorname{Re}c$, $\Delta S_{xp}=\operatorname{Im}c$. The corresponding direct computation integrates $\omega K_q(\omega)[n(|A|^2+|B|^2)+2\operatorname{Re}(cAB^*)]$, with $A=\int F e^{i\omega\eta}d\eta$ and $B=\int F^* e^{i\omega\eta}d\eta$. The vacuum contribution is $\int\omega K_q|B|^2 d\omega$. These formulas share a prediction and use independent numerical integrations.

We evaluated all 24 native momentum nodes, weighted the result by the archived $q^2 dq/(2\pi^2)$ rule, and divided by $a_f=1.23570620$ to obtain physical energy per initial coordinate volume at proper time 2.10, per $y^2$. Every row uses the same expanding quiet wall and wide switch. The four labels refer to the archived initial covariances.

| State | Excess from continuum response, $\mu=2$ | Excess from continuum response, $\mu=4$ | Excess from direct pair amplitudes, $\omega\le400$ |
|---|---:|---:|---:|
| Retain | −0.589546814 | −0.589546814 | −0.589546814 |
| Erase | −0.432813222 | −0.432813222 | −0.432813222 |
| Reflect | −0.276079629 | −0.276079629 | −0.276079629 |
| Kick | −0.437255261 | −0.437255261 | −0.437255261 |

The direct absolute vacuum energy is $23.10012154$, giving a retained total of $22.51057473$. In this preparation the coherent contribution lowers the energy relative to the matched vacuum. That sign is compatible with positive absolute radiation energy. Across all four states the largest discrepancy between the infinite-cutoff retarded and direct amplitude results is $4.2\times10^{-11}$. The largest change from doubling the subtraction scale is $8.9\times10^{-12}$. The direct excess varies by less than $2.5\times10^{-9}$ on increasing the energy cutoff from 160 to 400. These results use the archived finite band $0<q<8$ for the parent momentum. They establish subtraction-scale independence of an integrated state-dependent daughter energy under the specified class of local subtractions. Partitioning the total budget among the parent, clock, expanding background and mass source is a separate problem.

Run `OPENBLAS_NUM_THREADS=1 python memory/state_difference.py` from the capsule, or supply the preceding archive directory with `--native`. [Numerical output](STATE_DIFFERENCE.json)
