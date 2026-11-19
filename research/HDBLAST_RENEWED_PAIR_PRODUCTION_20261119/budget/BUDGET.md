# Who supplies the energy for the second pulse?

In the quiet-wall model with the original static compensation and comoving pair-energy cutoff 120, the late pulse adds **13.64794** units to the energy carried by the daughters. At the same readout the parent's energy correction decreases by **3.92214**. The scalar field that controls the gate does **9.97121** units of net work through the changing coupling and its local compensation. The parent's changing mass and the expansion account for the remaining signed work. These are coefficients of $y^2$ in proper energy per normalized comoving volume, on the prescribed background. [All measured components](BUDGET.json) · [Source and work calculation](budget.py)

## The matter Hamiltonian and the work equation

The calculation reuses the canonical covariance equations of the finite-bath study. For each retained parent momentum $q$, the bilinear interaction in the proper-time Hamiltonian is proportional to $y F(\phi)\chi\sum_jg_j Q_j$. The positive static completion adds $y^2 a F(\phi)^2 C_q\chi^2/2$, with $C_q=\sum_jg_j^2/\omega_j^2$ from the actual included oscillators. The bare calculation sets the completion to zero. The parent, pair bath and their correlations start at the same fork for every pulse history. The bath starts in its vacuum. The same couplings determine both its driven response and its quantum noise. [Prior Hamiltonian and source reconstruction](https://github.com/maldonado-research/HDblast/tree/49756c922760880a9f5da0e4da6f2d8cce36738a/research/HDBLAST_QUANTUM_MEMORY_20261119/source)

The coupling is either the original gate $A$, the later gate $B$, or their sum. Each gate is $\sin^4$ within its measured scalar-field support and zero outside. The inverse fine-quiet scalar clock calibrates gate $A$ at proper times 1.35–2.05 and gate $B$ at 2.055–2.145. Every other history uses those same scalar-field supports. The curvature and first derivative of both gates vanish at their endpoints. The evolution follows the parent, the correlations with each oscillator, and the emitted energy continuously from the original fork. The readout is shell time $s=H_0\tau=2.16$ on every wall.

For the order-$y^2$ energy density $\rho$, pressure $p$, scalar derivative $j$ and expansion rate $h=H/H_0$, the matter continuity equation becomes

$$
\frac{d(a^3\rho)}{ds}
=a^3\left[j\frac{d\phi}{ds}-3h p\right],
\qquad
j=j_{\rm parent}+j_{\rm int}+j_{\rm comp}.
$$

The parent term accounts for the changing parent mass. The interaction and static terms contain the derivative of the gate with respect to the physical scalar field. Together their integral gives signed gate work. The derivation differentiates the same Hamiltonian for the energy, pressure and scalar source. In particular, the interaction pressure vanishes and the static term has $p_{\rm comp}=-\rho_{\rm comp}/3$. Positive scalar work transfers energy into this matter sector. The prescribed-trajectory calculation does not specify the total energy available to the scalar or follow its depletion by a large transfer.

## The added energy and the changing force

Each entry compares $A+B$ with $A$ on the same background. The histories coincide until the second gate starts, and the interaction and compensation energies vanish again before the common readout. At cutoff 120 the original low band and the upper band each use 64 positive physical-frequency nodes for every parent mode.

| Wall, with static compensation | Daughter energy change | Parent energy change | Gate work | Mass work | Expansion work |
|---|---:|---:|---:|---:|---:|
| Fine quiet | +13.64794 | −3.92214 | +9.97121 | −0.05133 | −0.19408 |
| Fine plus | +13.66050 | −3.81516 | +10.08922 | −0.05006 | −0.19382 |
| Fine minus | +13.68378 | −4.38862 | +9.55317 | −0.05844 | −0.19957 |
| Fine far | +13.67071 | −3.92812 | +9.98988 | −0.05189 | −0.19540 |
| Coarse quiet | +13.65107 | −3.92538 | +9.97270 | −0.05182 | −0.19519 |

On fine quiet, the daughter energy increases from 22.27951 to 35.92745. Subtracting the isolated $B$ energy from that increment leaves an interference contribution of +0.002888. On the fine minus wall this contribution is −0.000500. These are small differences of positive totals. The direct $A-B$ evolution on quiet and minus also checks the sign reversal. The identity $T(A-B)=2T(A)+2T(B)-T(A+B)$ holds for all the order-$y^2$ density, pressure, scalar force and radiation histories, including their reported derivatives, to $1.1\times10^{-14}$ of each observable's maximum amplitude. [Sign and resolution checks](CONTROLS.json)

During the quiet late gate, the positive part of its signed work integrates to 13.56988 and the returned energy to 3.59867. The largest cumulative transfer from gate onset is 13.55137, while the net transfer is 9.97121. These values track energy exchanged with the scalar through the specified coupling and local term. They do not price the microscopic preparation of the time-dependent interaction.

In the bare alternative the same emitted daughter energy accompanies a parent energy change of −14.76710, gate work −0.74089, mass work −0.17850 and expansion work −0.19976. At this order the extra local potential changes the driven parent but leaves the leading daughter production unchanged. A particle count by itself therefore cannot determine which degree of freedom paid for the energy in this model.

## How much does the cutoff matter?

After the cutoff-120 experiment showed strong emission from the short pulse, we extended only the fine quiet *matter* calculation. Each new panel uses 64 positive frequency nodes. The static completion is recomputed for each actual spectral measure. The bare variant uses no local potential at every cutoff. The following values are again the changes induced by adding $B$, in proper energy units at $s=2.16$.

| Cutoff | Daughter energy change | Static parent change | Static gate work | Bare parent change | Bare gate work |
|---|---:|---:|---:|---:|---:|
| 120 | +13.64794 | −3.92214 | +9.97121 | −14.76710 | −0.74089 |
| 240 | +42.82915 | −21.11172 | +22.69989 | −65.00717 | −20.65796 |
| 480 | +44.12469 | −25.75967 | +19.56875 | −201.96820 | −154.48192 |

At 480 the positive static gate work is a residual of −154.48192 from the interaction and +174.05067 from the local potential. These changing terms make the energy attribution sensitive to the finite Hamiltonian. Doubling the time-step count from 4096 to 8192 changes the quiet 480 pair increment by $3.6\times10^{-7}$ and the static gate work by $5.1\times10^{-7}$ in absolute units. Extending the frequency band produces much larger changes. A physical continuum interpretation requires a specified renormalized model. The native-gravity experiment in this release uses cutoff 120.

Across all five primary backgrounds and both matter prescriptions, the largest pointwise continuity residual among the tested profiles is $4.0\times10^{-12}$ and the largest integrated proper-energy discrepancy is $2.4\times10^{-7}$ in absolute units. The separate Fourier-amplitude calculation agrees with the 120 Hamiltonian energies of $A$, $B$ and $A+B$ on all four fine walls within $3.3\times10^{-7}$. It agrees with all three quiet cutoff increments within $4.8\times10^{-7}$. [Independent emission calculation](../emission/EMISSION.md) · [Extended-cutoff budgets](CONTROLS.json)

## Replay

From the study root, run the following commands with the historical files present in their original repository paths. The temporary tables include the derivatives required by the gravity solver. The second command generates the independent negative-sign and resolution controls and the higher-cutoff quiet budgets.

~~~sh
OPENBLAS_NUM_THREADS=1 python -IB budget/budget.py --repo ../.. --tables /tmp/hdblast-repump-tables --jobs 2 --output /tmp/BUDGET.json
OPENBLAS_NUM_THREADS=1 python -IB budget/controls.py --repo ../.. --tables /tmp/hdblast-repump-tables --work /tmp/hdblast-repump-controls --jobs 2 --output /tmp/CONTROLS.json
~~~
