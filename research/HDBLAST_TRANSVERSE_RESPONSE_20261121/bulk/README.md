# Could the incoming bulk wave supply the missing sideways force after particle birth?

At the moment the first charged daughters appear, suppose the wall has zero additional transverse field. We let their currents act with the most favorable possible sign and add a transverse incoming wave from the five-dimensional side. How strong must that incoming wave be for the transverse momentum distribution to agree with the thermal prediction within one percentage point at every cutoff?

At the 1% whole-work marker, the required incoming amplitude for an arbitrary waveform has a time-and-space RMS of at least **0.0423** in the dedicated interaction and **0.0255** in the shared interaction. The actual inward axial wave reaching the wall during these two short windows has RMS **0.000144**. If the extra incoming component is bounded at each time by a fixed multiple of the actual axial waveform, the lower bounds on that multiple are approximately **323 and 196**, respectively. These are conditions on the residual wave arriving *after* particle birth. The earlier original inward wave reached a peak of about 8.50 in these units. Transverse fields stored before the first daughters form are subject to a separate initial-field test. [Numerical data](INCOMING.json) · [Initial transverse fields and density caps](../global/CLUMPING.json)

| At 1% of whole-run absolute electrical work | Dedicated | Shared |
|---|---:|---:|
| Lower bound on peak incoming spatial RMS | 0.0469 | 0.0283 |
| Lower bound on spacetime RMS for an arbitrary waveform | 0.0423 | 0.0255 |
| Measured axial incoming RMS in the birth window | 0.0001442 | 0.0001435 |
| Lower bound on pointwise multiple of the measured axial waveform | 323 | 196 |
| Same multiple with the coarser native time step | 320 | 193 |
| Compulsory RMS impulse at the strongest waveform witness, in daughter mass units | 521.7 | 657.5 |
| Upper bound supplied by freely chosen daughter currents at that witness | 235.9 | 414.0 |
| Upper bound from an extra wave equal to the actual axial tail there | 0.885 | 1.244 |

The spacetime RMS bounds refer to the full interval from the first birth-left knot to the marker. The impulse numbers are RMS over all daughters present, with correlations and collisions granted favorable choices subject to the specified density cap. A disturbance redirected from the original wave must satisfy the same amplitude requirement if the original axial wall history is to be preserved. Redirection can change the axial field, the births and the density. The true jointly driven evolution then needs to be recalculated.

## Boundary source and measurement

The finite-collar solver evolves the wall variable $K X$ from the incoming characteristic $i$, outgoing $o=X-i$, local guide field $B=2C/a$, and daughter current $J$. Its exact axial equation in wall time is

$$
\frac{d}{dt}(a^2 KX) = 4a^2 i-2a^2 X -2K aB-a^2J .
$$

The original simulator is rerun with an instrument that directly records the incoming radial characteristic $i$ at the wall. The separate electric, magnetic, current and expansion histories reproduce its signed integral via the other three terms. The full native and halved-time canonical-field residuals are below $4\times10^{-8}$ of peak canonical field. The output also records the independent reconstruction error relative to the absolute *incident signal within each birth window*. Since incoming and reflected characteristics both carry energy, $\int i^2dt$ used to compare the waveform tail with the earlier pulse is a measure of incident characteristic strength. The net energy flux into the wall depends on both inward and outward waves. The late large axial wall field and the contemporaneous small incoming wave can therefore coexist. [Capture and checks](incoming.py)

## Bound on the incoming transverse characteristic

We apply the axial source normalization $+4a^2 i$ to an additional transverse characteristic. This is a local Maxwell boundary model for an inhomogeneous tangential disturbance on the recorded geometry. It retains positive wave energy and the curl exchange between wall electric and magnetic fields. Computing an actual five-dimensional perturbed wave solution satisfying the full boundary and gravitational constraints is a separate problem.

Use the volume-normalized spatial norm and write $\mathcal F(t)=a^2\|(\delta\mathbf E,\delta\mathbf B)\|_2$. Let $I(t)=\|\delta\mathbf i_{\rm in}\|_2$ in the added spatial modes. The gauge coefficient $K=2.54$ is constant throughout these birth windows. Integrated over the closed wall, the Maxwell curl terms cancel in the field energy derivative. Dropping the damping contribution enlarges the response. The local density cap $C_{\!n}=1$ used in the particle count and current bound gives

$$
\mathcal F(t)\le\frac{Q}{K}\int_{t_0}^t\frac{D(s)}{a(s)}\,ds+
\frac4K\int_{t_0}^t a(s)^2 I(s)\,ds .
$$

Here $D$ is the accepted comoving daughter rest and $Q$ the benchmark charge-to-mass ratio. The first term places every available particle current in the added modes with its most favorable sign. The bound on $I$ is additional to that generous current allocation. The modeled source and receiving current use the same field normalization. The measured reference axial component supplies its time profile and amplitude for comparison.

For daughters born after a cutoff $c$, define $g(s)=a(s)^{-1}\sqrt{D(s)/D(T)}$. The spatial RMS of particle-sampled force is at most $\sqrt2$ times the electric-plus-magnetic norm under the stated cap. The contribution of the boundary source to the population RMS impulse is bounded above by

$$
\int_{t_0}^T k_c(u) I(u)\,du,\qquad
k_c(u)=\frac{4\sqrt2 Q}{K a(T)}a(u)^2\int_{\max(c,u)}^Tg(s)\,ds.
$$

For each cutoff $c$ we subtract the maximum current-supplied impulse from the required thermal-quantile impulse and retain the positive difference $L_c$. A viable source requires $\int k_c I\,du\ge L_c$. Hölder’s inequality gives a necessary peak $\max_c L_c/\int k_c$, and Cauchy–Schwarz gives a necessary spacetime RMS $\max_c L_c/\sqrt{(T-t_0)\int k_c^2}$. If $I(t)\le\eta|i_{\rm axial}(t)|$, then $\eta\ge\max_c L_c/\int k_c |i_{\rm axial}|$. Positive cell bounds enlarge the kernels and freely credit the entire terminal bracketing cell to the source, while discarding that cell’s births from the required population. The finite thermal-quantile cost is a lower sum. An explicit dense kernel checks the accumulated-sum implementation. [Particle-count and self-current inequalities](../global/clumping.py)

At 5% of whole-run electrical work, the allowed self-current contribution alone is already high enough that this inequality yields a zero lower bound on extra incoming waves for a one-percentage-point test. Establishing that self currents actually generate sufficient impulses there requires a dynamical calculation. Density concentrations, earlier transverse fields, extra particle interactions or a changed production history alter the assumptions.
