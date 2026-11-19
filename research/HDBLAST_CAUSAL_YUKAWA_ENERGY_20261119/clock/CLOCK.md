# Can the rolling shell scalar control the interaction?

Two different physical prescriptions can give exactly the same pair spectrum on the quiet shell. Their predictions separate when incident pulses change the scalar history. This exploratory comparison uses one calibration per numerical grid and applies it to each of the three pulse histories. [Study design recorded before the emission comparison](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-6998719638)

## Constructing a coupling on the quiet shell

Write $s=H_0\tau$ and let $f(s)$ be the original sine-to-the-fourth pulse on $[1.35,2.05]$. On the quiet shell the scalar $\phi_Q(s)$ increases strictly throughout the recorded history. Define the function

$$F(u)=f\!\left(\phi_Q^{-1}(u)\right) \qquad \text{for }\phi_Q(1.35)\le u\le\phi_Q(2.05),$$

and set $F(u)=0$ elsewhere. On the fine grid the onset, peak and turn-off occur at field values 0.338838, 0.771166 and 0.999332. The values at the two boundaries and the first three derivatives of the sine-to-the-fourth pulse vanish. The inverse across the sampled trajectory uses its native cubic Hermite interpolant. We check the minimum derivative on every cubic span by evaluating its endpoints and the interior extrema of its quadratic derivative. All eight full archived trajectories increase strictly, through $s\simeq2.168$, beyond the common readout at 2.10. A different future trajectory could return to the active field interval.

With a proper-time coupling the operator on each trajectory has the factor $y f(s)$. With the scalar coupling it has $yF(\phi_b(s))$. By construction these agree pointwise on the quiet trajectory. The quiet spectrum alone cannot identify what controls the coupling. Scalar-dependent decay couplings have precedents in modulated reheating. Here the form of $F$ is reconstructed from one HDblast history. [Dvali, Gruzinov and Zaldarriaga, Eqs. (2)–(4)](https://arxiv.org/abs/astro-ph/0303591)

## Outgoing energies on three independent pulse histories

The initial retained covariance, mass law, proper-time readout and parent momentum band are common to the two prescriptions. Each simulated numerator is the pulse-minus-quiet outgoing energy. The denominator is the quiet proper-time energy at the same cutoff and resolution. At $\Lambda=120$ the quiet fine-wall energy is $22.51056933$ in units of $y^2H_0^4$ per comoving volume for both prescriptions.

| Fine native wall | Proper-time coupling, % | Scalar coupling, % | Change, percentage points |
|---|---:|---:|---:|
| Nearby positive pulse | −0.147278 | −0.365832 | −0.218554 |
| Nearby negative pulse | +1.009110 | +1.603178 | +0.594068 |
| Distant compensated pulse | +0.129570 | +0.029501 | −0.100069 |

The coarse-wall changes are −0.218612, +0.594436 and −0.100060 percentage points. Raising the cutoff to 240 gives fine-wall changes of −0.217730, +0.594989 and −0.099235, respectively. At 480 the corresponding values are −0.217006, +0.595766 and −0.098526. Each sign persists across the six cutoff and wall-resolution combinations. The outgoing pulse contrasts distinguish these two fixed couplings on this set of simulated histories. Each cutoff defines a finite band of daughter energy. On the fine pulsed histories the scalar-clock energies still move by 0.000160–0.000175 in model units from cutoff 240 to 480. The quoted pulse contrasts therefore refer to their specified finite bands. [All energies, gates and resolution checks](CLOCK.json)

## Which field supplies the work?

For a prescribed scalar coupling on a given history, the chain rule gives $\dot f=F'(\phi_b)\dot\phi_b$. Let $B=\sum_jg_j x_j$ and $S_q=\sum_jg_j^2/\omega_j^2$. With the optional static counterterm, the interaction and counterterm contribute a signed scalar current, per unit proper volume, of

$$
 j_{F,q}=\frac{1}{a^3}\left[yF'(\phi_b)\langle\chi B\rangle+y^2aS_q F(\phi_b)F'(\phi_b)\langle\chi^2\rangle\right].
$$

The bare model omits the second term. Integrating the current with $a^3\dot\phi_b$ gives the part of the explicit proper-time work previously assigned to changing the coupling when the profiles coincide. The explicit scale-factor derivative of the counterterm contributes separately to metric work. In a dynamic source model this extra current must enter the scalar's junction equation and the trajectories must evolve again. Here both couplings are evaluated on the archived walls and only the outgoing radiation is recomputed. The two-clock comparison does not itself determine how much energy a physical bulk field loses.

## Numerical reproduction

The script reads the retained covariance and Hermite-interpolated scale factor and scalar from the same native fork and wall data as the budget solver. It evolves vacuum-normalized canonical scalar modes and conformal time together in proper time using DOP853. For each frequency, time integration evaluates both phase branches $\int ds\,f(s)\chi_q(s)e^{\pm i\omega\eta(s)}$ with composite Simpson quadrature, where $ds=a\,d\eta$. Gauss–Legendre quadrature integrates the inclusive Dirac measure over $q<\omega<\Lambda$ and weights by $\omega/a(2.10)$ per pair. Each cut uses two time and frequency resolutions. The largest energy change upon quadrature refinement across the stored cases is $4.8\times10^{-8}$. At $\Lambda=120$ the proper-time energies agree with the separately integrated covariance budgets within $5.1\times10^{-8}$ across the two grids and each arm. The quiet field-gated energy reproduces the proper-time value at all three cutoffs. Independent Linux runners differed by about $5\times10^{-7}$ percentage points on the small far-arm contrast. The replay compares clock contrasts to an absolute tolerance of $3\times10^{-6}$ percentage points, accounting for subtraction and the conversion from energy to percent.

```sh
OPENBLAS_NUM_THREADS=1 python -I -B clock/clock_test.py /path/to/previous/native --check
```
